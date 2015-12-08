/* Fitz Hugh Nagumo neuron model with periodic BCs on a 2D flat surface. ICs: stable state everywhere apart from
 * a small rectangle centred on the surface.
 *
 * The model equations are:
 *
 * u' = 3u - u^3 - v + D*Laplacian(u)
 * v' = epsilon*(u - beta)
 *
 * where u is the activator variable, v is the inhibitor variable.
 *

       /-/--\
     (@~@) ' )/\
 ___/--  ,` '\  >
(oo)__ _  ,` '}/
 ^^___/ '   ' \  _
       \ ' , , }/ \
        (,     )   }
        | , ' , \_/
        /   ,  ` \  __
  ,     ( ,  , `  )/  \
__)----/ ' (  ) ' |    }
==-___( , / '`/` , \__/
 '  , { _^ ',) ,  ' } __
  __)---  , _>  ,' ,\/  \
  ==-______/&        |   }            _/-\_
   '  { ' , ', ' '  , \_/            / ' , \_
      | ,       ' `    )          __| ' _ `  |
      |   ,(' ,')  ' ,  }-\       \/ ' ( ) ` \
       (` (  ,  ,)      \  }      {','/   \ ' }
       | ,( `   ,) '  ,` \/    __{ , {     {`,|
       \  ( , ', )  '  '  |    \/   /      (  )
     ___\( `,   {  ,  ' ' )/\__{ ', }       `'
   =-     (, ; (__  ,  ,  '` '   ,,'
   =_____-(   )   \__ `  '  , `  /
       ___\ ' ;)     `---_____,-'
     =-    '___/
     =_____---

 */


// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "arkode/arkode.h"            			// prototypes for ARKode fcts., consts.
#include "nvector/nvector_parallel.h" 			// parallel N_Vector types, fcts., macros
#include "arkode/arkode_pcg.h"        			// prototype for ARKPcg solver
#include "sundials/sundials_types.h"  			// def. of type 'realtype'
#include <sundials/sundials_math.h>   			// math macros
#include "mpi.h"                      			// MPI header file
#include <boost/property_tree/ptree.hpp>		//
#include <boost/property_tree/ini_parser.hpp>
using namespace std;

// accessor macro between (x,y) location and 1D NVector array
#define IDX(x,y) (NVARS*(x) + NVARS*(y)*udata->nxl)

// Constants
#define PI RCONST(3.1415926535897932)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)
#define THREE RCONST(3.0)
#define SURFACEWIDTH RCONST(20.0) 	// Width of the flat surface

// System parameters
#define EPSILON 0.36			// Time scale separation of the two variables

// Number of variables in the system
#define NVARS 2

// Initialise boundary values
double XMIN, XMAX, YMIN, YMAX;

// Initialise parameters, found in ini file
double DIFF = 0.0;					// Diffusion parameter - default is 0.12
double BETA = 0.0;					// Bifurcation parameter - system is oscillatory for BETA < 1, stable for BETA > 1
double SURFACELENGTH = 0.0;	 		// Length of the flat surface, usually 80 or 40 to match with the corresponding torus
double WAVELENGTH = 0.0;			// Initial wave segment length as a percentage of total length of torus (phi)
double WAVEWIDTH = 0.0;				// Initial wave segment width as a percentage of total width of torus (theta)
int WAVEINSIDE =  0;				// Bool/int for whether the initial wave is centered on the inside of the torus (true=1) or outside (false=0)
int OUTPUT_TIMESTEP = 0; 			// Number of timesteps to output to file
double TBOUNDARY = 0.0;				// Time to turn off the absorbing boundary at phi = 0 (to eliminate backwards travelling waves) - set to 0 for no absorbing boundary
double TFINAL = 0.0;				// Time to run simulation
int NX = 0;							// Mesh size in theta direction
double BETAMIN = 0;					// Minimum beta value when varyBeta = 1
double BETAMAX = 0;					// Maximum beta value when varyBeta = 1
int INCLUDEALLVARS = 0;				// Bool/int for whether we write all variables to file (true=1) or only the main activator variable u (false=0)
int VARYBETA = 0;					// Bool/int for whether to vary beta over the surface of the torus (true=1) or keep it constant (false=0)

// user data structure
typedef struct {
	long int nx;          // global number of x grid points
	long int ny;          // global number of y grid points
	long int is;          // global x indices of this subdomain
	long int ie;
	long int js;          // global y indices of this subdomain
	long int je;
	long int nxl;         // local number of x grid points
	long int nyl;         // local number of y grid points
	realtype dx;          // x-directional mesh spacing
	realtype dy;          // y-directional mesh spacing
	realtype Diff;        // diffusion coefficient for equation 1
	MPI_Comm comm;        // communicator object
	int rank;             // MPI process ID/rank
	int nprocs;           // total number of MPI processes
	realtype *Erecv;      // receive buffers for neighbour exchange
	realtype *Wrecv;
	realtype *Nrecv;
	realtype *Srecv;
	realtype *Esend;      // send buffers for neighbour exchange
	realtype *Wsend;
	realtype *Nsend;
	realtype *Ssend;
} UserData;

// User-supplied Function Called by the Solver

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// Private functions

//    checks function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

//    sets default values into UserData structure
static int InitUserData(UserData *udata);

//    sets up parallel decomposition
static int SetupDecomp(UserData *udata);

//    performs neighbor exchange
static int Exchange(N_Vector y, UserData *udata);

//    frees memory allocated within UserData
static int FreeUserData(UserData *udata);

// printf("[%d]---%s:%d\n", rank, __FILE__, __LINE__);

// Main Program
int main(int argc, char* argv[])
{
	// Output error if user does not specify parameter ini file
	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <Config file path>";
		exit(EXIT_FAILURE);
	}

	// Obtain model parameters from ini file
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(argv[1], pt);
	DIFF = pt.get<double>("Parameters.diffusion");
	BETA = pt.get<double>("Parameters.beta");
	SURFACELENGTH = pt.get<double>("Parameters.majorCirc");		// This is the length of the flat surface
	WAVELENGTH = pt.get<double>("Parameters.waveLength");
	WAVEWIDTH = pt.get<double>("Parameters.waveWidth");
	OUTPUT_TIMESTEP = pt.get<int>("Parameters.outputTimestep");
	TBOUNDARY = pt.get<double>("Parameters.tBoundary");
	TFINAL = pt.get<double>("Parameters.tFinal");
	NX = pt.get<int>("Parameters.thetaMesh");				// This is the mesh of the X direction
	BETAMIN = pt.get<double>("Parameters.betaMin");
	BETAMAX = pt.get<double>("Parameters.betaMax");
	INCLUDEALLVARS = pt.get<int>("System.includeAllVars");
	VARYBETA = pt.get<int>("System.varyBeta");

	XMIN = 0.0;				                // grid boundaries in theta
	XMAX = SURFACEWIDTH - XMIN;
	YMIN = 0.0;			    	    		// grid boundaries in phi
	YMAX = SURFACELENGTH - YMIN;

	// general problem parameters
	realtype T0 = RCONST(0.0);   		// initial time
	realtype Tf = TFINAL;     			// final time
	int Nt = OUTPUT_TIMESTEP;     		// total number of output times
	long int lengthWidthRatio = SURFACELENGTH/SURFACEWIDTH;
	long int nx = NX;             		// spatial mesh size
	long int ny = NX*lengthWidthRatio;
	int WaveInside = WAVEINSIDE;
	realtype xx, yy;					// real x,y values
	realtype Diff = DIFF;				// Diffusion coefficient for eq 1
	realtype rtol = 1.e-5;       		// relative and absolute tolerances
	realtype atol = 1.e-10;
	realtype WaveLength = (YMAX-YMIN)*WAVELENGTH;
	realtype WaveWidth = (XMAX-XMIN)*WAVEWIDTH;
	realtype WaveMidpoint;
	realtype WaveXMIN;
	realtype WaveXMAX;

	UserData *udata = NULL;
	realtype *data;
	realtype *ydata;
	long int N, Ntot, i, j;

	// general problem variables
	int flag;                      // reusable error-checking flag
	int rank;                      // MPI process ID
	N_Vector y = NULL;             // empty vector for storing solution
	void *arkode_mem = NULL;       // empty ARKode memory structure

	// initialize MPI
	flag = MPI_Init(&argc, &argv);
	if (check_flag(&flag, "MPI_Init", 1)) return 1;
	flag = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

	// allocate and fill udata structure
	udata = new UserData;
	flag = InitUserData(udata);
	if (check_flag(&flag, "InitUserData", 1)) return 1;

	udata->rank = rank;
	udata->nx = nx;
	udata->ny = ny;
	udata->Diff = Diff;
	udata->dx = (XMAX-XMIN)/(1.0*nx-1.0);   // x mesh spacing
	udata->dy = (YMAX-YMIN)/(1.0*ny-1.0);   // y mesh spacing

	// Set up parallel decomposition **********************************************************************
	flag = SetupDecomp(udata);
	if (check_flag(&flag, "SetupDecomp", 1)) return 1;

	// Find stable state of ODE model dependent on beta
	// The model is simple so the solution can be found analytically:
	double Us, Vs;	// Stable states of U, V respectively
	Us = -BETA;
	Vs = BETA*BETA*BETA - 3*BETA;

	// Initial problem output
	bool outproc = (udata->rank == 0);
	if (outproc) {
		cout << "\n2D FHN model PDE problem on a torus:\n";
		cout << "   nprocs = " << udata->nprocs << "\n";
		cout << "   nx = " << udata->nx << "\n";
		cout << "   ny = " << udata->ny << "\n";
		cout << "   nxl = " << udata->nxl << "\n";
		cout << "   nyl = " << udata->nyl << "\n";
		cout << "   Diff = " << udata->Diff << "\n";
		cout << "   Tfinal = " << TFINAL << "\n";
		cout << "   Output timesteps = " << OUTPUT_TIMESTEP << "\n";
		cout << "   Surface length = " << SURFACELENGTH << "\n";
		cout << "   Absorbing boundary turn off time = " << TBOUNDARY << "\n";
		cout << "   Wavelength = " << WAVELENGTH << "\%\n";
		cout << "   Wavewidth = " << WAVEWIDTH << "\%\n";
		cout << "   rtol = " << rtol << "\n";
		cout << "   atol = " << atol << "\n";
		cout << "   Include all variables in output = " << INCLUDEALLVARS << "\n";
		if (VARYBETA == 0)
		{
			cout << "   Beta = " << BETA << "\n";
			cout << "   Stable state values: U = " << Us << ", V = " << Vs << "\n\n";
		}
		else
		{
			cout << "   Beta varied over torus\n\n";
		}
	}

	// Initialize data structures
	N = NVARS*(udata->nxl)*(udata->nyl);
	Ntot = NVARS*nx*ny;
	y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution
	if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;

	WaveMidpoint = SURFACEWIDTH/2.0;
	WaveXMIN  = WaveMidpoint - WaveWidth/2.0;
	WaveXMAX  = WaveMidpoint + WaveWidth/2.0;

	// Set initial conditions - these are model dependent.
	ydata = N_VGetArrayPointer(y);
	for (j=0; j<udata->nyl; j++)
	{
		yy = YMIN + (udata->js+j)*(udata->dy);						// Actual x values

		for (i=0; i<udata->nxl; i++)
		{
			xx = XMIN + (udata->is+i)*(udata->dx);					// Actual x values

//			if (VARYBETA == 1)
//			{
//				ydata[IDX(i,j)] = 1;
//				ydata[IDX(i,j) + 1] = 1;
//			}
//			else
			{
				// Set initial wave segment
				if ( xx >= WaveXMIN && xx <= WaveXMAX && yy >= WaveLength && yy <= (2.0*WaveLength) )
				{
					// Set perturbed wave segment to higher initial values
					ydata[IDX(i,j)] = Us + 2;				// u
					ydata[IDX(i,j) + 1] = Vs + 1.5;			// v
				}
				else
				{
					// Set rest of area to stable u,v
					ydata[IDX(i,j)] = Us;					// u
					ydata[IDX(i,j) + 1] = Vs;				// v

				}
			}
		}
	}

	arkode_mem = ARKodeCreate();                       // Create the solver memory
	if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;

	/* Call ARKodeInit to initialize the integrator memory and specify the
		 right-hand side function in y'=f(t,y), the inital time T0, and
		 the initial dependent variable vector y.  */
	flag = ARKodeInit(arkode_mem, f, NULL, T0, y);
	if (check_flag(&flag, "ARKodeInit", 1)) return 1;

	flag = ARKodeSStolerances(arkode_mem, rtol, atol);      // Specify tolerances
	if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

	// Set routines
	flag = ARKodeSetUserData(arkode_mem, (void *) udata);
	if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;

	flag = ARKodeSetMaxNumSteps(arkode_mem, 200000);
	if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return(1);

	// Each processor outputs subdomain information
	char outname[100];
	sprintf(outname, "FHNmodel_flat_subdomain.%03i.txt", udata->rank);
	FILE *UFID = fopen(outname,"w");
	fprintf(UFID, "%li  %li  %li  %li  %li  %li %f %f %f\n",
			udata->nx, udata->ny, udata->is, udata->ie, udata->js, udata->je, XMIN, XMAX, TFINAL);
	fclose(UFID);

	ydata = N_VGetArrayPointer(y);

	sprintf(outname, "FHNmodel_flat_u.%03i.txt", udata->rank);
	UFID = fopen(outname, "w");

	sprintf(outname, "FHNmodel_flat_v.%03i.txt", udata->rank);
	FILE *UFID2 = fopen(outname, "w");


	// Write initial conditions to files, one for each variable in each subdomain
	for (j=0; j<udata->nyl; j++)
	{
		for (i=0; i<udata->nxl; i++)
		{
			fprintf(UFID," %.16e", ydata[IDX(i,j)]);

			if (INCLUDEALLVARS == 1)
			{
				fprintf(UFID2," %.16e", ydata[IDX(i,j) + 1]);
			}
		}
	}
	fprintf(UFID,"\n");

	if (INCLUDEALLVARS == 1)
	{
	fprintf(UFID2,"\n");
	}


	/* Main time-stepping loop: calls ARKode to perform the integration, then
		 prints results.  Stops when the final time has been reached */
	realtype t = T0;
	realtype dTout = (Tf-T0)/Nt;
	realtype tout = T0+dTout;

	int iout;
	for (iout=0; iout<Nt; iout++)
	{

		flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);         // call integrator
		if (check_flag(&flag, "ARKode", 1)) break;

		if (flag >= 0)
		{                                            // successful solve: update output time
			tout += dTout;
			tout = (tout > Tf) ? Tf : tout;
		} else
		{                                                    // unsuccessful solve: break
			if (outproc)
				cerr << "Solver failure, stopping integration\n";
			break;
		}

		// output results to disk
		for (j=0; j<udata->nyl; j++)
		{
			for (i=0; i<udata->nxl; i++)
			{
				fprintf(UFID," %.16e", ydata[IDX(i,j)]);

				if (INCLUDEALLVARS == 1)
				{
				fprintf(UFID2," %.16e", ydata[IDX(i,j) + 1]);
				}
			}
		}
		fprintf(UFID,"\n");

		if (INCLUDEALLVARS == 1)
		{
		fprintf(UFID2,"\n");
		}

		// Output progress
		if (outproc)
		{
			if (iout > 0)
			{
				// Rewrite previous line
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			}
			printf("   %3d /%3d done", iout+1, Nt);
			fflush(stdout);
		}
	}
	if (outproc)  cout << "   ----------------------\n";
	fclose(UFID);

	if (INCLUDEALLVARS == 1)
	{
	fclose(UFID2);
	}

	// Clean up and return with successful completion
	N_VDestroy_Parallel(y);        // Free vectors
	FreeUserData(udata);         // Free user data
	delete udata;
	ARKodeFree(&arkode_mem);     // Free integrator memory
	flag = MPI_Finalize();       // Finalize MPI

	/* ************************************************************************************** */

	return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y). ydot = f(t,y)
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_VConst(0.0, ydot);                           // Initialize ydot to zero
	UserData *udata = (UserData *) user_data;      // access problem data
	long int nxl = udata->nxl;                     // set variable shortcuts
	long int nyl = udata->nyl;
	realtype xx, yy;
	realtype Diff = udata->Diff;
	realtype dx = udata->dx;
	realtype dy = udata->dy;
	realtype *yarray, *ydotarray;
	yarray = NV_DATA_P(y);           // access data arrays
	ydotarray = NV_DATA_P(ydot);

	// Exchange boundary data with neighbors
	int ierr = Exchange(y, udata);
	if (check_flag(&ierr, "Exchange", 1)) return -1;

	// Add diffusion term

	realtype cu1 = Diff/dx/dx;			// D/delx^2
	realtype cu2 = Diff/dy/dy;			// D/dely^2
	realtype cu3 = -TWO*(cu1 + cu2);
	long int i, j;
	for (j=1; j<nyl-1; j++)
	{
		for (i=1; i<nxl-1; i++)
		{
			// Fill in diffusion for u variable
			ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
							 + cu2*(yarray[IDX(i,j-1)] + yarray[IDX(i,j+1)])
							 + cu3*yarray[IDX(i,j)];
		}
	}
	// iterate over subdomain boundaries
	// West face
	i=0;
	for (j=1; j<nyl-1; j++)
	{
		ydotarray[IDX(i,j)] = cu1*(udata->Wrecv[NVARS*j]   + yarray[IDX(i+1,j)])
						 + cu2*(yarray[IDX(i,j-1)] + yarray[IDX(i,j+1)])
						 + cu3*yarray[IDX(i,j)];
	}
	// East face
	i=nxl-1;

	for (j=1; j<nyl-1; j++)
	{
		ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
						 + cu2*(yarray[IDX(i,j-1)] + yarray[IDX(i,j+1)])
						 + cu3*yarray[IDX(i,j)];
	}
	// South face: absorbing boundary at phi = 0, if js = 0 and time < TBOUNDARY, so that no backwards travelling waves occur
	j=0;
	for (i=1; i<nxl-1; i++)
	{
		ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
						 + cu2*(udata->Srecv[NVARS*i]   + yarray[IDX(i,j+1)])
						 + cu3*yarray[IDX(i,j)];
	}

	// North face
	j=nyl-1;
	for (i=1; i<nxl-1; i++)
	{
		ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
							 + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
							 + cu3*yarray[IDX(i,j)];
	}

	// South-West corner
	i = 0;
	j = 0;
	ydotarray[IDX(i,j)] = cu1*(udata->Wrecv[NVARS*j] + yarray[IDX(i+1,j)])
						   + cu2*(udata->Srecv[NVARS*i] + yarray[IDX(i,j+1)])
						   + cu3*yarray[IDX(i,j)];

	// North-West corner
	i = 0;
	j = nyl-1;
	ydotarray[IDX(i,j)] = cu1*(udata->Wrecv[NVARS*j]   + yarray[IDX(i+1,j)])
						   + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
						   + cu3*yarray[IDX(i,j)];

	// South-East corner
	i = nxl-1;
	j = 0;
	ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
						   + cu2*(udata->Srecv[NVARS*i]   + yarray[IDX(i,j+1)])
						   + cu3*yarray[IDX(i,j)];


	// North-East corner
	i = nxl-1;
	j = nyl-1;
	ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
					   + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
					   + cu3*yarray[IDX(i,j)];

	realtype b;

	// Add other terms in equations
	for (j=0; j<nyl; j++)
	{
		yy = YMIN + (udata->js+j)*(udata->dy);

		if (VARYBETA == 0)
		{
			b = BETA;
		}
		else
		{
			b = BETAMIN + yy*(BETAMAX - BETAMIN)/(YMAX - YMIN);
		}

		for (i=0; i<nxl; i++)
		{
			realtype u = yarray[IDX(i,j)];
			realtype v = yarray[IDX(i,j)+1];

			// If we are on the north or south boundary of the entire domain and t<TBOUNDARY, set u_t = 0 to simulate
			// Dirichlet boundary conditions with values equal to the initial conditions

			// North boundary of the domain
			if (udata->je == udata->ny-1 && t<TBOUNDARY && j == nyl-1)
			{
				ydotarray[IDX(i,j)] = 0; 	// u
				ydotarray[IDX(i,j)+1] = 0;  // v
			}
			// South boundary of the domain
			else if (udata->js == 0 && t<TBOUNDARY && j == 0)
			{
				ydotarray[IDX(i,j)] = 0; 	// u
				ydotarray[IDX(i,j)+1] = 0;  // v
			}
			else
			{
				// u variable: du/dt = 3u - u^3 - v + Diff
				ydotarray[IDX(i,j)] += 3.0*u - (u*u*u) - v;

				// v variable: dv/dt = eps(u + beta)
				ydotarray[IDX(i,j)+1] += EPSILON*(u + b);
			}
		}
	}

	return 0;                                      // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
 */
static int check_flag(void *flagvalue, const string funcname, int opt)
{
	int *errflag;

	// Check if SUNDIALS function returned NULL pointer - no memory allocated
	if (opt == 0 && flagvalue == NULL) {
		cerr << "\nSUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
		return 1; }

	// Check if flag < 0
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
			return 1;
		}
	}

	// Check if function returned NULL pointer - no memory allocated
	else if (opt == 2 && flagvalue == NULL) {
		cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
		return 1; }

	return 0;
}

// Set up parallel decomposition
static int SetupDecomp(UserData *udata)
{
	// check that this has not been called before
	if (udata->Erecv != NULL || udata->Wrecv != NULL ||
			udata->Srecv != NULL || udata->Nrecv != NULL) {
		cerr << "SetupDecomp warning: parallel decomposition already set up\n";
		return 1;
	}

	// get suggested parallel decomposition
	int ierr, dims[] = {0, 0};
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &(udata->nprocs));
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Comm_size = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Dims_create(udata->nprocs, 2, dims);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Dims_create = " << ierr << "\n";
		return -1;
	}

	// set up 2D Cartesian communicator
	int periods[] = {1,1};		// {1,1} for periodic BCs, {0,0} for non-periodic
	ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(udata->comm));
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Cart_create = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Comm_rank(udata->comm, &(udata->rank));
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Comm_rank = " << ierr << "\n";
		return -1;
	}

	// determine local extents
	int coords[2];
	ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Cart_get = " << ierr << "\n";
		return -1;
	}
	udata->is = (udata->nx)*(coords[0])/(dims[0]);
	udata->ie = (udata->nx)*(coords[0]+1)/(dims[0])-1;
	udata->js = (udata->ny)*(coords[1])/(dims[1]);
	udata->je = (udata->ny)*(coords[1]+1)/(dims[1])-1;
	udata->nxl = (udata->ie)-(udata->is)+1;
	udata->nyl = (udata->je)-(udata->js)+1;

	// determine if I have neighbors, and allocate exchange buffers

	udata->Wrecv = new realtype[NVARS*(udata->nyl)];
	udata->Wsend = new realtype[NVARS*(udata->nyl)];

	udata->Erecv = new realtype[NVARS*(udata->nyl)];
	udata->Esend = new realtype[NVARS*(udata->nyl)];

	udata->Srecv = new realtype[NVARS*(udata->nxl)];
	udata->Ssend = new realtype[NVARS*(udata->nxl)];

	udata->Nrecv = new realtype[NVARS*(udata->nxl)];
	udata->Nsend = new realtype[NVARS*(udata->nxl)];

	return 0;     // return with success flag
}

// Perform neighbor exchange
static int Exchange(N_Vector y, UserData *udata)
{
	// local variables
	MPI_Request reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
	MPI_Status stat;
	int ierr, i, ipW=-1, ipE=-1, ipS=-1, ipN=-1;			// ip* are the ranks of neighbours
	int coords[2], dims[2], periods[2];
	int nyl = udata->nyl;
	int nxl = udata->nxl;

	// access data array
	realtype *Y = N_VGetArrayPointer(y);
	if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

	// MPI equivalent of realtype type
#if defined(SUNDIALS_SINGLE_PRECISION)
#define REALTYPE_MPI_TYPE MPI_FLOAT
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define REALTYPE_MPI_TYPE MPI_DOUBLE
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define REALTYPE_MPI_TYPE MPI_LONG_DOUBLE
#endif

	// MPI neighborhood information - get ranks of neighbours (output ip*)
	ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Cart_get = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Cart_shift(udata->comm, 0, 1, &ipE, &ipW);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Cart_shift = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Cart_shift(udata->comm, 1, 1, &ipS, &ipN);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Cart_shift = " << ierr << "\n";
		return -1;
	}

	// Test to see if neighbours are correct
	//printf("Rank %d has neighbours N:%d, S:%d, E:%d, W:%d and coordinates (%d, %d)\n", udata->rank, ipN, ipS, ipE, ipW, coords[0], coords[1]);

	/* Open Irecv buffers
	 *  Wrecv: place for received data to go
	 *  NVARS*nyl: length of data
	 *  ipW: rank of process to receive from
	 */
	ierr = MPI_Irecv(udata->Wrecv, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipW,
			MPI_ANY_TAG, udata->comm, &reqRW);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Irecv = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Irecv(udata->Erecv, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipE,
			MPI_ANY_TAG, udata->comm, &reqRE);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Irecv = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Irecv(udata->Srecv, NVARS*(udata->nxl), REALTYPE_MPI_TYPE, ipS,
			MPI_ANY_TAG, udata->comm, &reqRS);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Irecv = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Irecv(udata->Nrecv, NVARS*(udata->nxl), REALTYPE_MPI_TYPE, ipN,
			MPI_ANY_TAG, udata->comm, &reqRN);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Irecv = " << ierr << "\n";
		return -1;
	}

	// Send data
	for (i=0; i<nyl; i++)
	{
		udata->Wsend[2*i] = Y[IDX(nxl-1,i)];
		udata->Wsend[2*i + 1] = Y[IDX(nxl-1,i)+1];
	}
	ierr = MPI_Isend(udata->Wsend, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipW, 0,
			udata->comm, &reqSW);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Isend = " << ierr << "\n";
		return -1;
	}

	for (i=0; i<nyl; i++)
	{
		udata->Esend[2*i] = Y[IDX(0,i)];
		udata->Esend[2*i + 1] = Y[IDX(0,i)+1];
	}
	ierr = MPI_Isend(udata->Esend, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipE, 1,
			udata->comm, &reqSE);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Isend = " << ierr << "\n";
		return -1;
	}

	for (i=0; i<nxl; i++)
	{
		udata->Ssend[2*i] = Y[IDX(i,nyl-1)];
		udata->Ssend[2*i + 1] = Y[IDX(i,nyl-1)+1];
	}
	ierr = MPI_Isend(udata->Ssend, NVARS*(udata->nxl), REALTYPE_MPI_TYPE, ipS, 2,
			udata->comm, &reqSS);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Isend = " << ierr << "\n";
		return -1;
	}

	for (i=0; i<nxl; i++)
	{
		udata->Nsend[2*i] = Y[IDX(i,0)];
		udata->Nsend[2*i + 1] = Y[IDX(i,0)+1];
	}
	ierr = MPI_Isend(udata->Nsend, NVARS*(udata->nxl), REALTYPE_MPI_TYPE, ipN, 3,
			udata->comm, &reqSN);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Isend = " << ierr << "\n";
		return -1;
	}


	// wait for messages to finish
	ierr = MPI_Wait(&reqRW, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Wait(&reqSW, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Wait(&reqRE, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Wait(&reqSE, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Wait(&reqRS, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Wait(&reqSS, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}

	ierr = MPI_Wait(&reqRN, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}
	ierr = MPI_Wait(&reqSN, &stat);
	if (ierr != MPI_SUCCESS) {
		cerr << "Error in MPI_Wait = " << ierr << "\n";
		return -1;
	}


	return 0;     // return with success flag
}

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
	udata->nx = 0;
	udata->ny = 0;
	udata->is = 0;
	udata->ie = 0;
	udata->js = 0;
	udata->je = 0;
	udata->nxl = 0;
	udata->nyl = 0;
	udata->dx = 0.0;
	udata->dy = 0.0;
	udata->Diff = 0.0;
	udata->comm = MPI_COMM_WORLD;
	udata->rank = 0;
	udata->nprocs = 0;
	udata->Erecv = NULL;
	udata->Wrecv = NULL;
	udata->Nrecv = NULL;
	udata->Srecv = NULL;
	udata->Esend = NULL;
	udata->Wsend = NULL;
	udata->Nsend = NULL;
	udata->Ssend = NULL;

	return 0;     // return with success flag
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
	// free exchange buffers
	if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
	if (udata->Wsend != NULL)  delete[] udata->Wsend;
	if (udata->Erecv != NULL)  delete[] udata->Erecv;
	if (udata->Esend != NULL)  delete[] udata->Esend;
	if (udata->Srecv != NULL)  delete[] udata->Srecv;
	if (udata->Ssend != NULL)  delete[] udata->Ssend;
	if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
	if (udata->Nsend != NULL)  delete[] udata->Nsend;

	return 0;     // return with success flag
}


//---- end of file ----




