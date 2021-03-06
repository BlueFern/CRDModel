/* Goldbeter Calcium Cell model with periodic BCs on a 2D torus.
 *
 * The model equations are:
 *
 * Z' = v0 + v1*BETA - v2 + v3 + kf*Y - k*Z + D*Laplacian(Z)
 * Y' = v2 + v3 + kf*Y
 *
 * where Z is the calcium concentration in the SMC cytosol, Y is the calcium concentration in the SMC stores.
 *
 * The Laplacian is different on a torus from a flat surface (google it)
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

// System parameters
#define v0 1.0
#define k 10.0
#define kf 1.0
#define v1 7.3
#define VM2 65.0
#define VM3 500.0
#define K2 1.0
#define KR 2.0
#define KA 0.9
#define m 2.0
#define n 2.0
#define p 4.0

// Number of variables in the system
#define NVARS 2

// Initialise boundary values
#define XMIN         RCONST(0.0)                // grid boundaries in theta
#define XMAX         RCONST(2.0*PI)
#define YMIN         RCONST(0.0)        		// grid boundaries in phi
#define YMAX         RCONST(2.0*PI)

// Initialise parameters, found in ini file
double DIFF = 0.0;					// Diffusion parameter - default is 0.12
double BETA = 0.0;					// Bifurcation parameter - system is oscillatory for BETA < 1, stable for BETA > 1
double MAJOR_CIRC = 0.0;		 	// Major circumference of the torus
double MINOR_CIRC = 0.0;		 	// Minor circumference of the torus
double WAVE_LENGTH = 0.0;			// Initial wave segment length as a percentage of total length of surface (phi)
double WAVE_WIDTH = 0.0;			// Initial wave segment width as a percentage of total width of surface (theta)
int OUTPUT_TIMESTEP = 0; 			// Number of timesteps to output to file
int WAVE_INSIDE = 0; 				// Int for whether the initial perturbation is centered on the inside (1) or outside (0)
double T_BOUNDARY = 0.0;			// Time to turn off the absorbing boundary at phi = 0 (to eliminate backwards travelling waves) - set to 0 for no absorbing boundary
double T_FINAL = 0.0;				// Time to run simulation
int X_MESH = 0;						// Mesh size in theta direction
double BETA_MIN = 0;				// Minimum beta value when it's spatially varied
double BETA_MAX = 0;				// Maximum beta value when it's spatially varied
int INCLUDE_ALL_VARS = 0;			// Int for whether we write all variables to file (1) or only the main variable Z (0)
int VARY_BETA = 0;					// Int for whether to vary beta over the surface of the torus (1) or keep it spatially constant (0)
int JUST_DIFFUSION = 0;				// Int for whether to have no reaction terms (1) to simulate the diffusion equation, or normal (0)
int IC_TYPE = 0;					// Int for the initial conditions when beta is spatially varied. 0: homogeneous ICs, 1: initial perturbation, 2: random ICs.

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
	double Zs;			  // Steady state of Z variable dependent on beta (Goldbeter model dependent)
	double Ys;			  // Steady state of Y variable dependent on beta (Goldbeter model dependent)
	realtype *Erecv;      // receive buffers for neighbour exchange
	realtype *Wrecv;
	realtype *Nrecv;
	realtype *Srecv;
	realtype *Esend;      // send buffers for neighbour exchange
	realtype *Wsend;
	realtype *Nsend;
	realtype *Ssend;
	realtype R;
	realtype r;
} UserData;

// User-supplied Function Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// Private functions shown later
//    checks function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

//    sets default values into UserData structure
static int init_user_data(UserData *udata);

//    sets up parallel decomposition
static int setup_decomp(UserData *udata);

//    performs neighbor exchange
static int exchange(N_Vector y, UserData *udata);

//    frees memory allocated within UserData
static int free_user_data(UserData *udata);

// Line useful for debugging
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
	MAJOR_CIRC = pt.get<double>("Parameters.surfaceLength");
	MINOR_CIRC = pt.get<double>("Parameters.surfaceWidth");
	WAVE_LENGTH = pt.get<double>("Parameters.waveLength");
	WAVE_WIDTH = pt.get<double>("Parameters.waveWidth");
	WAVE_INSIDE = pt.get<int>("Parameters.waveInside");
	OUTPUT_TIMESTEP = pt.get<int>("Parameters.outputTimestep");
	T_BOUNDARY = pt.get<double>("Parameters.tBoundary");
	T_FINAL = pt.get<double>("Parameters.tFinal");
	X_MESH = pt.get<int>("Parameters.xMesh");
	INCLUDE_ALL_VARS = pt.get<int>("System.includeAllVars");
	VARY_BETA = pt.get<int>("System.varyBeta");
	JUST_DIFFUSION = pt.get<int>("System.justDiffusion");

	// Time variables
	time_t start_t = 0;
	time_t end_t = 0;
	double total_t = 0;
	double eta = 0;

	// Record the time at the start of program
	time(&start_t);

	// general problem parameters
	realtype T0 = RCONST(0.0);   		// initial time
	realtype Tf = T_FINAL;     			// final time
	int Nt = OUTPUT_TIMESTEP;     		// total number of output times
	realtype r = MINOR_CIRC / (2.0*PI); 	// Minor radius
	realtype R = MAJOR_CIRC / (2.0*PI); 	// Major radius
	realtype radiusRatio = R / r;
	long int nx = X_MESH;             		// spatial mesh size
	long int ny = X_MESH * (radiusRatio);
	int WaveInside = WAVE_INSIDE;
	realtype xx, yy;					// real x,y values
	realtype Diff = DIFF;				// Diffusion coefficient for eq 1
	realtype rtol = 1.e-5;       		// relative and absolute tolerances
	realtype atol = 1.e-10;
	realtype wave_length = (YMAX - YMIN) * WAVE_LENGTH;
	realtype wave_width = (XMAX - XMIN) * WAVE_WIDTH;
	realtype wave_midpoint;
	realtype wave_XMIN;
	realtype wave_XMAX;

	UserData *udata = NULL;
	realtype *data;
	realtype *ydata;
	long int N, Ntot, i, j;

	// general problem variables
	int flag;                      // reusable error-checking flag
	int rank;                      // MPI process ID
	N_Vector y = NULL;             // empty vector for storing solution
	void *arkode_mem = NULL;       // empty ARKode memory structure

	// initialize MPI **********************************************************************
	flag = MPI_Init(&argc, &argv);
	if (check_flag(&flag, "MPI_Init", 1)) return 1;
	flag = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

	// allocate and fill udata structure
	udata = new UserData;
	flag = init_user_data(udata);
	if (check_flag(&flag, "init_user_data", 1)) return 1;

	udata->R = R;
	udata->r = r;
	udata->rank = rank;
	udata->nx = nx;
	udata->ny = ny;
	udata->Diff = Diff;
	udata->dx = (XMAX - XMIN) / (1.0 * nx - 1.0);   // theta mesh spacing
	udata->dy = (YMAX - YMIN) / (1.0 * ny - 1.0);   // phi mesh spacing

	// Set up parallel decomposition
	flag = setup_decomp(udata);
	if (check_flag(&flag, "setup_decomp", 1)) return 1;

	// Find stable state of ODE model dependent on beta: run python script to solve the ODE problem
	std::string command = "SolveGoldbeterODE.py " + pt.get<std::string>("Parameters.beta");
	FILE* in = popen(command.c_str(), "r");
	double Zs, Ys;
	// Import the last two values of the solution as the stable state of the model
	int error = fscanf(in, "[%lf] [%lf]", &Zs, &Ys);
	// Put into udata structure so it can be used elsewhere
	udata->Zs = Zs;
	udata->Ys = Ys;

	// Initial problem std output
	bool outproc = (udata->rank == 0); 	// only outputs if the process is rank 0
	if (outproc)
	{
		cout << "\n Goldbeter model PDE problem on a torus:\n";
		cout << "   nprocs = " << udata->nprocs << "\n";
		cout << "   nx = " << udata->nx << "\n";
		cout << "   ny = " << udata->ny << "\n";
		cout << "   nxl = " << udata->nxl << "\n";
		cout << "   nyl = " << udata->nyl << "\n";
		cout << "   Diff = " << udata->Diff << "\n";
		cout << "   Tfinal = " << T_FINAL << "\n";
		cout << "   Output timesteps = " << OUTPUT_TIMESTEP << "\n";
		cout << "   Major circumference = " << MAJOR_CIRC << "\n";
		cout << "   Minor circumference = " << MINOR_CIRC<< "\n";
		cout << "   Wavelength = " << WAVE_LENGTH*100 << "\%\n";
		cout << "   Wavewidth = " << WAVE_WIDTH*100 << "\%\n";
		cout << "   rtol = " << rtol << "\n";
		cout << "   atol = " << atol << "\n";

		if (JUST_DIFFUSION == 1)
		{
			cout << "   Diffusion Only\n\n";
		}
		else
		{
			cout << "   Include all variables in output = " << INCLUDE_ALL_VARS << "\n";
			cout << "   Absorbing boundary turn off time = " << T_BOUNDARY << "\n";
			if (VARY_BETA == 0)
			{
				cout << "   Beta = " << BETA << "\n";
				cout << "   Stable state values: Z = " << Zs << ", Y = " << Ys << "\n\n";
			}
			else if (VARY_BETA == 1)
			{
				cout << "   Beta varied over torus\n";
				if (IC_TYPE == 0) {cout << "   Homogeneous ICs\n\n";}
				if (IC_TYPE == 1) {cout << "   ICs: initial perturbation\n\n";}
				if (IC_TYPE == 2) {cout << "   Random ICs\n\n";}
			}
		}
	}

	// Initialize data structures
	N = NVARS*(udata->nxl)*(udata->nyl);
	Ntot = NVARS*nx*ny;
	y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution
	if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;

	// Check where the initial perturbation is centered and adjust wave midpoint and the limits of the segment accordingly
	if (WaveInside == 1)
	{
		wave_midpoint = PI;
		wave_XMIN  = wave_midpoint - wave_width/2.0;
		wave_XMAX  = wave_midpoint + wave_width/2.0;
	}
	else if (WaveInside == 0)
	{
		wave_midpoint = 0.0;
		wave_XMIN  = wave_midpoint - wave_width/2.0 + (XMAX-XMIN);	// Add (XMAX-XMIN) to convert from negative to positive number in the range [XMIN,XMAX]
		wave_XMAX  = wave_midpoint + wave_width/2.0;
	}
	else
	{
		printf("WaveInside must be 0 or 1");
	}

	// Set initial conditions - these are model dependent and rely on WAVE_LENGTH, WAVE_WIDTH, VARY_BETA and IC_TYPE
	ydata = N_VGetArrayPointer(y);
	for (j=0; j<udata->nyl; j++)
	{
		yy = YMIN + (udata->js+j)*(udata->dy);						// Actual phi values

		for (i=0; i<udata->nxl; i++)
		{
			xx = XMIN + (udata->is+i)*(udata->dx);					// Actual theta values

			// Beta is constant over the entire surface
			if (VARY_BETA == 0)
			{
				// Initial perturbation is centered on the inside
				if (WaveInside == 1)
				{
					// Set initial wave segment
					if ( xx >= wave_XMIN && xx <= wave_XMAX && yy >= 1.0*wave_length && yy <= 2.0*wave_length )
					{
						ydata[IDX(i,j)] = Zs + 1;
						ydata[IDX(i,j) + 1] = Ys + 1;
					}
					else
					{
						// Set rest of area to stable
						ydata[IDX(i,j)] = Zs;
						ydata[IDX(i,j) + 1] = Ys;
					}
				}
				// Initial perturbation is centered on the outside
				else if (WaveInside == 0)
				{
					// Set initial wave segment
					if ( (xx >= wave_XMIN || xx <= wave_XMAX) && yy >= wave_length && yy <= (2.0*wave_length) )
					{
						ydata[IDX(i,j)] = Zs + 1;
						ydata[IDX(i,j) + 1] = Ys + 1;
					}
					else
					{
						// Set rest of area to stable
						ydata[IDX(i,j)] = Zs;
						ydata[IDX(i,j) + 1] = Ys;

					}
				}
			}

			// Beta is linearly varied over the surface with respect to phi
			else if (VARY_BETA == 1)
			{

				// Homogeneous ICs
				if (IC_TYPE == 0)
				{
					ydata[IDX(i,j)] = 0.4;
					ydata[IDX(i,j)+1] = 1.6;
				}

				// Initial perturbation with higher concentrations in some area of the domain
				if (IC_TYPE == 1)
				{
					if ( xx >= wave_XMIN && xx <= wave_XMAX && yy >= 2.0*wave_length && yy <= (3.0*wave_length) )
					{
						// Set this area to higher initial concentration
						ydata[IDX(i,j)] = 1.4;
						ydata[IDX(i,j) + 1] = 2.6;
					}
					else
					{
						// Set rest of area lower
						ydata[IDX(i,j)] = 0.4;
						ydata[IDX(i,j) + 1] = 1.6;
					}
				}

				// Set random ICs for each mesh point within [0,1.4]
				if (IC_TYPE == 2)
				{
					ydata[IDX(i,j)] = (float)rand() / RAND_MAX * 1.4;
					ydata[IDX(i,j) + 1] = (float)rand() / RAND_MAX * 1.4;
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

	// Each processor outputs subdomain information to a txt file to be used by other scripts
	char outname[100];
	sprintf(outname, "GoldbeterModel_torus_subdomain.%03i.txt", udata->rank);
	FILE *UFID = fopen(outname,"w");
	fprintf(UFID, "%li  %li  %li  %li  %li  %li %f %f %f\n",
			udata->nx, udata->ny, udata->is, udata->ie, udata->js, udata->je, XMIN, XMAX, T_FINAL);
	fclose(UFID);

	ydata = N_VGetArrayPointer(y);

	// Create files for each subdomain and each variable
	sprintf(outname, "GoldbeterModel_torus_Z.%03i.txt", udata->rank);
	UFID = fopen(outname, "w");

	sprintf(outname, "GoldbeterModel_torus_Y.%03i.txt", udata->rank);
	FILE *UFID2 = fopen(outname, "w");


	// Write initial conditions to files, one for each variable in each subdomain
	for (j=0; j<udata->nyl; j++)
	{
		for (i=0; i<udata->nxl; i++)
		{
			fprintf(UFID," %.16e", ydata[IDX(i,j)]);
			if (INCLUDE_ALL_VARS == 1) {fprintf(UFID2," %.16e", ydata[IDX(i,j) + 1]);}
		}
	}
	fprintf(UFID,"\n");
	if (INCLUDE_ALL_VARS == 1) {fprintf(UFID2,"\n");}


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
				if (INCLUDE_ALL_VARS == 1) {fprintf(UFID2," %.16e", ydata[IDX(i,j) + 1]);}
			}
		}
		fprintf(UFID,"\n");
		if (INCLUDE_ALL_VARS == 1) {fprintf(UFID2,"\n");}


		// Time taken since beginning
		time(&end_t);

		// Update total time since beginning and reset start_t
		total_t += difftime(end_t, start_t);
		start_t = end_t;

		// Estimated time remaining based on number of iterations left and total time taken so far
		eta = ( Nt - (iout+1) ) * ( total_t / (iout+1) );

		// Output progress in terms of percentage remaining, time elapsed and estimated time remaining
		if (outproc)
		{
			if (iout > 0)
			{
				// Rewrite previous line (fix eventually to be nicer): \b means backspace
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			}
			printf("   %3d \% | %3d min %2d sec elapsed | %3d min %2d sec remaining", 100*(iout+1)/Nt, (int)(total_t/60), ((int)total_t % 60), (int)(eta/60), ((int)eta % 60));
			fflush(stdout);
		}
	}
	if (outproc)  cout << "\n   ----------------------\n";

	fclose(UFID);
	if (INCLUDE_ALL_VARS == 1) {fclose(UFID2);}

	// Clean up and return with successful completion
	N_VDestroy_Parallel(y);        // Free vectors
	free_user_data(udata);         // Free user data
	delete udata;
	ARKodeFree(&arkode_mem);     // Free integrator memory
	flag = MPI_Finalize();       // Finalize MPI

	// End of MPI processes **************************************************************************************

	return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y) where y' = f(t,y)
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_VConst(0.0, ydot);                           // Initialize y' to zero
	UserData *udata = (UserData *) user_data;      // access problem data
	long int nxl = udata->nxl;                     // set variable shortcuts
	long int nyl = udata->nyl;
	realtype xx, yy;
	realtype Diff = udata->Diff;
	realtype dx = udata->dx;
	realtype dy = udata->dy;
	realtype R = udata->R;
	realtype r = udata->r;
	realtype *yarray, *ydotarray;
	yarray = NV_DATA_P(y);           // access data arrays
	ydotarray = NV_DATA_P(ydot);

	// exchange boundary data with neighbors
	int ierr = exchange(y, udata);
	if (check_flag(&ierr, "exchange", 1)) return -1;


	// Add diffusion term, approximated by finite difference methods

	long int i, j;
	for (j=1; j<nyl-1; j++)
	{
		for (i=1; i<nxl-1; i++)
		{
			xx = XMIN + (udata->is+i)*(dx);

			// Fill in diffusion for Z variable
			// Note that diffusion is different depending on the coordinate system - this is for a torus
			ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - yarray[IDX(i-1,j)]) )/(2*dx)
							+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
							+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);
		}
	}

	// iterate over subdomain boundaries

	// West face
	i=0;
	xx = XMIN + (udata->is+i)*(dx);
	for (j=1; j<nyl-1; j++)
	{
		ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - udata->Wrecv[NVARS*j]) )/(2*dx)
							+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + udata->Wrecv[NVARS*j]))/(dx*dx)
							+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);
	}
	// East face
	i=nxl-1;
	xx = XMIN + (udata->is+i)*(dx);
	for (j=1; j<nyl-1; j++)
	{
		ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(udata->Erecv[NVARS*j] - yarray[IDX(i-1,j)]) )/(2*dx)
							+ Diff*( (1/(r*r))* (udata->Erecv[NVARS*j] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
							+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);

	}
	// South face
	j=0;
	for (i=1; i<nxl-1; i++)
	{
		xx = XMIN + (udata->is+i)*(dx);

		ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - yarray[IDX(i-1,j)]) )/(2*dx)
						+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + udata->Srecv[NVARS*i]))/(dy*dy);
	}

	// North face
	j=nyl-1;
	for (i=1; i<nxl-1; i++)
	{
		xx = XMIN + (udata->is+i)*(dx);

		ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - yarray[IDX(i-1,j)]) )/(2*dx)
						+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (udata->Nrecv[NVARS*i] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);
	}

	// South-West corner
	i = 0;
	j = 0;
	xx = XMIN + (udata->is+i)*(dx);
	ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - udata->Wrecv[NVARS*j]) )/(2*dx)
						+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + udata->Wrecv[NVARS*j]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + udata->Srecv[NVARS*i]))/(dy*dy);

	// North-West corner
	i = 0;
	j = nyl-1;
	xx = XMIN + (udata->is+i)*(dx);
	ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(yarray[IDX(i+1,j)] - udata->Wrecv[NVARS*j]) )/(2*dx)
						+ Diff*( (1/(r*r))* (yarray[IDX(i+1,j)] - 2*yarray[IDX(i,j)] + udata->Wrecv[NVARS*j]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (udata->Nrecv[NVARS*i] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);

	// South-East corner
	i = nxl-1;
	j = 0;
	xx = XMIN + (udata->is+i)*(dx);
	ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(udata->Erecv[NVARS*j] - yarray[IDX(i-1,j)]) )/(2*dx)
						+ Diff*( (1/(r*r))* (udata->Erecv[NVARS*j] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (yarray[IDX(i,j+1)] - 2*yarray[IDX(i,j)] + udata->Srecv[NVARS*i]))/(dy*dy);


	// North-East corner
	i = nxl-1;
	j = nyl-1;
	xx = XMIN + (udata->is+i)*(dx);
	ydotarray[IDX(i,j)] = Diff*( (-sin(xx)/(r*(R+r*cos(xx))))*(udata->Erecv[NVARS*j] - yarray[IDX(i-1,j)]) )/(2*dx)
						+ Diff*( (1/(r*r))* (udata->Erecv[NVARS*j] - 2*yarray[IDX(i,j)] + yarray[IDX(i-1,j)]))/(dx*dx)
						+ Diff*( (1/(((R+r*cos(xx)))*((R+r*cos(xx)))))* (udata->Nrecv[NVARS*i] - 2*yarray[IDX(i,j)] + yarray[IDX(i,j-1)]))/(dy*dy);

	// Parameter beta (b)
	realtype b;

	// Add terms other than diffusion (model dependent), in this case
	// Z' += v0 + v1*b - v2 + v3 + kf*Y - k*Z
	// Y' += v2 + v3 + kf*Y

	if (JUST_DIFFUSION == 0) // If JUST_DIFFUSION = 1 then omit this step leaving only diffusion terms
	{
		for (j=0; j<nyl; j++)
		{
			yy = YMIN + (udata->js+j)*(udata->dy);

			// Set b to either spatially constant BETA if VARY_BETA = 0, or varying linearly with phi
			if (VARY_BETA == 0)
			{
				b = BETA;
			}
			else if (VARY_BETA == 1)
			{
				// Beta varying linearly from BETA_MIN to BETA_MAX
				b = BETA_MIN + yy * (BETA_MAX - BETA_MIN) / (YMAX - YMIN);
			}


			for (i=0; i<nxl; i++)
			{
				xx = XMIN + (udata->is+i)*(udata->dx);

				realtype Z = yarray[IDX(i,j)];
				realtype Y = yarray[IDX(i,j)+1];

				// Algebraic equations
				realtype v2 = VM2 * pow(Z,n) / ( pow(K2,n) + pow(Z,n) );
				realtype v3 = VM3 * pow(Y,m) * pow(Z,p) / ( (pow(KR,m) + pow(Y,m)) * (pow(KA,p) + pow(Z,p)) );

				// If we are on the north or south boundary of the entire domain and t<T_BOUNDARY, set u_t = 0 to simulate
				// Dirichlet boundary conditions with values equal to the initial conditions

				// North boundary of the domain
				if (udata->je == udata->ny-1 && t<T_BOUNDARY && j == nyl-1)
				{
					ydotarray[IDX(i,j)] = 0; 	// Z
					ydotarray[IDX(i,j)+1] = 0;  // Y
				}
				// South boundary of the domain
				else if (udata->js == 0 && t<T_BOUNDARY && j == 0)
				{
					ydotarray[IDX(i,j)] = 0; 	// Z
					ydotarray[IDX(i,j)+1] = 0;  // Y
				}
				else
				{
					// Add rest of the RHS of the ODEs
					ydotarray[IDX(i,j)] += v0 + v1*b - v2 + v3 + kf*Y - k*Z;  // Z
					ydotarray[IDX(i,j)+1] += v2 - v3 - kf*Y;				  // Y
				}

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
static int setup_decomp(UserData *udata)
{
	// check that this has not been called before
	if (udata->Erecv != NULL || udata->Wrecv != NULL ||
			udata->Srecv != NULL || udata->Nrecv != NULL) {
		cerr << "setup_decomp warning: parallel decomposition already set up\n";
		return 1;
	}

	// get suggested parallel decomposition based on mesh size, outputs dimensions of each subdomain (dims)
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
static int exchange(N_Vector y, UserData *udata)
{
	// local variables
	MPI_Request reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
	MPI_Status stat;
	int ierr, i, ipW=-1, ipE=-1, ipS=-1, ipN=-1;			// ip* are the ranks of neighbours defined later by MPI_Cart_shift
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

	// Send West side
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

	// Send East side
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

	// Send South side
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

	// Send North side
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
static int init_user_data(UserData *udata)
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
	udata->R = 0.0;
	udata->r = 0.0;

	return 0;     // return with success flag
}

// Free memory allocated within Userdata
static int free_user_data(UserData *udata)
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




