/* FHN neuron model with periodic BCs on a Cartesian surface. ICs: stable state everywhere apart from
 * small rectangle centred at y=0

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
#include "arkode/arkode.h"            // prototypes for ARKode fcts., consts.
#include "nvector/nvector_parallel.h" // parallel N_Vector types, fcts., macros
#include "arkode/arkode_pcg.h"        // prototype for ARKPcg solver
#include "sundials/sundials_types.h"  // def. of type 'realtype'
#include <sundials/sundials_math.h>   // math macros
#include "mpi.h"                      // MPI header file

using namespace std;

// accessor macros between (x,y) location and 1D NVector array
#define IDX(x,y) (NVARS*(x) + NVARS*(y)*nxl)
#define PI RCONST(3.1415926535897932)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)
#define THREE RCONST(3.0)
#define NVARS 2

// Equation parameters
#define BETA 1.4
#define EPSILON 0.36
#define DIFF 0.12

/* TMAX is the end time for the absorbing boundary at the bottom, set to 0 for no absorbing boundary */
#define TMAX 1000.0
//#define BTHRES 1.344
//#define B0 1.32
//#define K 1458

#define XMIN         RCONST(0.0)                /* grid boundaries in x  */
#define XMAX         RCONST(20.0)

#define YMIN         RCONST(0.0)        		/* grid boundaries in y  */
#define YMAX         RCONST(80.0)

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
	  realtype Diff;        // diffusion coefficient for eq 1
	  realtype Diff2;		// diffusion coefficient for eq 2
	  N_Vector h;           // heat source vector
	  N_Vector d;           // inverse of Jacobian diagonal
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

// User-supplied Functions Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
//static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
//		booleantype *jcurPtr, realtype gamma, void *user_data,
//		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
//static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r,
//		N_Vector z, realtype gamma, realtype6 delta, int lr,
//		void *user_data, N_Vector tmp);

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

// Main Program
int main(int argc, char* argv[]) {

	  // general problem parameters
	  realtype T0 = RCONST(0.0);   	// initial time
	  realtype Tf = RCONST(1000.0);    // final time
	  int Nt = 100;                 	// total number of output times
	  long int nx = 100;            // spatial mesh size
	  long int ny = 400;
	  long int offset;
	  realtype xx, yy;				// real x,y values
	  realtype Diff = DIFF;			// Diffusion coefficient for eq 1
	  realtype Diff2 = 0;			// Diffusion coefficient for eq 2 (0 for FHN model)
	  realtype rtol = 1.e-5;       	// relative and absolute tolerances
	  realtype atol = 1.e-10;

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
	  udata->Diff2 = Diff2;
	  udata->dx = (XMAX-XMIN)/(1.0*nx-1.0);   // x mesh spacing       *** Change xlim and ylim here (2*pi)
	  udata->dy = (YMAX-YMIN)/(1.0*ny-1.0);   // y mesh spacing

	  // Set up parallel decomposition **********************************************************************
	  flag = SetupDecomp(udata);
	  if (check_flag(&flag, "SetupDecomp", 1)) return 1;

	  // Initial problem output
	  bool outproc = (udata->rank == 0);
	  if (outproc) {
		cout << "\n2D FHN model PDE problem:\n";
		cout << "   nprocs = " << udata->nprocs << "\n";
		cout << "   nx = " << udata->nx << "\n";
		cout << "   ny = " << udata->ny << "\n";
		cout << "   Diff = " << udata->Diff << "\n";
		cout << "   Diff2 = " << udata->Diff2 << "\n";
		cout << "   rtol = " << rtol << "\n";
		cout << "   atol = " << atol << "\n";
		cout << "   nxl (proc 0) = " << udata->nxl << "\n";
		cout << "   is (proc 0) = " << udata->is << "\n";
		cout << "   ie (proc 0) = " << udata->ie << "\n";
		cout << "   js (proc 0) = " << udata->js << "\n";
		cout << "   je (proc 0) = " << udata->je << "\n";
		cout << "   nyl (proc 0) = " << udata->nyl << "\n\n";
	  }

	  // Initialize data structures
	  N = NVARS*(udata->nxl)*(udata->nyl);
	  Ntot = NVARS*nx*ny;
	  y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution
	  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;

      // Set initial conditions
	  ydata = N_VGetArrayPointer(y);
	  for (j=0; j<udata->nyl; j++)
	  {
		for (i=0; i<udata->nxl; i++)
		{
			xx = (udata->is+i)*(udata->dx);					// Actual x and y values
	  	  	yy = (udata->js+j)*(udata->dy);
			if ( xx>(XMIN + XMAX/4.0) && xx<(XMAX - XMAX/4.0) && yy>(YMIN + YMAX/8.0) && yy<(YMIN + YMAX/4.0) )	// Set initial wave segment
			{
				ydata[offset] = RCONST(-BETA+2);								// u
				ydata[offset + 1] = RCONST(BETA*BETA*BETA - 3*BETA + 1.5);		// v
			}
			else
			{
				ydata[offset] = RCONST(-BETA);							// u	// Set rest of area to stable u,v
				ydata[offset + 1] = RCONST(BETA*BETA*BETA - 3*BETA);	// v

				// Set all area to perturbed state
//				ydata[offset] = RCONST(-BETA+2);								// u
//				ydata[offset + 1] = RCONST(BETA*BETA*BETA - 3*BETA + 1.5);		// v
			}
			offset += NVARS;
		}
	  }
	  udata->h = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for heat source
	  if (check_flag((void *) udata->h, "N_VNew_Parallel", 0)) return 1;
	  udata->d = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for Jacobian diagonal
	  if (check_flag((void *) udata->d, "N_VNew_Parallel", 0)) return 1;

	  arkode_mem = ARKodeCreate();                       // Create the solver memory
	  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;

	  // fill in the heat source array - move to f later (zero at the moment, maybe get rid of altogether)
	  data = N_VGetArrayPointer(udata->h);
	  offset = 0;
	  for (j=0; j<udata->nyl; j++)
		for (i=0; i<udata->nxl; i++)
		  data[offset] = RCONST(0.0);
	  	  data[offset+1] = RCONST(0.0);
	  	  offset += NVARS;

	  /* Call ARKodeInit to initialize the integrator memory and specify the
		 right-hand side function in y'=f(t,y), the inital time T0, and
		 the initial dependent variable vector y.  */
	  // f, NULL - explicit, switch for implicit
	  flag = ARKodeInit(arkode_mem, f, NULL, T0, y);
	  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

	  flag = ARKodeSStolerances(arkode_mem, rtol, atol);      // Specify tolerances
	  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

	  // Set routines
	  flag = ARKodeSetUserData(arkode_mem, (void *) udata);   // Pass heat source array udata to user functions
	  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;

	  flag = ARKodeSetMaxNumSteps(arkode_mem, 200000);
	  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return(1);


//	  flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);      // Update solver convergence coeff.
//	  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
//
//	  // Linear solver specification
//	  flag = ARKPcg(arkode_mem, 1, 20);                           // Specify the PCG solver
//	  if (check_flag(&flag, "ARKPcg", 1)) return 1;
//	  flag = ARKSpilsSetPreconditioner(arkode_mem, PSet, PSol);   // Specify the Preconditioner
//
//	  if (check_flag(&flag, "ARKSpilsSetPreconditioner", 1)) return 1;
//
//	  // Specify linearly implicit RHS, with non-time-dependent preconditioner
//	  flag = ARKodeSetLinear(arkode_mem, 0);
//	  if (check_flag(&flag, "ARKodeSetLinear", 1)) return 1;

	  // Each processor outputs subdomain information
	  char outname[100];
	  sprintf(outname, "FHNmodel_subdomain.%03i.txt", udata->rank);
	  FILE *UFID = fopen(outname,"w");
	  fprintf(UFID, "%li  %li  %li  %li  %li  %li\n",
		  udata->nx, udata->ny, udata->is, udata->ie, udata->js, udata->je);
	  fclose(UFID);

	  ydata = N_VGetArrayPointer(y);

	  sprintf(outname, "FHNmodel_u.%03i.txt", udata->rank);
	  UFID = fopen(outname, "w");										// w is write
	  sprintf(outname, "FHNmodel_v.%03i.txt", udata->rank);
	  FILE *UFID2 = fopen(outname, "w");

	  // Write initial conditions to files, one for each variable in each subdomain
	  offset = 0;
	  for (j=0; j<udata->nyl; j++)
	  {
		for (i=0; i<udata->nxl; i++)
		{
			  fprintf(UFID," %.16e", ydata[offset]);
			  fprintf(UFID2," %.16e", ydata[offset+1]);
			  offset += NVARS;
		}
	  }
	  fprintf(UFID,"\n");
	  fprintf(UFID2,"\n");





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
		  offset = 0;
		  for (j=0; j<udata->nyl; j++)
		  {
			for (i=0; i<udata->nxl; i++)
			{
				  fprintf(UFID," %.16e", ydata[offset]);
				  fprintf(UFID2," %.16e", ydata[offset+1]);
				  offset += NVARS;
			}
		  }
			fprintf(UFID,"\n");
			fprintf(UFID2,"\n");
	  }
	  if (outproc)  cout << "   ----------------------\n";
	  fclose(UFID);
	  fclose(UFID2);

	  // Clean up and return with successful completion
	  N_VDestroy_Parallel(y);        // Free vectors
	  N_VDestroy_Parallel(udata->h);
	  N_VDestroy_Parallel(udata->d);
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
	  long int offset = 0;
	  realtype xx, yy;
	  realtype Diff = udata->Diff;
	  realtype Diff2 = udata->Diff2;
	  realtype dx = udata->dx;
	  realtype dy = udata->dy;
	  realtype *yarray, *ydotarray;
	  yarray = NV_DATA_P(y);           // access data arrays
	  ydotarray = NV_DATA_P(ydot);

	  // Exchange boundary data with neighbors
	  int ierr = Exchange(y, udata);
	  if (check_flag(&ierr, "Exchange", 1)) return -1;

	  // iterate over subdomain *interior*, computing approximation to RHS
	  realtype cu1 = Diff/dx/dx;			// D/delx^2
	  realtype cu2 = Diff/dy/dy;
	  realtype cu3 = -TWO*(cu1 + cu2);
	  realtype cv1 = Diff2/dx/dx;			// D2/delx^2
	  realtype cv2 = Diff2/dy/dy;
	  realtype cv3 = -TWO*(cv1 + cv2);
	  long int i, j;
	  for (j=1; j<nyl-1; j++)                        // diffusive terms
	  {
		for (i=1; i<nxl-1; i++)
		{
			// Fill in diffusion for u variable
			ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
							 + cu2*(yarray[IDX(i,j-1)] + yarray[IDX(i,j+1)])
							 + cu3*yarray[IDX(i,j)];
	  	  	// Diffusion for other variable
			ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + yarray[IDX(i+1,j)+1])
							 + cv2*(yarray[IDX(i,j-1)+1] + yarray[IDX(i,j+1)+1])
							 + cv3*yarray[IDX(i,j)+1];

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
			ydotarray[IDX(i,j)+1] = cv1*(udata->Wrecv[NVARS*j+1]   + yarray[IDX(i+1,j)+1])
							 + cv2*(yarray[IDX(i,j-1)+1] + yarray[IDX(i,j+1)+1])
							 + cv3*yarray[IDX(i,j)+1];
		}
		// East face
		i=nxl-1;
		for (j=1; j<nyl-1; j++)
		{
			ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
							 + cu2*(yarray[IDX(i,j-1)] + yarray[IDX(i,j+1)])
							 + cu3*yarray[IDX(i,j)];
			ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + udata->Erecv[NVARS*j+1])
							 + cv2*(yarray[IDX(i,j-1)+1] + yarray[IDX(i,j+1)+1])
							 + cv3*yarray[IDX(i,j)+1];
		}
	  // South face
		j=0;
		if (udata->js == 0 && t<TMAX)
		{
			for (i=1; i<nxl-1; i++)
			{
				ydotarray[IDX(i,j)] = 0;
			}
		}
		else
		{
			for (i=1; i<nxl-1; i++)
			{
				xx = (udata->is+i)*(dx);
				yy = (udata->js+j)*(dy);

				ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
								 + cu2*(udata->Srecv[NVARS*i]   + yarray[IDX(i,j+1)])
								 + cu3*yarray[IDX(i,j)];
				ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + yarray[IDX(i+1,j)+1])
								 + cv2*(udata->Srecv[NVARS*i+1]   + yarray[IDX(i,j+1)+1])
								 + cv3*yarray[IDX(i,j)+1];
			}
		}

	  // North face
		j=nyl-1;
		for (i=1; i<nxl-1; i++)
		{
			ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + yarray[IDX(i+1,j)])
							 + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
							 + cu3*yarray[IDX(i,j)];
			ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + yarray[IDX(i+1,j)+1])
							 + cv2*(yarray[IDX(i,j-1)+1] + udata->Nrecv[NVARS*i+1])
							 + cv3*yarray[IDX(i,j)+1];
		}
	  // South-West corner
		i = 0;
		j = 0;
		if (udata->js == 0 && t<TMAX)
		{
			ydotarray[IDX(i,j)] = 0;
		}
		else
		{
			ydotarray[IDX(i,j)] = cu1*(udata->Wrecv[NVARS*j] + yarray[IDX(i+1,j)])
							   + cu2*(udata->Srecv[NVARS*i] + yarray[IDX(i,j+1)])
							   + cu3*yarray[IDX(i,j)];
			ydotarray[IDX(i,j)+1] = cv1*(udata->Wrecv[NVARS*j+1] + yarray[IDX(i+1,j)+1])
							   + cv2*(udata->Srecv[NVARS*i+1] + yarray[IDX(i,j+1)+1])
							   + cv3*yarray[IDX(i,j)+1];
		}
	  // North-West corner
		i = 0;
		j = nyl-1;
		ydotarray[IDX(i,j)] = cu1*(udata->Wrecv[NVARS*j]   + yarray[IDX(i+1,j)])
						   + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
						   + cu3*yarray[IDX(i,j)];
		ydotarray[IDX(i,j)+1] = cv1*(udata->Wrecv[NVARS*j+1]   + yarray[IDX(i+1,j)+1])
						   + cv2*(yarray[IDX(i,j-1)+1] + udata->Nrecv[NVARS*i+1])
						   + cv3*yarray[IDX(i,j)+1];
	  // South-East corner
		i = nxl-1;
		j = 0;
		if (udata->js == 0 && t<TMAX)
		{
			ydotarray[IDX(i,j)] = 0;
		}
		else
		{
			ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
							   + cu2*(udata->Srecv[NVARS*i]   + yarray[IDX(i,j+1)])
							   + cu3*yarray[IDX(i,j)];
			ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + udata->Erecv[NVARS*j+1])
							   + cv2*(udata->Srecv[NVARS*i+1]   + yarray[IDX(i,j+1)+1])
							   + cv3*yarray[IDX(i,j)+1];
		}
	  // North-East corner
		i = nxl-1;
		j = nyl-1;
		ydotarray[IDX(i,j)] = cu1*(yarray[IDX(i-1,j)] + udata->Erecv[NVARS*j])
						   + cu2*(yarray[IDX(i,j-1)] + udata->Nrecv[NVARS*i])
						   + cu3*yarray[IDX(i,j)];
		ydotarray[IDX(i,j)+1] = cv1*(yarray[IDX(i-1,j)+1] + udata->Erecv[NVARS*j+1])
						   + cv2*(yarray[IDX(i,j-1)+1] + udata->Nrecv[NVARS*i+1])
						   + cv3*yarray[IDX(i,j)+1];

//		// Wavesize = integral of u over the grid (where u > 0)
//		realtype S = 0.0;
//		for (j=0; j<nyl; j++)
//		{
//			for (i=0; i<nxl; i++)
//			{
//				realtype u = yarray[IDX(i,j)];
//
//				if (u > 0)
//				{
//					S += u*dx*dy;
//				}
//			}
//		}



		// Add other terms in equations
		for (j=0; j<nyl; j++)
		{
	  	  	yy = (udata->js+j)*(dy);

			for (i=0; i<nxl; i++)
			{
				xx = (udata->is+i)*(dx);					// Actual x and y values

				realtype u = yarray[IDX(i,j)];
				realtype v = yarray[IDX(i,j)+1];

//				realtype BETAt = B0 + K*S;
		  	  	realtype BETAt = 0.5 + (2.0-0.5)*yy/(YMAX-YMIN);		// Vary beta between __ and __ over y axis

				// u variable: du/dt = 3u - u^3 - v + Diff
		  	    ydotarray[IDX(i,j)] += THREE*u - u*u*u - v;

		  	    // v variable: dv/dt = eps(u + beta)
				ydotarray[IDX(i,j)+1] += EPSILON*(u + BETA);
			}
		}

	  return 0;                                      // Return with success
}

 // Preconditioner setup routine (fills inverse of Jacobian diagonal)
//static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
//		booleantype *jcurPtr, realtype gamma, void *user_data,
//		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
//{
//	  UserData *udata = (UserData *) user_data;      // variable shortcuts
//	  realtype Diff = udata->Diff;
//	  realtype dx = udata->dx;
//	  realtype dy = udata->dy;
//	  realtype *diag = N_VGetArrayPointer(tmp1);  // access data arrays
//	  if (check_flag((void *) diag, "N_VGetArrayPointer", 0)) return -1;
//
//	  // set all entries of tmp1 to the diagonal values of interior
//	  // (since boundary RHS is 0, set boundary diagonals to the same)
//	  realtype c = ONE + gamma*TWO*(Diff/dx/dx + Diff/dy/dy);
//	  N_VConst(c, tmp1);
//	  N_VInv(tmp1, udata->d);      // set d to inverse of diagonal
//	  return 0;                    // Return with success
//}
//
//// Preconditioner solve routine
//static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r,
//		N_Vector z, realtype gamma, realtype delta, int lr,
//		void *user_data, N_Vector tmp)
//{
//	  UserData *udata = (UserData *) user_data;  // access user_data structure
//	  N_VProd(r, udata->d, z);                   // perform Jacobi iteration
//	  return 0;                                  // Return with success
//}

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
	  int periods[] = {1,1};		// 1 for periodic BCs
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
	  int ierr, i, ipW=-1, ipE=-1, ipS=-1, ipN=-1;			// ip* are ranks of neighbours
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
		long int offset = 0;
		for (i=0; i<nyl; i++)
		{
			udata->Wsend[offset] = Y[IDX(nxl-1,i)];		// Fill data to send with Y values on West side
			udata->Wsend[offset+1] = Y[IDX(nxl-1,i)+1];
			offset += NVARS;
		}
		ierr = MPI_Isend(udata->Wsend, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipW, 0,
				  udata->comm, &reqSW);
		if (ierr != MPI_SUCCESS) {
		  cerr << "Error in MPI_Isend = " << ierr << "\n";
		  return -1;
		}

		offset = 0;
		for (i=0; i<nyl; i++)
		{
			udata->Esend[offset] = Y[IDX(0,i)];
			udata->Esend[offset+1] = Y[IDX(0,i)+1];
			offset += NVARS;
		}
		ierr = MPI_Isend(udata->Esend, NVARS*(udata->nyl), REALTYPE_MPI_TYPE, ipE, 1,
				  udata->comm, &reqSE);
		if (ierr != MPI_SUCCESS) {
		  cerr << "Error in MPI_Isend = " << ierr << "\n";
		  return -1;
		}

		offset = 0;
		for (i=0; i<nxl; i++)
		{
			udata->Ssend[offset] = Y[IDX(i,nyl-1)];
			udata->Ssend[offset+1] = Y[IDX(i,nyl-1)+1];
			offset += NVARS;
		}
		ierr = MPI_Isend(udata->Ssend, NVARS*(udata->nxl), REALTYPE_MPI_TYPE, ipS, 2,
				  udata->comm, &reqSS);
		if (ierr != MPI_SUCCESS) {
		  cerr << "Error in MPI_Isend = " << ierr << "\n";
		  return -1;
		}

		offset = 0;
		for (i=0; i<nxl; i++)
		{
			udata->Nsend[offset] = Y[IDX(i,0)];
			udata->Nsend[offset+1] = Y[IDX(i,0)+1];
			offset += NVARS;
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
	  udata->Diff2 = 0.0;
	  udata->h = NULL;
	  udata->d = NULL;
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




