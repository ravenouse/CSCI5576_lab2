//  ================================================================================
//  ||                                                                            ||
//  ||              fd                                                            ||
//  ||              -------------------------------------------                   ||
//  ||              F I N I T E   D I F F E R E N C E                             ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              -------------------------------------------                   ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020-2022 Scott Runnels                                ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================

#include "mpi.h"
#include "fd.h"
#include "mpiInfo.h"


//  ==
//  ||
//  ||   M E S H  +  F I E L D   V A R I A B L E   N U M B E R I N G   S Y S T E M
//  ||
//  ==
//
//  The mesh for a 3x3 mesh is shown below. 
//
//  Some key variables
//  ------------------------------------
//
//  nRealx = number of **physical** nodes in the x-direction
//  nRealy =    "   "       "         "    "  "  y-    "
//
//  In the figure below,
//
//  nRealx = 3   nRealy = 3
//
//  Field variable/natural node numbering
//  -------------------------------------
//
//  The ID numbers of the field variables, i.e., variable we're solving for, is shown
//  in () and [] in the figure.  Note:
//
//  [...] are Physical nodes --> those that lie in the physical doman for this processor
//  (...) are Ghost    nodes --> those that lie outside the physical doman for this processor
//
//  i-j node numbering
//  -------------------------------------
//
//  In addition to natural node numbering, used in the linear system solver and in the
//  field variable, nodes are numbered using an i-j logical numbering scheme, also as
//  shown in the figure below.  
//
//  Physical nodes --> range from 1,1 to nRealx,nRealy
//  Ghost    nodes --> are the candy-coating, with lower left at 0,0 and upper right at nRealx+1,nRealy+1
//
//    (21)----(22)----(23)----(24)-----(25)  <--- nRealx+1 , nRealy+1
//     |       |       |       |        |   
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (16)----[17]----[18]----[19]-----(20)
//     |    1,3|    2,3|    3,3|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (11)----[12]----[13]----[14]-----(15)
//     |    1,2|    2,2|    3,2|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (6)-----[7]-----[8]-----[9]------(10) 
//     |    1,1|    2,1|    3,1|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (1)-----(2)-----(3)-----(4)------(5)
//  0,0                              4,0
//


//  ==
//  ||
//  ||      C L A S S:   L A P L A C I A N O N G R I D 
//  ||
//  == 

class LaplacianOnGrid
{

public:

  double x0, x1, y0, y1;
  VD x,y;
  int nRealx    , nRealy   , nField;
  double dx, dy;
  VDD A ;    VD  phi ;  VD  b ;
  int myPE;

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  LaplacianOnGrid(double _x0   , double _x1,  double _y0, double _y1 , int ncell_x , int ncell_y, mpiInfo &myMPI )
  {

    x0 = _x0;    x1 = _x1;
    y0 = _y0;    y1 = _y1;
    
    nRealx     = ncell_x + 1; 
    nRealy     = ncell_y + 1; 
    nField  = (nRealx+2)*(nRealy+2);  
    dx     = (x1-x0)/ncell_x;
    dy     = (y1-y0)/ncell_y;

    // Allocate memory -- Note that the node numbers and field variable numbers must
    // match.  So even though we only will be caring about real nodes, their node
    // numbers are naturally ordered, so must be of size nField.

    x.resize(nField+1); y.resize(nField+1);
		       
    for ( int i = 1 ; i <= nRealx ; ++i )
      for ( int j = 1 ; j <= nRealy ; ++j )
	{
	  int p = pid(i,j);
	  x[p] = x0 + (i-1)*dx;
	  y[p] = y0 + (j-1)*dy;
	}

    A.resize(nField+1 ); rLOOP A[r].resize(nField+1);
    b.resize(nField+1 );
    phi.resize(nField+1);

    myPE = myMPI.myPE;
  }

  //  ==
  //  ||
  //  ||  Form Linear System Ax = b
  //  ||
  //  ==

  void FormLS(mpiInfo &myMPI)
  {
    
    rLOOP cLOOP A[r][c] = 0.;  // Initialize linear system
    rLOOP b[r] = 0.;
    rLOOP phi[r] = myPE;
    
    double dx2 = dx*dx;           // Form matrix entries for the interior grid points
    double dy2 = dy*dy;           // Form matrix entries for the interior grid points
    
    iLOOP      // A*phi = b has extra rows and columns to accommodate nodes outside the physical domain.
      jLOOP    // That is why we can loop all the from 1 to nRealx and 1 to nRealy.
      {
	int p = pid(i,j);
	A[ p ][ pid( i  , j  ) ] = -2./dx2 - 2./dy2;
	A[ p ][ pid( i+1, j  ) ] =  1./dx2;
	A[ p ][ pid( i-1, j  ) ] =  1./dx2;
	A[ p ][ pid( i  , j+1) ] =  1./dy2;
	A[ p ][ pid( i  , j-1) ] =  1./dy2;
      }

    // Populate arrays, one for each side, containing the Dirichlet BCs

    VD phiL, phiR, phiT, phiB; phiL.resize(nRealy+2);  phiR.resize(nRealy+2); phiB.resize(nRealx+2); phiT.resize(nRealx+2);

    for ( int i = 0 ; i <= nRealx + 1 ; ++i ) { phiB[i] = 1.; phiT[i] = -1.; }
    for ( int j = 0 ; j <= nRealy + 1 ; ++j ) { phiL[j] = 1.; phiR[j] = -1.; }

    // Put a diagonal on each ghost cell, regardless of PE.  The actual phi values used here
    // will not be used in serial cases and will get overwritten during parallel communication.

    ApplyBCs(   0     , nRealx+1 ,        0 ,        0 , phiB );
    ApplyBCs(   0     , nRealx+1 , nRealy+1 , nRealy+1 , phiT );
    ApplyBCs(   0     ,        0 ,        0 , nRealy+1 , phiL );
    ApplyBCs(nRealx+1 , nRealx+1 ,        0 , nRealy+1 , phiR );
    
    // Now, apply actual BCs
    
    if ( myMPI.jPE  == 0             ) ApplyBCs(  0      ,   nRealx+1 ,   0      , 0          , phiB);
    if ( myMPI.jPE  == myMPI.nPEy-1  ) ApplyBCs(  0      ,   nRealx+1 ,   nRealy , nRealy     , phiT);
    if ( myMPI.iPE  == 0             ) ApplyBCs(  0      ,   0        ,   0      , nRealy+1   , phiL);
    if ( myMPI.iPE  == myMPI.nPEx-1  ) ApplyBCs(  nRealx ,   nRealx   ,   0      , nRealy+1   , phiR);
      
  }

  void ApplyBCs(int iMin, int iMax , int jMin, int jMax, VD &phiValues)
  {
    int count = 0;
    for ( int i = iMin ; i <= iMax ; ++i )
      for ( int j = jMin ; j <= jMax ; ++j )
	{
	  int p = pid(i,j);
	  cLOOP A[p][c] = 0.;
	  A[ p ] [ p ] = 1. ;  b[ p ] =  phiValues[count];
	  ++count;
	}

    return;
  }

  //  ==
  //  ||
  //  ||  Utility routines
  //  ||
  //  ==

  int pid(int i,int j) { return (i+1) + (j)*(nRealx+2); }  // Given i-j, return point ID.  Here i-j is the physical grid.
                                                           // The row/col numbers must include the candy-coating.
  #include "plotter.h"
  #include "gauss_seidel.h"
};



//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==


int main(int argc, char *argv[])
{
  
   mpiInfo myMPI;
   MPI_Init     (&argc         , &argv       );
   MPI_Comm_size(MPI_COMM_WORLD, &myMPI.numPE);
   MPI_Comm_rank(MPI_COMM_WORLD,&myMPI.myPE  );


  int nPEx, nPEy, nCellx, nCelly;

    if ( myMPI.myPE == 0 )
      {
       cout << "\n";
       cout << "---------------------------------------------\n";
       cout << "\n";
       cout << " F I N I T E   D I F F E R E N C E         \n";
       cout << " D E M O   C O D E                           \n";
       cout << "\n";
       cout << " Running on " << myMPI.numPE << " processors \n";
       cout << "\n";
       cout << "---------------------------------------------\n";
       cout << "\n";
      }
  
   for (int count =  0 ; count < argc; ++count)
     {
       if ( !strcmp(argv[count],"-nPEx"  ) ) nPEx   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nPEy"  ) ) nPEy   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCellx") ) nCellx = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCelly") ) nCelly = atoi(argv[count+1]);
     }

   if ( myMPI.myPE == 0 )
     {
       cout << endl;
       cout << "Input Summary: " << endl;
       cout << "--------------------------- " << endl;
       cout << "No. PE   in  x-direction: " << nPEx   << endl;
       cout << "             y-direction: " << nPEy   << endl;
       cout << "No. Cells in x-direction: " << nCellx << endl;
       cout << "             y-direction: " << nCelly << endl;
       cout << endl;
     }

   myMPI.GridDecomposition(nPEx,nPEy,nCellx,nCelly);

   double x0, x1;
   double y0, y1;
   double length = 1.;

   x0 =  length*myMPI.iPE;   x1 = x0 + length;
   y0 =  length*myMPI.jPE;   y1 = y0 + length;

   cout << "myPE: " << myMPI.myPE << " x0 = " << x0 << endl;
   cout << "myPE: " << myMPI.myPE << " x1 = " << x1 << endl;
   cout << "myPE: " << myMPI.myPE << " y0 = " << y0 << endl;
   cout << "myPE: " << myMPI.myPE << " y1 = " << y1 << endl;
   
   LaplacianOnGrid F(x0,x1,y0,y1,nCellx,nCelly, myMPI );

   MPI_Barrier(MPI_COMM_WORLD);

   if ( myMPI.myPE == 0 )
     {
       cout << "                                          \n";
       cout << "Form Linear System:                       \n";
       cout << "------------------------------------------\n";
       cout << "                                          \n";
     }

   F.FormLS(myMPI);

   MPI_Barrier(MPI_COMM_WORLD);

   if ( myMPI.myPE == 0 )
     {
       cout << "                                          \n";
       cout << "Solve Linear System:\n";
       cout << "------------------------------------------\n";
       cout << "                                          \n";
     }

   F.GS_or_Jacobi(500, F.b , F.phi , myMPI , 2);

   MPI_Barrier(MPI_COMM_WORLD);
  
   if ( myMPI.myPE == 0 )
     {
       cout << "                                          \n";
       cout << "Plot Results:\n";
       cout << "------------------------------------------\n";
       cout << "                                          \n";
     }

   F.plot("F_phi",F.phi,myMPI);

   MPI_Finalize();
   return 0;

}
