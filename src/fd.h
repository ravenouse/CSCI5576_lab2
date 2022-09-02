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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include "stdio.h"
#include "math.h"
#include "string.h"

using std :: string;
using std :: vector;
using std :: stringstream;
using std :: cout;
using std :: endl;

#define rLOOP  for ( int r = 1 ; r <= nField    ; ++r )
#define cLOOP  for ( int c = 1 ; c <= nField    ; ++c )
#define iLOOP  for ( int i = 1 ; i <= nRealx    ; ++i )
#define jLOOP  for ( int j = 1 ; j <= nRealy    ; ++j )
#define sLOOP  for ( int s = 0 ; s <= nRealx+1  ; ++s )
#define tLOOP  for ( int t = 0 ; t <= nRealy+1  ; ++t )


typedef  vector<double>          VD;
typedef  vector<vector<double> > VDD;

