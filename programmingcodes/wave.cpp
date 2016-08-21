//
//  main_wave.cpp
//  1d wave_equation: u_tt - c^2 u_xx=0 (c:const.)
//
//  Created by 難波時永 on 5/28/16.
//  Copyright (c) 2016 難波時永. All rights reserved.
//

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

# include "fd1d_wave.hpp"


void main()
{
	
}

//compute lambda
double fd1d_wave_alpha ( int x_num, double x1, double x2, int t_num, double t1, double t2, double c )
{
	double delta_t;
	double delta_x;
	double lambda;

	delta_t = ( t2 - t1 ) / ( double ) ( t_num - 1 );
	delta_x = ( x2 - x1 ) / ( double ) ( x_num - 1 );
	alpha = c * delta_t / delta_x;

	cout << "\n";
	cout << "Stability condition LAMBDA = C * DT / DX = " << lambda << "\n";

	if ( 1.0 < r8_abs ( alpha ) ){
    	cerr << "\n";
    	cerr << "FD1D_WAVE_ALPHA - Warning!\n";
    	cerr << "The stability condition |ALPHA| <= 1 fails.\n";
    	cerr << "Computed results are liable to be inaccurate.\n";
  	}

  	return alpha;
}

//first step: solve wave eq.
double *fd1d_wave_start ( int x_num, double x_vec[], double t, double t_delta, double alpha,
						double u_x1 ( double t ), double u_x2 ( double t ),
						double *ut_t1 ( int x_num, double x_vec[] ), double u1[] )
{
  	int j;
  	double *u2;
  	double *ut;

  	ut = ut_t1 ( x_num, x_vec );

  	u2 = new double[x_num];

  	u2[0] = u_x1 ( t );

  	for ( j = 1; j < x_num - 1; j++ )
  	{
    	u2[j] =           alpha * alpha   * u1[j+1] / 2.0
        	    + ( 1.0 - alpha * alpha ) * u1[j] 
            	+         alpha * alpha   * u1[j-1] / 2.0
            	+         t_delta         * ut[j];
  	}

  	u2[x_num-1] = u_x2 ( t );

  	delete [] ut;

  	return u2;
}

//computes steps
double *fd1d_wave_step ( int x_num, double t, double alpha, 
  double u_x1 ( double t ), double u_x2 ( double t ), double u1[], double u2[] )
{
	int j;
  	double *u3;

 	u3 = new double[x_num];

  	u3[0] = u_x1 ( t );

  	for ( j = 1; j < x_num - 1; j++ ){
    	u3[j] =                 alpha * alpha   * u2[j+1] 
        	  	+ 2.0 * ( 1.0 - alpha * alpha ) * u2[j]
            	+               alpha * alpha   * u2[j-1] 
            	-                                 u1[j];
  	}

  	u3[x_num-1] = u_x2 ( t );

  	return u3;
}

//evaluate a piecewise linear spline
double *piecewise_linear ( int nd, double xd[], double yd[], int nv, double xv[] )
{
  	int id;
  	int iv;
  	double *yv;

  	yv = new double[nv];

  	for ( iv = 0; iv < nv; iv++ ){
    	if ( xv[iv] < xd[0] ){
     		yv[iv] = yd[0];
    	}
    	else if ( xd[nd-1] < xv[iv] ){
   			yv[iv] = yd[nd-1];
    	}
    	else{
    		for ( id = 1; id < nd; id++ ){
        		if ( xv[iv] < xd[id] ){
            	yv[iv] = ( ( xd[id] - xv[iv]            ) * yd[id-1] 
                	   + (          xv[iv] - xd[id-1] ) * yd[id] ) 
                   	   / ( xd[id]          - xd[id-1] );
          		break;
   	        	}
        	}
    	}
  	}
	return yv;
}

//return the absolute vale of an R8
double r8_abs ( double x )
{
	double value;

  	if ( 0.0 <= x ){
 	   value = + x;
  	}
  	else{
    	value = - x;
  	}
  	return value;
}

//writes an R8MAT file
void r8mat_write ( string output_filename, int m, int n, double table[] )
{
	int i;
  	int j;
  	ofstream output;
//
//  Open the file.
//
  	output.open ( output_filename.c_str ( ) );

  	if ( !output ){
    	cerr << "\n";
    	cerr << "R8MAT_WRITE - Fatal error!\n";
    	cerr << "Could not open the output file.\n";
    	exit ( 1 );
  	}
//
//  Write the data.
//
	for ( j = 0; j < n; j++ ){
		for ( i = 0; i < m; i++ ){
	    	output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    	}
    	output << "\n";
	}
//
//  Close the file.
//
  	output.close ( );

  	return;
}

//creates a vector of linearly spaced values
double *r8vec_linspace_new ( int n, double a_first, double a_last )
{
	double *a;
	int i;

  	a = new double[n];

  	if ( n == 1 ){
    	a[0] = ( a_first + a_last ) / 2.0;
    }
  	else{
    	for ( i = 0; i < n; i++ ){
		    a[i] = ( ( double ) ( n - 1 - i ) * a_first 
        	     + ( double ) (         i ) * a_last ) 
            	 / ( double ) ( n - 1     );
   		}
  	}
  	return a;
}

//prints the current YMDHMS date as a time stamp
void timestamp ( )
{
	# define TIME_SIZE 40

  	static char time_buffer[TIME_SIZE];
  	const struct std::tm *tm_ptr;
  	size_t len;
  	std::time_t now;

  	now = std::time ( NULL );
  	tm_ptr = std::localtime ( &now );

  	len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  	std::cout << time_buffer << "\n";

  	return;
	# undef TIME_SIZE
}


