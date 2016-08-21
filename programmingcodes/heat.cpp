//
//  heat.cpp
//  non-homogeneous heat equation with Dirichlet or
//  Neumann boundary conditions: u_t - D*Delta u = f
//
//  Created by 難波時永 on 5/28/16.
//  Copyright (c) 2016 難波時永. All rights reserved.
//

#include <iostream>
#include <cstdio>
#include <math.h>
using namespace std;

#define nmax 500 //number of unknown functions

int main()
{
	int i, n;
	double D = 1.0;			//diffusion coefficient
	int N = 15;				//number of partitions !! N+1 < nmax !!
	double h = 1.0/(N+1);
	double lam = 0.4;		//tau*D/(h*h);	
	double tau = lam*h*h/D;	//0.4*h*h/D
	double tmax = 0.5;
	double m = (int)(tmax/tau);
	double x[nmax],u[nmax],v[nmax];
	double t;

	for(i=0; i<=N+1; i++){
		x[i] = i*h;
		//u[i] = phi(x[i]);
		u[i] = sin(x[i]); //!!compatibility condition!!
	}
	t = 0.0;

	//printf("N:");
	cout << "N:" << N << " " << "h:" << h << " " << "tau:" << tau
	 << " " << "lam:" << lam << " " << "tmax:" << tmax << endl;
	for(i=0; i<=N+1; i++){
		//printf("%lf %lf %lf\n", x[i],t,u[i]);
		cout << x[i] << " " << t << " " << u[i] << endl;
	}
	//printf("\n\n");
	cout << endl;
	cout << endl;z2

	for(n=0; n<=m; n++){
		t = n*tau;
		v[0] = 0.0;
		for(i=0; i<=N; i++){
			v[i] = (1.0 - 2.0*lam)*u[i] + lam*(u[i-1] + u[i+1]);
		}
		v[N+1] = 0.0;
		for(i=0; i<=N+1; i++){
			u[i] = v[i];
		}
		for(i=0; i<=N+1; i++){
			//printf("%lf %lf %lf\n",x[i],t,u[i]);
			cout << x[i] << " " << t << " " << u[i] << endl;
		}
		//prinf("\n\n");
		cout << endl;
		cout << endl;
	}

	return(0);
}

/*doube boundary(*u, char)
{
	if(char == 1){
		u[0] = 0.0;
		u[N+1] = 0.0;
	}else{
		u[0] = u[1];
		u[N+1] = u[N];
	}
}*/