#include <math.h>
#include <complex>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

double ** FT(int ** y, int inverse, int M, int N);
double ** FT1(double ** y, int inverse, int M, int N);
int main()
{
	return 0;
}
double ** FT1(double ** y, int inverse, int M, int N)
{
	double ** freq=new double*[M];
	for(int i=0;i<M;i++)
	{
		freq = new double[N];
	}

	double inv=-1.0;
	double c=M*N;
	if (inverse==1)
	{
		inv=1.0;
	}
	double com=inv*2.0;
	double pi = acos(-1);
	//Matriz con M filas y N columnas	
	double ** trans= new double*[M];
	for(int i=0;i<M;i++)
	{
		trans[i]=new double[N];
	}
	complex<double> * e_1 = new complex<double> [M];
	complex<double> * e_2 = new complex<double> [N];
	
	for(int k=0;k<M;k++)
	{
		for(int m=0;m<M;m++)
		{
			double f=com * pi * (m*k/M);
			e_1[k]=polar(1.0,f);
		}	
	}
	for(int k=0;k<N;k++)
	{
		for(int m=0;m<N;m++)
		{
			double f=com * pi * (m*k/N);
			e_2[k]=polar(1.0,f);
		}	
	}
	complex<double> ** e = new complex<double> * [M];
	for(int i =0;i<M;i++)
	{
		e[i]=new complex<double>[N];
	}	
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			e[j][i]= e_2[i] * e_1[j]/c;
		}
	}
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			freq[i][j]= 1.0 / arg(e[i][j]);
		}
	}
	complex<double> ** res = new complex<double> * [M];
	for(int i =0;i<M;i++)
	{
		res[i]=new complex<double>[N];
	}
	for(int i=0;i<M;i++)
	{
		complex<double> t;
		for(int j=0;j<N;j++)
		{
			res[i][j]= (double)y[i][j] * e[i][j];		
		}
		t=0.0;
		for(int k=0;k<M;k++)
		{	
			for(int l=0;l<N;l++)
			{
				trans[k][l]= abs( trans[k][l] + res[k][l] );
			}
		}
	}
	return trans;
}
