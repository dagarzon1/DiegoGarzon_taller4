#include <iostream>
#include <png.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <stdio.h>
#include <tgmath.h>
#include <math.h>
#include <complex>


using std::cout;
using std::endl;

void FT(double y[], double trans[]);
int main()
{
cout<<"Hola"<<endl;
return 0;
}
void FT(double* y, double* trans)
{
	double fun=*y;
	int n= (sizeof(y)/sizeof(y[0]));
	std::complex<double> res[n];
	std::complex<double> com=(0.0,-2.0);
	double resp[n];
	double pi = std::acos(-1);
	for(int i=0;i<n;i++)
	{	
		for(int j=0;j<n;j++)
		{
			std::complex<double> fo = ( com * pi * (double)i * (double)j ) / (double)n;
			res[i] += y[j]* exp(fo);
		}
		resp[i]=std::abs(res[i]);
	}
	*trans = resp;	
}
