#include<iostream>
#include<fstream>
#include<string.h>
#include<sstream>
#include<stdlib.h>
#include<complex>
#include<math.h>

using namespace std;

int main(int argc,char ** argv)
{
	ifstream f;
	string line;
	f.open(argv[1]);
	int lin_max=-1;
	double pi = -2.0 * acos(-1);
	while(!f.eof())
	{
		getline(f,line);
		lin_max = lin_max +1;
	}
	f.close();
	f.open(argv[1]);
	double * t=new double[lin_max];
	double *f_t=new double[lin_max];
	string d1=" ";
	const char* d=d1.c_str();
	for(int i=0;i<lin_max;i++)
	{
		getline(f,line);
		char* chr = strtok( strdup(line.c_str()) , d );
		t[i]=atof(chr);
		
		char* chr2 = strtok( strdup(line.c_str()) , d );
		f_t[i]=atof(chr2);
	}
	double * x = new double[lin_max];
	double h=( t[lin_max-1]-t[0] ) / (double) (lin_max-1); 
	for(int i=0;i<lin_max;i++)
	{
		x[i]=t[0]+ (h*i) ;
	}
	double * l= new double[lin_max];
	for(int i=0;i<lin_max;i++)
	{
		double mul=1.0;
		for(int j=0;j<lin_max;j++)
		{
			if( ( t[i] - t[j] ) != 0.0 )
			{
			mul= mul * ( x[i] - t[j] ) / ( t[i] - t[j] );
			}
		}
		l[i]=mul;
	}
	double * y = new double[lin_max];
	for(int i=0;i<lin_max;i++)	
	{
		y[i]=f_t[i] * l[i];
	}
	complex<double> * trans =new complex<double>[lin_max];
	double * freq =new double[lin_max];
	for(int i=0;i<lin_max;i++)
	{	
		for(int j=0;j<lin_max;j++)
		{	
			double theta= pi*j*i/(double)lin_max;
			trans[i]= trans[i] + ( polar(y[j],theta) ); 
		}
		trans[i] = trans[i] / (double) lin_max; 
		freq[i] = 1.0 / arg(trans[i]);
	}
	f.close();
	ofstream f_nuevo("transformada.txt");
	for(int i=0;i<lin_max;i++)
	{
		f_nuevo << freq[i] << " " << real( trans[i] ) << " " << imag( trans[i] ) << " "<<endl;
	}
}

