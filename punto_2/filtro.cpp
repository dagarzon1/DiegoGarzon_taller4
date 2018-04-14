#include <math.h>
#include <complex>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

using namespace std;

double ** FT(int ** y, int inverse, int M, int N, double ** freq);
int ** read_image(char *filename,int * height, int * width, int *b, int *co);
void write_image(char *filen, int h, int w, int bit, int color, double ** y);
double ** FT1(double ** y, int inverse, int M, int N, double ** freq);
double ** filter(double ** y, double ** freq, int M, int N, string ftr);
int main(int argc,char ** argv)
{
	//n es columnas
	int n,bit,color;
	//m es filas
	int m;
	int ** mat = read_image(argv[1],&m,&n,&bit,&color);
	double ** f= new double*[M];
	for(int i=0;i<M;i++)
	{
		f[i]=new double[N];
	}
	double ** imgft = FT(mat,0,m,n,f);
	double ** imgfft = filter(imgft,f,m,n,atof(argv[2]));
	double ** r = FT1(imgffr,1,m,n,f);
	
	
	
	return 0;
}
void write_image(char *filen, int h, int w, int bit, int color, double ** y)
{
	FILE *f=fopen(filen, "wb");
	png_structp ptr;
	png_infop info;
	ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	info=png_create_info_struct(ptr);
	png_init_io(ptr,f);
	png_set_IHDR(ptr,info,w,h,bit,color,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
	png_write_info(ptr,info);
	png_bytepp img=(png_bytepp) y;
	png_write_image(ptr,img);
	png_write_end(ptr,NULL); 
	fclose(f);	
}
double ** FT1(int ** y, int inverse, int M, int N, double ** freq)
{
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
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)
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
double ** FT(int ** y, int inverse, int M, int N, double ** freq)
{

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
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)
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
				trans[k][l]= abs( trans[k][l] + res[k][l]);
			}
		}}
	return trans;
}
double ** filter(double ** y, double ** freq, int M, int N, string ftr)
{
	int cut=2500;
	int w=200;
	double pi = acos(-1);
	if(strcmp(ftr,"bajo")==0)
	{
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<N;j++)
			{
				if( freq[i][j] < ( cut-w ) )
				{
					y[i][j]=1;
				}
				else if( freq[i][j] > ( cut+w ) )
				{
					y[i][j]=0;
				}
				else
				{
					y[i][j]= 0.5 * ( 1.0 - sin( pi * ( f[i][j] - c ) ) / ( 2.0 * w ) )
				}
			}
		}
	}
	if(strcmp(ftr,"alto")==0)
	{
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<N;j++)
			{
				if( freq[i][j] < ( cut-w ) )
				{
					y[i][j]=0;
				}
				else if( freq[i][j] > ( cut+w ) )
				{
					y[i][j]=1;
				}
				else
				{
					y[i][j]= 0.5 * ( 1.0 - sin( pi * ( f[i][j] - c ) ) / ( 2.0 * w ) )
				}
			}
		}
	}
	return y;
}

