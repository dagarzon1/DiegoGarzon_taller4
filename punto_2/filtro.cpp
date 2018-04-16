#include <math.h>
#include <complex>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

using namespace std;

double ** FT(double ** y, int inverse, int M, int N, double ** freq);
double ** read_image(char *filename,int * height, int * width, int *b, int *co);
void write_image(char *filen, int h, int w, int bit, int color, double ** y);
double ** FT1(double ** y, int inverse, int M, int N, double ** freq);
double ** filter(double ** y, double ** freq, int M, int N, int ftr);
double ** conv(double ** x , double ** y, int m, int n);
int main(int argc,char ** argv)
{
	//n es columnas
	int n,bit,color;
	//m es filas
	int m;
	double ** mat = read_image(argv[1],&m,&n,&bit,&color);
	double ** f= new double*[m];
	for(int i=0;i<m;i++)
	{
		f[i]=new double[n];
	}
	double ** imgft = FT(mat,0,m,n,f);
	if(strcmp(argv[2],"alto")==0)
	{
		double ** fil = filter(imgft,f,m,n,1);
		double ** S= conv( mat , fil , m , n );
		double ** r = FT(S,1,m,n,f);
		write_image("altas.png", m,n,bit, color, r);
	}
	if(strcmp(argv[2],"bajo")==0)
	{
		double ** fil = filter(imgft,f,m,n,0);
		double ** S= conv( mat , fil , m , n );
		double ** r = FT(S,1,m,n,f);
		write_image("bajas.png", m,n,bit, color, r);
	}
	return 0;
}
double ** read_image(char *filename, int * height, int *width,int *b, int *co)
{	
	FILE * f;
	int bit, color, bytes;
	png_uint_32 w,h;
	png_structp ptr;
	png_infop info;
	png_bytepp filas;
	f= fopen(filename,"r");
	ptr=png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	info=png_create_info_struct(ptr);
	png_init_io(ptr,f);
	png_read_png(ptr,info,0,0);
	png_get_IHDR(ptr,info, &w, &h, &bit, &color, NULL, NULL, NULL);
	filas=png_get_rows(ptr, info);
	double ** img=new double*[h];
	for(int i=0;i<h;i++)
	{
		img[i]=new double[w];
	}
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			img[i][j]=(double) (int) filas[i][j];
		}
	}
	*height=h;
	*width=w;
	*b=bit;
	*co=color;
	fclose(f);
	return img;	
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
double ** FT(double ** y, int inverse, int M, int N, double ** freq)
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
			res[i][j]= y[i][j] * e[i][j];		
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
double ** filter(double ** y, double ** freq, int M, int N, int ftr)
{
	int cut=1000;
	int w=50;
	double pi = acos(-1);
	double ** fil=new double*[M];
	for(int i=0;i<M;i++)
	{
		fil[i]=new double[N];
	}
	if(ftr==0)
	{
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<N;j++)
			{
				if( freq[i][j] < ( cut-w ) )
				{
					fil[i][j]=1;
				}
				else if( freq[i][j] > ( cut+w ) )
				{
					fil[i][j]=0;
				}
				else
				{
					fil[i][j]= 0.5 * ( 1.0 - sin( pi * ( freq[i][j] - cut ) ) / ( 2.0 * w ) );
				}
			}
		}
	}
	if(ftr==1)
	{
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<N;j++)
			{
				if( freq[i][j] < ( cut-w ) )
				{
					fil[i][j]=0;
				}
				else if( freq[i][j] > ( cut+w ) )
				{
					fil[i][j]=1;
				}
				else
				{
					fil[i][j]= 0.5 * ( 1.0 - sin( pi * ( freq[i][j] - cut ) ) / ( 2.0 * w ) );
				}
			}
		}
	}
	return fil;
}
double ** conv(double ** x , double ** y, int m, int n)
{
	double ** S=new double * [m];
	for(int i=0;i<m;i++)
	{
		S[i]=new double[n];
	}
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			S[i][j] = x[i][j] * y[i][j];
		}
	}
	return S;
}

