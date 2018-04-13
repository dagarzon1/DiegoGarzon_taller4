#include <math.h>
#include <complex>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

double ** FT(int ** y, int inverse, int M, int N);
int ** read_image(char *filename,int * height, int * width, int *b, int *co);
double ** GS(double s, int m, int n);
void imp_m(double ** y,int n, int m);
double ** conv(double ** x , double ** y, int m, int n);
double ** FT1(double ** y, int inverse, int M, int N);
int main()
{
//n es columnas
int n,bit,color;
//m es filas
int m;
int ** mat = read_image("prueba.png",&m,&n,&bit,&color);
double ** G = GS(1,m,n);
double ** f = FT(mat,0,m,n);
double ** f_1 = FT1(G,0,m,n);
double ** c = conv(f,f_1,m,n);
double ** r = FT1(c,1,m,n);

imp_m(r,n,m);

return 0;
}
double ** FT1(double ** y, int inverse, int M, int N)
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
	complex<double> ** e = new complex<double> * [N];
	for(int i =0;i<N;i++)
	{
		e[i]=new complex<double>[M];
	}
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			e[j][i]= e_2[i] * e_1[j]/c;
		}
	}
	for(int k=0;k<M;k++)
	{	
			for(int l=0;l<N;l++)
			{
				complex<double> t=0.0;
				for(int m=0;m<M;m++)
				{
					for(int n=0;n<N;n++)
						{
							t= t + ( (double) y[m][n] * e[m][n] );
						}
				}
				trans[k][l]=abs( t );
			}
	}
	return trans;
}
double ** FT(int ** y, int inverse, int M, int N)
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
	complex<double> ** e = new complex<double> * [N];
	for(int i =0;i<N;i++)
	{
		e[i]=new complex<double>[M];
	}
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			e[j][i]= e_2[i] * e_1[j]/c;
		}
	}
	for(int k=0;k<M;k++)
	{	
			for(int l=0;l<N;l++)
			{
				complex<double> t=0.0;
				for(int m=0;m<M;m++)
				{
					for(int n=0;n<N;n++)
						{
							t= t + ( (double) y[m][n] * e[m][n] );
						}
				}
				trans[k][l]=abs( t );
			}
	}
	return trans;
}
int ** read_image(char *filename, int * height, int *width,int *b, int *co)
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
	int ** img=new int*[h];
	for(int i=0;i<h;i++)
	{
		img[i]=new int[w];
	}
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			img[i][j]=filas[i][j];
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
	pg_init_io(ptr,f);
	png_set_IHDR(ptr,info,w,h,bit,color,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
	png_write_info(ptr,info);
	int ** img=new int*[h];
	for(int i=0;i<h;i++)
	{
		img[i]=new int[w];
	}
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			img[i][j]=(int) y[i][j];
		}
	}
	png_write_image(ptr,img);
	png_write_end(ptr,NULL); 	
}
double ** GS(double s, int m, int n)
{
	double pi=3.14159265359;
	double c=1.0 / ( s * pow( 2.0 * pi , 0.5 ) );
	double * G1=new double[m];
	double * G2=new double[n];
	double X1=0.0;
	double X2=0.0;
	double inicio=-20;
	double S1=40.0/m;
	double S2=40.0/n;
	for(int i=0;i<m;i++)
	{	
		X1=inicio+( i * S1 );
		G1[i]=c * exp( - pow( X1 , 2 ) / ( 2.0 * pow( s , 2 ) ) );
	}
	for(int j=0;j<n;j++)
	{
		X2=inicio+( j * S2 );
		G2[j]=c * exp( - pow( X2 , 2 ) / ( 2.0 * pow( s , 2 ) ) );
	}
	double ** GS=new double * [m];
	for(int i=0;i<m;i++)
	{
		GS[i]=new double[n];
	}
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			GS[j][i]= G2[i] * G1[j];
		}
	}
	return GS;
}
void imp_m(double ** y,int n,int m)
{
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout<<y[i][j]<<" ";
		}
		cout<<endl;
	}
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



