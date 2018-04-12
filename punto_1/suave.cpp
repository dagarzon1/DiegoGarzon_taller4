#include <math.h>
#include <complex>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

double ** FT(int ** y, int inverse, int M, int N);
int ** read_image(char *filename,int * height, int * width);
double ** GS(double s, int m, int n);
void imp_m(double ** y,int n, int m);
int main()
{
//n es filas
int n;
//m es columnas
int m;
int ** mat = read_image("prueba.png",&m,&n);

//double ** f = FT(mat,0,m,n);

double ** G = GS(5,m,n);

//imp_m(f,n,m);
//imp_m(G,n,m);

return 0;
}
double ** FT(int ** y, int inverse, int M, int N)
{
	//Matriz con N filas y M columnas
	double ** trans= new double*[N];
	for(int i=0;i<N;i++)
	{
		trans[i]=new double[M];
	}

	double inv=-1.0;
	if (inverse==1)
	{
		inv=1.0;
	}

	double com=inv*2.0;
	double pi = acos(-1);
	int k,l;
	for(k=0;k<N;k++)
	{	
			for(l=0;;l<M,l++)
			{
				complex<double> t=(0.0,0.0);
				int n,m;
				for(m=0;m<N;m++)
				{
					for(n=0;n<M;n++)
					{
						double fo = com * pi * ((m*k/M)+(n*l/N));
						t = t + ( (double) y[n][m] * polar(1.0,fo));
					}
				}
				double r=abs(t / ((double)N *(double)M ) );
				trans[k][l]=r;
			}
	}
	return trans;
}
int ** read_image(char *filename, int * height, int *width)
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
	fclose(f);
	return img;	
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
	double S1=40/m;
	double S2=40/n;
	int i,j;
	for(i=0,j=0;i<m || j<n;j++,i++)
	{	
		X1=inicio+( i * S1 );
		X2=inicio+( j * S2 );
		G1[i]=c * exp( - pow( X1 , 2 ) / ( 2.0 * pow( s , 2 ) ) );
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
			cout<<GS[j][i]<<endl;
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
			cout<<y[i][j]<<" "<<j<<" "<<i<<endl;
		}
	}
}



