#include <math.h>
#include <complex>
#include <vector>
#include <png.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

vector< vector<double> > FT(vector< vector <double> >& y, int inverse);
vector< vector<int> > readimage(char *filename);
int main()
{

return 0;
}
vector< vector<double> > FT(vector< vector<double> >& y, int inverse)
{
	int M= y.size();
	int N= y[0].size();
	vector< vector< complex<double> > > res(M, vector< complex<double> >(N,0) );
	vector< vector<double> > trans(M, vector < double > (N,0) );
	double inv=-1.0;
	if (inverse==1)
	{
		inv=1.0;
	}
	double com=inv*2.0;
	double pi = acos(-1);
	int k,l;
	for(k=0,l=0;k<M || l<N;k++,l++)
	{	
			complex<double> t=(0.0,0.0);
			int n,m;
			for(m=0,n=0;m<M || n<N;n++,m++)
			{
				double fo = com * pi * ((m*k/M)+(n*l/N));
				t = t + (y[m][n] * polar(1.0,fo));
			}
			res[k][l]=t;
			double r=abs(res[k][l] / ((double)N *(double)M ) );
			trans[k][l]=r;
	}
	return trans;
}
vector< vector<int> > read_image(char *filename)
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
	bytes=png_get_rowbytes(ptr,info);
	vector< vector<int> > img(h, vector < int > (w,0) );
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			img[i][j]=filas[i][j];
		}
	}
	fclose(f);
	return img;
	
}

