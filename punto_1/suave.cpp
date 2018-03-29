#include <iostream>
#include <string.h>
#include <math.h>
#include <complex>
#include <vector>
#include <png.h>


using namespace std;

vector<vector<double> > FT(vector< vector <double> >& y, int inverse);
int main()
{

return 0;
}
vector<vector<double> > FT(vector< vector<double> >& y, int inverse)
{
	int M= y.size();
	int N= y[0].size();
	vector< vector< complex<double> > > res(M, vector< complex<double> >(N,0) );
	vector< vector<double> > trans(M, vector < double > (N,0) );
	double inv=1.0;
	if (inverse==1)
	{
		inv=-1.0;
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
