#include "common.h"
#include "matrix.h"
#include "io.h"

int n=10000;
int m=2000;

int a[300][300][300];

void lemon()
{
	int width=int(pow(n,1.0/3.0));
	n=width*width*width;
	
	int m = width*width*(width-1)*3-(n-1);

	int now=0;
	int flag=1, flag2=1;
	rep(k,0,width-1)
	{
		if (flag2)
			rep(i,0,width-1)
			{
				if (flag)
					repd(j,width-1,0)
						a[k][i][j]=now, now++;
				else
					rep(j,0,width-1)
						a[k][i][j]=now, now++;
				flag=1-flag;
			}
		else
			repd(i,width-1,0)
			{
				if (flag)
					repd(j,width-1,0)
						a[k][i][j]=now, now++;
				else
					rep(j,0,width-1)
						a[k][i][j]=now, now++;
				flag=1-flag;
			}
		flag2=1-flag2;
	}
	
	Mat A(n,n), T(n,n);
	rep(i,1,n-1)
	{
		double r=double(rand())/double(RAND_MAX)*10+1;
		r=1.0/r;
		A.entryAddValue(i-1,i,r); A.entryAddValue(i,i-1,r); 
		T.entryAddValue(i-1,i,r); T.entryAddValue(i,i-1,r); 
	}
	
	rep(i,0,width-1)
		rep(j,0,width-1)
			rep(k,0,width-2)
			{
				if (abs(a[i][j][k]-a[i][j][k+1])!=1)
				{
					double r=double(rand())/double(RAND_MAX)*10+1;
					r=1.0/r;
					A.entryAddValue(a[i][j][k],a[i][j][k+1],r); 
					A.entryAddValue(a[i][j][k+1],a[i][j][k],r); 
				}
				if (abs(a[i][k][j]-a[i][k+1][j])!=1)
				{
					double r=double(rand())/double(RAND_MAX)*10+1;
					r=1.0/r;
					A.entryAddValue(a[i][k][j],a[i][k+1][j],r);
					A.entryAddValue(a[i][k+1][j],a[i][k][j],r);
				}
				if (abs(a[k][i][j]-a[k+1][i][j])!=1)
				{
					double r=double(rand())/double(RAND_MAX)*10+1;
					r=1.0/r;
					A.entryAddValue(a[k][i][j],a[k+1][i][j],r);
					A.entryAddValue(a[k+1][i][j],a[k][i][j],r);
				}
			}
	
	IO::saveMM(A,"graph.mtx");
	IO::saveMM(T,"tree.mtx");
	Vec x(n);
	rep(i,0,n-1) x[i]=(double(rand())/double(RAND_MAX)-0.5)*10;
	Mat L=IO::readMMA("graph.mtx");
	Vec b=A*x;
	IO::saveMMVec(x,"x.vec");
	IO::saveMMVec(b,"b.vec");
}

int main(int argc, char *argv[])
{
	#ifndef ONLINE_JUDGE
		freopen("","r",stdin);
	#endif
	lemon();
	return 0;
}

