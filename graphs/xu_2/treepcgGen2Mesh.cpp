#include "common.h"
#include "matrix.h"
#include "io.h"

int n=10000;
int m=2000;

int a[1000][1000];

void lemon()
{
	int width=int(sqrt(n));
	n=width*width;
	
	int m = width*(width-1)*2-(n-1);

	int now=0;
	rep(i,0,width-1) a[0][i]=now, now++;
	rep(i,1,width-1) a[i][width-1]=now, now++;
	int flag=1;
	repd(i,width-1,1)
	{
		if (flag)
			repd(j,width-2,0)
				a[i][j]=now, now++;
		else
			rep(j,0,width-2)
				a[i][j]=now, now++;
		flag=1-flag;
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
		rep(j,0,width-2)
			if (abs(a[i][j]-a[i][j+1])!=1)
			{
				double r=double(rand())/double(RAND_MAX)*10+1;
				r=1.0/r;
				A.entryAddValue(a[i][j],a[i][j+1],r); 
				A.entryAddValue(a[i][j+1],a[i][j],r); 
			}
	
	rep(i,0,width-1)
		rep(j,0,width-2)
			if (abs(a[j][i]-a[j+1][i])!=1)
			{
				double r=double(rand())/double(RAND_MAX)*10+1;
				r=1.0/r;
				A.entryAddValue(a[j][i],a[j+1][i],r); 
				A.entryAddValue(a[j+1][i],a[j][i],r); 
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

