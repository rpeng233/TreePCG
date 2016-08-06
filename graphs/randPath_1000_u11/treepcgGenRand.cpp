#include "common.h"
#include "matrix.h"
#include "io.h"

int n=1000;
int m=2000;

void lemon()
{
	Mat A(n,n), T(n,n);
	rep(i,1,n-1)
	{
		double r=double(rand())/double(RAND_MAX)*10+1;
		r=1.0/r;
		A.entryAddValue(i-1,i,r); A.entryAddValue(i,i-1,r); 
		T.entryAddValue(i-1,i,r); T.entryAddValue(i,i-1,r); 
	}
	rep(i,1,m-n+1)
	{
		double r=double(rand())/double(RAND_MAX)*10+1;
		r=1.0/r;
		int x,y;
		while (1)
		{
			x=rand()%n; y=rand()%n;
			if (x!=y) break;
		}
		A.entryAddValue(x,y,r); A.entryAddValue(y,x,r);
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

