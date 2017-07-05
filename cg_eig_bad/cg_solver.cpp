// Copyright 2016 Haoran Xu

#include "../haoran_code/io.h"
#include "../haoran_code/graph.h"
#include "../haoran_code/cg.h"
#include "../haoran_code/jacobiprecon.h"


int main(int argc, char *argv[]) {
    string dir = argv[1];
    if (dir[dir.length()-1] != '/') dir+="/";
    Mat A = IO::readMMA(dir+"mat.mtx");
//    Mat tree = IO::readMML(dir+"tree.mtx");
    
    // Vec x(A.n);
    // rep(i,0,A.n-1) x[i]=(double(rand())/double(RAND_MAX)-0.5)*10;
    // Vec b=A*x;
    Vec x = IO::readMMVec(dir+"x.vec");
    Vec b = A*x;
    AbstractSolver S = PCG(A);

    int flag; FLOAT relres; int iter; vector<FLOAT> resvec;
    clock_t t_start = clock();
    tie(x, flag, relres, iter, resvec) = S.solve(b, 1e-6);
    clock_t t_end = clock();
    double tcost = static_cast<double>(t_end - t_start) /
        static_cast<double>(CLOCKS_PER_SEC);
    printf("n = %d iter = %d time = %.3lfs relres = %.16lf\n",
        A.n, iter, tcost, printFloat((A*x-b).norm()/b.norm()));
    return 0;
}
