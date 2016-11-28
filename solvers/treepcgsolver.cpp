// Copyright 2016 Haoran Xu

#include "../haorancode/io.h"
#include "../haorancode/graph.h"
#include "../haorancode/cg.h"
#include "../haorancode/jacobiprecon.h"
#include "../haorancode/treeprecon.h"

int main(int argc, char *argv[]) {
    string dir = argv[1];
    if (dir[dir.length()-1] != '/') dir+="/";
    Mat A = IO::readMMA(dir+"graph.mtx");
    Mat tree = IO::readMML(dir+"tree.mtx");
    GraphSP g = IO::specifyTree(IO::convertLaplacianMatrixToGraph(A), tree);
    // Vec x(A.n);
    // rep(i,0,A.n-1) x[i]=(double(rand())/double(RAND_MAX)-0.5)*10;
    // Vec b=A*x;
    Vec x = IO::readMMVec(dir+"x.vec");
    Vec b = A*x;
    AbstractSolver S = PCG(A, TreePreconditioner(g));

    printf("Avgstr = %.16lf\n",
        printFloat(StretchCalculator::calculateTotalStretch(g) / g.o.size()));
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
