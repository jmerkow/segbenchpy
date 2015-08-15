
#include <string.h>

#include "Matrix.hh"
#include "csa.hh"
#include "match.hh"

static const double maxDistDefault = 0.05;
static const double outlierCostDefault = 1000;

void correspondPixels(const double* bmap1, const int rows1, const int cols1, 
                      const double* bmap2, const int rows2, const int cols2,
                      double* match1, int m1, int n1, 
                      double* match2, int m2, int n2,
                      double &cost, double& oc,
                      double maxDist=0.005)
{

    // todo if rows/cols not the same
    // do the computation
    double outlierCost = outlierCostDefault;
    const double idiag = sqrt( rows1*rows2 + cols1*cols2 );
    oc = outlierCost*maxDist*idiag;
    // Matrix mat1, mat2;
    cost = matchEdgeMaps2D(
        bmap1, rows1, cols1, 
        bmap2, rows2, cols2,
        match1, m1, n1,
        match2, m2, n2,
        maxDist*idiag, oc, 6);
    
    // set output arguments
    // memcpy(match1,mat1.data(),mat1.numel()*sizeof(double));
    // memcpy(match2,mat2.data(),mat2.numel()*sizeof(double));
}

void correspondVoxels(const double* bmap1, const int rows1, const int cols1, const int aisles1,
                      const double* bmap2, const int rows2, const int cols2, const int aisles2,
                      double* match1, int m1, int n1, int k1, 
                      double* match2, int m2, int n2, int k2,
                      double &cost, double& oc,
                      double maxDist=0.005, int degree=6)
{

    // todo if rows/cols not the same
    // do the computation
    double outlierCost = outlierCostDefault;
    const double idiag = sqrt( rows1*rows2 + cols1*cols2 +aisles1*aisles2);
    oc = outlierCost*maxDist*idiag;
    // Matrix mat1, mat2;
    cost = matchEdgeMaps3D(
        bmap1, rows1, cols1, aisles1,
        bmap2, rows2, cols2, aisles2,
        match1, m1, n1, k1,
        match2, m2, n2, k2,
        maxDist*idiag, oc, degree);
    
    // set output arguments
    // memcpy(match1,mat1.data(),mat1.numel()*sizeof(double));
    // memcpy(match2,mat2.data(),mat2.numel()*sizeof(double));
}


