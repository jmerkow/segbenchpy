
#include <string.h>

#include "Matrix.hh"
#include "csa.hh"
#include "match.hh"

static const double maxDistDefault = 0.05;
static const double outlierCostDefault = 100;

void correspondPixels(double* bmap1, const int rows1, const int cols1, 
                      double* bmap2, const int rows2, const int cols2,
                      double* match1, int m1, int n1, 
                      double* match2, int m2, int n2,
                      double maxDist,
                      double &cost, double& oc)
{

    // todo if rows/cols not the same
    // do the computation
    double outlierCost = outlierCostDefault;
    const double idiag = sqrt( rows1*rows2 + cols1*cols2 );
    oc = outlierCost*maxDist*idiag;
    Matrix mat1, mat2;
    cost = matchEdgeMaps(
        Matrix(rows1,cols1,bmap1), Matrix(rows2,cols2,bmap2),
        maxDist*idiag, oc,
        mat1, mat2);
    
    // set output arguments
    memcpy(match1,mat1.data(),mat1.numel()*sizeof(double));
    memcpy(match2,mat2.data(),mat2.numel()*sizeof(double));
}




