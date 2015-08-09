
#ifndef __match_hh__
#define __match_hh__

class Matrix;

// returns the cost of the assignment
double matchEdgeMaps2D(const double* bmap1, const int height, const int width, 
                     const double* bmap2, const int height2, const int width2,
                     double* m1, int m1h, int m1w, 
                     double* m2, int m2h, int m2w,
                     double maxDist, double outlierCost, int degree=6);
double matchEdgeMaps3D(const double* bmap1, const int height, const int width, const int depth, 
                     const double* bmap2, const int height2, const int width2, const int depth2,
                     double* m1, int m1h, int m1w, int m1z,
                     double* m2, int m2h, int m2w, int m2z,
                     double maxDist, double outlierCost, int degree=6);

double
matchEdgeMaps (
    Matrix& bmap1, Matrix& bmap2,
    double maxDist, double outlierCost,
    Matrix& m1, Matrix& m2);

#endif // __match_hh__
