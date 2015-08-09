%module benchmark_swig

%{
    #define SWIG_FILE_WITH_INIT
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}



%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* bmap1, const int rows1, const int cols1), 
                                                (double* bmap2, const int rows2, const int cols2)};

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* bmap1, const int height, const int width),
                                                (const double* bmap2, const int height2, const int width2)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* match1, int m1, int n1),(double* match2, int m2, int n2)};

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* m1, int m1h, int m1w),(double* m2, int m2h, int m2w)};


%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) 
                                                {(const double* bmap1, const int height, const int width, const int depth),
                                                (const double* bmap2, const int height2, const int width2, const int depth2)};
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) 
                                                {(double* m1, int m1h, int m1w, int m1z),
                                                 (double* m2, int m2h, int m2w, int m2z) };


%include "typemaps.i"
%apply double& OUTPUT {double& cost};
%apply double& OUTPUT {double& oc};


%inline %{
void correspondPixels(double* bmap1, const int rows1, const int cols1, 
                      double* bmap2, const int rows2, const int cols2,
                      double* match1, int m1, int n1, 
                      double* match2, int m2, int n2,
                      double maxDist,
                      double &cost, double& oc);
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
%}