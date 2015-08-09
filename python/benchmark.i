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
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* match1, int m1, int n1)}; 
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* match2, int m2, int n2)};

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
%}