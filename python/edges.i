%module edges_swig

%{
    #define SWIG_FILE_WITH_INIT
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}



%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double *E0, int h, int w), (const double *O, int ho, int wo)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *E, int h2, int w2)};

%rename(_edgeNms2d) edgeNms2d;
%inline %{
void edgeNms2d(const double *E0, int h, int w,
               const double *O, int ho, int wo,
               double *E, int h2, int w2,
               int r, int s, float m);
%}