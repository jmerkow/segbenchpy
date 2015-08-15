%module edges_swig

%{
    #define SWIG_FILE_WITH_INIT
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}



%apply (double* IN_ARRAY2, int DIM1, int DIM2) {
(const double *E0, int h, int w), (const double *O, int ho, int wo)};

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *E, int h2, int w2)};

%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
(const double *I, int w, int h, int d),
(const double *E0, int h, int w, int d),
(const double *dX, int hdx, int wdx, int ddx),
(const double *dY, int hdy, int wdy, int ddy),
(const double *dZ, int hdz, int wdz, int ddz)};

%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
(double *E, int h2, int w2, int d2)};

%rename(_edgeNms2d) edgeNms2d;
%rename(_interp3) interp3;
%rename(_edgeNms3d) edgeNms3d;
%inline %{
void edgeNms2d(const double *E0, int h, int w,
               const double *O, int ho, int wo,
               double *E, int h2, int w2,
               int r, int s, float m);
double interp3( const double *I, int w, int h, int d, 
    double y, double x, double z );

void edgeNms3d(
    const double *E0, int h, int w, int d,
    const double *dX, int hdx, int wdx, int ddx,
    const double *dY, int hdy, int wdy, int ddy,
    const double *dZ, int hdz, int wdz, int ddz,
    double *E, int h2, int w2, int d2,
    float r, int s, float m);
%}