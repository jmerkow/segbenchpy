%module csa_swig
%include <std_string.i>
%include <std_map.i>
%include <stdint.i>
%include "exception.i"

%{
    #define SWIG_FILE_WITH_INIT
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* g, const int m, const int three)};


%include "std_vector.i"
namespace std {
  %template(VecDouble) vector<double>;
  %template(VecVecdouble) vector< vector<double> >;
}
%rename(_csaAssign) csaAssign;
%inline %{
    std::vector< std::vector<double> > csaAssign (const int n, 
                const double* g, const int m, const int three);
%}

