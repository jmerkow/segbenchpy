%module csa
%include <std_vector.i>
%include <std_string.i>
%include <std_map.i>
%include <stdint.i>
%include "exception.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "csa.hh"
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
  import_array();
%}

%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int m, int k, const int* graph)};
//%apply (int* IN_ARRAY1, int DIM1) {(int* graph, int len1)};

%include "typemaps.i"
%apply int* INPUT {const int* graph};
%apply int& OUTPUT {int& a};
%apply int& OUTPUT {int& b};
%apply int& OUTPUT {int& cost};


%include "csa.hh"