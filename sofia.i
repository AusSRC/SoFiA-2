/* sofia.i - swig wrapper for sofia code */
%module sofia
%{
#define SWIG_FILE_WITH_INIT
int mainline(double* dataPtr, int datasize);
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (float* INPLACE_ARRAY1, int DIM1) {(double* dataPtr, int datasize)}
int mainline(double* dataPtr,int datasize);
