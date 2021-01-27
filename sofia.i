/* sofia.i - swig wrapper interface declarations for sofia library code */
%module sofia

%{
  #define SWIG_FILE_WITH_INIT
  #include "sofia.h"
  #include <setjmp.h>
  static __thread int infunc = 0;
  static __thread jmp_buf buf;

  static void exitFromC(int code, void *data) {
    if (!infunc) return;
    (void)data;
    longjmp(buf,code);
  }
%}
%include "numpy.i"
%init %{
  import_array();
  on_exit(exitFromC, NULL);

%}
%apply (float* INPLACE_ARRAY_FLAT, int DIM_FLAT) {(float* dataPtr, int datasize)}
%apply (int ** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int ** channels, int * chanlen)}
%apply (float ** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float ** moment0, int * mom0len)}
%apply (float ** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float ** moment1, int * mom1len)}
%apply (float ** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float ** moment2, int * mom2len)}
%apply (signed char ** ARGOUTVIEWM_ARRAY1, int* DIM1) {(signed char ** catlog, int * catlen)}

%exception {
  infunc = 1;
  int err = 0;
  if (!(err=setjmp(buf))) {
    $action
  }
  else {
    // Raise exception, code=err
    if (err == 2) PyErr_SetString(PyExc_StopIteration,"Null Pointer error");
    else if (err == 3) PyErr_SetString(PyExc_MemoryError,"Memory allocation error");
    else if (err == 4) PyErr_SetString(PyExc_IndexError,"Range error");
    else if (err == 5) PyErr_SetString(PyExc_IOError,"File access error");
    else if (err == 6) PyErr_SetString(PyExc_OverflowError,"Integer overflow error");
    else if (err == 7) PyErr_SetString(PyExc_TypeError,"User input error - check parameters");
    else if (err == 8) PyErr_SetString(PyExc_SystemExit,"No sources found!");
    else PyErr_SetString(PyExc_RuntimeError,"General Error");
    infunc = 0;
    on_exit(exitFromC, NULL);
    SWIG_fail;
  }
  infunc = 0;
}


%include "sofia.h"
/*
%inline %{
	int sofia_mainline(float* dataPtr, int datasize, char *headerPtr, int headersize, char *path_to_par, int parsize, int ** channels, int* chanlen, float **moment0, int* mom0len, float **moment1, int* mom1len, float **moment2, int* mom2len, signed char **catlog, int* catlen) {
		return mainline(dataPtr, datasize, headerPtr, headersize, path_to_par, parsize, channels,chanlen,moment0,mom0len,moment1,mom1len,moment2,mom2len,catlog,catlen);
	}
%}
*/

