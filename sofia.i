/* sofia.i - swig wrapper for sofia code */
%module sofia

%{
#define SWIG_FILE_WITH_INIT
void mainline(double* dataPtr, int datasize, char* headerPtr, int headersize, char *path_to_par, int parsize);
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (double* IN_ARRAY1, int DIM1) {(double* dataPtr, int datasize)}
/* %apply (signed char* IN_ARRAY1, int DIM1) {(char* headerPtr, int headersize)}
 %apply (char *STRING, int LENGTH) {(char *headerPtr, int headersize)} */

%inline %{
	void sofia_mainline(double* dataPtr, int datasize, char *headerPtr, int headersize, char *path_to_par, int parsize) {
		return mainline(dataPtr, datasize, headerPtr, headersize, path_to_par, parsize);
	}
%}

