/* sofia.i - swig wrapper for sofia library code */
%module sofia

%{
#define SWIG_FILE_WITH_INIT
void mainline(float* dataPtr, int datasize, char* headerPtr, int headersize, char *path_to_par, int parsize);
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (float* IN_ARRAY1, int DIM1) {(float* dataPtr, int datasize)}

%inline %{
	void sofia_mainline(float* dataPtr, int datasize, char *headerPtr, int headersize, char *path_to_par, int parsize) {
		return mainline(dataPtr, datasize, headerPtr, headersize, path_to_par, parsize);
	}
%}

