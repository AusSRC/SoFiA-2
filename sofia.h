#pragma once

int sofia_mainline(float* dataPtr, int datasize, char* headerPtr, int headersize, char *path_to_par, int parsize, int ** channels, int* chanlen, float **moment0, int* mom0len, float **moment1, int* mom1len, float **moment2, int* mom2len, signed char **catlog, int* catlen);
 
