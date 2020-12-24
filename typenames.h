#include <stdio.h>
#include <stddef.h>
#include <stdint.h>


#define typename(x) _Generic((x),        /* Get the name of a type */             \
                                                                                  \
        _Bool: "_Bool",                  unsigned char: "unsigned char",          \
         char: "char",                     signed char: "signed char",            \
    short int: "short int",         unsigned short int: "unsigned short int",     \
          int: "int",                     unsigned int: "unsigned int",           \
     long int: "long int",           unsigned long int: "unsigned long int",      \
long long int: "long long int", unsigned long long int: "unsigned long long int", \
        float: "float",                         double: "double",                 \
  long double: "long double",                   char *: "pointer to char",        \
       void *: "pointer to void",                int *: "pointer to int",         \
	 double *: "pointer to double",      long double *: "pointer to long double", \
      void **: "pointer to void *",             int **: "pointer to int *",       \
	double **: "pointer to double *",   long double **: "pointer to long double *", \
      char **: "pointer to char *",            default: "other")
