/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (common.h) - Source Finding Application                  ///
/// Copyright (C) 2020 Tobias Westmeier                                  ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// Address:  Tobias Westmeier                                           ///
///           ICRAR M468                                                 ///
///           The University of Western Australia                        ///
///           35 Stirling Highway                                        ///
///           Crawley WA 6009                                            ///
///           Australia                                                  ///
///                                                                      ///
/// E-mail:   tobias.westmeier [at] uwa.edu.au                           ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// This program is free software: you can redistribute it and/or modify ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// This program is distributed in the hope that it will be useful,      ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
/// ____________________________________________________________________ ///
///                                                                      ///

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdbool.h>
#include <time.h>

// SoFiA version number
#define SOFIA_VERSION "2.2.1"
#define SOFIA_VERSION_FULL "SoFiA 2.2.1"
#define SOFIA_CREATION_DATE "18-Nov-2020"

// Define value of pi
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

// Define value of conversion factor between MAD and standard deviation
// of normal distribution, calculated as 1.0 / scipy.stats.norm.ppf(3.0 / 4.0)
#ifndef MAD_TO_STD
#define MAD_TO_STD 1.482602218505602
#endif

// Define value of 1 / sqrt(2 * pi)
#ifndef INV_SQRT_TWO_PI
#define INV_SQRT_TWO_PI 0.3989422804014327
#endif

// Check for NaN
#ifndef IS_NAN
#define IS_NAN(x) ((x) != (x))
#endif
#ifndef IS_NOT_NAN
#define IS_NOT_NAN(x) ((x) == (x))
#endif
#ifndef FILTER_NAN
#define FILTER_NAN(x) ((x) == (x) ? (x) : 0)
#endif

// Check if odd
#ifndef IS_ODD
#define IS_ODD(x) ((x) & 1)
#endif
#ifndef IS_EVEN
#define IS_EVEN(x) (!((x) & 1))
#endif

// Define memory allocation modes
#define MALLOC 0
#define CALLOC 1

// Define maximum RMS measurement sample size
// NOTE: This is chosen to be a prime number to reduce the risk of
//       obtaining a stride that is a multiple of the x-axis size.
#define NOISE_SAMPLE_SIZE 999983

// Define size of kB, MB and GB
#define KILOBYTE       1024
#define MEGABYTE    1048576
#define GIGABYTE 1073741824

// Define object-oriented terminology
#define CLASS struct
#define PUBLIC extern
#define PRIVATE static

// Define error codes
#define ERR_SUCCESS      0
#define ERR_FAILURE      1
#define ERR_NULL_PTR     2
#define ERR_MEM_ALLOC    3
#define ERR_INDEX_RANGE  4
#define ERR_FILE_ACCESS  5
#define ERR_INT_OVERFLOW 6
#define ERR_USER_INPUT   7
#define ERR_NO_SRC_FOUND 8

// ----------------------------------------------------------------- //
// OUTTYPE                                                           //
// ----------------------------------------------------------------- //
// This enum defines the output type use for results data            //
// The default is to write a FITS format file to the filesystem.     //
// Other types are:                                                  //
//                                                                   //
//		"NONE" - no output data will be saved (for testing only)     //
//		"BUFF" - output data will be written to memory buffers       //
//		"OUTTYPERANGE" - should never equal this - keep this as the  //
//					   last element of the enum.                     //
// ----------------------------------------------------------------- //
typedef enum {FILES, NONE, BUFF, OUTTYPERANGE} OUTTYPE;

// Generic compile time check; should result in a compiler error if
// condition is false due to attempt to create array of negative size.
// NOTE: This does not actually create a physical array, but merely
//       defines a new type.
#define COMPILE_TIME_CHECK(condition, message) typedef char message[(condition) ? 1 : -1]

// Check condition and exit if not met
void ensure(const bool condition, const int errorCode, const char *format, ...);
void check_null(const void *ptr);

// Print info and warning messages
void message(const char *format, ...);
void message_verb(const bool verbosity, const char *format, ...);
void status(const char *format, ...);
void warning(const char *format, ...);
void warning_verb(const bool verbosity, const char *format, ...);

// Display progress bar and time stamp
void progress_bar(const char *text, const size_t progress, const size_t maximum);
void timestamp(const time_t start, const clock_t start_clock);

// Memory allocation
void *memory(const int mode, const size_t n_blocks, const size_t block_size);
void *memory_realloc(void *ptr, const size_t n_blocks, const size_t block_size);

// String functions
char *trim_string(char *str);
void int_to_str(char *str, const size_t size, const long int value);

// Swap two values
void swap(double *val1, double *val2);

// Plotting aids
double auto_tick(const double range, const size_t n);
void write_eps_header(FILE *fp, const char *title, const char *creator, const char *bbox);
void write_eps_footer(FILE *fp);

// Byte-order functions
bool is_little_endian(void);
void swap_byte_order(char *word, const size_t size);

#endif
