/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (DataCube.h) - Source Finding Application                ///
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

#ifndef DATACUBE_H
#define DATACUBE_H

#include <stdlib.h>
#include <stdbool.h>

#include "common.h"
#include "String.h"
#include "Stack.h"
#include "Array_dbl.h"
#include "Array_siz.h"
#include "Map.h"
#include "Catalog.h"
#include "LinkerPar.h"
#include "Header.h"
#include "WCS.h"

#define DESTROY  false
#define PRESERVE true
typedef enum {NOISE_STAT_STD, NOISE_STAT_MAD, NOISE_STAT_GAUSS} noise_stat;

// ----------------------------------------------------------------- //
// DATATYPE                                                          //
// ----------------------------------------------------------------- //
// This enum defines the source to use for the data to read into    //
// the cube. The default is to read a FITS format file to retrieve   //
// data and metadata (in the header). Other types are:               //
//                                                                   //
//		"MEM" - use a pointer to a flat array of floats and a        //
//              pointer to a char array of key/value pairs.          //
//		"SOURCETYPERANGE" - should never equal this - keep this as the    //
//					   last element of the enum.                     //
// ----------------------------------------------------------------- //
typedef enum {FITS, MEM, SOURCETYPERANGE} SOURCETYPE;

// ----------------------------------------------------------------- //
// Class 'DataCube'                                                  //
// ----------------------------------------------------------------- //
// The purpose of this class is to handle up to three-dimensional    //
// astronomical data cubes. The class is intended for reading and    //
// manipulating FITS data cubes by providing methods for loading and //
// saving FITS files and manipulating the header and data units of a //
// FITS file. Currently, only single-HDU files are supported.        //
// ----------------------------------------------------------------- //

typedef CLASS DataCube DataCube;

// Constructor and destructor
PUBLIC DataCube  *DataCube_new              (const bool verbosity);
PUBLIC DataCube  *DataCube_copy             (const DataCube *source);
PUBLIC DataCube  *DataCube_blank            (const size_t nx, const size_t ny, const size_t nz, const int type, const bool verbosity);
PUBLIC void       DataCube_delete           (DataCube *self);

// Public methods
// Loading data into the cube
PUBLIC void       DataCube_load             (DataCube *self, char *dataSrc, const Array_siz *region, SOURCETYPE source, char *header);
PUBLIC void       DataCube_readFITS         (DataCube *self, const char *filename, const Array_siz *region);
PUBLIC void       DataCube_readMEM          (DataCube *self, float *dataPtr, int datasize, char *headerPtr, int headersize);
// Saving to FITS format or memory
PUBLIC void       DataCube_save             (const DataCube *self, const char *filename, const bool overwrite, const bool preserve, const OUTTYPE OUTPUTS);
PUBLIC size_t     DataCube_writeMemFloat    (const DataCube *self, float **memPtr, char * reallocOK);
PUBLIC size_t     DataCube_writeMemInt      (const DataCube *self, int **memPtr, char * reallocOK);

// Getting basic information
PUBLIC size_t     DataCube_get_size         (const DataCube *self);
PUBLIC size_t     DataCube_get_axis_size    (const DataCube *self, const size_t axis);

// Wrappers for relevant methods of class Header
PUBLIC long int   DataCube_gethd_int        (const DataCube *self, const char *key);
PUBLIC double     DataCube_gethd_flt        (const DataCube *self, const char *key);
PUBLIC bool       DataCube_gethd_bool       (const DataCube *self, const char *key);
PUBLIC int        DataCube_gethd_str        (const DataCube *self, const char *key, char *value);
PUBLIC String    *DataCube_gethd_string     (const DataCube *self, const char *key);
PUBLIC int        DataCube_puthd_int        (DataCube *self, const char *key, const long int value);
PUBLIC int        DataCube_puthd_flt        (DataCube *self, const char *key, const double value);
PUBLIC int        DataCube_puthd_bool       (DataCube *self, const char *key, const bool value);
PUBLIC int        DataCube_puthd_str        (DataCube *self, const char *key, const char *value);
PUBLIC size_t     DataCube_chkhd            (const DataCube *self, const char *key);
PUBLIC bool       DataCube_cmphd            (const DataCube *self, const char *key, const char *value, const size_t n);
PUBLIC int        DataCube_delhd            (DataCube *self, const char *key);
PUBLIC void       DataCube_copy_wcs         (const DataCube *source, DataCube *target);

// Extract data values
PUBLIC double     DataCube_get_data_flt     (const DataCube *self, const size_t x, const size_t y, const size_t z);
PUBLIC long int   DataCube_get_data_int     (const DataCube *self, const size_t x, const size_t y, const size_t z);

// Manipulate data values
PUBLIC void       DataCube_set_data_flt     (DataCube *self, const size_t x, const size_t y, const size_t z, const double value);
PUBLIC void       DataCube_set_data_int     (DataCube *self, const size_t x, const size_t y, const size_t z, const long int value);
PUBLIC void       DataCube_add_data_flt     (DataCube *self, const size_t x, const size_t y, const size_t z, const double value);
PUBLIC void       DataCube_add_data_int     (DataCube *self, const size_t x, const size_t y, const size_t z, const long int value);
PUBLIC void       DataCube_fill_flt         (DataCube *self, const double value);

// Arithmetic operations
PUBLIC void       DataCube_divide           (DataCube *self, const DataCube *divisor);
PUBLIC void       DataCube_apply_weights    (DataCube *self, const DataCube *weights);
PUBLIC void       DataCube_multiply_const   (DataCube *self, const double factor);
PUBLIC void       DataCube_add_const        (DataCube *self, const double summand);

// Statistical measurements
PUBLIC double     DataCube_stat_std         (const DataCube *self, const double value, const size_t cadence, const int range);
PUBLIC double     DataCube_stat_mad         (const DataCube *self, const double value, const size_t cadence, const int range);
PUBLIC double     DataCube_stat_gauss       (const DataCube *self, const size_t cadence, const int range);

// Noise scaling
PUBLIC void       DataCube_scale_noise_spec (const DataCube *self, const noise_stat statistic, const int range);
PUBLIC DataCube  *DataCube_scale_noise_local(DataCube *self, const noise_stat statistic, const int range, size_t window_spat, size_t window_spec, size_t grid_spat, size_t grid_spec, const bool interpolate);

// Spatial and spectral smoothing
PUBLIC void       DataCube_boxcar_filter    (DataCube *self, size_t radius);
PUBLIC void       DataCube_gaussian_filter  (DataCube *self, const double sigma);

// Continuum subtraction and ripple removal
PUBLIC void       DataCube_contsub          (DataCube *self, unsigned int order, size_t shift, const size_t padding, double threshold);
PUBLIC void       DataCube_spatial_filter   (DataCube *self, const int statistic, const size_t window_spat, const size_t radius_spec);

// Masking
PUBLIC void       DataCube_mask             (const DataCube *self, DataCube *maskCube, const double threshold);
PUBLIC void       DataCube_mask_8           (const DataCube *self, DataCube *maskCube, const double threshold, const uint8_t value);
PUBLIC void       DataCube_set_masked       (DataCube *self, const DataCube *maskCube, const double value);
PUBLIC void       DataCube_set_masked_8     (DataCube *self, const DataCube *maskCube, const double value);
PUBLIC void       DataCube_reset_mask_32    (DataCube *self, const int32_t value);
PUBLIC void       DataCube_filter_mask_32   (DataCube *self, const Map *filter);
PUBLIC size_t     DataCube_copy_mask_8_32   (DataCube *self, const DataCube *source, const int32_t value);
PUBLIC void       DataCube_dilate_mask_xy   (const DataCube *self, DataCube *mask, Catalog *cat, const size_t iter_max, const double threshold);
PUBLIC void       DataCube_dilate_mask_z    (const DataCube *self, DataCube *mask, Catalog *cat, const size_t iter_max, const double threshold);
PUBLIC DataCube  *DataCube_2d_mask          (const DataCube *self);

// Flagging
PUBLIC void       DataCube_flag_regions     (DataCube *self, const Array_siz *region);
PUBLIC void       DataCube_copy_blanked     (DataCube *self, const DataCube *source);
PUBLIC void       DataCube_autoflag         (const DataCube *self, const double threshold, const unsigned int mode, Array_siz *region);
PUBLIC size_t     DataCube_flag_infinity    (const DataCube *self, Array_siz *region);

// Source finding
PUBLIC void       DataCube_run_scfind       (const DataCube *self, DataCube *maskCube, const Array_dbl *kernels_spat, const Array_siz *kernels_spec, const double threshold, const double maskScaleXY, const noise_stat method, const int range, const int scaleNoise, const noise_stat snStatistic, const int snRange, const size_t snWindowXY, const size_t snWindowZ, const size_t snGridXY, const size_t snGridZ, const bool snInterpol, const time_t start_time, const clock_t start_clock);
PUBLIC void       DataCube_run_threshold    (const DataCube *self, DataCube *maskCube, const bool absolute, double threshold, const noise_stat method, const int range);

// Linking
PUBLIC LinkerPar *DataCube_run_linker       (const DataCube *self, DataCube *mask, const size_t radius_x, const size_t radius_y, const size_t radius_z, const size_t min_size_x, const size_t min_size_y, const size_t min_size_z, const size_t max_size_x, const size_t max_size_y, const size_t max_size_z, const bool positivity, const double rms);

// Parameterisation
PUBLIC void       DataCube_parameterise     (const DataCube *self, const DataCube *mask, Catalog *cat, bool use_wcs, bool physical, const char *prefix);

// Create moment maps and cubelets
PUBLIC void       DataCube_create_moments   (const DataCube *self, const DataCube *mask, DataCube **mom0, DataCube **mom1, DataCube **mom2, DataCube **chan, const char *obj_name, bool use_wcs, const bool positive);
PUBLIC void       DataCube_create_cubelets  (const DataCube *self, const DataCube *mask, const Catalog *cat, const char *basename, const bool overwrite, bool use_wcs, bool physical, const size_t margin,const OUTTYPE OUTPUTS);

// WCS
PUBLIC WCS       *DataCube_extract_wcs      (const DataCube *self);

// Private methods
PRIVATE inline size_t DataCube_get_index       (const DataCube *self, const size_t x, const size_t y, const size_t z);
PRIVATE        void   DataCube_get_xyz         (const DataCube *self, const size_t index, size_t *x, size_t *y, size_t *z);
PRIVATE        void   DataCube_process_stack   (const DataCube *self, DataCube *mask, Stack *stack, const size_t radius_x, const size_t radius_y, const size_t radius_z, const int32_t label, LinkerPar *lpar, const double rms);
PRIVATE        void   DataCube_grow_mask_xy    (const DataCube *self, DataCube *mask, const long int src_id, const size_t radius, const long int mask_value, double *f_sum, double *f_min, double *f_max, size_t *n_pix, long int *flag, size_t *x_min, size_t *x_max, size_t *y_min, size_t *y_max, const size_t z_min, const size_t z_max);
PRIVATE        double DataCube_get_beam_area   (const DataCube *self);
PRIVATE        void   DataCube_get_wcs_info    (const DataCube *self, String **unit_flux_dens, String **unit_flux, String **label_lon, String **label_lat, String **label_spec, String **ucd_lon, String **ucd_lat, String **ucd_spec, String **unit_lon, String **unit_lat, String **unit_spec, double *beam_area, double *chan_size);
PRIVATE        void   DataCube_create_src_name (const DataCube *self, String **source_name, const char *prefix, const double longitude, const double latitude, const String *label_lon);
PRIVATE        void   DataCube_swap_byte_order (const DataCube *self);
PRIVATE        size_t DataCube_countNonZeros   (const DataCube *self);
PRIVATE        size_t DataCube_sparseFloatCsr  (const DataCube *self, float **dataPtr);
PRIVATE        size_t DataCube_sparseIntCsr    (const DataCube *self, int **dataPtr);

// TEST
PUBLIC void DataCube_continuum_flagging(DataCube *self, const char *filename, const int coord_system, const long int radius);

#endif
