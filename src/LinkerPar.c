/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (LinkerPar.c) - Source Finding Application               ///
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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "LinkerPar.h"
#include "String.h"

// Set to 1 if measurement of flux-weighted centroid required
#define MEASURE_CENTROID_POSITION 0



// ----------------------------------------------------------------- //
// Declaration of properties of class LinkerPar                      //
// ----------------------------------------------------------------- //

CLASS LinkerPar
{
	size_t  size;
	int     verbosity;
	size_t *label;
	size_t *n_pix;
	size_t *x_min;
	size_t *x_max;
	size_t *y_min;
	size_t *y_max;
	size_t *z_min;
	size_t *z_max;
	#if MEASURE_CENTROID_POSITION
	double *x_ctr;
	double *y_ctr;
	double *z_ctr;
	#endif
	double *f_min;
	double *f_max;
	double *f_sum;
	double *rel;
	unsigned char *flags;
};



// ----------------------------------------------------------------- //
// Standard constructor                                              //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   No arguments.                                                   //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Pointer to newly created LinkerPar object.                      //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Standard constructor. Will create a new and empty LinkerPar ob- //
//   ject and return a pointer to the newly created object. No me-   //
//   mory will be allocated other than for the object itself. Note   //
//   that the destructor will need to be called explicitly once the  //
//   object is no longer required to release any memory allocated    //
//   during the lifetime of the object.                              //
// ----------------------------------------------------------------- //

PUBLIC LinkerPar *LinkerPar_new(const bool verbosity)
{
	LinkerPar *self = (LinkerPar *)memory(MALLOC, 1, sizeof(LinkerPar));
	
	self->verbosity = verbosity;
	self->size = 0;
	
	self->label = NULL;
	self->n_pix = NULL;
	self->x_min = NULL;
	self->x_max = NULL;
	self->y_min = NULL;
	self->y_max = NULL;
	self->z_min = NULL;
	self->z_max = NULL;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr = NULL;
	self->y_ctr = NULL;
	self->z_ctr = NULL;
	#endif
	self->f_min = NULL;
	self->f_max = NULL;
	self->f_sum = NULL;
	self->rel   = NULL;
	self->flags = NULL;
	
	return self;
}



// ----------------------------------------------------------------- //
// Destructor                                                        //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Destructor. Note that the destructor must be called explicitly  //
//   if the object is no longer required. This will release the me-  //
//   mory occupied by the object.                                    //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_delete(LinkerPar *self)
{
	if(self != NULL)
	{
		self->size = 0;
		LinkerPar_reallocate_memory(self);
		free(self);
	}
	
	return;
}



// ----------------------------------------------------------------- //
// Return number of sources in LinkerPar object                      //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Number of sources stored in LinkerPar object.                   //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the size of the LinkerPar object,   //
//   i.e. the number of sources it currently contains.               //
// ----------------------------------------------------------------- //

PUBLIC size_t LinkerPar_get_size(const LinkerPar *self)
{
	return self == NULL ? 0 : self->size;
}



// ----------------------------------------------------------------- //
// Insert a new object at the end of the current list                //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Label of the new object.                         //
//   (3) x        - x-position of the new object.                    //
//   (4) y        - y-position of the new object.                    //
//   (5) z        - z-position of the new object.                    //
//   (6) flux     - Flux value of the new object.                    //
//   (7) flag     - Flag values to be set; 1 = spatial edge;         //
//                  2 = spectral edge; 4 = blanked pixels; 8 = other //
//                  sources.                                         //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for adding a new object to the end of the current //
//   list pointed to by 'self'. The label will be set to 'label',    //
//   the number of pixels to 1, and the (x, y, z) position will be   //
//   used as the initial x_min, x_max, etc. The memory allocation of //
//   the object will automatically be expanded if necessary.         //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_push(LinkerPar *self, const size_t label, const size_t x, const size_t y, const size_t z, const double flux, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	
	// Increment size counter
	++self->size;
	
	// Allocate additional memory
	LinkerPar_reallocate_memory(self);
	
	// Insert new element at end
	self->label[self->size - 1] = label;
	self->n_pix[self->size - 1] = 1;
	self->x_min[self->size - 1] = x;
	self->x_max[self->size - 1] = x;
	self->y_min[self->size - 1] = y;
	self->y_max[self->size - 1] = y;
	self->z_min[self->size - 1] = z;
	self->z_max[self->size - 1] = z;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr[self->size - 1] = flux * x;
	self->y_ctr[self->size - 1] = flux * y;
	self->z_ctr[self->size - 1] = flux * z;
	#endif
	self->f_min[self->size - 1] = flux;
	self->f_max[self->size - 1] = flux;
	self->f_sum[self->size - 1] = flux;
	self->rel  [self->size - 1] = 0.0;  // NOTE: Must be 0 (default for neg. sources), as only pos. sources will be updated later!
	self->flags[self->size - 1] = flag;
	
	return;
}



// ----------------------------------------------------------------- //
// Remove last object from list                                      //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for removing the most recently added object from  //
//   the list. The process will be terminated if the original list   //
//   is empty.                                                       //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_pop(LinkerPar *self)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_FAILURE, "Failed to pop element from empty LinkerPar object.");
	
	// Decrement size
	--self->size;
	
	// Reallocate memory
	LinkerPar_reallocate_memory(self);
	
	return;
}



// ----------------------------------------------------------------- //
// Add another pixel to the last object in the list                  //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) x        - x-position of the new pixel.                     //
//   (3) y        - y-position of the new pixel.                     //
//   (4) z        - z-position of the new pixel.                     //
//   (5) flux     - Flux value of the new pixel.                     //
//   (6) flag     - Flag values to be set; 1 = spatial edge;         //
//                  2 = spectral edge; 4 = blanked pixels; 8 = other //
//                  sources.                                         //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for adding another pixel to the last object in    //
//   the current linker list. The object's x_min, x_max, y_min, etc. //
//   values will be checked against the newly added pixel and up-    //
//   dated if necessary. The programme will terminate if the list is //
//   found to be empty.                                              //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_update(LinkerPar *self, const size_t x, const size_t y, const size_t z, const double flux, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_USER_INPUT, "Failed to update LinkerPar object; list is currently empty.");
	
	// Get index
	const size_t index = self->size - 1;
	
	++self->n_pix[index];
	if(x < self->x_min[index]) self->x_min[index] = x;
	if(x > self->x_max[index]) self->x_max[index] = x;
	if(y < self->y_min[index]) self->y_min[index] = y;
	if(y > self->y_max[index]) self->y_max[index] = y;
	if(z < self->z_min[index]) self->z_min[index] = z;
	if(z > self->z_max[index]) self->z_max[index] = z;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr[index] += flux * x;
	self->y_ctr[index] += flux * y;
	self->z_ctr[index] += flux * z;
	#endif
	if(flux > self->f_max[index]) self->f_max[index] = flux;
	if(flux < self->f_min[index]) self->f_min[index] = flux;
	self->f_sum[index] += flux;
	self->flags[index] |= flag;
	
	return;
}

// Same, but only the flag will get updated

PUBLIC void LinkerPar_update_flag(LinkerPar *self, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_USER_INPUT, "Failed to update LinkerPar object; list is currently empty.");
	
	// Set flag
	self->flags[self->size - 1] |= flag;
	
	return;
}



// ----------------------------------------------------------------- //
// Get the size of an object in x, y or z                            //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Index of the object to be retrieved.             //
//   (3) axis     - Axis for which size should be returned; 0 = x,   //
//                  1 = y, and 2 = z.                                //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Size of the object in pixels along the specified axis.          //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the size of the specified object    //
//   along the specified axis. The programme will terminate if the   //
//   axis or label are out of range.                                 //
// ----------------------------------------------------------------- //

PUBLIC size_t LinkerPar_get_obj_size(const LinkerPar *self, const size_t label, const int axis)
{
	// Sanity checks
	check_null(self);
	ensure(axis >= 0 && axis <= 2, ERR_USER_INPUT, "Invalid axis selection (%d) in LinkerPar object.", axis);
	
	// Determine index
	const size_t index = LinkerPar_get_index(self, label);
	
	if(axis == 0) return self->x_max[index] - self->x_min[index] + 1;
	if(axis == 1) return self->y_max[index] - self->y_min[index] + 1;
	return self->z_max[index] - self->z_min[index] + 1;
}



// ----------------------------------------------------------------- //
// Get the number of pixels of an object                             //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Index of the object to be retrieved.             //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Number of pixels of the specified object.                       //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the number of pixels that have been //
//   recorded for the specified object. The programme will terminate //
//   if the label is out of range.                                   //
// ----------------------------------------------------------------- //

PUBLIC size_t LinkerPar_get_npix(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->n_pix[LinkerPar_get_index(self, label)];
}



// ----------------------------------------------------------------- //
// Get the total flux of an object                                   //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Label of the object to be retrieved.             //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Total flux of the specified object.                             //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the total flux of the object speci- //
//   fied by 'label'. The programme will terminate if the label is   //
//   out of range.                                                   //
// ----------------------------------------------------------------- //

PUBLIC double LinkerPar_get_flux(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->f_sum[LinkerPar_get_index(self, label)];
}



// ----------------------------------------------------------------- //
// Get the reliability of an object                                  //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Label of the object to be retrieved.             //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Reliability of the specified object.                            //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the reliability of the object spe-  //
//   cified by 'label'. The programme will terminate if the label is //
//   out of range.                                                   //
// ----------------------------------------------------------------- //

PUBLIC double LinkerPar_get_rel(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->rel[LinkerPar_get_index(self, label)];
}



// ----------------------------------------------------------------- //
// Get the label of an object by index                               //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) index    - Index of the object to be retrieved.             //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Reliability of the specified object.                            //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for returning the label of the object with the    //
//   specified index. The programme will terminate if the index is   //
//   out of range.                                                   //
// ----------------------------------------------------------------- //

PUBLIC size_t LinkerPar_get_label(const LinkerPar *self, const size_t index)
{
	// Sanity checks
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Index out of range. Cannot retrieve label.");
	
	return self->label[index];
}



// ----------------------------------------------------------------- //
// Get bounding box of an object                                     //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Label of the object for which the bounding box   //
//                  is to be retrieved.                              //
//   (3)-(8) x_min, x_max, y_min, y_max, z_min, z_max                //
//                - Pointer for holding the bounding box values.     //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for retrieving the bounding box of the object     //
//   with the specified label. The process will be terminated if the //
//   label does not exist.                                           //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_get_bbox(const LinkerPar *self, const size_t label, size_t *x_min, size_t *x_max, size_t *y_min, size_t *y_max, size_t *z_min, size_t *z_max)
{
	// Sanity checks
	check_null(self);
	
	// Determine index
	size_t index = LinkerPar_get_index(self, label);
	
	// Copy bounding box values
	*x_min = self->x_min[index];
	*x_max = self->x_max[index];
	*y_min = self->y_min[index];
	*y_max = self->y_max[index];
	*z_min = self->z_min[index];
	*z_max = self->z_max[index];
	
	return;
}



// ----------------------------------------------------------------- //
// Create source catalogue from LinkerPar object                     //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self      - Object self-reference.                          //
//   (2) filter    - Map object containing old and new label pairs   //
//                   of only those sources considered as reliable.   //
//   (3) flux_unit - String containing the flux unit of the data.    //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Pointer to newly created catalogue.                             //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for generating a source catalogue from the speci- //
//   fied LinkerPar object. A pointer to the newly created catalogue //
//   will be returned. Note that the user will assume ownership of   //
//   the catalogue and will be responsible for explicitly calling    //
//   the destructor if the catalogue is no longer required.          //
//   Unreliable sources can be excluded from the catalogue by provi- //
//   ding a non-empty filter object that contains the old and new    //
//   labels of all sources deemed reliable. The old labels will be   //
//   replaced by the new ones in the source ID column of the cata-   //
//   logue. If filter is NULL or empty, all sources will be copied   //
//   into the catalogue without filtering.                           //
// ----------------------------------------------------------------- //

PUBLIC Catalog *LinkerPar_make_catalog(const LinkerPar *self, const Map *filter, const char *flux_unit)
{
	// Sanity checks
	check_null(self);
	check_null(flux_unit);
	
	// Check if reliability filtering requested
	const bool remove_unreliable = filter != NULL && Map_get_size(filter);
	
	// Create an empty source catalogue
	Catalog *cat = Catalog_new();
	
	// Create string for holding identifier
	String *identifier = String_new("");
	
	// Loop over all LinkerPar entries
	for(size_t i = 0; i < self->size; ++i)
	{
		const size_t new_label = remove_unreliable && Map_key_exists(filter, self->label[i]) ? Map_get_value(filter, self->label[i]) : self->label[i];
		
		if(!remove_unreliable || Map_key_exists(filter, self->label[i]))
		{
			// Create a new source
			Source *src = Source_new(self->verbosity);
			
			// Set the identifier to the current label
			String_set_int(identifier, "%zu", self->label[i]);
			Source_set_identifier(src, String_get(identifier));
			
			// Set other parameters
			Source_add_par_int(src, "id",    new_label,                       "",        "meta.id");
			#if MEASURE_CENTROID_POSITION
			Source_add_par_flt(src, "x",     self->x_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.x");
			Source_add_par_flt(src, "y",     self->y_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.y");
			Source_add_par_flt(src, "z",     self->z_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.z");
			#else
			// Insert placeholder if no centroid requested
			Source_add_par_flt(src, "x",     0.0,                             "pix",     "pos.cartesian.x");
			Source_add_par_flt(src, "y",     0.0,                             "pix",     "pos.cartesian.y");
			Source_add_par_flt(src, "z",     0.0,                             "pix",     "pos.cartesian.z");
			#endif
			Source_add_par_int(src, "x_min", self->x_min[i],                  "pix",     "pos.cartesian.x;stat.min");
			Source_add_par_int(src, "x_max", self->x_max[i],                  "pix",     "pos.cartesian.x;stat.max");
			Source_add_par_int(src, "y_min", self->y_min[i],                  "pix",     "pos.cartesian.y;stat.min");
			Source_add_par_int(src, "y_max", self->y_max[i],                  "pix",     "pos.cartesian.y;stat.max");
			Source_add_par_int(src, "z_min", self->z_min[i],                  "pix",     "pos.cartesian.z;stat.min");
			Source_add_par_int(src, "z_max", self->z_max[i],                  "pix",     "pos.cartesian.z;stat.max");
			Source_add_par_int(src, "n_pix", self->n_pix[i],                  "",        "meta.number;instr.pixel");
			Source_add_par_flt(src, "f_min", self->f_min[i],                  flux_unit, "phot.flux.density;stat.min");
			Source_add_par_flt(src, "f_max", self->f_max[i],                  flux_unit, "phot.flux.density;stat.max");
			Source_add_par_flt(src, "f_sum", self->f_sum[i],                  flux_unit, "phot.flux");
			Source_add_par_flt(src, "rel",   self->rel  [i],                  "",        "stat.probability");
			Source_add_par_int(src, "flag",  self->flags[i],                  "",        "meta.code.qual");
			
			// Add source to catalogue
			Catalog_add_source(cat, src);
		}
	}
	
	// Clean up
	String_delete(identifier);
	
	// Return catalogue
	return cat;
}



// ----------------------------------------------------------------- //
// Print some basic information about the LinkerPar object           //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for printing some basic information on the size   //
//   and memory usage of the LinkerPar object pointed to by 'self'.  //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_print_info(const LinkerPar *self)
{
	// Sanity checks
	check_null(self);
	
	// Calculate memory usage
	#if MEASURE_CENTROID_POSITION
	const double memory_usage = (double)(self->size * (8 * sizeof(size_t) + 7 * sizeof(double) + 1 * sizeof(char)));
	#else
	const double memory_usage = (double)(self->size * (8 * sizeof(size_t) + 4 * sizeof(double) + 1 * sizeof(char)));
	#endif
	
	// Print size and memory information
	message("Linker status:");
	message(" - No. of objects:  %zu", self->size);
	if(memory_usage < MEGABYTE) message(" - Memory usage:    %.2f kB\n", memory_usage / KILOBYTE);
	else message(" - Memory usage:    %.2f MB\n", memory_usage / MEGABYTE);
	
	return;
}



// ----------------------------------------------------------------- //
// Return the index of the element identified by label               //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) label    - Label of the element the index of which is to be //
//                  returned.                                        //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Index of the element identified by 'label'.                     //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Private method for finding and returning the index of the ele-  //
//   ment with the specified label. The process will be terminated   //
//   if the requested label does not exist.                          //
// ----------------------------------------------------------------- //

PRIVATE size_t LinkerPar_get_index(const LinkerPar *self, const size_t label)
{
	size_t index = 0;
	while(index < self->size && self->label[index] != label) ++index;
	ensure(self->size && self->label[index] == label, ERR_USER_INPUT, "Label not found.");
	return index;
}



// ----------------------------------------------------------------- //
// Reallocate memory for LinkerPar object                            //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Private method for reallocating the memory requirements of the  //
//   specified LinkerPar object, e.g. as necessitated by a change in //
//   size. If the new size is 0, all memory will be de-allocated and //
//   the pointers will be set to NULL.                               //
// ----------------------------------------------------------------- //

PRIVATE void LinkerPar_reallocate_memory(LinkerPar *self)
{
	if(self->size)
	{
		// Reallocate memory
		self->label = (size_t *)memory_realloc(self->label, self->size, sizeof(size_t));
		self->n_pix = (size_t *)memory_realloc(self->n_pix, self->size, sizeof(size_t));
		self->x_min = (size_t *)memory_realloc(self->x_min, self->size, sizeof(size_t));
		self->x_max = (size_t *)memory_realloc(self->x_max, self->size, sizeof(size_t));
		self->y_min = (size_t *)memory_realloc(self->y_min, self->size, sizeof(size_t));
		self->y_max = (size_t *)memory_realloc(self->y_max, self->size, sizeof(size_t));
		self->z_min = (size_t *)memory_realloc(self->z_min, self->size, sizeof(size_t));
		self->z_max = (size_t *)memory_realloc(self->z_max, self->size, sizeof(size_t));
		#if MEASURE_CENTROID_POSITION
		self->x_ctr = (double *)memory_realloc(self->x_ctr, self->size, sizeof(double));
		self->y_ctr = (double *)memory_realloc(self->y_ctr, self->size, sizeof(double));
		self->z_ctr = (double *)memory_realloc(self->z_ctr, self->size, sizeof(double));
		#endif
		self->f_min = (double *)memory_realloc(self->f_min, self->size, sizeof(double));
		self->f_max = (double *)memory_realloc(self->f_max, self->size, sizeof(double));
		self->f_sum = (double *)memory_realloc(self->f_sum, self->size, sizeof(double));
		self->rel   = (double *)memory_realloc(self->rel,   self->size, sizeof(double));
		self->flags = (unsigned char *)memory_realloc(self->flags, self->size, sizeof(unsigned char));
	}
	else
	{
		free(self->label);
		free(self->n_pix);
		free(self->x_min);
		free(self->x_max);
		free(self->y_min);
		free(self->y_max);
		free(self->z_min);
		free(self->z_max);
		#if MEASURE_CENTROID_POSITION
		free(self->x_ctr);
		free(self->y_ctr);
		free(self->z_ctr);
		#endif
		free(self->f_min);
		free(self->f_max);
		free(self->f_sum);
		free(self->rel);
		free(self->flags);
		
		self->label = NULL;
		self->n_pix = NULL;
		self->x_min = NULL;
		self->x_max = NULL;
		self->y_min = NULL;
		self->y_max = NULL;
		self->z_min = NULL;
		self->z_max = NULL;
		#if MEASURE_CENTROID_POSITION
		self->x_ctr = NULL;
		self->y_ctr = NULL;
		self->z_ctr = NULL;
		#endif
		self->f_min = NULL;
		self->f_max = NULL;
		self->f_sum = NULL;
		self->rel   = NULL;
		self->flags = NULL;
	}
	
	return;
}



// ----------------------------------------------------------------- //
// Determine reliability of detections                               //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self         - Object self-reference.                       //
//   (2) scale_kernel - The size of the convolution kernel used in   //
//                      determining the density of positive and ne-  //
//                      gative detections in parameter space will be //
//                      scaled by this factor. If set to 1, the ori- //
//                      ginal covariance matrix derived from the     //
//                      distribution of negative sources is used.    //
//   (3) fmin         - Value of the fmin parameter, where fmin =    //
//                      sum / sqrt(N).                               //
//   (4) rel_cat      - Table of pixel coordinates on the sky. All   //
//                      negative detections with bounding boxes in-  //
//                      cluding those positions will be removed be-  //
//                      fore reliability calculation. NULL can be    //
//                      used to disable this feature.                //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Covariance matrix from the negative detections.                 //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for measuring the reliability of all the sources  //
//   in the specified LinkerPar object. This will set the rel pro-   //
//   perty of the object, but not yet filter out unreliable sources. //
//   Reliability measurement works by comparing the density of posi- //
//   tive and negative detections in an N-dimensional parameter      //
//   space. For this purpose, the covariance matrix of the distribu- //
//   tion of negative sources in parameter space is first calcula-   //
//   ted. The covariance matrix is assumed to describe the multi-    //
//   variate normal distribution of the Gaussian noise of the data.  //
//   Next, the sum of the probability density functions of all posi- //
//   tive and negative sources is evaluated at the location of each  //
//   positive detection (multivariate Gaussian kernel density esti-  //
//   mation). From this, the reliability is estimated as             //
//                                                                   //
//      R = (P - N) / N,                                             //
//                                                                   //
//   where P is the sum of the PDFs of the positive sources, and N   //
//   is the sum of the PDFs of the negative sources. If N > P, R is  //
//   set to 0 to ensure that the resulting reliability is always in  //
//   the range of 0 to 1. Note that the reliability will only be de- //
//   termined for positive sources above the fmin threshold, where   //
//   fmin is the summed flux divided by the square root of the num-  //
//   ber of pixels contributing to a source.                         //
//   In order to be able to exclude certain negative artefacts from  //
//   affecting the reliability calculation, the user has the option  //
//   of specifying a table of (x, y) pixel positions using the para- //
//   meter rel_cat. All negative detections the (x, y) bounding box  //
//   of which contains one of those positions will be excluded from  //
//   the reliability calculation. rel_cat must contain exactly two   //
//   columns (x and y in pixels). If set to NULL, this feature will  //
//   be disabled altogether.                                         //
// ----------------------------------------------------------------- //

PUBLIC Matrix *LinkerPar_reliability(LinkerPar *self, const double scale_kernel, const double fmin, const Table *rel_cat)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_NO_SRC_FOUND, "No sources left after linking. Cannot proceed.");
	
	// Dimensionality of parameter space
	const int dim = 3;
	
	// Define a few parameters
	size_t n_neg = 0;
	size_t n_pos = 0;
	size_t counter_neg = 0;
	size_t counter_pos = 0;
	const size_t threshold_warning = 50;
	const double log_fmin_squared = 2.0 * log10(fmin);
	
	// Determine number of positive and negative detections
	for(size_t i = self->size; i--;)
	{
		if(self->f_sum[i] < 0.0) ++n_neg;
		else if(self->f_sum[i] > 0.0) ++n_pos;
	}
	
	ensure(n_neg, ERR_FAILURE, "No negative sources found. Cannot proceed.");
	ensure(n_pos, ERR_FAILURE, "No positive sources found. Cannot proceed.");
	message("Found %zu positive and %zu negative sources.", n_pos, n_neg);
	if(n_neg < threshold_warning) warning("Only %zu negative detections found.\n         Reliability calculation may not be accurate.", n_neg);
	
	// Extract relevant parameters
	double *par_pos = (double *)memory(MALLOC, dim * n_pos, sizeof(double));
	size_t *idx_pos = (size_t *)memory(MALLOC, n_pos, sizeof(size_t));
	double *par_neg = (double *)memory(MALLOC, dim * n_neg, sizeof(double));
	size_t *idx_neg = (size_t *)memory(MALLOC, n_neg, sizeof(size_t));
	
	for(size_t i = self->size; i--;)
	{
		if(self->f_sum[i] < 0.0)
		{
			bool include_source = true;
			
			if(rel_cat != NULL)
			{
				for(size_t row = 0; row < Table_rows(rel_cat); ++row)
				{
					const double cat_x = Table_get(rel_cat, row, 0);
					const double cat_y = Table_get(rel_cat, row, 1);
					
					if(cat_x >= (double)(self->x_min[i]) && cat_x <= (double)(self->x_max[i]) && cat_y >= (double)(self->y_min[i]) && cat_y <= (double)(self->y_max[i]))
					{
						include_source = false;
						break;
					}
				}
			}
			
			if(include_source)
			{
				ensure(self->f_min[i] < 0.0, ERR_FAILURE, "Non-negative minimum assigned to source with negative flux!");
				par_neg[dim * counter_neg + 0] = log10(-self->f_min[i]);
				par_neg[dim * counter_neg + 1] = log10(-self->f_sum[i]);
				par_neg[dim * counter_neg + 2] = log10(-self->f_sum[i] / self->n_pix[i]);
				idx_neg[counter_neg] = i;
				++counter_neg;
			}
		}
		else if(self->f_sum[i] > 0.0)
		{
			ensure(self->f_max[i] > 0.0, ERR_FAILURE, "Non-positive maximum assigned to source with positive flux!");
			par_pos[dim * counter_pos + 0] = log10(self->f_max[i]);
			par_pos[dim * counter_pos + 1] = log10(self->f_sum[i]);
			par_pos[dim * counter_pos + 2] = log10(self->f_sum[i] / self->n_pix[i]);
			idx_pos[counter_pos] = i;
			++counter_pos;
		}
	}
	
	// Adjust array sizes if necessary (some negative sources may have been removed)
	if(counter_neg < n_neg)
	{
		message("Excluding %zu out of %zu negative sources from reliability analysis.", n_neg - counter_neg, n_neg);
		n_neg = counter_neg;
		ensure(n_neg, ERR_FAILURE, "No negative sources found. Cannot proceed.");
		if(n_neg < threshold_warning) warning("Only %zu negative detections found.\n         Reliability calculation may not be accurate.", n_neg);
		par_neg = (double *)memory_realloc(par_neg, dim * n_neg, sizeof(double));
		idx_neg = (size_t *)memory_realloc(idx_neg, n_neg, sizeof(size_t));
	}
	else message("Retaining all negative detections.");
	
	// Determine covariance matrix from negative detections
	Matrix *covar = Matrix_new(dim, dim);
	double mean[dim];
	
	// Calculate mean values first
	for(size_t i = dim; i--;)
	{
		mean[i] = 0.0;
		for(size_t j = 0; j < n_neg; ++j) mean[i] += par_neg[dim * j + i];
		mean[i] /= n_neg;
	}
	
	// Then calculate the covariance matrix
	for(size_t i = dim; i--;)
	{
		for(size_t j = dim; j--;)
		{
			for(size_t k = 0; k < dim * n_neg; k += dim)
			{
				Matrix_add_value(covar, i, j, (par_neg[k + i] - mean[i]) * (par_neg[k + j] - mean[j]));
			}
			Matrix_mul_value(covar, i, j, scale_kernel * scale_kernel / n_neg);  // NOTE: Variance = sigma^2, hence scale_kernel^2 here.
		}
	}
	
	// Invert covariance matrix + sanity check
	Matrix *covar_inv = Matrix_invert(covar);
	ensure(covar_inv != NULL, ERR_FAILURE, "Covariance matrix is not invertible; cannot measure reliability.\n       Ensure that there are enough negative detections.");
	
	// Inverse of the square root of |2 * pi * covar| = (2 pi)^n |covar|
	// This is the scale factor needed to calculate the PDF of the multivariate normal distribution later on.
	const double scal_fact = 1.0 / sqrt(Matrix_det(covar, 2.0 * M_PI));
	// const double scal_fact = 1.0;
	// NOTE: This can be set to 1, as we don’t really care about the correct
	//       normalisation of the Gaussian kernel, so we might as well normalise
	//       the amplitude rather than the integral. The normalisation factor 
	//       does matter for the Skellam plot generation further down, though.
	
	
	// Create Skellam array
	/*Array_dbl *skellam = Array_dbl_new(n_neg);
	
	// Loop over all negative sources to derive Skellam distribution
	#pragma omp parallel
	{
		Matrix *vector = Matrix_new(dim, 1);
		
		#pragma omp for schedule(static)
		for(size_t i = 0; i < n_neg; ++i)
		{
			double p1 = par_neg[dim * i];
			double p2 = par_neg[dim * i + 1];
			double p3 = par_neg[dim * i + 2];
			
			// Multivariate kernel density estimation for negative detections
			double pdf_neg_sum = 0.0;
			
			for(double *ptr = par_neg + n_neg * dim; ptr > par_neg;)
			{
				// Set up relative position vector
				Matrix_set_value_nocheck(vector, 2, 0, *(--ptr) - p3);
				Matrix_set_value_nocheck(vector, 1, 0, *(--ptr) - p2);
				Matrix_set_value_nocheck(vector, 0, 0, *(--ptr) - p1);
				
				pdf_neg_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
			}
			
			// Multivariate kernel density estimation for positive detections
			double pdf_pos_sum = 0.0;
			
			for(double *ptr = par_pos + n_pos * dim; ptr > par_pos;)
			{
				// Set up relative position vector
				Matrix_set_value_nocheck(vector, 2, 0, *(--ptr) - p3);
				Matrix_set_value_nocheck(vector, 1, 0, *(--ptr) - p2);
				Matrix_set_value_nocheck(vector, 0, 0, *(--ptr) - p1);
				
				pdf_pos_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
			}
			
			// Determine Skellam parameter
			Array_dbl_set(skellam, i, (pdf_pos_sum - pdf_neg_sum) / sqrt(scal_fact * pdf_pos_sum + pdf_neg_sum));
		}
		
		Matrix_delete(vector);
	}
	
	// Create Skellam plot
	LinkerPar_skellam_plot(skellam, "skellam_plot_test.eps", true);
	
	// Delete Skellam array again
	Array_dbl_delete(skellam);*/
	
	
	// Loop over all positive detections to measure their reliability
	const size_t cadence = (n_pos / 100) ? n_pos / 100 : 1;  // Only needed for progress bar
	size_t progress = 0;
	message("");
	
	#pragma omp parallel
	{
		Matrix *vector = Matrix_new(dim, 1);
		
		#pragma omp for schedule(static)
		for(size_t i = 0; i < n_pos; ++i)
		{
			#pragma omp critical
			if(++progress % cadence == 0 || progress == n_pos) progress_bar("Progress: ", progress, n_pos);
			
			const double p2 = par_pos[dim * i + 1];
			const double p3 = par_pos[dim * i + 2];
			
			// Only process sources above fmin
			if(p2 + p3 > log_fmin_squared)
			{
				const double p1 = par_pos[dim * i];
				
				// Multivariate kernel density estimation for negative detections
				double pdf_neg_sum = 0.0;
				
				for(double *ptr = par_neg + n_neg * dim; ptr > par_neg;)
				{
					// Set up relative position vector
					Matrix_set_value_nocheck(vector, 2, 0, *(--ptr) - p3);
					Matrix_set_value_nocheck(vector, 1, 0, *(--ptr) - p2);
					Matrix_set_value_nocheck(vector, 0, 0, *(--ptr) - p1);
					
					pdf_neg_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
				}
				
				// Multivariate kernel density estimation for positive detections
				double pdf_pos_sum = 0.0;
				
				for(double *ptr = par_pos + n_pos * dim; ptr > par_pos;)
				{
					// Set up relative position vector
					Matrix_set_value_nocheck(vector, 2, 0, *(--ptr) - p3);
					Matrix_set_value_nocheck(vector, 1, 0, *(--ptr) - p2);
					Matrix_set_value_nocheck(vector, 0, 0, *(--ptr) - p1);
					
					pdf_pos_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
				}
				
				// Determine reliability
				self->rel[idx_pos[i]] = pdf_pos_sum > pdf_neg_sum ? (pdf_pos_sum - pdf_neg_sum) / pdf_pos_sum : 0.0;
			}
		}
		
		Matrix_delete(vector);
	}
	
	// Release memory again
	Matrix_delete(covar_inv);
	free(par_pos);
	free(par_neg);
	free(idx_pos);
	free(idx_neg);
	
	return covar;
}



// ----------------------------------------------------------------- //
// Create reliability diagnostic plots                               //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self         - Object self-reference.                       //
//   (2) threshold    - Reliability threshold.                       //
//   (3) fmin         - Threshold for SNR filtering.                 //
//   (4) covar        - Covariance matrix.                           //
//   (5) filename     - Name of the output EPS file.                 //
//   (6) overwrite    - If true, overwrite output file, otherwise do //
//                      not overwrite.                               //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for creating a diagnostic plot for the reliabili- //
//   ty measurement. The method will generate an Encapsulated Post-  //
//   Script (EPS) file that shows the distribution of negative and   //
//   positive sources in 2-D projections of the parameter space and  //
//   highlight the ones that do or don't fulfil the reliability      //
//   threshold or fmin requirements. In addition, the location of    //
//   fmin will be plotted, and the Gaussian smoothing kernel used in //
//   the reliability measurement will be overplotted as an ellipse   //
//   based on the provided covariance matrix.                        //
//   Note that if overwrite is set the false and the output file al- //
//   ready exists, the process will be terminated.                   //
// ----------------------------------------------------------------- //

PUBLIC void LinkerPar_rel_plots(const LinkerPar *self, const double threshold, const double fmin, const Matrix *covar, const char *filename, const bool overwrite)
{
	// Sanity checks
	check_null(self);
	if(self->size == 0)
	{
		warning("No sources found; cannot generate reliability plots.");
		return;
	}
	ensure(filename != NULL && strlen(filename), ERR_USER_INPUT, "Empty file name for reliability plot provided.");
	
	// Some settings
	const size_t plot_size_x = 300;  // pt
	const size_t plot_size_y = 300;  // pt
	const size_t plot_offset_y = 50;  // pt
	
	const char *colour_neg = "1 0.4 0.4";
	const char *colour_pos = "0.4 0.4 1";
	const char *colour_rel = "0 0 0";
	const char *colour_kernel = "0.8 0.8 0.8";
	const char *colour_fmin = "0.5 0.5 0.5";
	const char *colour_axes = "0 0 0";
	
	const char *par_space_x[3] = {"log\\(peak / rms\\)", "log\\(peak / rms\\)", "log\\(sum / rms\\)"};
	const char *par_space_y[3] = {"log\\(sum / rms\\)", "log\\(mean / rms\\)", "log\\(mean / rms\\)"};
	
	// Create arrays for parameters
	double *data_x = (double *)memory(MALLOC, self->size, sizeof(double));
	double *data_y = (double *)memory(MALLOC, self->size, sizeof(double));
	
	// Open PS file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", filename);
	
	message("Creating postscript file: %s", strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/') + 1);
	
	// Print PS header
	write_eps_header(fp, "SoFiA Reliability Plots", SOFIA_VERSION_FULL, "0 10 1060 360");
	
	for(int n = 0; n < 3; ++n)
	{
		// Read values and determine plotting range
		size_t plot_offset_x = 50 + n * (plot_size_x + 50);
		
		double data_min_x =  999.9;
		double data_max_x = -999.9;
		double data_min_y =  999.9;
		double data_max_y = -999.9;
		
		double radius_maj, radius_min, pa;
		if(n == 0) Matrix_err_ellipse(covar, 0, 1, &radius_maj, &radius_min, &pa);
		else if(n == 1) Matrix_err_ellipse(covar, 0, 2, &radius_maj, &radius_min, &pa);
		else if(n == 2) Matrix_err_ellipse(covar, 1, 2, &radius_maj, &radius_min, &pa);
		
		for(size_t i = 0; i < self->size; ++i)
		{
			// Extract relevant parameters
			switch(n)
			{
				case 0:
					if(self->f_sum[i] < 0.0)
					{
						data_x[i] = log10(-self->f_min[i]);
						data_y[i] = log10(-self->f_sum[i]);
					}
					else
					{
						data_x[i] = log10(self->f_max[i]);
						data_y[i] = log10(self->f_sum[i]);
					}
					break;
				case 1:
					if(self->f_sum[i] < 0.0)
					{
						data_x[i] = log10(-self->f_min[i]);
						data_y[i] = log10(-self->f_sum[i] / self->n_pix[i]);
					}
					else
					{
						data_x[i] = log10(self->f_max[i]);
						data_y[i] = log10(self->f_sum[i] / self->n_pix[i]);
					}
					break;
				case 2:
					if(self->f_sum[i] < 0.0)
					{
						data_x[i] = log10(-self->f_sum[i]);
						data_y[i] = log10(-self->f_sum[i] / self->n_pix[i]);
					}
					else
					{
						data_x[i] = log10(self->f_sum[i]);
						data_y[i] = log10(self->f_sum[i] / self->n_pix[i]);
					}
					break;
			}
			
			if(data_min_x > data_x[i]) data_min_x = data_x[i];
			if(data_max_x < data_x[i]) data_max_x = data_x[i];
			if(data_min_y > data_y[i]) data_min_y = data_y[i];
			if(data_max_y < data_y[i]) data_max_y = data_y[i];
		}
		
		double data_range_x = data_max_x - data_min_x;
		double data_range_y = data_max_y - data_min_y;
		
		// Add a little bit of margin
		data_min_x -= 0.05 * data_range_x;
		data_max_x += 0.05 * data_range_x;
		data_min_y -= 0.05 * data_range_y;
		data_max_y += 0.05 * data_range_y;
		
		// Ensure that ranges are equal in x and y
		/*if(data_range_x > data_range_y)
		{
			data_min_y -= 0.5 * (data_range_x - data_range_y);
			data_max_y += 0.5 * (data_range_x - data_range_y);
			data_range_y = data_max_y - data_min_y;
		}
		else if(data_range_x < data_range_y)
		{
			data_min_x -= 0.5 * (data_range_y - data_range_x);
			data_max_x += 0.5 * (data_range_y - data_range_x);
			data_range_x = data_max_x - data_min_x;
		}*/
		
		// Determine optimal tick mark increments
		const double tick_inc_x = auto_tick(data_max_x - data_min_x, 4);
		const double tick_inc_y = auto_tick(data_max_y - data_min_y, 4);
		
		// Determine the mean of negative sources
		double mean_x = 0.0;
		double mean_y = 0.0;
		size_t counter = 0;
		
		for(size_t i = 0; i < self->size; ++i)
		{
			if(self->f_sum[i] < 0.0)
			{
				mean_x += data_x[i];
				mean_y += data_y[i];
				++counter;
			}
		}
		
		mean_x /= counter;
		mean_y /= counter;
		
		const double centre_x = (mean_x - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		const double centre_y = (mean_y - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		const double radius_x = radius_maj * plot_size_x / (data_max_x - data_min_x);
		const double radius_y = radius_min * plot_size_x / (data_max_x - data_min_x);
		
		// Determine scale factor of kernel ellipse
		const double scale_factor = data_range_x / data_range_y;
		
		// Plot negative sources
		fprintf(fp, "%s rgb\n", colour_neg);
		fprintf(fp, "0.5 lw\n");
		fprintf(fp, "np\n");
		
		for(size_t i = self->size; i--;)
		{
			if(self->f_sum[i] < 0.0)
			{
				const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
				const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				
				fprintf(fp, "%.1f %.1f 1 0 360 af\n", plot_x, plot_y);
			}
		}
		
		// Plot unreliable positive sources
		fprintf(fp, "%s rgb\n", colour_pos);
		
		for(size_t i = self->size; i--;)
		{
			if(self->f_sum[i] > 0.0 && self->rel[i] < threshold)
			{
				const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
				const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				
				fprintf(fp, "%.1f %.1f 1 0 360 af\n", plot_x, plot_y);
			}
		}
		
		// Plot reliable positive sources
		fprintf(fp, "%s rgb\n", colour_rel);
		
		for(size_t i = self->size; i--;)
		{
			if(self->f_sum[i] > 0.0 && self->rel[i] >= threshold)
			{
				const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
				const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				
				if(self->f_sum[i] / sqrt(self->n_pix[i]) > fmin) fprintf(fp, "%.1f %.1f 2 0 360 af\n", plot_x, plot_y);
				else fprintf(fp, "%.1f %.1f 2 0 360 as\n", plot_x, plot_y);
			}
		}
		
		// Plot kernel ellipse
		fprintf(fp, "gsave\n");
		fprintf(fp, "%s rgb\n", colour_kernel);
		fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
		fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, radius_x, radius_y, 180.0 * pa / M_PI, scale_factor);
		fprintf(fp, "[2 2] 0 setdash\n");
		fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, 2.0 * radius_x, 2.0 * radius_y, 180.0 * pa / M_PI, scale_factor);
		fprintf(fp, "[0.5 1.5] 0 setdash\n");
		fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, 3.0 * radius_x, 3.0 * radius_y, 180.0 * pa / M_PI, scale_factor);
		fprintf(fp, "grestore\n");
		
		// Plot fmin line if possible
		if(n == 2)
		{
			double plot_x = plot_offset_x;
			double plot_y = (2.0 * log10(fmin) - data_min_x - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
			
			fprintf(fp, "gsave\n");
			fprintf(fp, "%s rgb\n", colour_fmin);
			fprintf(fp, "[3 3] 0 setdash\n");
			fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
			fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
			
			plot_x = plot_offset_x + plot_size_x;
			plot_y = (2.0 * log10(fmin) - data_max_x - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
			
			fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
			fprintf(fp, "grestore\n");
		}
		
		// Plot frame
		fprintf(fp, "%s rgb\n", colour_axes);
		fprintf(fp, "[] 0 setdash\n");
		fprintf(fp, "np\n");
		fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
		fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
		fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
		fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
		fprintf(fp, "cp s\n");
		
		// Plot tick marks
		for(double tm = ceil(data_min_x / tick_inc_x) * tick_inc_x; tm <= data_max_x; tm += tick_inc_x)
		{
			if(fabs(tm) < 0.001) tm = 0.0;
			const double plot_x = (tm - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
			fprintf(fp, "np %.2f %zu m %.2f %zu l s\n", plot_x, plot_offset_y, plot_x, plot_offset_y + 5);
			fprintf(fp, "np %.2f %zu m (%.1f) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_x, plot_offset_y - 14, tm);
		}
		
		for(double tm = ceil(data_min_y / tick_inc_y) * tick_inc_y; tm <= data_max_y; tm += tick_inc_y)
		{
			if(fabs(tm) < 0.001) tm = 0.0;
			const double plot_y = (tm - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
			fprintf(fp, "np %zu %.2f m %zu %.2f l s\n", plot_offset_x, plot_y, plot_offset_x + 5, plot_y);
			fprintf(fp, "np %zu %.2f m (%.1f) dup stringwidth pop neg 0 rmoveto show\n", plot_offset_x - 4, plot_y - 4.0, tm);
		}
		
		// Print labels
		fprintf(fp, "np %zu 20 m (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_offset_x + plot_size_x / 2, par_space_x[n]);
		
		fprintf(fp, "np %zu %zu m gsave 90 rotate (%s) dup stringwidth pop 2 div neg 0 rmoveto show grestore\n", plot_offset_x - 34, plot_offset_y + plot_size_y / 2, par_space_y[n]);
	}
	
	// Print EPS footer
	write_eps_footer(fp);
	
	// Close output file
	fclose(fp);
	
	// Clean up
	free(data_x);
	free(data_y);
	
	return;
}



// ----------------------------------------------------------------- //
// Create Skellam diagnostic plot                                    //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) skellam      - Array of (P - N) / sqrt(P + N) values.       //
//   (2) filename     - Output file name for plot.                   //
//   (3) overwrite    - If true, overwrite existing file.            //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Private function for generating a Skellam diagnostic plot show- //
//   ing the cumulative distribution of values in the Array called   //
//   'skellam' which should contain values of (P - N) / sqrt(P + N)  //
//   generated by the reliability module. The resulting plot will be //
//   written to an EPS file with the specified file name.            //
// ----------------------------------------------------------------- //

PRIVATE void LinkerPar_skellam_plot(Array_dbl *skellam, const char *filename, const bool overwrite)
{
	// Sanity checks
	check_null(skellam);
	const size_t size = Array_dbl_get_size(skellam);
	ensure(size, ERR_USER_INPUT, "Failed to create Skellam plot; no valid data found.");
	
	// Open output file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", filename);
	
	message("Creating postscript file: %s", strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/') + 1);
	
	// Print PS header
	write_eps_header(fp, "SoFiA Skellam Plot", SOFIA_VERSION_FULL, "0 0 480 360");
	
	// Sort the Skellam array
	Array_dbl_sort(skellam);
	//double data_min_x = Array_dbl_get(skellam, 0);
	//double data_max_x = Array_dbl_get(skellam, size - 1);
	//if(data_min_x < -10.0) data_min_x = -10.0;
	//else if(data_min_x > -4.0) data_min_x = -4.0;
	//if(data_max_x > 10.0) data_max_x = 10.0;
	//else if(data_max_x < 4.0) data_max_x = 4.0;
	//if(data_max_x < -data_min_x) data_max_x = -data_min_x;
	//else data_min_x = -data_max_x;
	const double data_min_x = -4.0;
	const double data_max_x = 4.0;
	const double data_min_y = 0.0;
	const double data_max_y = 1.0;
	
	// Plot geometry
	const size_t plot_size_x = 410;
	const size_t plot_size_y = 310;
	const size_t plot_offset_x = 50;
	const size_t plot_offset_y = 40;
	
	// Determine optimal tick mark increments
	const double tick_inc_x = auto_tick(data_max_x - data_min_x, 5);
	const double tick_inc_y = auto_tick(data_max_y - data_min_y, 5);
	
	// Colours
	const char *colour_data = "1 0 0";
	const char *colour_erf  = "0.3 0.3 0.3";
	const char *colour_zero = "0.7 0.7 0.7";
	const char *colour_axes = "0 0 0";
	
	// Labels
	const char *label_x = "\\(P - N\\) / sqrt\\(P + N\\)";
	const char *label_y = "Cumulative fraction";
	
	// Set clip path
	fprintf(fp, "gsave\n");
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
	fprintf(fp, "cp clip\n");
	
	// Plot vertical line at 0
	fprintf(fp, "np\n");
	fprintf(fp, "%.1f %.1f m\n", -data_min_x * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, -data_min_y * plot_size_y / (data_max_y - data_min_y) + plot_offset_y);
	fprintf(fp, "%.1f %.1f l\n", -data_min_x * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, (1.0 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y);
	fprintf(fp, "%s rgb\n", colour_zero);
	fprintf(fp, "[5 3] 0 setdash\n");
	fprintf(fp, "s\n");
	
	// Plot error function for Gaussian with sigma = 1
	fprintf(fp, "np\n");
	for(int i = -100; i <= 100; ++i)
	{
		const double x = 0.05 * i;
		const double plot_x = (x - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		const double plot_y = (0.5 - 0.5 * erf(-x / sqrt(2.0)) - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		
		fprintf(fp, "%.1f %.1f %s\n", plot_x, plot_y, i == -100 ? "m" : "l");
	}
	fprintf(fp, "%s rgb\n", colour_erf);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	
	// Plot data
	fprintf(fp, "np\n");
	for(size_t i = 0; i < size; ++i)
	{
		const double plot_x = (Array_dbl_get(skellam, i) - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		const double plot_y = (((double)(i) / (double)(size - 1)) - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		
		fprintf(fp, "%.1f %.1f %s\n", plot_x, plot_y, i ? "l" : "m");
	}
	fprintf(fp, "%s rgb\n", colour_data);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "grestore\n");
	
	// Plot frame
	fprintf(fp, "%s rgb\n", colour_axes);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
	fprintf(fp, "cp s\n");
	
	// Plot tick marks
	for(double tm = ceil(data_min_x / tick_inc_x) * tick_inc_x; tm <= data_max_x; tm += tick_inc_x)
	{
		if(fabs(tm) < 0.001) tm = 0.0;
		const double plot_x = (tm - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		fprintf(fp, "np %.2f %zu m %.2f %zu l s\n", plot_x, plot_offset_y, plot_x, plot_offset_y + 5);
		fprintf(fp, "np %.2f %zu m (%.1f) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_x, plot_offset_y - 14, tm);
	}
	
	for(double tm = ceil(data_min_y / tick_inc_y) * tick_inc_y; tm <= data_max_y; tm += tick_inc_y)
	{
		if(fabs(tm) < 0.001) tm = 0.0;
		const double plot_y = (tm - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		fprintf(fp, "np %zu %.2f m %zu %.2f l s\n", plot_offset_x, plot_y, plot_offset_x + 5, plot_y);
		fprintf(fp, "np %zu %.2f m (%.1f) dup stringwidth pop neg 0 rmoveto show\n", plot_offset_x - 4, plot_y - 4.0, tm);
	}
	
	// Print labels
	fprintf(fp, "np %zu 10 m (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_offset_x + plot_size_x / 2, label_x);
	fprintf(fp, "np %zu %zu m gsave 90 rotate (%s) dup stringwidth pop 2 div neg 0 rmoveto show grestore\n", plot_offset_x - 34, plot_offset_y + plot_size_y / 2, label_y);
	
	// Plot legend
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 20, plot_offset_y + plot_size_y - 20);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + 40, plot_offset_y + plot_size_y - 20);
	fprintf(fp, "%s rgb\n", colour_data);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 46, plot_offset_y + plot_size_y - 24);
	fprintf(fp, "(Data) show\n");
	
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 20, plot_offset_y + plot_size_y - 35);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + 40, plot_offset_y + plot_size_y - 35);
	fprintf(fp, "%s rgb\n", colour_erf);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 46, plot_offset_y + plot_size_y - 39);
	fprintf(fp, "(Gaussian \\(sigma = 1\\)) show\n");
	
	// Print EPS footer
	write_eps_footer(fp);
	
	// Close output file
	fclose(fp);
	
	return;
}
