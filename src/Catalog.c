/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (Catalog.c) - Source Finding Application                 ///
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
#include <string.h>
#include <math.h>
#include <time.h>

#include "Catalog.h"



// ----------------------------------------------------------------- //
// Declaration of properties of class Catalog                        //
// ----------------------------------------------------------------- //

CLASS Catalog
{
	size_t size;
	Source **sources;
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
//   Pointer to newly created Catalog object.                        //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Standard constructor. Will create a new and empty Catalog ob-   //
//   ject and return a pointer to the newly created object. No memo- //
//   ry will be allocated other than for the object itself. Note     //
//   that the destructor will need to be called explicitly once the  //
//   object is no longer required to release any memory allocated    //
//   during the lifetime of the object.                              //
// ----------------------------------------------------------------- //

PUBLIC Catalog *Catalog_new(void)
{
	// Allocate memory for new catalog
	Catalog *self = (Catalog *)memory(MALLOC, 1, sizeof(Catalog));
	
	// Initialise properties
	self->size = 0;
	self->sources = NULL;
	
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
//   NOTE that the destructor will explicitly call the destructor on //
//   all source objects stored in the catalogue. Hence, deleting a   //
//   catalogue will automatically delete all sources associated with //
//   that catalogue!                                                 //
// ----------------------------------------------------------------- //

PUBLIC void Catalog_delete(Catalog *self)
{
	if(self != NULL)
	{
		if(self->sources != NULL)
		{
			// Call the destructor on individual sources first
			Source **ptr = self->sources + self->size;
			while(ptr --> self->sources) Source_delete(*ptr);
			
			// Then de-allocate memory for pointers to those sources
			free(self->sources);
		}
		
		// Lastly, de-allocate memory for catalog object
		free(self);
	}
	
	return;
}



// ----------------------------------------------------------------- //
// Add a new source to a catalogue                                   //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) src      - Pointer to the source to be added                //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for adding a new source to the specified cata-    //
//   logue. Note that the function does not check if a source with   //
//   the same name already exists; a new source will always be added //
//   to the existing source list.                                    //
// ----------------------------------------------------------------- //

PUBLIC void Catalog_add_source(Catalog *self, Source *src)
{
	// Sanity checks
	check_null(self);
	check_null(src);
	ensure(!Catalog_source_exists(self, src, NULL), ERR_USER_INPUT, "Source \'%s\' is already in catalogue.", Source_get_identifier(src));
	
	Catalog_append_memory(self);
	*(self->sources + self->size - 1) = src;
	
	return;
}



// ----------------------------------------------------------------- //
// Get source index                                                  //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) src      - Pointer to the source to be checked              //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Returns the index (i.e. row number) of the source within the    //
//   catalogue if the source was found. Otherwise, SIZE_MAX will be  //
//   returned.                                                       //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for checking if the specified source is included  //
//   in the catalogue. If so, the function will return the row num-  //
//   ber of the source in the catalogue (starting with 0). If the    //
//   source is not found, the function will return SIZE_MAX.         //
// ----------------------------------------------------------------- //

PUBLIC size_t Catalog_get_index(const Catalog *self, const Source *src)
{
	// Sanity checks
	check_null(self);
	check_null(src);
	
	for(size_t i = 0; i < self->size; ++i)
	{
		if(self->sources[i] == src) return i;
	}
	
	return SIZE_MAX;
}



// ----------------------------------------------------------------- //
// Check if source exists in catalogue                               //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self     - Object self-reference.                           //
//   (2) src      - Pointer to the source to be checked              //
//   (3) index    - Pointer to index variable that will be set to    //
//                  the catalogue index of the source.               //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Returns true if the source is included in the catalogue and     //
//   false otherwise.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for checking if the specified source is included  //
//   in the catalogue. If so, the function will return true, other-  //
//   wise false. If the source is found, the variable 'index' will   //
//   be set to the catalogue index of the source. Otherwise, it will //
//   be left untouched. If no index is required, a NULL pointer can  //
//   instead be provided.
// ----------------------------------------------------------------- //

PUBLIC bool Catalog_source_exists(const Catalog *self, const Source *src, size_t *index)
{
	// Sanity checks
	check_null(self);
	check_null(src);
	
	for(size_t i = 0; i < self->size; ++i)
	{
		if(self->sources[i] == src)
		{
			if(index != NULL) *index = i;
			return true;
		}
	}
	
	return false;
}



// ----------------------------------------------------------------- //
// Retrieve a source from the catalogue by index                     //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self       - Object self-reference.                         //
//   (2) index      - Index of requested source.                     //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Returns a pointer to the requested source.                      //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for extracting a specific source from the cata-   //
//   logue by its index. A pointer to the source will be returned.   //
//   NOTE that the returned pointer must not be freed or deleted, as //
//   it is still owned by the Catalog object.                        //
// ----------------------------------------------------------------- //

PUBLIC Source *Catalog_get_source(const Catalog *self, const size_t index)
{
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Catalogue index out of range.");
	
	return self->sources[index];
}




// ----------------------------------------------------------------- //
// Return size of catalogue                                          //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self - Object self-reference.                               //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Returns the current size of the catalogue pointed to by 'self'. //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Returns the size of the catalogue, i.e. the number of sources   //
//   it contains. For empty catalogues a value of zero will be re-   //
//   turned.                                                         //
// ----------------------------------------------------------------- //

PUBLIC size_t Catalog_get_size(const Catalog *self)
{
	return self->size;
}



// ----------------------------------------------------------------- //
// Save catalogue to file                                            //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self      - Object self-reference.                          //
//   (2) filename  - Full path to the output file.                   //
//   (3) format    - Output format; can be CATALOG_FORMAT_ASCII for  //
//                   plain text ASCII files, CATALOG_FORMAT_XML for  //
//                   VOTable format or CATALOG_FORMAT_SQL for SQL    //
//                   table format.                                   //
//   (4) overwrite - Overwrite existing file (true) or not (false)?  //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for saving the current catalogue under the speci- //
//   fied name in the specified file format. The file name will be   //
//   relative to the process execution directory unless the full     //
//   path to the output directory is specified. Available formats    //
//   are plain text ASCII, VOTable XML format and SQL format.        //
// ----------------------------------------------------------------- //

PUBLIC void Catalog_save(const Catalog *self, const char *filename, const file_format format, const bool overwrite, const OUTTYPE OUTPUTS)
{
	// Sanity checks
	if(OUTPUTS == NONE) return;  // No output required

	check_null(self);
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "File name is empty.");
	
	if(!self->size)
	{
		warning("Failed to save catalogue; no sources found.");
		return;
	}
	
	// Open output file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", filename);
	
	// Some initial definitions
	const char char_comment = '#';
	const char char_nocomment = ' ';
	
	// Get current date and time
	char current_time_string[64];
	time_t current_time = time(NULL);
	strftime(current_time_string, 64, "%a, %d %b %Y, %H:%M:%S", localtime(&current_time));
	
	// Get first source to extract parameter names and units
	Source *src0 = self->sources[0];
	
	if(format == CATALOG_FORMAT_XML)
	{
		const char *data_type_names[2] = {"long", "double"};
		const char *indentation[7] = {"", "\t", "\t\t", "\t\t\t", "\t\t\t\t", "\t\t\t\t\t", "\t\t\t\t\t\t"}; // Better readability
		//const char *indentation[7] = {"", "", "", "", "", "", ""}; // Smaller file size
		
		// Write XML catalogue (VOTable)
		fprintf(fp, "%s<?xml version=\"1.0\" ?>\n", indentation[0]);
		fprintf(fp, "%s<VOTABLE version=\"1.3\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.ivoa.net/xml/VOTable/v1.3\">\n", indentation[0]);
		fprintf(fp, "%s<RESOURCE>\n", indentation[1]);
		fprintf(fp, "%s<DESCRIPTION>Source catalogue created by the Source Finding Application (SoFiA %s)</DESCRIPTION>\n", indentation[2], SOFIA_VERSION);
		fprintf(fp, "%s<PARAM name=\"Creator\" datatype=\"char\" arraysize=\"*\" value=\"SoFiA %s\" ucd=\"meta.id;meta.software\"/>\n", indentation[2], SOFIA_VERSION);
		fprintf(fp, "%s<PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\" value=\"%s\" ucd=\"time.creation\"/>\n", indentation[2], current_time_string);
		//fprintf(fp, "%s<COOSYS ID=\"wcs\" system=\"ICRS\" equinox=\"J2000\"/>\n", indentation[2]);
		// WARNING: COOSYS needs to be sorted out; see http://www.ivoa.net/documents/VOTable/ for documentation
		fprintf(fp, "%s<TABLE ID=\"SoFiA_source_catalogue\" name=\"SoFiA source catalogue\">\n", indentation[2]);
		
		// Column descriptors
		fprintf(fp, "%s<FIELD arraysize=\"32\" datatype=\"char\" name=\"name\" unit=\"\" ucd=\"meta.id\"/>\n", indentation[3]);
		for(size_t j = 0; j < Source_get_num_par(src0); ++j)
		{
			fprintf(fp, "%s<FIELD datatype=\"%s\" name=\"%s\" unit=\"%s\" ucd=\"%s\"/>\n", indentation[3], data_type_names[Source_get_type(src0, j)], Source_get_name(src0, j), Source_get_unit(src0, j), Source_get_ucd(src0, j));
		}
		
		// Start of data table
		fprintf(fp, "%s<DATA>\n", indentation[3]);
		fprintf(fp, "%s<TABLEDATA>\n", indentation[4]);
		
		// Data rows
		for(size_t i = 0; i < self->size; ++i)
		{
			Source *src = self->sources[i];
			fprintf(fp, "%s<TR>\n", indentation[5]);
			
			fprintf(fp, "%s<TD>%s</TD>\n", indentation[6], Source_get_identifier(src));
			
			for(size_t j = 0; j < Source_get_num_par(src); ++j)
			{
				if(Source_get_type(src, j) == SOURCE_TYPE_INT)
				{
					// Integer value
					const long int value = Source_get_par_int(src, j);
					fprintf(fp, "%s<TD>%ld</TD>\n", indentation[6], value);
				}
				else
				{
					// Floating-point value
					const double value = Source_get_par_flt(src, j);
					fprintf(fp, "%s<TD>%.15e</TD>\n", indentation[6], value);
				}
			}
			
			fprintf(fp, "%s</TR>\n", indentation[5]);
		}
		
		// End of data table
		fprintf(fp, "%s</TABLEDATA>\n", indentation[4]);
		fprintf(fp, "%s</DATA>\n", indentation[3]);
		
		// Finalise XML file
		fprintf(fp, "%s</TABLE>\n", indentation[2]);
		fprintf(fp, "%s</RESOURCE>\n", indentation[1]);
		fprintf(fp, "%s</VOTABLE>\n", indentation[0]);
	}
	else if(format == CATALOG_FORMAT_SQL)
	{
		// Write SQL catalogue
		const char *catalog_name = "SoFiA-Catalogue";
		
		fprintf(fp, "-- SoFiA source catalogue\n-- Creator: %s\n-- Time:    %s\n\n", SOFIA_VERSION_FULL, current_time_string);
		fprintf(fp, "SET SQL_MODE = \"NO_AUTO_VALUE_ON_ZERO\";\n\n");
		fprintf(fp, "CREATE TABLE IF NOT EXISTS `%s` (\n", catalog_name);
		fprintf(fp, "\t`name` VARCHAR(255) NOT NULL,\n");
		
		for(size_t j = 0; j < Source_get_num_par(src0); ++j)
		{
			if(Source_get_type(src0, j) == SOURCE_TYPE_INT) fprintf(fp, "\t`%s` INTEGER NOT NULL,\n", Source_get_name(src0, j));
			else fprintf(fp, "\t`%s` DOUBLE PRECISION NOT NULL,\n", Source_get_name(src0, j));
		}
		
		fprintf(fp, "\tPRIMARY KEY (`id`),\n\tKEY (`id`)\n) COMMENT=\'SoFiA source catalogue; created with SoFiA version %s\';\n\n", SOFIA_VERSION);
		fprintf(fp, "INSERT INTO `SoFiA-Catalogue` (`name`, ");
		
		for(size_t j = 0; j < Source_get_num_par(src0); ++j)
		{
			fprintf(fp, "`%s`", Source_get_name(src0, j));
			if(j + 1 < Source_get_num_par(src0)) fprintf(fp, ", ");
			else fprintf(fp, ") VALUES\n");
		}
		
		// Loop over all sources to write parameters
		for(size_t i = 0; i < self->size; ++i)
		{
			fprintf(fp, "(");
			
			Source *src = self->sources[i];
			
			fprintf(fp, "\'%s\', ", Source_get_identifier(src));
			
			for(size_t j = 0; j < Source_get_num_par(src); ++j)
			{
				if(Source_get_type(src, j) == SOURCE_TYPE_INT) fprintf(fp, "%ld", Source_get_par_int(src, j));
				else fprintf(fp, "%.15e", Source_get_par_flt(src, j));
				if(j + 1 < Source_get_num_par(src)) fprintf(fp, ", ");
			}
			
			if(i + 1 < self->size) fprintf(fp, "),\n");
			else fprintf(fp, ");\n");
		}
	}
	else
	{
		// Write ASCII catalogue
		fprintf(fp, "# SoFiA source catalogue\n# Creator: %s\n# Time:    %s\n#\n", SOFIA_VERSION_FULL, current_time_string);
		fprintf(fp, "# Header rows:\n#   1 = column number\n#   2 = parameter name\n#   3 = parameter unit\n%c\n%c", char_comment, char_comment);
		
		fprintf(fp, "%*d", 2 * CATALOG_COLUMN_WIDTH, 1);
		for(size_t j = 0; j < Source_get_num_par(src0); ++j) fprintf(fp, "%*zu", CATALOG_COLUMN_WIDTH, j + 2);
		fprintf(fp, "\n%c", char_comment);
		
		fprintf(fp, "%*s", 2 * CATALOG_COLUMN_WIDTH, "name");
		for(size_t j = 0; j < Source_get_num_par(src0); ++j) fprintf(fp, "%*s", CATALOG_COLUMN_WIDTH, Source_get_name(src0, j));
		fprintf(fp, "\n%c", char_comment);
		
		fprintf(fp, "%*s", 2 * CATALOG_COLUMN_WIDTH, " ");
		for(size_t j = 0; j < Source_get_num_par(src0); ++j) fprintf(fp, "%*s", CATALOG_COLUMN_WIDTH, Source_get_unit(src0, j));
		fprintf(fp, "\n\n");
		
		// Loop over all sources to write parameters
		for(size_t i = 0; i < self->size; ++i)
		{
			Source *src = self->sources[i];
			
			fprintf(fp, "%c", char_nocomment);
			fprintf(fp, "%*s", 2 * CATALOG_COLUMN_WIDTH, Source_get_identifier(src));
			
			for(size_t j = 0; j < Source_get_num_par(src); ++j)
			{
				if(Source_get_type(src, j) == SOURCE_TYPE_INT)
				{
					// Integer value
					const long int value = Source_get_par_int(src, j);
					fprintf(fp, "%*ld", CATALOG_COLUMN_WIDTH, value);
				}
				else
				{
					// Floating-point value
					const double value = Source_get_par_flt(src, j);
					if(value != 0.0 && (fabs(value) >= 1.0e+4 || fabs(value) < 1.0e-3)) fprintf(fp, "%*.5e", CATALOG_COLUMN_WIDTH, value);
					else fprintf(fp, "%*.6f", CATALOG_COLUMN_WIDTH, value);
				}
			}
			
			fprintf(fp, "\n");
		}
	}
	
	fclose(fp);
	
	return;
}



// ----------------------------------------------------------------- //
// Reallocate memory for one additional source                       //
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
//   Private method for allocating additional memory for one more    //
//   source in the specified catalogue. Note that this will not cre- //
//   ate a new source yet, but just allocate the memory needed to    //
//   append a source at the end of the catalogue. The function       //
//   should be called from public member functions that will add     //
//   sources to a catalogue prior to inserting the new source.       //
// ----------------------------------------------------------------- //

PRIVATE void Catalog_append_memory(Catalog *self)
{
	self->sources = (Source **)memory_realloc(self->sources, ++(self->size), sizeof(Source *));
	return;
}
