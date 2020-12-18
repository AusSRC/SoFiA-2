/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (Table.c) - Source Finding Application                   ///
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
#include <string.h>
#include <ctype.h>
#include "Table.h"



// ----------------------------------------------------------------- //
// Declaration of properties of class Table                          //
// ----------------------------------------------------------------- //

CLASS Table
{
	size_t cols;
	size_t rows;
	double *data;
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
//   Pointer to newly created Table object.                          //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Standard constructor. Will create a new and empty Table object. //
//   Note that the destructor will have to be called explicitly once //
//   the object is no longer required to release its memory again.   //
//   NOTE: The standard constructor has been made private to disable //
//   the creation of empty Table objects, as there are no methods to //
//   change the size of a Table object. The Table_from_file() con-   //
//   structor must instead be used to directly generate Tables from  //
//   tabulated data stored in a text file.                           //
// ----------------------------------------------------------------- //

PRIVATE Table *Table_new(void)
{
	Table *self = (Table *)memory(MALLOC, 1, sizeof(Table));
	
	self->cols = 0;
	self->rows = 0;
	self->data = NULL;
	
	return self;
}



// ----------------------------------------------------------------- //
// Alternative constructor for creating tables from text file        //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) filename   - Path to the file from which to read data.      //
//   (2) delimiters - List of delimiting characters to be used to    //
//                    separate columns in the data file.             //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Pointer to newly created Table object.                          //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Alternative constructor. Will create a new table object from    //
//   tabulated data read from the specified text file. The list of   //
//   delimiters to be used to separate data columns can be specified //
//   using the delimiters string (e.g. " \t" will use space and tab  //
//   characters as delimiters). Consecutive delimiting characters    //
//   will always be merged even when of mixed type. A new object     //
//   filled with the data from the file will be returned.            //
//   If no valid data are found in the specified file, then an empty //
//   Table object with 0 rows and columns will be returned. The size //
//   of the table will be defined by the first data row in the file. //
//   All subsequent rows must have the same number of columns. If    //
//   fewer columns are encountered, the method will terminate with   //
//   an error message. If more columns are encountered, then any ex- //
//   cess columns will be silently ignored. If an entry is encoun-   //
//   tered that does not constitute a floating-point number, then    //
//   the resulting table entry will be set to 0.                     //
//   The data file can contain empty lines and lines beginning with  //
//   a comment character (#); these will be ignored.                 //
// ----------------------------------------------------------------- //

PUBLIC Table *Table_from_file(const char *filename, const char *delimiters)
{
	// Sanity checks
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "Empty file name provided.");
	check_null(delimiters);
	
	// Allocate memory for a single line
	char *line = (char *)memory(MALLOC, TABLE_MAX_LINE_SIZE, sizeof(char));
	size_t counter = 0;
	
	// Try to open file
	FILE *fp = fopen(filename, "r");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open input file: %s.", filename);
	
	// Read beginning of file to establish number of columns
	while(fgets(line, TABLE_MAX_LINE_SIZE, fp))
	{
		// Trim line and check for comments and empty lines
		char *trimmed = trim_string(line);
		if(strlen(trimmed) == 0 || !isalnum(trimmed[0])) continue;
		
		char *entry = strtok(line, delimiters);
		while(entry != NULL)
		{
			++counter;
			entry = strtok(NULL, delimiters);
		}
		
		if(counter > 0) break;
	}
	
	// Return empty table if no valid data found
	if(counter == 0)
	{
		warning("No valid data found in file %s.\n         Returning empty table.", filename);
		fclose(fp);
		return Table_new();
	}
	
	// Otherwise reset file pointer
	rewind(fp);
	
	// Create a new table
	Table *self = Table_new();
	self->cols  = counter;
	
	while(fgets(line, TABLE_MAX_LINE_SIZE, fp))
	{
		// Trim line and check for comments and empty lines
		char *trimmed = trim_string(line);
		if(strlen(trimmed) == 0 || !isalnum(trimmed[0])) continue;
		
		// Add row
		self->rows += 1;
		self->data = memory_realloc(self->data, self->rows * self->cols, sizeof(double));
		
		char *entry = strtok(trimmed, delimiters);
		for(size_t i = 0; i < self->cols; ++i)
		{
			ensure(entry != NULL, ERR_USER_INPUT, "Inconsistent number of data columns in file %s.\n       %zu columns expected, but only %zu columns found in data row %zu.", filename, self->cols, i, self->rows);
			self->data[TABLE_INDEX(self->rows - 1, i)] = strtod(entry, NULL);
			entry = strtok(NULL, delimiters);
		}
	}
	
	// Close file again and clean up
	fclose(fp);
	free(line);
	
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
//   if the object is no longer required. This will release all me-  //
//   mory occupied by the object.                                    //
// ----------------------------------------------------------------- //

PUBLIC void Table_delete(Table *self)
{
	if(self != NULL) free(self->data);
	free(self);
	return;
}



// ----------------------------------------------------------------- //
// Return the number of table rows or columns                        //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self       - Object self-reference.                         //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Number of table rows.                                           //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public methods for retrieving the number of rows or columns in  //
//   the specified Table object.                                     //
// ----------------------------------------------------------------- //

PUBLIC size_t Table_rows(const Table *self)
{
	return self == NULL ? 0 : self->rows;
}

PUBLIC size_t Table_cols(const Table *self)
{
	return self == NULL ? 0 : self->cols;
}



// ----------------------------------------------------------------- //
// Get table value at specified row and column                       //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self       - Object self-reference.                         //
//   (2) row        - Requested row.                                 //
//   (3) column     - Requested column.                              //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   Table value at specified row and column.                        //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for retrieving the table entry at the specified   //
//   row and column of the Table object. The method will terminate   //
//   with an error message if the requested row or column is out of  //
//   range.                                                          //
// ----------------------------------------------------------------- //

PUBLIC double Table_get(const Table *self, const size_t row, const size_t col)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Requested Table column or row out of range.");
	
	return self->data[TABLE_INDEX(row, col)];
}



// ----------------------------------------------------------------- //
// Set table value at specified row and column                       //
// ----------------------------------------------------------------- //
// Arguments:                                                        //
//                                                                   //
//   (1) self       - Object self-reference.                         //
//   (2) row        - Requested row.                                 //
//   (3) column     - Requested column.                              //
//   (4) value      - Value to write into the table.                 //
//                                                                   //
// Return value:                                                     //
//                                                                   //
//   No return value.                                                //
//                                                                   //
// Description:                                                      //
//                                                                   //
//   Public method for writing the specified value into the speci-   //
//   fied row and column of the Table object, thus overwriting any   //
//   existing value. The method will terminate with an error message //
//   if the requested row or column is out of range.                 //
// ----------------------------------------------------------------- //

PUBLIC void Table_set(Table *self, const size_t row, const size_t col, const double value)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Requested Table column or row out of range.");
	
	self->data[TABLE_INDEX(row, col)] = value;
	return;
}
