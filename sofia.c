/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.2.1 (sofia.c) - Source Finding Application                   ///
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
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>

// WARNING: The following is compiler-specific and
//          not part of the C standard library.
#ifdef _OPENMP
	#include <omp.h>
#endif

// WARNING: The following will only work on POSIX-compliant
//          systems, but is needed for mkdir().
#include <errno.h>
#include <sys/stat.h>

//#include "sofia.h"
#include "src/common.h"
#include "src/Table.h"
#include "src/Path.h"
#include "src/Array_dbl.h"
#include "src/Array_siz.h"
#include "src/Map.h"
#include "src/Matrix.h"
#include "src/Parameter.h"
#include "src/Catalog.h"
#include "src/WCS.h"
#include "src/DataCube.h"
#include "src/LinkerPar.h"


// ----------------------------------------------------------------- //
// This file contains the actual SoFiA pipeline that will read user  //
// parameters and data files, call the requested processing modules  //
// and write out catalogues and images.                              //
// ----------------------------------------------------------------- //

//int mainline(char *par_file, char *argv2, char *argv3)
//void mainline(const sofia_args *arg)
void mainline(double *dataPtr, char *headerPtr, int datasize, int headersize)
{

	SOURCETYPE DATATYPE = MEM;
	char *path_to_par = "sofia.par";
/*
	char *headerPtr = arg->headerPtr;
	double *dataPtr = arg->dataPtr;
	int datasize = arg->datasize;
	int headersize = arg->headersize;
*/
	// If headerPtr is NULL, we've been called in FITS mode
	if (headerPtr == NULL) {
		DATATYPE = FITS;
		path_to_par = (char *)dataPtr;
	}

	// ---------------------------- //
	// Record starting time         //
	// ---------------------------- //
	
	const time_t start_time = time(NULL);
	const clock_t start_clock = clock();
	
	
	
	// ---------------------------- //
	// A few global definitions     //
	// ---------------------------- //
	
	const char *noise_stat_name[] = {"standard deviation", "median absolute deviation", "Gaussian fit to flux histogram"};
	const char *average_stat_name[] = {"mean", "median"};
	const char *flux_range_name[] = {"negative", "full", "positive"};
	double global_rms = 1.0;
	#ifdef _OPENMP
		const int n_cpu_cores = omp_get_num_procs();
	#endif
	
	
	
	// ---------------------------- //
	// Print basic information      //
	// ---------------------------- //
	
	status("Pipeline started");
	message("Using:    Source Finding Application (SoFiA)");
	message("Version:  %s (%s)", SOFIA_VERSION, SOFIA_CREATION_DATE);
	#ifdef _OPENMP
		message("CPU:      %d %s available", n_cpu_cores, n_cpu_cores == 1 ? "core" : "cores");
	#else
		message("CPU:      OpenMP disabled");
	#endif
	message("Time:     %s", ctime(&start_time));
	
	
	// ---------------------------- //
	// Set default parameters       //
	// ---------------------------- //
	
	status("Loading parameter settings");
	
	message("Activating SoFiA default parameter settings.");
	
	Parameter *par = Parameter_new(false);
	Parameter_default(par);
	
	
	
	// ---------------------------- //
	// Load user parameters         //
	// ---------------------------- //
	
	//message("Loading user parameter file: \'%s\'.\n", argv[1]);
	//Parameter_load(par, argv[1], PARAMETER_UPDATE);
	message("Loading user parameter file: \'%s\'.\n", path_to_par);
	Parameter_load(par, path_to_par, PARAMETER_UPDATE);
	
	
	
	// ---------------------------- //
	// Set number of threads        //
	// ---------------------------- //
	
	#ifdef _OPENMP
		const int n_threads = Parameter_get_int(par, "pipeline.threads");
		if(n_threads > 0 && n_threads < n_cpu_cores)
		{
			omp_set_num_threads(n_threads);
			message("Using %d out of %d available CPU cores.\n", n_threads, n_cpu_cores);
		}
		else
		{
			omp_set_num_threads(n_cpu_cores);
			message("Using all %d available CPU cores.\n", n_cpu_cores);
		}
	#else
		warning("Multi-threading is currently disabled. To enable it, please re-\n         install SoFiA with the '-fopenmp' option.");
	#endif
	
	
	
	// ---------------------------- //
	// Extract important settings   //
	// ---------------------------- //
	
	const bool verbosity         = Parameter_get_bool(par, "pipeline.verbose");
	const bool use_region        = strlen(Parameter_get_str(par, "input.region")) ? true : false;
	const bool use_gain          = strlen(Parameter_get_str(par, "input.gain"))   ? true : false;
	const bool use_noise         = strlen(Parameter_get_str(par, "input.noise"))  ? true : false;
	const bool use_weights       = strlen(Parameter_get_str(par, "input.weights"))? true : false;
	const bool use_mask          = strlen(Parameter_get_str(par, "input.mask"))   ? true : false;
	const bool use_invert        = Parameter_get_bool(par, "input.invert");
	      bool use_flagging      = strlen(Parameter_get_str(par, "flag.region"))  ? true : false;
	const bool use_flagging_cat  = strlen(Parameter_get_str(par, "flag.catalog")) ? true : false;
	const bool autoflag_log      = Parameter_get_bool(par, "flag.log");
	const bool use_cont_sub      = Parameter_get_bool(par, "contsub.enable");
	const bool use_noise_scaling = Parameter_get_bool(par, "scaleNoise.enable");
	const bool use_sc_scaling    = Parameter_get_bool(par, "scaleNoise.scfind");
	const bool use_spat_filter   = Parameter_get_bool(par, "spatFilter.enable");
	const bool use_scfind        = Parameter_get_bool(par, "scfind.enable");
	const bool use_threshold     = Parameter_get_bool(par, "threshold.enable");
	const bool keep_negative     = Parameter_get_bool(par, "linker.keepNegative");
	const bool use_reliability   = Parameter_get_bool(par, "reliability.enable");
	const bool use_rel_plot      = Parameter_get_bool(par, "reliability.plot");
	const bool use_rel_cat       = strlen(Parameter_get_str(par, "reliability.catalog")) ? true : false;
	const bool use_mask_dilation = Parameter_get_bool(par, "dilation.enable");
	const bool use_parameteriser = Parameter_get_bool(par, "parameter.enable");
	const bool use_wcs           = Parameter_get_bool(par, "parameter.wcs");
	const bool use_physical      = Parameter_get_bool(par, "parameter.physical");
	const bool use_pos_offset    = Parameter_get_bool(par, "parameter.offset");
	
	const bool write_ascii       = Parameter_get_bool(par, "output.writeCatASCII");
	const bool write_xml         = Parameter_get_bool(par, "output.writeCatXML");
	const bool write_sql         = Parameter_get_bool(par, "output.writeCatSQL");
	const bool write_noise       = Parameter_get_bool(par, "output.writeNoise");
	const bool write_filtered    = Parameter_get_bool(par, "output.writeFiltered");
	const bool write_mask        = Parameter_get_bool(par, "output.writeMask");
	const bool write_mask2d      = Parameter_get_bool(par, "output.writeMask2d");
	const bool write_rawmask     = Parameter_get_bool(par, "output.writeRawMask");
	const bool write_moments     = Parameter_get_bool(par, "output.writeMoments");
	const bool write_cubelets    = Parameter_get_bool(par, "output.writeCubelets");
	const bool overwrite         = Parameter_get_bool(par, "output.overwrite");
	
	const double rel_threshold   = Parameter_get_flt(par, "reliability.threshold");
	const double rel_fmin        = Parameter_get_flt(par, "reliability.fmin");
	
	// For defining the type of data source - file, memory etc
	if (strcmp(Parameter_get_str(par,"input.source"), "FITS") == 0) DATATYPE = FITS;
	else if (strcmp(Parameter_get_str(par,"input.source"), "MEM") == 0) DATATYPE = MEM;

	unsigned int autoflag_mode = 0;
	if     (strcmp(Parameter_get_str(par, "flag.auto"), "channels") == 0) autoflag_mode = 1;
	else if(strcmp(Parameter_get_str(par, "flag.auto"), "pixels")   == 0) autoflag_mode = 2;
	else if(strcmp(Parameter_get_str(par, "flag.auto"), "true")     == 0) autoflag_mode = 3;
	
	// Statistic and range for noise measurement in
	// noise scaling, S+C finder and threshold finder
	noise_stat sn_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "scaleNoise.statistic"), "mad") == 0) sn_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "scaleNoise.statistic"), "gauss") == 0) sn_statistic = NOISE_STAT_GAUSS;
	
	int sn_range = 0;
	if(strcmp(Parameter_get_str(par, "scaleNoise.fluxRange"), "negative") == 0) sn_range = -1;
	else if(strcmp(Parameter_get_str(par, "scaleNoise.fluxRange"), "positive") == 0) sn_range = 1;
	
	noise_stat sc_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "scfind.statistic"), "mad") == 0) sc_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "scfind.statistic"), "gauss") == 0) sc_statistic = NOISE_STAT_GAUSS;
	
	int sc_range = 0;
	if(strcmp(Parameter_get_str(par, "scfind.fluxRange"), "negative") == 0) sc_range = -1;
	else if(strcmp(Parameter_get_str(par, "scfind.fluxRange"), "positive") == 0) sc_range = 1;
	
	noise_stat tf_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "threshold.statistic"), "mad") == 0) tf_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "threshold.statistic"), "gauss") == 0) tf_statistic = NOISE_STAT_GAUSS;
	
	int tf_range = 0;
	if(strcmp(Parameter_get_str(par, "threshold.fluxRange"), "negative") == 0) tf_range = -1;
	else if(strcmp(Parameter_get_str(par, "threshold.fluxRange"), "positive") == 0) tf_range = 1;
	
	// Spatial filter statistic
	int spat_filter_statistic = 0;
	if(strcmp(Parameter_get_str(par, "spatFilter.statistic"), "median") == 0) spat_filter_statistic = 1;
	
	// Noise and weights sanity check
	if(use_noise && use_weights) warning("Applying both a weights cube and a noise cube.");
	
	// Negative detections sanity check
	ensure(!keep_negative || !use_reliability, ERR_USER_INPUT, "With the reliability filter enabled, negative detections would always\n       be discarded irrespective of the value of linker.keepNegative! Please\n       set linker.keepNegative = false or disable reliability filtering.");
	
	
	
	// ---------------------------- //
	// Define file names and paths  //
	// ---------------------------- //
	
	const char *base_dir  = Parameter_get_str(par, "output.directory");
	const char *base_name = Parameter_get_str(par, "output.filename");
	
	// Set up input paths
	Path *path_data_in = Path_new();
	Path_set(path_data_in, Parameter_get_str(par, "input.data"));
	
	Path *path_gain_in = Path_new();
	if(use_gain) Path_set(path_gain_in, Parameter_get_str(par, "input.gain"));
	
	Path *path_noise_in = Path_new();
	if(use_noise) Path_set(path_noise_in, Parameter_get_str(par, "input.noise"));
	
	Path *path_weights_in = Path_new();
	if(use_weights) Path_set(path_weights_in, Parameter_get_str(par, "input.weights"));
	
	Path *path_mask_in = Path_new();
	if(use_mask) Path_set(path_mask_in, Parameter_get_str(par, "input.mask"));
	
	// Choose appropriate output file and directory names depending on user input
	String *output_file_name = String_new(strlen(base_name) ? base_name : Path_get_file(path_data_in));
	String *output_dir_name  = String_new(strlen(base_dir) ? base_dir : (strlen(Path_get_dir(path_data_in)) ? Path_get_dir(path_data_in) : "."));
	
	// Ensure that output file name ends with .fits or similar
	String *check_mime_type = String_new("");
	String_to_lower(String_set_delim(check_mime_type, String_get(output_file_name), '.', false, false));
	if(!String_compare(check_mime_type, "fits") && !String_compare(check_mime_type, "fit")) String_append(output_file_name, ".fits");
	String_delete(check_mime_type);
	
	// Set up output paths
	Path *path_cat_ascii = Path_new();
	Path *path_cat_xml   = Path_new();
	Path *path_cat_sql   = Path_new();
	Path *path_noise_out = Path_new();
	Path *path_filtered  = Path_new();
	Path *path_mask_out  = Path_new();
	Path *path_mask_2d   = Path_new();
	Path *path_mask_raw  = Path_new();
	Path *path_mom0      = Path_new();
	Path *path_mom1      = Path_new();
	Path *path_mom2      = Path_new();
	Path *path_chan      = Path_new();
	Path *path_rel_plot  = Path_new();
	Path *path_skel_plot = Path_new();
	Path *path_flag      = Path_new();
	Path *path_cubelets  = Path_new();
	
	// Set up output directory names
	Path_set_dir(path_cat_ascii, String_get(output_dir_name));
	Path_set_dir(path_cat_xml,   String_get(output_dir_name));
	Path_set_dir(path_cat_sql,   String_get(output_dir_name));
	Path_set_dir(path_noise_out, String_get(output_dir_name));
	Path_set_dir(path_filtered,  String_get(output_dir_name));
	Path_set_dir(path_mask_out,  String_get(output_dir_name));
	Path_set_dir(path_mask_2d,   String_get(output_dir_name));
	Path_set_dir(path_mask_raw,  String_get(output_dir_name));
	Path_set_dir(path_mom0,      String_get(output_dir_name));
	Path_set_dir(path_mom1,      String_get(output_dir_name));
	Path_set_dir(path_mom2,      String_get(output_dir_name));
	Path_set_dir(path_chan,      String_get(output_dir_name));
	Path_set_dir(path_rel_plot,  String_get(output_dir_name));
	Path_set_dir(path_skel_plot, String_get(output_dir_name));
	Path_set_dir(path_flag,      String_get(output_dir_name));
	Path_set_dir(path_cubelets,  String_get(output_dir_name));
	
	// Set up output file names
	Path_set_file_from_template(path_cat_ascii,  String_get(output_file_name), "_cat",      ".txt");
	Path_set_file_from_template(path_cat_xml,    String_get(output_file_name), "_cat",      ".xml");
	Path_set_file_from_template(path_cat_sql,    String_get(output_file_name), "_cat",      ".sql");
	Path_set_file_from_template(path_noise_out,  String_get(output_file_name), "_noise",    ".fits");
	Path_set_file_from_template(path_filtered,   String_get(output_file_name), "_filtered", ".fits");
	Path_set_file_from_template(path_mask_out,   String_get(output_file_name), "_mask",     ".fits");
	Path_set_file_from_template(path_mask_2d,    String_get(output_file_name), "_mask-2d",  ".fits");
	Path_set_file_from_template(path_mask_raw,   String_get(output_file_name), "_mask-raw", ".fits");
	Path_set_file_from_template(path_mom0,       String_get(output_file_name), "_mom0",     ".fits");
	Path_set_file_from_template(path_mom1,       String_get(output_file_name), "_mom1",     ".fits");
	Path_set_file_from_template(path_mom2,       String_get(output_file_name), "_mom2",     ".fits");
	Path_set_file_from_template(path_chan,       String_get(output_file_name), "_chan",     ".fits");
	Path_set_file_from_template(path_rel_plot,   String_get(output_file_name), "_rel",      ".eps");
	Path_set_file_from_template(path_skel_plot,  String_get(output_file_name), "_skellam",  ".eps");
	Path_set_file_from_template(path_flag,       String_get(output_file_name), "_flags",    ".log");
	
	// Set up cubelet directory and file base name
	Path_append_dir_from_template(path_cubelets, String_get(output_file_name), "_cubelets");
	Path_set_file_from_template(path_cubelets,   String_get(output_file_name), "", "");
	
	// Delete temporary strings again
	String_delete(output_file_name);
	String_delete(output_dir_name);
	
	
	
	// ---------------------------- //
	// Check output settings        //
	// ---------------------------- //
	
	// Try to create cubelet directory
	if(write_cubelets)
	{
		errno = 0;
		mkdir(Path_get_dir(path_cubelets), 0755);
		ensure(errno == 0 || errno == EEXIST, ERR_FILE_ACCESS, "Failed to create cubelet directory; please check write permissions.");
	}
	
	// Check overwrite conditions
	if(!overwrite)
	{
		if(write_cubelets) {
			ensure(errno != EEXIST, ERR_FILE_ACCESS,
				"Cubelet directory already exists. Please delete the directory\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_ascii) {
			ensure(!Path_file_is_readable(path_cat_ascii), ERR_FILE_ACCESS,
				"ASCII catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_xml) {
			ensure(!Path_file_is_readable(path_cat_xml), ERR_FILE_ACCESS,
				"XML catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_sql) {
			ensure(!Path_file_is_readable(path_cat_sql), ERR_FILE_ACCESS,
				"SQL catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_noise) {
			ensure(!Path_file_is_readable(path_noise_out), ERR_FILE_ACCESS,
				"Noise cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_filtered) {
			ensure(!Path_file_is_readable(path_filtered), ERR_FILE_ACCESS,
				"Filtered cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_mask) {
			ensure(!Path_file_is_readable(path_mask_out), ERR_FILE_ACCESS,
				"Mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_mask2d) {
			ensure(!Path_file_is_readable(path_mask_2d), ERR_FILE_ACCESS,
				"2-D mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_rawmask) {
			ensure(!Path_file_is_readable(path_mask_raw), ERR_FILE_ACCESS,
				"Raw mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_moments) {
			ensure(!Path_file_is_readable(path_mom0) && !Path_file_is_readable(path_mom1) && !Path_file_is_readable(path_mom2), ERR_FILE_ACCESS,
				"Moment maps already exist. Please delete the files\n"
				"       or set \'output.overwrite = true\'.");
			ensure(!Path_file_is_readable(path_chan), ERR_FILE_ACCESS,
				"Channel map already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(use_reliability && use_rel_plot) {
			ensure(!Path_file_is_readable(path_rel_plot), ERR_FILE_ACCESS,
				"Reliability plot already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
			/*ensure(!Path_file_is_readable(path_skel_plot), ERR_FILE_ACCESS,
				"Skellam plot already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");*/
		}
		if(autoflag_log) {
			ensure(!Path_file_is_readable(path_flag), ERR_FILE_ACCESS,
				"Flagging log file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
	}
	
	
	
	// ---------------------------- //
	// Load data cube               //
	// ---------------------------- //
	// Set up region if required
	Array_siz *region = use_region ? Array_siz_new_str(Parameter_get_str(par, "input.region")) : NULL;
	
	// Set up flagging region if required
	Array_siz *flag_regions = use_flagging ? Array_siz_new_str(Parameter_get_str(par, "flag.region")) : Array_siz_new(0);
	
	// Load data cube
	status("Loading data cube");
	DataCube *dataCube = DataCube_new(verbosity);
	if (DATATYPE == FITS) {
		status("    - FITS mode");
		DataCube_readFITS(dataCube, Path_get(path_data_in), region);
	}
	else {
		status("    - MEM mode");
		DataCube_readMEM(dataCube, dataPtr,datasize,headerPtr,headersize);
	}
	// Search for values of infinity and append affected pixels to flagging region
	// (Yes, some data cubes do contain those!)
	if(DataCube_flag_infinity(dataCube, flag_regions)) use_flagging = true;
	
	// Apply flags if required
	if(use_flagging) DataCube_flag_regions(dataCube, flag_regions);
	
	// Invert cube if requested
	if(use_invert)
	{
		message("Inverting data cube");
		DataCube_multiply_const(dataCube, -1.0);
	}
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Apply flagging catalogue     //
	// ---------------------------- //
	
	if(use_flagging_cat)
	{
		status("Loading and applying flagging catalogue");
		message("Catalogue file:   %s", Parameter_get_str(par, "flag.catalog"));
		message("Flagging radius:  %ld", Parameter_get_int(par, "flag.radius"));
		DataCube_continuum_flagging(dataCube, Parameter_get_str(par, "flag.catalog"), 1, Parameter_get_int(par, "flag.radius"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load and apply noise cube    //
	// ---------------------------- //
	
	if(use_noise)
	{
		status("Loading and applying noise cube");
		DataCube *noiseCube = DataCube_new(verbosity);
		DataCube_readFITS(noiseCube, Path_get(path_noise_in), region);
		
		// Divide data by noise cube
		DataCube_divide(dataCube, noiseCube);
		
		// Delete noise cube again
		DataCube_delete(noiseCube);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load and apply weights cube  //
	// ---------------------------- //
	
	if(use_weights)
	{
		status("Loading and applying weights cube");
		DataCube *weightsCube = DataCube_new(verbosity);
		DataCube_readFITS(weightsCube, Path_get(path_weights_in), region);
		
		// Multiply data by square root of weights cube
		DataCube_apply_weights(dataCube, weightsCube);
		
		// Delete weights cube again
		DataCube_delete(weightsCube);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Continuum subtraction        //
	// ---------------------------- //
	
	if(use_cont_sub)
	{
		status("Continuum subtraction");
		message("Subtracting residual continuum emission.");
		message("- Polynomial order:  %ld",   Parameter_get_int(par, "contsub.order"));
		message("- Clip threshold:    %.1f",  Parameter_get_flt(par, "contsub.threshold"));
		message("- Shift:             %ld",   Parameter_get_int(par, "contsub.shift"));
		message("- Padding:           %ld\n", Parameter_get_int(par, "contsub.padding"));
		
		DataCube_contsub(dataCube, Parameter_get_int(par, "contsub.order"), Parameter_get_int(par, "contsub.shift"), Parameter_get_int(par, "contsub.padding"), Parameter_get_flt(par, "contsub.threshold"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Scale data by noise level    //
	// ---------------------------- //
	
	if(use_noise_scaling)
	{
		status("Scaling data by noise");
		
		if(strcmp(Parameter_get_str(par, "scaleNoise.mode"), "local") == 0)
		{
			// Local noise scaling
			message("Correcting for local noise variations.");
			message("- Noise statistic:  %s", noise_stat_name[sn_statistic]);
			message("- Flux range:       %s\n", flux_range_name[sn_range + 1]);
			
			DataCube *noiseCube = DataCube_scale_noise_local(
				dataCube,
				sn_statistic,
				sn_range,
				Parameter_get_int(par, "scaleNoise.windowXY"),
				Parameter_get_int(par, "scaleNoise.windowZ"),
				Parameter_get_int(par, "scaleNoise.gridXY"),
				Parameter_get_int(par, "scaleNoise.gridZ"),
				Parameter_get_bool(par, "scaleNoise.interpolate")
			);
			
			if(write_noise)
			{
				// Apply flags to noise cube
				if(use_flagging) DataCube_flag_regions(noiseCube, flag_regions);
				DataCube_save(noiseCube, Path_get(path_noise_out), overwrite, DESTROY);
			}
			DataCube_delete(noiseCube);
		}
		else
		{
			// Global noise scaling along spectral axis
			message("Correcting for noise variations along spectral axis.");
			message("- Noise statistic:  %s", noise_stat_name[sn_statistic]);
			message("- Flux range:       %s\n", flux_range_name[sn_range + 1]);
			DataCube_scale_noise_spec(dataCube, sn_statistic, sn_range);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Automatic data flagging      //
	// ---------------------------- //
	
	if(autoflag_mode)
	{
		status("Auto-flagging");
		
		// Set up auto-flagging if requested
		Array_siz *autoflag_regions = Array_siz_new(0);
		DataCube_autoflag(dataCube, Parameter_get_flt(par, "flag.threshold"), autoflag_mode, autoflag_regions);
		
		const size_t size = Array_siz_get_size(autoflag_regions);
		
		// Apply flags if necessary
		if(size)
		{
			DataCube_flag_regions(dataCube, autoflag_regions);  // Apply auto-flagging regions
			Array_siz_cat(flag_regions, autoflag_regions);      // Append auto-flagging regions to general flagging regions
			use_flagging = true;                                // Update flagging switch
		}
		else message("No flagging required.");
		
		// Write auto-flags to log file if requested
		if(size && autoflag_log)
		{
			// Try to open output file
			FILE *fp;
			if(overwrite) fp = fopen(Path_get(path_flag), "wb");
			else fp = fopen(Path_get(path_flag), "wxb");
			
			// If successful...
			if(fp != NULL)
			{	
				// ...write out flags...
				message("Writing log file:     %s", Path_get_file(path_flag));
				fprintf(fp, "# Auto-flagging log file\n");
				fprintf(fp, "# Creator: %s\n#\n", SOFIA_VERSION_FULL);
				fprintf(fp, "# Flagging codes:\n");
				fprintf(fp, "#   C z            =  spectral channel (z)\n");
				fprintf(fp, "#   P x y          =  spatial pixel (x, y)\n");
				fprintf(fp, "#   R x1 x2 y1 y2  =  spatial region (x1:x2, y1:y2)\n");
				fprintf(fp, "# Note that coordinates will be relative to subregion\n");
				fprintf(fp, "# unless parameter.offset was set to true.\n\n");
				
				for(size_t i = 0; i < size; i += 6)
				{
					const size_t x_min = Array_siz_get(autoflag_regions, i);
					const size_t x_max = Array_siz_get(autoflag_regions, i + 1);
					const size_t y_min = Array_siz_get(autoflag_regions, i + 2);
					const size_t y_max = Array_siz_get(autoflag_regions, i + 3);
					const size_t z_min = Array_siz_get(autoflag_regions, i + 4);
					const size_t z_max = Array_siz_get(autoflag_regions, i + 5);
					
					// NOTE: Subregion offset will be added if requested (use_pos_offset == true).
					if(z_min == z_max)
					{
						fprintf(fp, "C %zu\n", z_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 4) : 0));
					}
					else if(x_min == x_max && y_min == y_max)
					{
						fprintf(fp, "P %zu %zu\n", x_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), y_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0));
					}
					else if(z_min == 0 && z_max == DataCube_get_axis_size(dataCube, 2) - 1)
					{
						fprintf(fp, "R %zu %zu %zu %zu\n", x_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), x_max + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), y_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0), y_max + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0));
					}
				}
				
				// ...and close output file again
				fclose(fp);
			}
			else warning("Failed to write flagging log file: %s", Path_get_file(path_flag));
		}
		
		// Clean up
		Array_siz_delete(autoflag_regions);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Spatial averaging filter     //
	// ---------------------------- //
	
	if(use_spat_filter)
	{
		status("Applying spatial filter");
		
		long int spat_filter_window = Parameter_get_int(par, "spatFilter.window");
		if(spat_filter_window < 30) warning("Adjusting window size to minimum of %zu.", spat_filter_window = 30);
		
		const long int spat_filter_kernel_size = Parameter_get_int(par, "spatFilter.boxcar");
		const long int spat_filter_kernel_radius = spat_filter_kernel_size ? spat_filter_kernel_size / 2 : 0;
		if(spat_filter_kernel_size && spat_filter_kernel_size % 2 == 0) warning("Forcing boxcar size to be odd.");
		
		message("Using the following parameters:");
		message("- Window size:   %ld x %ld", spat_filter_window, spat_filter_window);
		message("- Statistic:     %s", average_stat_name[spat_filter_statistic]);
		message("- Boxcar width:  %ld\n", 2 * spat_filter_kernel_radius + 1);
		
		DataCube_spatial_filter(dataCube, spat_filter_statistic, spat_filter_window, spat_filter_kernel_radius);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Write filtered cube          //
	// ---------------------------- //
	
	if(write_filtered && (use_region || use_flagging || use_flagging_cat || use_cont_sub || use_noise || use_weights || use_noise_scaling || use_spat_filter))  // ALERT: Add conditions here as needed.
	{
		status("Writing filtered cube");
		DataCube_save(dataCube, Path_get(path_filtered), overwrite, PRESERVE);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Measure global noise level   //
	// ---------------------------- //
	
	// NOTE: This is necessary so the linker and reliability module can
	//       divide all flux values by the RMS later on.
	// NOTE: This is currently being applied even when a noise cube has 
	//       been applied before or noise scaling is enabled.
	//       This is needed, as other algorithms, such as continuum sub-
	//       traction, might alter the noise level.
	
	status("Measuring global noise level");
	
	size_t cadence = DataCube_get_size(dataCube) / NOISE_SAMPLE_SIZE;          // Stride for noise calculation
	if(cadence < 2) cadence = 1;
	else if(cadence % DataCube_get_axis_size(dataCube, 0) == 0) cadence -= 1;  // Ensure stride is not equal to multiple of x-axis size
	
	global_rms = MAD_TO_STD * DataCube_stat_mad(dataCube, 0.0, cadence, -1);
	message("Global RMS:  %.3e  (using stride of %zu)", global_rms, cadence);
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Run source finder            //
	// ---------------------------- //
	
	// Terminate if no source finder is run, but no input mask is provided either
	ensure(use_scfind || use_threshold || use_mask, ERR_USER_INPUT, "No mask provided and no source finder selected. Cannot proceed.");
	
	// Create temporary 8-bit mask to hold source finding output
	DataCube *maskCubeTmp = DataCube_blank(DataCube_get_axis_size(dataCube, 0), DataCube_get_axis_size(dataCube, 1), DataCube_get_axis_size(dataCube, 2), 8, verbosity);
	DataCube_copy_wcs(dataCube, maskCubeTmp);
	DataCube_puthd_str(maskCubeTmp, "BUNIT", " ");
	
	// S+C finder
	if(use_scfind)
	{
		status("Running S+C finder");
		message("Using the following parameters:");
		message("- Kernels");
		message("  - spatial:        %s", Parameter_get_str(par, "scfind.kernelsXY"));
		message("  - spectral:       %s", Parameter_get_str(par, "scfind.kernelsZ"));
		message("- Flux threshold:   %s * rms", Parameter_get_str(par, "scfind.threshold"));
		message("- Noise statistic:  %s", noise_stat_name[sc_statistic]);
		message("- Flux range:       %s\n", flux_range_name[sc_range + 1]);
		
		// Extract and sort kernel sizes to ensure that smallest kernel comes first
		Array_dbl *kernels_spat = Array_dbl_new_str(Parameter_get_str(par, "scfind.kernelsXY"));
		Array_siz *kernels_spec = Array_siz_new_str(Parameter_get_str(par, "scfind.kernelsZ"));
		Array_dbl_sort(kernels_spat);
		Array_siz_sort(kernels_spec);
		
		// Sanity checks
		for(size_t i = 0; i < Array_dbl_get_size(kernels_spat); ++i)
		{
			const double ks = Array_dbl_get(kernels_spat, i);
			ensure(ks >= 0.0 && ks < DataCube_get_axis_size(dataCube, 0) && ks < DataCube_get_axis_size(dataCube, 1), ERR_USER_INPUT, "Illegal spatial kernel size encountered.");
			if(ks > 0.0 && ks < 3.0) warning("Spatial kernel sizes of < 3 cannot be accurately modelled.");
		}
		for(size_t i = 0; i < Array_siz_get_size(kernels_spec); ++i)
		{
			const size_t ks = Array_siz_get(kernels_spec, i);
			ensure(ks < DataCube_get_axis_size(dataCube, 2), ERR_USER_INPUT, "Illegal spectral kernel size encountered.");
			if(ks != 0 && ks % 2 == 0) warning("Spectral kernel size of %zu is even, will be treated as %zu!", ks, ks + 1);
			else if(ks == 1) warning("Spectral kernel size of 1 found, will be treated as 0!");
		}
		if(Array_dbl_get(kernels_spat, 0) > 0.0) warning("Including spatial kernel size of 0 is strongly advised.");
		if(Array_siz_get(kernels_spec, 0) > 0) warning("Including spectral kernel size of 0 is strongly advised.");
		
		// Run S+C finder to obtain mask
		DataCube_run_scfind(
			dataCube,
			maskCubeTmp,
			kernels_spat,
			kernels_spec,
			Parameter_get_flt(par, "scfind.threshold"),
			Parameter_get_flt(par, "scfind.replacement"),
			sc_statistic,
			sc_range,
			(use_noise_scaling && use_sc_scaling) ? (strcmp(Parameter_get_str(par, "scaleNoise.mode"), "local") == 0 ? 2 : 1) : 0,
			sn_statistic,
			sn_range,
			Parameter_get_int(par, "scaleNoise.windowXY"),
			Parameter_get_int(par, "scaleNoise.windowZ"),
			Parameter_get_int(par, "scaleNoise.gridXY"),
			Parameter_get_int(par, "scaleNoise.gridZ"),
			Parameter_get_bool(par, "scaleNoise.interpolate"),
			start_time,
			start_clock
		);
		
		// Clean up
		Array_dbl_delete(kernels_spat);
		Array_siz_delete(kernels_spec);
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCubeTmp, flag_regions);
	}
	
	// Threshold finder
	if(use_threshold)
	{
		// Determine mode
		const bool absolute = (strcmp(Parameter_get_str(par, "threshold.mode"), "absolute") == 0);
		
		status("Running threshold finder");
		message("Using the following parameters:");
		message("- Mode:             %s", absolute ? "absolute" : "relative");
		message("- Flux threshold:   %s%s", Parameter_get_str(par, "threshold.threshold"), absolute ? "" : " * rms");
		if(!absolute)
		{
			message("- Noise statistic:  %s", noise_stat_name[tf_statistic]);
			message("- Flux range:       %s", flux_range_name[tf_range + 1]);
		}
		
		// Run threshold finder
		DataCube_run_threshold(
			dataCube,
			maskCubeTmp,
			absolute,
			Parameter_get_flt(par, "threshold.threshold"),
			tf_statistic,
			tf_range
		);
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCubeTmp, flag_regions);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load mask cube if specified  //
	// ---------------------------- //
	
	DataCube *maskCube = NULL;
	
	if(use_mask)
	{
		// Load mask cube
		status("Loading mask cube");
		maskCube = DataCube_new(verbosity);
		DataCube_readFITS(maskCube, Path_get(path_mask_in), region);
		
		// Ensure that mask has the right type and size
		ensure(DataCube_gethd_int(maskCube, "BITPIX") == 32, ERR_USER_INPUT, "Mask cube must be of 32-bit integer type.");
		ensure(
			DataCube_gethd_int(maskCube, "NAXIS1") == DataCube_gethd_int(dataCube, "NAXIS1") &&
			DataCube_gethd_int(maskCube, "NAXIS2") == DataCube_gethd_int(dataCube, "NAXIS2") &&
			DataCube_gethd_int(maskCube, "NAXIS3") == DataCube_gethd_int(dataCube, "NAXIS3"),
			ERR_USER_INPUT, "Data cube and mask cube have different sizes."
		);
		
		// Set all masked pixels to -1
		DataCube_reset_mask_32(maskCube, -1);
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCube, flag_regions);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	else
	{
		// Else create an empty mask cube
		maskCube = DataCube_blank(DataCube_get_axis_size(dataCube, 0), DataCube_get_axis_size(dataCube, 1), DataCube_get_axis_size(dataCube, 2), 32, verbosity);
		
		// Copy WCS header elements from data cube to mask cube
		DataCube_copy_wcs(dataCube, maskCube);
		
		// Set BUNIT keyword of mask cube
		DataCube_puthd_str(maskCube, "BUNIT", " ");
	}
	
	
	
	// ---------------------------- //
	// Sort out masks               //
	// ---------------------------- //
	
	// Copy SF mask before linking
	const size_t n_pix_det = DataCube_copy_mask_8_32(maskCube, maskCubeTmp, -1);
	message("%zu pixels detected (%.3f%%).\n", n_pix_det, 100.0 * (double)(n_pix_det) / (double)(DataCube_get_size(maskCube)));
	
	// Write raw binary mask if requested
	if(write_rawmask)
	{
		status("Writing raw binary mask");
		DataCube_save(maskCubeTmp, Path_get(path_mask_raw), overwrite, DESTROY);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	// Delete temporary SF mask again
	DataCube_delete(maskCubeTmp);
	
	
	
	// ---------------------------- //
	// Run linker                   //
	// ---------------------------- //
	
	status("Running Linker");
	
	const bool remove_neg_src = !use_reliability && !keep_negative;  // ALERT: Add conditions here as needed.
	
	LinkerPar *lpar = DataCube_run_linker(
		dataCube,
		maskCube,
		Parameter_get_int(par, "linker.radiusXY"),
		Parameter_get_int(par, "linker.radiusXY"),
		Parameter_get_int(par, "linker.radiusZ"),
		Parameter_get_int(par, "linker.minSizeXY"),
		Parameter_get_int(par, "linker.minSizeXY"),
		Parameter_get_int(par, "linker.minSizeZ"),
		Parameter_get_int(par, "linker.maxSizeXY"),
		Parameter_get_int(par, "linker.maxSizeXY"),
		Parameter_get_int(par, "linker.maxSizeZ"),
		remove_neg_src,
		global_rms
	);
	
	// Print time
	timestamp(start_time, start_clock);
	
	// Terminate pipeline if no sources left after linking
	ensure(LinkerPar_get_size(lpar), ERR_NO_SRC_FOUND, "No sources left after linking. Terminating pipeline.");
	
	
	
	// ---------------------------- //
	// Run reliability filter       //
	// ---------------------------- //
	
	Map *rel_filter = Map_new();  // Empty container for storing old and new labels of reliable sources
	
	if(use_reliability)
	{
		status("Measuring reliability");
		
		// Check if catalogue supplied
		Table *rel_cat = NULL;
		if(use_rel_cat)
		{
			message("Reading in reliability catalogue.");
			
			// Read catalogue into table
			rel_cat = Table_from_file(Parameter_get_str(par, "reliability.catalog"), " \t,|");
			
			if(Table_rows(rel_cat) == 0 || Table_cols(rel_cat) != 2)
			{
				warning("Reliability catalogue non-compliant; must contain 2 data columns.\n         Catalogue file will be ignored.");
				Table_delete(rel_cat);
				rel_cat = NULL;
			}
			else
			{
				message("Extracting %zu position%s from catalogue.", Table_rows(rel_cat), Table_rows(rel_cat) > 1 ? "s" : "");
				
				// Extract WCS information
				WCS *wcs = DataCube_extract_wcs(dataCube);
				
				if(wcs != NULL)
				{
					// Loop over all rows and convert WCS to pixels
					for(size_t row = 0; row < Table_rows(rel_cat); ++row)
					{
						double lon = -1e+30;
						double lat = -1e+30;
						WCS_convertToPixel(wcs, Table_get(rel_cat, row, 0), Table_get(rel_cat, row, 1), 0.0, &lon, &lat, NULL);
						Table_set(rel_cat, row, 0, lon);
						Table_set(rel_cat, row, 1, lat);
					}
				}
				else
				{
					warning("WCS conversion failed; cannot apply reliability catalogue.");
					Table_delete(rel_cat);
					rel_cat = NULL;
				}
				
				// Delete WCS object again
				WCS_delete(wcs);
			}
		}
		
		// Calculate reliability values
		Matrix *covar = LinkerPar_reliability(lpar, Parameter_get_flt(par, "reliability.scaleKernel"), rel_fmin, rel_cat);
		
		// Create plots if requested
		if(use_rel_plot) LinkerPar_rel_plots(lpar, rel_threshold, rel_fmin, covar, Path_get(path_rel_plot), overwrite);
		
		// Delete covariance matrix and catalogue table again
		Matrix_delete(covar);
		Table_delete(rel_cat);
		
		// Set up relabelling filter by recording old and new label pairs of reliable sources
		size_t new_label = 1;
		
		for(size_t i = 0; i < LinkerPar_get_size(lpar); ++i)
		{
			const size_t old_label = LinkerPar_get_label(lpar, i);
			
			// Keep source if reliability > threshold and fmin parameter satisfied
			if(LinkerPar_get_rel(lpar, old_label) >= rel_threshold && LinkerPar_get_flux(lpar, old_label) / sqrt(LinkerPar_get_npix(lpar, old_label)) > rel_fmin) Map_push(rel_filter, old_label, new_label++);
		}
		
		// Check if any reliable sources left
		ensure(Map_get_size(rel_filter), ERR_NO_SRC_FOUND, "No reliable sources found. Terminating pipeline.");
		message("%zu reliable sources found.", Map_get_size(rel_filter));
		
		// Apply filter to mask cube, so unreliable sources are removed
		// and reliable ones relabelled in consecutive order
		DataCube_filter_mask_32(maskCube, rel_filter);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create initial catalogue     //
	// ---------------------------- //
	
	// Extract flux unit from header
	String *unit_flux = String_trim(DataCube_gethd_string(dataCube, "BUNIT"));
	if(!String_size(unit_flux))
	{
		warning("No flux unit (\'BUNIT\') defined in header.");
		String_set(unit_flux, "???");
	}
	
	// Generate catalogue of reliable sources from linker output
	Catalog *catalog = LinkerPar_make_catalog(lpar, rel_filter, String_get(unit_flux));
	
	// Delete linker parameters, reliability filter and flux unit string, as they are no longer needed
	LinkerPar_delete(lpar);
	Map_delete(rel_filter);
	String_delete(unit_flux);
	
	// Terminate if catalogue is empty
	ensure(Catalog_get_size(catalog), ERR_NO_SRC_FOUND, "No reliable sources found. Terminating pipeline.");
	
	
	
	// ---------------------------- //
	// Mask dilation if requested   //
	// ---------------------------- //
	// NOTE: It is not yet clear if mask dilation should happen in the noise-normalised data cube
	//       or the original data cube. Some more though will need to go into this...
	
	if(use_mask_dilation)
	{
		status("Mask dilation");
		
		message("Spectral dilation");
		DataCube_dilate_mask_z(dataCube, maskCube, catalog, Parameter_get_int(par, "dilation.iterationsZ"), Parameter_get_flt(par, "dilation.threshold"));
		
		message("Spatial dilation");
		DataCube_dilate_mask_xy(dataCube, maskCube, catalog, Parameter_get_int(par, "dilation.iterationsXY"), Parameter_get_flt(par, "dilation.threshold"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	// ---------------------------- //
	// Reload data cube if required //
	// ---------------------------- //
	
	if(use_noise || use_weights || use_noise_scaling)  // ALERT: Add conditions here as needed.
	{
		status("Reloading data cube for parameterisation");
		DataCube_readFITS(dataCube, Path_get(path_data_in), region);
		
		// Apply flags if required
		if(use_flagging) DataCube_flag_regions(dataCube, flag_regions);
		
		// Apply flagging catalogue if required		
		if(use_flagging_cat) DataCube_continuum_flagging(dataCube, Parameter_get_str(par, "flag.catalog"), 1, Parameter_get_int(par, "flag.radius"));
		
		// Invert cube if requested
		if(use_invert)
		{
			message("Inverting data cube");
			DataCube_multiply_const(dataCube, -1.0);
		}
		
		// Apply gain cube if provided
		if(use_gain)
		{
			status("Loading and applying gain cube");
			DataCube *gainCube = DataCube_new(verbosity);
			DataCube_readFITS(gainCube, Path_get(path_gain_in), region);
			
			// Divide by gain cube
			DataCube_divide(dataCube, gainCube);
			
			// Delete gain cube again
			DataCube_delete(gainCube);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Parameterise sources         //
	// ---------------------------- //
	
	if(use_parameteriser)
	{
		status("Measuring source parameters");
		DataCube_parameterise(dataCube, maskCube, catalog, use_wcs, use_physical, Parameter_get_str(par, "parameter.prefix"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create and save cubelets     //
	// ---------------------------- //
	
	if(write_cubelets)
	{
		status("Creating cubelets");
		DataCube_create_cubelets(dataCube, maskCube, catalog, Path_get(path_cubelets), overwrite, use_wcs, use_physical, Parameter_get_int(par, "output.marginCubelets"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create and save moment maps  //
	// ---------------------------- //
	
	if(write_moments)
	{
		status("Creating moment maps");
		
		// Generate moment maps
		DataCube *mom0 = NULL;
		DataCube *mom1 = NULL;
		DataCube *mom2 = NULL;
		DataCube *chan = NULL;
		DataCube_create_moments(dataCube, maskCube, &mom0, &mom1, &mom2, &chan, NULL, use_wcs, true);
		
		// Save moment maps to disk
		if(mom0 != NULL) DataCube_save(mom0, Path_get(path_mom0), overwrite, DESTROY);
		if(mom1 != NULL) DataCube_save(mom1, Path_get(path_mom1), overwrite, DESTROY);
		if(mom2 != NULL) DataCube_save(mom2, Path_get(path_mom2), overwrite, DESTROY);
		if(chan != NULL) DataCube_save(chan, Path_get(path_chan), overwrite, DESTROY);
		
		// Delete moment maps again
		DataCube_delete(mom0);
		DataCube_delete(mom1);
		DataCube_delete(mom2);
		DataCube_delete(chan);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Save mask cube               //
	// ---------------------------- //
	
	if(write_mask || write_mask2d)
	{
		status("Writing mask cube");
		
		// Create and save projected 2-D mask image
		if(write_mask2d)
		{
			DataCube *maskImage = DataCube_2d_mask(maskCube);
			DataCube_save(maskImage, Path_get(path_mask_2d), overwrite, DESTROY);
			DataCube_delete(maskImage);
		}
		
		// Write 3-D mask cube
		if(write_mask) DataCube_save(maskCube, Path_get(path_mask_out), overwrite, DESTROY);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Save catalogue(s)            //
	// ---------------------------- //
	
	if(write_ascii || write_xml || write_sql)
	{
		status("Writing source catalogue");
		
		// Correct x, y and z for subregion offset if requested
		// WARNING: This will alter the original x, y and z positions!
		if(use_region && use_pos_offset)
		{
			for(size_t i = Catalog_get_size(catalog); i--;)
			{
				Source *src = Catalog_get_source(catalog, i);
				Source_offset_xyz(src, Array_siz_get(region, 0), Array_siz_get(region, 2), Array_siz_get(region, 4));
			}
		}
		
		if(write_ascii)
		{
			message("Writing ASCII file:   %s", Path_get_file(path_cat_ascii));
			Catalog_save(catalog, Path_get(path_cat_ascii), CATALOG_FORMAT_ASCII, overwrite);
		}
		
		if(write_xml)
		{
			message("Writing VOTable file: %s", Path_get_file(path_cat_xml));
			Catalog_save(catalog, Path_get(path_cat_xml), CATALOG_FORMAT_XML, overwrite);
		}
		
		if(write_sql)
		{
			message("Writing SQL file:     %s", Path_get_file(path_cat_sql));
			Catalog_save(catalog, Path_get(path_cat_sql), CATALOG_FORMAT_SQL, overwrite);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Clean up and exit            //
	// ---------------------------- //
	
	// Delete data cube and mask cube
	DataCube_delete(maskCube);
	DataCube_delete(dataCube);
	
	// Delete sub-cube region
	Array_siz_delete(region);
	
	// Delete flagging regions
	Array_siz_delete(flag_regions);
	
	// Delete input parameters
	Parameter_delete(par);
	
	// Delete file paths
	Path_delete(path_data_in);
	Path_delete(path_gain_in);
	Path_delete(path_noise_in);
	Path_delete(path_weights_in);
	Path_delete(path_mask_in);
	Path_delete(path_cat_ascii);
	Path_delete(path_cat_xml);
	Path_delete(path_cat_sql);
	Path_delete(path_mask_out);
	Path_delete(path_mask_2d);
	Path_delete(path_mask_raw);
	Path_delete(path_noise_out);
	Path_delete(path_filtered);
	Path_delete(path_mom0);
	Path_delete(path_mom1);
	Path_delete(path_mom2);
	Path_delete(path_chan);
	Path_delete(path_rel_plot);
	Path_delete(path_skel_plot);
	Path_delete(path_flag);
	Path_delete(path_cubelets);
	
	// Delete source catalogue
	Catalog_delete(catalog);
	
	// Print status message
	status("Pipeline finished.");
	
//  return ERR_SUCCESS;
}


int main(int argc, char **argv)
{
	int n =1;
	ensure(argc == 2, ERR_USER_INPUT, "Unexpected number of command line arguments.\nUsage: %s <parameter_file>", argv[0]);
/*
	sofia_args *arg;
	arg->dataPtr = argv[1];
	arg->datasize = 0;
	arg->headerPtr = NULL;
	arg->headersize = 0;
	mainline(arg);
*/
	mainline(NULL,argv[1],1,1);
	return 0;
}
