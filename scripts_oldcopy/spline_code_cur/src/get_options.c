/***************************************************************************//**
 * @file
 * @brief Functions for reading modis_smooth_t program options
 * 
 * Functions in this file handle reading runtime options from the command line 
 * and a specially formatted configuration file.
 ******************************************************************************/
 
#include "modis_smooth_t.h"
#include <sys/types.h> // stat
#include <sys/stat.h> // stat()
#include <getopt.h> // getopt_long() and friends
#include <string.h> // strtok(), strdup()

// local marcos
#define OPT_HELP 0          /**<Option label used by getopt_long() */
#define OPT_VERBOSE 1       /**<Option label used by getopt_long() */
#define OPT_CONFIG 2        /**<Option label used by getopt_long() */
#define OPT_TILE 3          /**<Option label used by getopt_long() */
#define OPT_SMOOTH 4         /**<Option label used by getopt_long() */
#define OPT_RESIDUAL 5      /**<Option label used by getopt_long() */
#define OPT_CHUNK 6         /**<Option label used by getopt_long() */
#define OPT_START_YEAR 7    /**<Option label used by getopt_long() */
#define OPT_NUM_YRS 8      /**<Option label used by getopt_long() */
#define OPT_STRIDE 9       /**<Option label used by getopt_long() */
#define OPT_OUT_STEM 10      /**<Option label used by getopt_long() */
#define OPT_OUT_DIR 11       /**<Option label used by getopt_long() */
#define OPT_NBAR_DIR 12     /**<Option label used by getopt_long() */
#define OPT_QA_DIR 13       /**<Option label used by getopt_long() */
#define OPT_CONFIG_OUT 14   /**<Option label used by getopt_long() */
#define OPT_PREFILTER 15    /**<Option label used by getopt_long() */
#define OPT_QUANTILE 16    /**<Option label used by getopt_long() */
#define OPT_CLIM_DIR 17	    /**<Option label used by getopt_long() */	
#define OPT_EXTRACT_PATH 18	/**<Option label used by getopt_long() */	
 
// local prototypes
void check_command_line(int, char**, struct option*);
void get_options_default(config*);
void get_options_config_file(int, char**, struct option*, config*);
void get_options_command_line(int, char**, struct option*, config*);
void parse_options(int, char**, struct option*, config*);
void check_options(config);
void print_options(config*, FILE*, char*);
void print_usage();
void csv_to_flags(const char*, int*);
void set_read_flags(int*, int*, int*);
void set_constants(config*);

/***************************************************************************//**
 * @brief Read program options from defaults, config file, and command line   
 ******************************************************************************/
void get_options(int argc, char **argv, config *param) {
	
	// define options flags, last field contains integer labels
	struct option opt[] = {
	{"help",       no_argument,       0, OPT_HELP},
	{"verbose",    no_argument,       0, OPT_VERBOSE},
	{"config",     required_argument, 0, OPT_CONFIG},
	{"tile_id",    required_argument, 0, OPT_TILE},
	{"smooth",     required_argument, 0, OPT_SMOOTH},
	{"residual",   required_argument, 0, OPT_RESIDUAL},
	{"chunk",      required_argument, 0, OPT_CHUNK},
	{"start_year", required_argument, 0, OPT_START_YEAR},
	{"num_yrs",    required_argument, 0, OPT_NUM_YRS},
	{"stride",     required_argument, 0, OPT_STRIDE},
	{"out_stem",   required_argument, 0, OPT_OUT_STEM},
	{"out_dir",    required_argument, 0, OPT_OUT_DIR},
	{"nbar_dir",   required_argument, 0, OPT_NBAR_DIR},
	{"qa_dir",     required_argument, 0, OPT_QA_DIR},
	{"config_out", required_argument, 0, OPT_CONFIG_OUT},  
	{"prefilter",  required_argument, 0, OPT_PREFILTER},  
	{"quantile",  required_argument, 0, OPT_QUANTILE},
	{"clim_dir",  required_argument, 0, OPT_CLIM_DIR},
	{"extract_path",  required_argument, 0, OPT_EXTRACT_PATH},
	{0, 0, 0, 0}
	};

    // populate constants
    set_constants(param);

	// check for sane command line options
	check_command_line(argc, argv, opt);

	// start with default values
	get_options_default(param);

	// overwrite with config file values, if provided
	get_options_config_file(argc, argv, opt, param);
	
	// overwrite with command line values, if provided
	get_options_command_line(argc, argv, opt, param);

	// get derived parameters
    sscanf(param->tile_id, "h%dv%d", &param->tile_column, &param->tile_row);
    csv_to_flags(param->smooth_str, param->smooth);
    csv_to_flags(param->prefilter_str, param->prefilter);
    csv_to_flags(param->residual_str, param->residual);
    set_read_flags(param->prefilter, param->smooth, param->read);

    // sanity check
    check_options(*param); 
    
    // write current configuration
    if (strcmp(param->config_out, "NONE") != 0) {
		FILE *fp;
		assert( (fp = fopen(param->config_out, "w")) != NULL );
		print_options(param, fp, "");
		assert( fclose(fp) != EOF );
	}
   
	// verbose output
	if (param->verbose == 1) {
		printf("Options:\n");
		print_options(param, stdout, "  ");
		printf("\n");
	}

	return;
}


/***************************************************************************//**
 * @brief Set constant values in program options struct
 ******************************************************************************/
void set_constants(config *pp) {

    sprintf(pp->nbar_sds_name[0], "%s", SDS_NAME_NBAR_BAND1);
    sprintf(pp->nbar_sds_name[1], "%s", SDS_NAME_NBAR_BAND2);
    sprintf(pp->nbar_sds_name[2], "%s", SDS_NAME_NBAR_BAND3);
    sprintf(pp->nbar_sds_name[3], "%s", SDS_NAME_NBAR_BAND4);
    sprintf(pp->nbar_sds_name[4], "%s", SDS_NAME_NBAR_BAND5);
    sprintf(pp->nbar_sds_name[5], "%s", SDS_NAME_NBAR_BAND6);
    sprintf(pp->nbar_sds_name[6], "%s", SDS_NAME_NBAR_BAND7);

    sprintf(pp->qa_sds_name[0], "%s", SDS_NAME_QA_BAND1);
    sprintf(pp->qa_sds_name[1], "%s", SDS_NAME_QA_BAND2);
    sprintf(pp->qa_sds_name[2], "%s", SDS_NAME_QA_BAND3);
    sprintf(pp->qa_sds_name[3], "%s", SDS_NAME_QA_BAND4);
    sprintf(pp->qa_sds_name[4], "%s", SDS_NAME_QA_BAND5);
    sprintf(pp->qa_sds_name[5], "%s", SDS_NAME_QA_BAND6);
    sprintf(pp->qa_sds_name[6], "%s", SDS_NAME_QA_BAND7);

    return;
}

/***************************************************************************//**
 * @brief Check for parse-able command line parameters
 ******************************************************************************/
void check_command_line(int argc, char **argv, struct option *opt) {
	
	int this;
	
	optind = 1; // reset getopt	
	while ((this = getopt_long(argc, argv, "", opt, NULL)) != -1) {
		switch (this) {
			case OPT_HELP:
				print_usage();			
			case '?':
				print_usage();	
		}		
	}
	if (optind < argc) {
		printf("Some command line parameters could not be parsed, exiting\n\n");
		print_usage();
	}	
	
	return;
}


/***************************************************************************//**
 * @brief Set program configuration parameters to default values
 ******************************************************************************/
void get_options_default(config *pp) {
	
	pp->verbose =        DEFAULT_VERBOSE;
	pp->chunk_num_row =  DEFAULT_CHUNK;
	pp->start_year =     DEFAULT_START_YEAR;
	pp->num_years =      DEFAULT_NUM_YRS;
	pp->stride =         DEFAULT_STRIDE;
	pp->quantile =	     DEFAULT_QUANTILE;

	sprintf(pp->tile_id,       "%s", DEFAULT_TILE_ID);
   	sprintf(pp->smooth_str,    "%s", DEFAULT_SMOOTH_STR);
   	sprintf(pp->prefilter_str, "%s", DEFAULT_PREFILTER_STR);
    	sprintf(pp->residual_str,  "%s", DEFAULT_RESID_STR);
	sprintf(pp->out_stem,      "%s", DEFAULT_OUT_STEM);
	sprintf(pp->out_dir,       "%s", DEFAULT_OUT_DIR);
	sprintf(pp->nbar_dir,      "%s", DEFAULT_NBAR_DIR);
	sprintf(pp->qa_dir,        "%s", DEFAULT_QA_DIR);
	sprintf(pp->clim_dir,    "%s", DEFAULT_CLIM_DIR);
	sprintf(pp->extract_path,    "%s", DEFAULT_EXTRACT_PATH);
	sprintf(pp->config_out,    "%s", DEFAULT_CONFIG_OUT);
	
	return;
}


/***************************************************************************//**
 * @brief Find configuration file and read program options flags therein
 ******************************************************************************/
void get_options_config_file(int argc, char **argv, struct option *opt, config *pp) {

	int this;
	char line[LONG_STR], *token; 
	FILE *fp;
	
	char cfg_file[LONG_STR] = "none";
	int cfg_argc = 1; // skip 0th element, just like command line
	char *cfg_argv[CONFIG_MAX_NUM_OPTS];

	// search for configuration file option on command line	
	optind = 1; // reset getopt	
	while ((this = getopt_long(argc, argv, "", opt, NULL)) != -1) {
		if (this == OPT_CONFIG) sprintf(cfg_file, "%s", optarg);
	}
	
	if (strcmp(cfg_file, "none") != 0) {
				
		// get argument count and values from configuration file
		assert( (fp = fopen(cfg_file, "r")) != NULL );
		while (fgets(line, sizeof(line), fp) != NULL) {
			if ( (token = strtok(line, CONFIG_DELIM)) != NULL ) { // only 1st strtok call refs string
				cfg_argv[cfg_argc] = strdup(token);
				cfg_argc++;
			}				
			while ( (token = strtok(NULL, CONFIG_DELIM)) != NULL ) {
				cfg_argv[cfg_argc] = strdup(token);
				cfg_argc++;
			}
		}
		assert( fclose(fp) != EOF );
		
		// parse options from argument count and argument value
		parse_options(cfg_argc, cfg_argv, opt, pp);

		// free memory allocated with strdup()
		int ii;
		for (ii=0; ii<cfg_argc; ii++) {
			free(cfg_argv[ii]);
		}
	}
	
	return;

}


/***************************************************************************//**   
 * @brief Read program options flags from command line 
 ******************************************************************************/
void get_options_command_line(int argc, char **argv, struct option *opt, config *pp) {
	
	parse_options(argc, argv, opt, pp);	
	return;
}


/***************************************************************************//** 
 * @brief Parse program configuration options flags from argc and argv
 ******************************************************************************/
void parse_options(int argc, char **argv, struct option *opt, config *pp)
{
	int this;
	
	optind = 1; // reset getopt	
	while ((this = getopt_long(argc, argv, "", opt, NULL)) != -1) {
		switch (this) {
			case OPT_VERBOSE:
				pp->verbose = 1;
				break;
			case OPT_TILE:	
				sprintf(pp->tile_id, "%s", optarg);
				break;
			case OPT_SMOOTH:
				sprintf(pp->smooth_str, "%s", optarg);
				break;
			case OPT_RESIDUAL:
				sprintf(pp->residual_str, "%s", optarg);
				break;
			case OPT_PREFILTER:
                		sprintf(pp->prefilter_str, "%s", optarg);
                		break;
			case OPT_CHUNK:
				pp->chunk_num_row = atoi(optarg);
				break;
			case OPT_START_YEAR:
				pp->start_year = atoi(optarg);
				break;
			case OPT_NUM_YRS:
				pp->num_years = atoi(optarg);
				break;
            		case OPT_STRIDE:
               			pp->stride = atoi(optarg);
                		break;
			case OPT_QUANTILE:
               			pp->quantile = atoi(optarg);
                		break;
			case OPT_OUT_STEM:
				sprintf(pp->out_stem, "%s", optarg);
				break;
			case OPT_OUT_DIR:
				sprintf(pp->out_dir, "%s", optarg);
				break;
			case OPT_NBAR_DIR:
				sprintf(pp->nbar_dir, "%s", optarg);
				break;
			case OPT_QA_DIR:
				sprintf(pp->qa_dir, "%s", optarg);
				break;
			case OPT_CLIM_DIR:
				sprintf(pp->clim_dir, "%s", optarg);
				break;
			case OPT_EXTRACT_PATH:
				sprintf(pp->extract_path, "%s", optarg);
				break;
			case OPT_CONFIG_OUT:
				sprintf(pp->config_out, "%s", optarg);
				break;

			case '?':
				print_usage();	
		}		
	}
	
	return;
}


/***************************************************************************//**
 * @brief Check for sane parameter values, exit with friendly message upon
 * failure
 ******************************************************************************/
void check_options(config pp) {
	
	struct stat info;
    int ii, num_layers;
	
    // valid tile row and column?
    if (pp.tile_column < MODIS_TILE_HMIN || pp.tile_column > MODIS_TILE_HMAX || \
            pp.tile_row < MODIS_TILE_VMIN || pp.tile_row > MODIS_TILE_VMAX) {
        printf("check_options: invalid tile row and column\n");
        exit(EXIT_FAILURE);
    }

    // at least one data layer is selected?
    num_layers = 0;
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        if (pp.smooth[ii] == 1) num_layers++;
    }
    if (num_layers == 0) {
        printf("check_options: 0 data layers selected\n");
        exit(EXIT_FAILURE);
    }

    // only request residuals for layers that are smoothed?
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        if (pp.residual[ii] == 1 && pp.smooth[ii] == 0) {
            printf("check options: cannot compute residuals for layer %d unless "
                   "it is selected for smoothing\n", ii+1);
            exit(EXIT_FAILURE);
        }
    }

    // only prefilter layers that are smoothed?
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        if (pp.prefilter[ii] == 1 && pp.smooth[ii] == 0) {
            printf("check options: cannot prefilter layer %d unless "
                   "it is selected for smoothing\n", ii+1);
            exit(EXIT_FAILURE);
        }
    }

    // prefiltering is only implemented for the EVI2 index
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        if ((pp.prefilter[ii] == 1) && (ii != EVI2)) {
            printf("check_options: prefiltering is only implemented for the EVI2 "
                   "index (layer %d)\n", EVI2+1);
            exit(EXIT_FAILURE);
        }
    }

    // valid chunk size?
    if (pp.chunk_num_row < 1 || pp.chunk_num_row > MODIS_NROW) {
        printf("check_options: invalid chunk size\n");
        exit(EXIT_FAILURE);
    }

    // valid directories?
    if (stat(pp.nbar_dir, &info) != 0) {
        printf("check_options: invalid NBAR input directory\n");
        exit(EXIT_FAILURE);
    }

    if (stat(pp.qa_dir, &info) != 0) {
        printf("check_options: invalid QA input directory\n");
        exit(EXIT_FAILURE);
    }

    if (stat(pp.out_dir, &info) != 0) {
        printf("check_options: invalid output directory\n");
        exit(EXIT_FAILURE);
    }

    if (stat(pp.clim_dir, &info) != 0) {
        printf("check_options: invalid clim directory\n");
        exit(EXIT_FAILURE);
    }

     // ## if the extract file is specified then check if the file exists. If not, throw error.
     if (strcmp(pp.extract_path, "NONE") != 0) {
	if (stat(pp.extract_path, &info) != 0) {
        	printf("check_options: invalid extract file given.\n");
        	exit(EXIT_FAILURE);
    	}
    }

    // valid flags?
    if (pp.verbose != 1 && pp.verbose != 0) {
        printf("check_options: verbose flag must be 1 or 0\n");
        exit(EXIT_FAILURE);
    }

    // valid number of years?
    if (pp.num_years < 3) {
        printf("check_options: invalid num_years, minimum is 3\n");
        exit(EXIT_FAILURE);
    }
    
    // valid stride?
    if (pp.stride < 1 || pp.stride > MODIS_PERIOD) {
        printf("check_options: invalid value for stride\n");
        exit(EXIT_FAILURE);
    }

    // valid quantile?
    if (pp.quantile != 1 && pp.quantile != 0) {
        printf("check_options: invalid value for quantile\n");
        exit(EXIT_FAILURE);
    }

    // valid quantile inputs? 
    // the program needs bands 1,2,4,5,6,7,evi
    if (pp.quantile == 1) {
        if ((pp.smooth[RED]   != 1) || (pp.smooth[NIR]   != 1) || 
            (pp.smooth[GREEN] != 1) || (pp.smooth[SWIR1] != 1) || 
            (pp.smooth[SWIR2] != 1) || (pp.smooth[SWIR3] != 1) || 
            (pp.smooth[EVI2]  != 1)) {

            printf("check_options: missing smooth layers required for quantile\n");
            printf("\trequires RED, NIR, GREEN, SWIR1, SWIR2, SWIR3, EVI2 layers "
                    "(%d,%d,%d,%d,%d,%d,%d)\n", RED+1, NIR+1, GREEN+1,
                    SWIR1+1, SWIR2+1, SWIR3+1, EVI2+1);
            exit(EXIT_FAILURE);
        }
    }

    // prefilter assumes NDSI is the lasy layer
    if (NDSI != (NUM_DATA_LAYERS-1)) {
        printf("check_options: bad macro values, NDSI must be the last data layer\n");
        exit(EXIT_FAILURE);
    }

	return; 
}

/***************************************************************************//**
 * @brief Parse layer and residual options, converting CSV string to array of
 * logical flags
 *
 * Used for both the layer and residual arrays. The latter may be set to 0 to
 * skip all data layers, so this behavior is included here.
 ******************************************************************************/
void csv_to_flags(const char *csv, int *flags) {

    char *token, str[SHORT_STR];
    int ii, layer_index; 

    // initialize layer flags to 0
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        flags[ii] = 0;
    }

    if (strcmp(csv, "0") == 0) {
        // keep all flags as 0 and exit
        return;
    }

    // copy string literals to var b/c strtok() modifies the string
    sprintf(str, "%s", csv); 

    // read first layer: only 1st strtok call refs string, must succeed once
    assert((token = strtok(str, ",")) != NULL); 
    layer_index = atoi(token)-1;
    assert(layer_index >= 0 && layer_index < NUM_DATA_LAYERS);
    flags[layer_index] = 1;

    // read remaining layers
    while ((token = strtok(NULL, ",")) != NULL) {
        layer_index = atoi(token)-1;
        assert(layer_index >= 0 && layer_index < NUM_DATA_LAYERS);
        flags[layer_index] = 1;
    }

    return;
}


/***************************************************************************//**
 * @brief Determine which data layers need to be read / computed
 * @param [in] prefilter 1D array of logical flags indicating if layers should be prefiltered
 * @param [in] smooth 1D array of logical flags for each layer
 * @param [out] read 1D array of logical flags for each layer
 *
 * Read flags are computed such that all the dependencies are met for
 * prefiltering and index computation
 ******************************************************************************/
void set_read_flags(int *prefilter, int *smooth, int *read) {

    int ii, any_prefilter;

    // all layers to be smoothed must be read
    for (ii = 0; ii <  NUM_DATA_LAYERS; ii++) {
        read[ii] = smooth[ii];
    }

    // prefiltering requires NDSI index
    any_prefilter = 0;
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
        if (prefilter[ii] == 1) {
            any_prefilter = 1;
            break;
        }
    }
    if (any_prefilter == 1) {
        read[NDSI] = 1;
    }

    // NDSI index requires green and swir2 bands
    if (read[NDSI] == 1) {
        read[GREEN] = 1;
        read[SWIR2] = 1;
    }
    
    // EVI2 index requires red and nir bands
    if (read[EVI2] == 1) {
        read[RED] = 1;
        read[NIR] = 1;
    }

    return;
}


/***************************************************************************//**
 * @brief Print formatted values of program options
 ******************************************************************************/
void print_options(config *pp, FILE *fp, char *prefix) {
	
	fprintf(fp, "%s--tile_id %s\n",     prefix, pp->tile_id);
	fprintf(fp, "%s--smooth %s\n",      prefix, pp->smooth_str);
	fprintf(fp, "%s--prefilter %s\n",   prefix, pp->prefilter_str);
	fprintf(fp, "%s--residual %s\n",    prefix, pp->residual_str);
	fprintf(fp, "%s--quantile %d\n",      prefix, pp->quantile);
	fprintf(fp, "%s--chunk %zu\n",      prefix, pp->chunk_num_row);
	fprintf(fp, "%s--start_year %d\n",  prefix, pp->start_year);
	fprintf(fp, "%s--num_yrs %d\n",     prefix, pp->num_years);
	fprintf(fp, "%s--stride %d\n",      prefix, pp->stride);
	fprintf(fp, "%s--out_stem %s\n",    prefix, pp->out_stem);
	fprintf(fp, "%s--out_dir %s\n",     prefix, pp->out_dir);
	fprintf(fp, "%s--nbar_dir %s\n",    prefix, pp->nbar_dir);
	fprintf(fp, "%s--qa_dir %s\n",      prefix, pp->qa_dir);
	fprintf(fp, "%s--clim_dir %s\n",      prefix, pp->clim_dir);
	fprintf(fp, "%s--extract_path %s\n",      prefix, pp->extract_path);
	fprintf(fp, "%s--config_out %s\n",  prefix, pp->config_out);
    	
	return;
} 


/***************************************************************************//**
 * @brief Print command-line usage help message and exit program
 ******************************************************************************/
void print_usage() 
{
	printf("\nUsage: modis_spline_t [OPTIONS]\n\n");
	printf("  --help\t\tPrint this help message\n\n");
	printf("  --config\t\tPath to configuration file\n\n");
	printf("  --tile_id ID\t\tTile id string, default = " DEFAULT_TILE_ID "\n\n");
	printf("  --smooth CSV\t\tData layers to smooth, default = %s\n", DEFAULT_SMOOTH_STR);
    printf("  \t\t\t\tred=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9\n\n");
	printf("  --prefilter CSV\tData layers to prefilter using NDSI, default = %s\n", DEFAULT_PREFILTER_STR);
    printf("  \t\t\t\tred=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9\n\n");
	printf("  --residual CSV\tData layers to save smooth-original residuals, default = %s\n", DEFAULT_RESID_STR);
    printf("  \t\t\t\tred=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9\n\n");
	printf("  --quantile \tFlag whether to include quantile outputs (0=off and 1=on), default = %d\n\n", DEFAULT_QUANTILE);
	printf("  --chunk NUM\t\tNumber of rows read and processed at once, default = %d\n\n", DEFAULT_CHUNK);
	printf("  --start_year YEAR\tFirst year to be processed, default = %d\n\n", DEFAULT_START_YEAR);
	printf("  --num_yrs NUM\t\tNumber of years to be processed, default = %d\n\n", DEFAULT_NUM_YRS);
	printf("  --stride NUM\t\tInterval in days between MODIS input files, others are ignored, default = %d\n\n", DEFAULT_STRIDE);
	printf("  --out_stem STR\tOutput filename stem (prefix), default = " DEFAULT_OUT_STEM "\n\n");
	printf("  --out_dir PATH\tDirectory to write output files, must exist, default = " DEFAULT_OUT_DIR "\n\n");
	printf("  --nbar_dir PATH\tDirectory containing NBAR data, default = " DEFAULT_NBAR_DIR "\n\n");
	printf("  --qa_dir PATH\tDirectory containing QA data, default = " DEFAULT_QA_DIR "\n\n");
	printf("  --clim_dir PATH\tDirectory containing WORLDCLIM TMIN data, default = %s\n\n",DEFAULT_CLIM_DIR);
	printf("  --extract_path PATH\tPath to file containing csv with pixels to be extracted. Should have four columns in order of tile_x, tile_y, pix_x, pix_y. Default = %s\n\n",DEFAULT_EXTRACT_PATH);
	printf("  --config_out PATH\tOutput file name for current configuration, default = " DEFAULT_CONFIG_OUT "\n\n");
	
	exit(EXIT_FAILURE);
	
	return;
}
