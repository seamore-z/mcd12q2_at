/***************************************************************************//**
 * @file
 * @brief Prototypes, macros, etc. for the MODIS_SMOOTH_T program
 ******************************************************************************/

#ifndef MODIS_SMOOTH_T_INCLUDED
#define MODIS_SMOOTH_T_INCLUDED

// common headers
#include <stdio.h>  // printf(), fopen(),  etc
#include <stdlib.h> // size_t, EXIT_SUCCESS, EXIT_FAILURE
#include <assert.h> // assert()
#include <omp.h>    // omp_get_num_threads, etc
#include "mfhdf.h"  // HDF API, int16 etc, TRUE etc

// macros
#define SHORT_STR 100 /**<Array size for "long" strings, arbitrary  */
#define LONG_STR 500  /**<Array size for "short" strings, arbitrary */

#define DEFAULT_VERBOSE 0			           /**<Default --verbose option    */
#define DEFAULT_TILE_ID "h10v02"	           /**<Default --tile option       */
#define DEFAULT_SMOOTH_STR "1,2,3,4,5,6,7,8,9" /**<Default --band option       */
#define DEFAULT_RESID_STR "0"                  /**<Default --residual option   */
#define DEFAULT_PREFILTER_STR "0"              /**<Default --prefilter option  */
#define DEFAULT_CHUNK 100			           /**<Default --chunk option      */
#define DEFAULT_START_YEAR 2001		           /**<Default --start_year option */
#define DEFAULT_NUM_YRS 3			           /**<Default --num_yrs option    */
#define DEFAULT_OUT_STEM "out"		           /**<Default --out_stem option   */
#define DEFAULT_OUT_DIR "./out"		           /**<Default --out_dir option    */
#define DEFAULT_NBAR_DIR "./nbar"              /**<Default --nbar_dir option   */
#define DEFAULT_QA_DIR "./qa"                  /**<Default --qa_dir option     */
#define DEFAULT_CONFIG_OUT "NONE"	           /**<Default --config_out option */
#define DEFAULT_STRIDE 8                       /**<Default --stride option     */
#define DEFAULT_QUANTILE 0                       /**<Default --quantile option     */
#define DEFAULT_CLIM_DIR "/projectnb/modislc/data/climate/worldclim/worldclim_data/tmin"   /**<Default --clim_dir  */
#define DEFAULT_EXTRACT_PATH "NONE"	           /**<Default --extract_path option */

#define CONFIG_MAX_NUM_OPTS 50	/**<Maximum number of delimited tokens in config file */
#define CONFIG_DELIM " \n" 		/**<Token delimiters in config file                   */

#define MODIS_TILE_HMIN 0	/**<Minimum column number in MODIS tile grid */
#define MODIS_TILE_HMAX 35	/**<Maximum column number in MODIS tile grid */
#define MODIS_TILE_VMIN 0	/**<Minimum row number in MODIS tile grid    */
#define MODIS_TILE_VMAX 17	/**<Maximum row number in MODIS tile grid    */
#define MODIS_NROW 2400		/**<Rows in a MODIS tile                     */
#define MODIS_NCOL 2400		/**<Columns in a MODIS tile                  */
#define MODIS_PERIOD 365	/**<Number of observations in MODIS year     */
#define MODIS_NUM_BANDS 7   /**<Number of bands in MODIS dataset         */

#define NUM_DATA_LAYERS 9 /**<Number of data layers used in this analysis, includes MODIS bands, indices, and other data */
#define RED 0             /**<Index of MODIS red band in data arrays */
#define NIR 1             /**<Index of MODIS nir band in data arrays */
#define BLUE 2            /**<Index of MODIS blue band in data arrays */
#define GREEN 3           /**<Index of MODIS green band in data arrays */
#define SWIR1 4           /**<Index of MODIS swir1 band in data arrays */
#define SWIR2 5           /**<Index of MODIS swir2 band in data arrays */
#define SWIR3 6           /**<Index of MODIS swir3 band in data arrays */
#define EVI2 7            /**<Index of MODIS evi2 index in data arrays */
#define NDSI 8            /**<Index of MODIS ndsi index in data arrays */

#define NUM_SPL_INDICES 4   /**<Number of extra indices to be computed from splined data        */
#define NUM_QUANT_BANDS 11   /**<Number of output features to have a quantile calculated       */
#define NUM_QUANTS 6   /**<Number of quantiles to be created for each band and index        */

#ifdef C5 // temporary: FILL_QA is set to 3 for C6 because preliminary C6 dataset is missing many QA files
# define FILL_QA 15     /**<Fill value for the QA flags */
#else
# define FILL_QA 4      /**<Fill value for the QA flags */
#endif
#define FILL_NBAR 32767 /**<Fill value for the spectral bands */

# define SDS_NAME_NBAR_BAND1 "Nadir_Reflectance_Band1"      /**<HDF SDS name for NBAR band 1, C5 & C6 */
# define SDS_NAME_NBAR_BAND2 "Nadir_Reflectance_Band2"      /**<HDF SDS name for NBAR band 2, C5 & C6 */
# define SDS_NAME_NBAR_BAND3 "Nadir_Reflectance_Band3"      /**<HDF SDS name for NBAR band 3, C5 & C6 */
# define SDS_NAME_NBAR_BAND4 "Nadir_Reflectance_Band4"      /**<HDF SDS name for NBAR band 4, C5 & C6 */
# define SDS_NAME_NBAR_BAND5 "Nadir_Reflectance_Band5"      /**<HDF SDS name for NBAR band 5, C5 & C6 */
# define SDS_NAME_NBAR_BAND6 "Nadir_Reflectance_Band6"      /**<HDF SDS name for NBAR band 6, C5 & C6 */
# define SDS_NAME_NBAR_BAND7 "Nadir_Reflectance_Band7"      /**<HDF SDS name for NBAR band 7, C5 & C6 */

#ifdef C5
# define SDS_NAME_QA_BAND1 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 1, C5 */
# define SDS_NAME_QA_BAND2 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 2, C5 */
# define SDS_NAME_QA_BAND3 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 3, C5 */
# define SDS_NAME_QA_BAND4 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 4, C5 */
# define SDS_NAME_QA_BAND5 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 5, C5 */
# define SDS_NAME_QA_BAND6 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 6, C5 */
# define SDS_NAME_QA_BAND7 "BRDF_Albedo_Band_Quality" /**<HDF SDS name for QA band 7, C5 */
#else
# define SDS_NAME_QA_BAND1 "BRDF_Albedo_Band_Quality_Band1" /**<HDF SDS name for QA band 1, C6 */
# define SDS_NAME_QA_BAND2 "BRDF_Albedo_Band_Quality_Band2" /**<HDF SDS name for QA band 2, C6 */
# define SDS_NAME_QA_BAND3 "BRDF_Albedo_Band_Quality_Band3" /**<HDF SDS name for QA band 3, C6 */
# define SDS_NAME_QA_BAND4 "BRDF_Albedo_Band_Quality_Band4" /**<HDF SDS name for QA band 4, C6 */
# define SDS_NAME_QA_BAND5 "BRDF_Albedo_Band_Quality_Band5" /**<HDF SDS name for QA band 5, C6 */
# define SDS_NAME_QA_BAND6 "BRDF_Albedo_Band_Quality_Band6" /**<HDF SDS name for QA band 6, C6 */
# define SDS_NAME_QA_BAND7 "BRDF_Albedo_Band_Quality_Band7" /**<HDF SDS name for QA band 7, C6 */
#endif

#define INPUT_STEM_NBAR "MCD43A4" /**<First token in NBAR data file names */
#define INPUT_STEM_QA "MCD43A2"   /**<First token in QA data file names   */

//  ## for output kind flags - used in prefilter and calc_quantile functions
#define EVI2_FLOOR_DEFAULT 1500  /**< Default value used for EVI2 quantile if too few points to compute */
#define KIND_NSNOW_OBS 0     /**< Flag value for non-snow, observed pixels */
#define KIND_NSNOW_INTERP 1  /**< Flag value for non-snow, interpolated pixels */
#define KIND_NSNOW_OBS_MIN 2 /**< Flag value for non-snow, observed, pixels truncated to the mininum value */
#define KIND_SNOW_OBS 3     /**< Flag value for snow, observed pixels */
#define KIND_SNOW_INTERP 4   /**< Flag value for snow, interpolated pixels */
#define KIND_SNOW_INTERP_OFF 5 /**< Flag value for non-snow, observed, pixels truncated to the mininum value */

// utility functions
static inline int max ( int a, int b ) { return a > b ? a : b; }
static inline int min ( int a, int b ) { return a < b ? a : b; }

// datatypes
/** @brief Program configuration variables */
typedef struct config 
{
    char    tile_id[SHORT_STR];	                      /**<Tile identification string, default = DEFAULT_TILE_ID          */
    int     tile_column; 		                      /**<Tile position, derived from tile_id                    */
    int     tile_row; 			                      /**<Tile position, derived from tile_id          */
    char    nbar_dir[LONG_STR];	                      /**<Directory containing MODIS data                      */
    char    qa_dir[LONG_STR];	                      /**<Directory containing MODIS data                            */  
    char    smooth_str[SHORT_STR];                    /**<Data layers to process, CSV, red=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9      */
    char    prefilter_str[SHORT_STR];                 /**<Data layers to prefilter using NDSI, CSV, red=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9      */
    char    residual_str[SHORT_STR];                  /**<Residuals to compute, CSV, red=1, nir=2, blue=3, green=4, swir1=5, swir2=6, swir3=7, evi2=8, ndsi=9      */
    int     read[NUM_DATA_LAYERS];                    /**<Data layers to read or compute, logical flags for red=[0], nir=[1], blue=[2], green=[3], swir1=[4], swir2=[5], swir3=[6], evi2=[7], ndsi=[8] */
    int     prefilter[NUM_DATA_LAYERS];           /**<Data layers to prefilter using NDSI values, logical flags for red=[0], nir=[1], blue=[2], green=[3], swir1=[4], swir2=[5], swir3=[6], evi2=[7], ndsi=[8] */
    int     smooth[NUM_DATA_LAYERS];                  /**<Data layers to smooth, logical flags for red=[0], nir=[1], blue=[2], green=[3], swir1=[4], swir2=[5], swir3=[6], evi2=[7], ndsi=[8] */
    int     residual[NUM_DATA_LAYERS];                /**<Residuals to compute, logical flags for red=[0], nir=[1], blue=[2], green=[3], swir1=[4], swir2=[5], swir3=[6], evi2=[7], ndsi=[8] */
    char    nbar_sds_name[MODIS_NUM_BANDS][LONG_STR]; /**<Subdataset name for NBAR band in the input HDF files                        */
    char    qa_sds_name[MODIS_NUM_BANDS][LONG_STR];   /**<Subdataset name for QA band in the input HDF files                            */
    char    out_dir[LONG_STR];	                      /**<Directory to write output files, must exist, default = DEFAULT_OUT_DIR           */
    char    out_stem[SHORT_STR];                      /**<Directory where WORLCLIM TMIN climatology exists filename stem/prefix   */ 
    char    clim_dir[LONG_STR];	                      /**<Directory to read climate files for prefiltering step, must exist, default = DEFAULT_CLIM_DIR  */     
    char    extract_path[LONG_STR];	                      /**<Full path to csv file with pixels to extract into text file, must exist, default = DEFAULT_EXTRACT_PATH  */      
    char    config_out[LONG_STR];                     /**<Path where current configuration should be saved, default = DEFAULT_CONFIG_OUT       */   
    int     stride;                                   /**<Interval in days between MODIS files to be read, skipped files use fill values        */
    int     verbose;			                      /**<Enable verbose output messages (1) or do not (0)                */  
    int     start_year;			                      /**<First year to be processed                                       */
    int     num_years;			                      /**<Number of years to be processed                                             */
    int     quantile;			                      /**<Flag for whether quantiles of bands will be an output                     */
    size_t  num_step;                                 /**<Number of timesteps in chunk */
    size_t  chunk_num_row;   	                      /**<Number of rows read and processed at once, default = DEFAULT_CHUNK               */   
    size_t  num_elem_step;                            /**<Elements in one timestep of chunk */
    size_t  num_elem_chunk;                           /**<Elements in all timesteps of chunk */
} config;

/** @brief HDF file stack information */
typedef struct {
	int num_files; /**<Number of file locations stored in this struct             */
	char **names;  /**<Array of file names (paths)                                */
	int *is_fill;  /**<Flag indicating if an expected file is fill (1) or not (0) */
} data_loc;

// shared prototypes
void get_options(int, char**, config*);                                    // get_options.c

void get_input_files(config, data_loc*, data_loc*);                        // get_input_files.c

FILE*** get_output_files(config, int);                                     // write_output.c

void close_loc(data_loc*);                                                 // write_output.c

void write_chunk(config, size_t, int16**, FILE***, int);                   // write_output.c

void close_output_files(config, int, FILE***);                             // write_output.c

void get_chunk(config, data_loc, data_loc, int, int16**, uint8**); // get_input_data.c

void apply_spline(config, int16**, uint8**, int16**);       // apply_spline.c

int sbart(double*, double*, double*, int, double*, int, double*, double*,
        double, int, double*, int*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int);        // sbart.c

float wallclock_timer(int);                                                // timer.c

void get_indices(config, int16**, uint8**); // get_index.c

void get_indices_spl(config, int16**, int16**); // get_index.c

void prefilter(config, int16**, uint8**,  int16**, int);                  // prefilter.c

void calc_all_quantiles(config, int16**,int16**, int16*, FILE**);	// calc_quantiles.c

FILE** get_output_files_quant(config);	// write_output.c
void close_output_files_quant(config, FILE**);	// write_output.c

#endif // MODIS_SMOOTH_T_INCLUDED
