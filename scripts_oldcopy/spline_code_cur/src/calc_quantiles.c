/***************************************************************************//**
 * @file
 * @brief Functions for computing quantiles on bands and indices
 ******************************************************************************/

#include "modis_smooth_t.h"
#include <math.h> // sqrt

// local macros
#define QUANTILES 0.1, 0.33, 0.5, 0.66, 0.9 /**< Quantiles to compute for each band, comma-separated values */

// ## local definitions
static void get_quants(config,int16*, int16*, int16*, int, int);
static int compare (const void*, const void*);

/***************************************************************************//**
 * @brief Compute all the quantiles for all the bands for each year and write out to disk
 *
 * @param [in] pp Current config parameters
 * @param [in] 2D data_chunk Will contain all the smoothed bands and EVI
 * @param [in] 2D index_chunk Will contain the indices from smoothed input bands
 * @param [in] 1D flag_chunk Will contain the snow flags from EVI2 prefilter
 * @param [out] quant_files File pointers to the outputs for each year
 ******************************************************************************/
void calc_all_quantiles(config pp, int16** data_chunk, int16** index_chunk, int16* flag_chunk, FILE** quant_files) {

    int16 *quant;	
    int year, ii, jj;
    size_t pixel, idx, total_quants;

    // ## we want NUM_QUANTS outputs for each band and one snow flag count across all bands
    total_quants = NUM_QUANT_BANDS*NUM_QUANTS + 1;
    assert((quant  = malloc(pp.num_elem_step*total_quants*sizeof(int16))) != NULL);

    //## loop through each band for each year to create the outputs for that year
    for(year = 1; year < (pp.num_years-1); year++) {
	// ## jj keeps track of the current output band
        jj = 0;
        // ## loop through splined band data
        for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
            // ## dont do blue band or NDSI
            if(ii == 2 || ii == 8) continue;
	    // ## call generic parallelized get_quants function and increment jj
            get_quants(pp, data_chunk[ii], flag_chunk, quant, year, jj);
            jj++;
        }
        // ## loop through spectral indices - call same function
        for (ii = 0; ii < NUM_SPL_INDICES; ii++) {
            get_quants(pp, index_chunk[ii], flag_chunk, quant, year, jj);
            jj++;
        }
	// ## at this point jj should equal NUM_QUANT_BANDS
        // ## write data for that year directly to disk
        for (pixel = 0; pixel < pp.num_elem_step; pixel++) {
            idx = pixel*total_quants;
            assert( fwrite((&quant[idx]), sizeof(int16), total_quants, quant_files[year-1]) == total_quants);
        }
    }

    free(quant);
}

/***************************************************************************//**
 * @brief Get quantiles across time for each pixel from a given splined output
 *
 * @param [in] pp Current config parameters
 * @param [in] input Input band for current chunk
 * @param [in] flag Input snow flag for current chunk
 * @param [out] output Output quantiles - will contain all quantiles for all bands
 * @param [in] cur_year Current year that is being processed
 * @param [in] cur_band Band location for quantiles to be saved into output array
 ******************************************************************************/
static void get_quants(config pp, int16* input, int16* flag, int16* output, int cur_year, int cur_band) {

    	size_t ii, nn;
	int16 *quant_pix;
	size_t in_loc, new_out_loc, snow_count,good_count;
	float cur_quant[NUM_QUANTS-1] = { QUANTILES }; // -1 because last "quantile" is standard deviation
	double sum, mean, var;
	short std;

	// ## use parallel processing to do efficient sorting
	#pragma omp parallel default(none) \
        shared(pp, input, flag, output, cur_year, cur_band, cur_quant) \
        private(in_loc,new_out_loc,quant_pix,sum,mean,var,std,ii,nn,good_count,snow_count)
	{
	// ## temporary array to hold the band values for one year
	quant_pix = malloc(sizeof(int16)*MODIS_PERIOD);	

	// ## loop through each pixel
	 #pragma omp for schedule(static)
	for( ii = 0; ii < pp.num_elem_step; ii++ ) {

		// ## this is the location within the input where the current year's values start for that pixel
		in_loc = ii*MODIS_PERIOD*pp.num_years + MODIS_PERIOD*cur_year;
		// ## output location needs to account for all bands and all pixels
		// ## organized for each pixel should contain NUM_QUANTS repeating metrics for all bands plus one snow count value
		new_out_loc = ii*((NUM_QUANTS*NUM_QUANT_BANDS)+1) + cur_band*NUM_QUANTS;
		
		// ## fill the quant_pix array and also compute the mean - assume no values are fill
		sum = 0;
		good_count = 0;
		for(nn = 0; nn < MODIS_PERIOD; nn++) {
			// ## only use non-snow values in the mean and var
			if(flag[in_loc + nn] < KIND_SNOW_OBS)
			{
				// ## only fill the array with non-snow values so we use count to index them				
				quant_pix[good_count] = input[in_loc + nn];
				sum += quant_pix[good_count];
				// ## the count will provide the number of snow-free values in that year
				good_count++;
			}
		}
		// ## only compute the metrics if there are more than 10 useful (non-snow) values to composite
		if(good_count > 10)
		{
			// ## calculate mean and use it to calculate standard deviation			
			mean = (double)sum/good_count;
			// ## loop through data again to compute the variance/standard deviation
			var = 0;
			good_count = 0;
			for(nn = 0; nn < MODIS_PERIOD; nn++) {
				// ## only use non-snow values in the mean and var
				if(flag[in_loc + nn] < KIND_SNOW_OBS)
				{
					// ## increment the sum of squared differences					
					var += ((quant_pix[good_count]-mean)*(quant_pix[good_count]-mean));
					good_count++;
				}
			}
			// ## get the standard deviation for that year 
			std = (int16)(sqrt(var/good_count));
			// ## sort the quant_pix array (from tuples 0 to count) from lowest to highest
			qsort(quant_pix, good_count, sizeof(int16), compare);
			// ## write outputs to output array
			for(nn = 0; nn < (NUM_QUANTS-1); nn++) {
   				output[new_out_loc+nn] = quant_pix[(size_t)((good_count-1) * cur_quant[nn]) ];
			}
			output[new_out_loc+(NUM_QUANTS-1)] = std;
		}
		// ## if there are no valid values we use snow contaminated values instead
		else
		{			
			// ## fill the quant_pix array and also compute the mean - assume no values are fill
			sum = 0;
			snow_count = 0;
			for(nn = 0; nn < MODIS_PERIOD; nn++) {
				// ## only use snow values in the mean and var
				if(flag[in_loc + nn] > KIND_NSNOW_OBS_MIN && flag[in_loc + nn] < 6)
				{
					// ## only fill the array with snow values so we use count to index them				
					quant_pix[snow_count] = input[in_loc + nn];
					sum += quant_pix[snow_count];
					// ## the count will provide the number of snow values in that year
					snow_count++;
				}
			}
			// ## for snow contaminated versions - use a lower threshold
			if(snow_count > 0)
			{
				// ## calculate mean and use it to calculate standard deviation			
				mean = (double)sum/snow_count;
				// ## loop through data again to compute the variance/standard deviation
				var = 0;
				snow_count = 0;
				for(nn = 0; nn < MODIS_PERIOD; nn++) {
					// ## only use snow values in the mean and var
					if(flag[in_loc + nn] > KIND_NSNOW_OBS_MIN && flag[in_loc + nn] < 6)
					{
						// ## increment the sum of squared differences					
						var += ((quant_pix[snow_count]-mean)*(quant_pix[snow_count]-mean));
						snow_count++;
					}
				}
				// ## get the standard deviation for that year 
				std = (int16)(sqrt(var/snow_count));
				// ## sort the quant_pix array (from tuples 0 to count) from lowest to highest
				qsort(quant_pix, snow_count, sizeof(int16), compare);
				// ## write outputs to output array
				for(nn = 0; nn < (NUM_QUANTS-1); nn++) {
	   				output[new_out_loc+nn] = quant_pix[(size_t)((snow_count-1) * cur_quant[nn]) ];
				}
				output[new_out_loc+(NUM_QUANTS-1)] = std;
			}
			// ## if we still have no good values we fill the whole array with fill
			else
			{
				for(nn = 0; nn < NUM_QUANTS; nn++) {
   					output[new_out_loc+nn] = FILL_NBAR;
				}
			}
		}
		// ## only do this once for the last band but we include the number of snow_flags as modis_period - good_count	 
		if(cur_band == (NUM_QUANT_BANDS-1))
		{	
			output[new_out_loc+NUM_QUANTS] = MODIS_PERIOD - good_count;
		}
	}
	free(quant_pix);

	} // ## end omp parallel

}

/***************************************************************************//**
 * @brief Compare two numbers for QSORT
 * @param [in] a First number to  be compared
 * @param [in] b Second number to be compared
 ******************************************************************************/
static int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
