/***************************************************************************//**
 * @file
 * @brief Functions to prefilter input data
 ******************************************************************************/

#include "modis_smooth_t.h"

// local macros
#define NDSI_THRESHOLD_MAX -3000     /**< Maximum NDSI for "good/non-snowy" pixels, used in prefiltering */
#define QUANTILE_MIN_NUM_VALID 3 /**< Minimum number of valid points needed to compute a quantile */
#define EVI2_FLOOR_QUANT 0.10    /**< Quantile used to compute the floor for EVI2 data */

// local prototypes
static void update_qa_evi2(size_t, int16*, uint8*);
static void prefilter_evi2(size_t, int, int, int16*, int16*, int16*, int, short**, int16*);
short* get_clim_dat(config, int, int);
static int16 masked_quantile(int16*, int16*, size_t, float, int16);
static int compare (const void*, const void*);
static void adjust_bit(short* data);
static void read_prev_spline(config, int, int16*);

/***************************************************************************//**
 * @brief Filters input data chunk (public)
 * @param [in] pp Program configuration options
 * @param [in, out] data 2D array, values to be filtered, if selected in pp.prefilter
 * @param [in, out] qa 2D array, QA flags corresponding to points in 'data'
 * @param [out] kind 2D array, flags recording the processing action taken at
 *  each point. See macro definitions for valid kind flags.
 *
 * Filters input data chunk, replacing "snow" pixels identified with the NDSI
 * index with a local fill value. Fill value may be a local quantile (1-year or
 * 3-year) or an alternative constant. See the code for the decision logic. 
 ******************************************************************************/
void prefilter(config pp, int16 **data, uint8 **qa, int16 **kind, int row) {
    short **tmin;	
    int m;
    int16 *prev_snow_fill;

    // run prefilter for EVI2 layer
    if (pp.prefilter[EVI2] == 1) {
        // ## will have 12 months of the year
        tmin = malloc(sizeof(short*) * 12);
        for(m=0;m<12;m++)
        {
            // ## each month has tmin values for that row chunk		
            tmin[m] = get_clim_dat(pp, row, m); 
        }
	    prev_snow_fill = malloc(sizeof(int16)*pp.num_elem_step);
	    read_prev_spline(pp, row, prev_snow_fill);
        prefilter_evi2(pp.num_elem_step, pp.num_years, pp.stride, data[EVI2], 
                data[NDSI], kind[EVI2], pp.tile_row,tmin,prev_snow_fill);
        update_qa_evi2(pp.num_elem_chunk, kind[EVI2], qa[EVI2]);
        // ## free the tmin array
        for(m=0;m<12;m++)
        {
            free(tmin[m]);
        }	
        free(tmin);
	free(prev_snow_fill);

    }

    // ## prefiltering is not implemented for other layers besides EVI2 
    // ## although the QA layer is used to filter out snow contimination

    return;
}

/***************************************************************************//**
 * @brief Reset QA flags using locally generated kind flags
 * @param [in] numel Number of elements in the input arrays
 * @param [in] kind Array of local kind flags
 * @param [in, out]  Array of QA flags
 *
 * Functions in the prefiltering routine make changes to the data arrays that
 * must be reflected in the smoothing spline weights. To make these changes,
 * the kind flags (which store the actions taken for each pixel) are used to
 * update the QA flags (which are used to set the spline wieghts).
 ******************************************************************************/
static void update_qa_evi2(size_t numel, int16 *kind, uint8 *qa) {

    size_t ii;

    for (ii = 0; ii < numel; ii++) {
        if (kind[ii] == KIND_NSNOW_OBS) {
            // do nothing
        } else if (kind[ii] == KIND_SNOW_OBS) {
            qa[ii] = 2;
        } else if (kind[ii] == KIND_NSNOW_INTERP) {
            qa[ii] = FILL_QA;
        } else if (kind[ii] == KIND_SNOW_INTERP) {
            qa[ii] = 3;
        } else if (kind[ii] == KIND_NSNOW_OBS_MIN) {
            qa[ii] = 2;
        } else if (kind[ii] == KIND_SNOW_INTERP_OFF) {
            qa[ii] = FILL_QA;
        } else { 
            printf("Invalid kind flag at index %zu, exiting.\n", ii);
            exit(EXIT_FAILURE);
        }
    }

    return;
}

/***************************************************************************//**
 * @brief Apply prefilter to the EVI2 layer (private)
 * @param [in] num_pix Number of pixels in chunk 
 * @param [in] num_years Number of years in data arrays
 * @param [in, out] evi2 Input layer data to be filtered
 * @param [in] ndsi NDSI index used in prefiltering
 * @param [out] kind Integer flags used to encode actions taken for each pixel
 * @param [in] tile_column the current latitude of the site
 * @param [in] tmin 2-D array for the monthly minimum temperatures
 ******************************************************************************/
static void prefilter_evi2(size_t num_pix, int num_years, int stride, int16 *evi2, 
      int16 *ndsi, int16 *kind, int tile_row, short **tmin, int16 *prev_snow_fill) {

	int16 yr_quant_1, yr_quant_min, yr_quant_max,*yr_floor, *day_floor,last_val,*day_ceil;
	// ## should be updated if make year dependable on leap year
	int months[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	int m, d, step, half_year, num_steps, yr, day, prev_step, q_start;
	int *day_mon,cur_mon,start_day, end_day, gap_len;
	size_t pix, pix_start, yr_start, yr_count, ii, prev_ii;
	// ## measured in degrees C x 10
	int tmin_thres = 0;

    	// compute constants
    	num_steps = num_years*MODIS_PERIOD;
    	half_year = MODIS_PERIOD/2;

    	// ## fill if ndsi is not available
        // ## as of 7/24/17 dont do this
    	//for(ii = 0; ii < num_pix*num_steps; ii++) {
        //	if (ndsi[ii] == FILL_NBAR) evi2[ii] = FILL_NBAR;
    	//}

    	// initialize kind flags 
    	for(ii = 0; ii < num_pix*num_steps; ii++) {
        	kind[ii] = (evi2[ii] == FILL_NBAR) ? KIND_NSNOW_INTERP : KIND_NSNOW_OBS;
    	}

    	// allocate local arrays
    	assert((yr_floor = malloc(num_years*sizeof(int16))) != NULL);
    	assert((day_floor = malloc(num_steps*sizeof(int16))) != NULL);
	    assert((day_ceil = malloc(num_steps*sizeof(int16))) != NULL);
    	assert((day_mon = malloc(MODIS_PERIOD*sizeof(int))) != NULL);

	// ## to fill day_mon array - this array contains a month value for each day
	// ## need this because we use a monthly climatology
    	start_day=0;
    	for(m=0;m<12;m++)
    	{
		    // ## figure out what the end day for that month will be        	
		    end_day = start_day + months[m];	
		    // ## loop through and fill the right set of days with that month
        	for(d=start_day;d<end_day;d++)
        	{
            		day_mon[d] = m;
        	}
		    // ## update the start day for the next month
        	start_day = end_day;
    	}

    	// loop over all pixels to do the prefiltering
    	for (pix = 0; pix < num_pix; pix++) {

        	// index of the first element for this pixel
        	pix_start = pix*num_steps;

		    // ## compute one minimum for whole three year period we are splining
		    yr_start = pix_start;
		    yr_count = MODIS_PERIOD*num_years;
		    // ## returns a value indicating a low quantile of the snow-free evi2 data for that pixel to be used as a snow-fill
		    yr_quant_1 = masked_quantile(evi2 + yr_start, ndsi + yr_start, yr_count, EVI2_FLOOR_QUANT, EVI2_FLOOR_DEFAULT); 
		    // ## also calculate the maximum quantile
		    yr_quant_max = masked_quantile(evi2 + yr_start, ndsi + yr_start, yr_count, (1-EVI2_FLOOR_QUANT), EVI2_FLOOR_DEFAULT); 

		    // ## now we assign this minimum quantile to each year in the series
		    // ## first look at previously splined values for the first year in the series
		    if(prev_snow_fill[pix]!=FILL_NBAR)
		    {
			    // ## depending on southern or northern hemisphere we either:
			    //  ## fill first two winters with previous values in northern hemisphere
			    // ## only fill first winter with previous value in southern hemisphere
			    q_start = (tile_row < 9) ? 2 : 1;	
			    for (yr = 0; yr < q_start; yr++) {
		     		yr_floor[yr] = prev_snow_fill[pix];
			    }

			    // ## check to make sure the current min is not a large increase relative to the previous year min 
			    // ## in that case ignore it
				// JMG July, 2019: if the current min is greater than 125% of the previous min then we assign the previous min?
				// that doesn't seem right...
			    yr_quant_min = ((double)yr_quant_1/(double)prev_snow_fill[pix] > 1.25) ? prev_snow_fill[pix] : yr_quant_1;
		
			    // ## fill the new years' minimums with the current quantile value
			    for (yr = q_start; yr < num_years; yr++) {
		     		yr_floor[yr] = yr_quant_min;
			    }
		    }
		    else
		    {		
			    // ## if the previous year spline doesnt exist we fill all years with the same value
			    for (yr = 0; yr < num_years; yr++) {
		     		yr_floor[yr] = yr_quant_1;
			    }
		    }  // ## end if prev_snow_fill

            //##  convert yearly floor values to daily floor values
            // ## ...first half year, use first year
            // ## only do this if the tile is in the northern hemisphere
            if(tile_row < 9) 
            {
			    // ## ... fitting day minimum values based on day of year
	      		for (step = 0; step < num_steps; step++) 
                {
				    // ## for each day in the 3 year period we determine the current year
                    		yr = step/MODIS_PERIOD;
                    		if ((step % MODIS_PERIOD) < half_year) 
                            {
                        		// ## first half of calendar year - fills with current year floor
                        		day_floor[step] = yr_floor[yr];
					            // ## also fill the ceiling value
					            day_ceil[step] = yr_quant_max;
                    		} 
                            else 
                            {                            
                        		// ## second half of calendar year - fills with next year's floor
                        		day_floor[step] = yr_floor[yr+1];
					            // ## also fill the ceiling value
					            day_ceil[step] = yr_quant_max;
                    		}
                    }  // ## end for step

               // ## if the tile is in the southern hemisphere then the year boundaries actually 
                // ## correspond to peak summer so it can be left as is
            	} 
                else 
                {
                    for (step = 0; step < num_steps; step++) 
                    {
                    	yr = step/MODIS_PERIOD;
                    	day_floor[step] = yr_floor[yr];
				        // ## also fill the ceiling value
				        day_ceil[step] = yr_quant_max;
                     }
            } // ## end if tile_v

        	// ## replace snowy pixels and truncate using daily floor values
		    // ## keep track of long gaps
		    last_val = 0;
		    gap_len = 0;
        	for (step = 0; step < num_steps; step++) 
            {
			    ii = pix_start+step;
			    // ## calculate the current year and the current julian day
			    yr = step/MODIS_PERIOD;
			    day = step - yr*MODIS_PERIOD;
			    // ## calculate the month belonging to that day
			    cur_mon = day_mon[day];
			    if (evi2[ii] != FILL_NBAR && tmin[cur_mon][pix] < tmin_thres && tmin[cur_mon][pix] != FILL_NBAR) 
                {
				    // ## if evi2 is good but the ndsi indicates snow - fill with snow
				    // ## ndsi will be fill whenever evi2 is fill already
				    // ## for cases where ndsi is corrupt we look at mean climatology of the minimum temperature for that month
                   	if (ndsi[ii] >= NDSI_THRESHOLD_MAX && ndsi[ii] != FILL_NBAR) 
                    {
					    evi2[ii] = day_floor[step];
					    kind[ii] = KIND_SNOW_OBS;				
                    } 
				    // ## if evi2 is good but the evi2 value is lower than the period minimum during a period of potential snow
                    // ## will fill with special min flag
                    if (ndsi[ii] < NDSI_THRESHOLD_MAX && evi2[ii] < day_floor[step] && ndsi[ii] != FILL_NBAR) 
                    { 
                        evi2[ii] = day_floor[step];
                        kind[ii] = KIND_NSNOW_OBS_MIN;
                    }
                 
				     // ## new special case 7/11/16 if the evi is higher than the normal ceiling but its during a winter month
				     // ## current observation is unreliable and probably snow
				     if (ndsi[ii] < NDSI_THRESHOLD_MAX && evi2[ii] >= day_ceil[step] && ndsi[ii] != FILL_NBAR) 
                     {
					     evi2[ii] = day_floor[step];
					     kind[ii] = KIND_SNOW_OBS;				
                     } 
                
				      // ## keep track of the last good kind value
				      last_val = kind[ii];
                 }  // ## end if not fill
			     else
			     {
				    // ## if fill we keep track of the gap length
				    gap_len++;
			     }

			    // ## special fill for 2 month gaps in the data record with no snow in between
			    // ## to fill that gap we need to ensure it is during a winter month - so we use the tmin climatology
			    if(gap_len>(MODIS_PERIOD/18) && //last_val == KIND_NSNOW_OBS_MIN &&
				    ((day % stride)==0) && (tmin[cur_mon][pix] < tmin_thres) && (tmin[cur_mon][pix] != FILL_NBAR)) 
                {
				    evi2[ii] = day_floor[step];
                    kind[ii] = KIND_SNOW_INTERP;
				    gap_len=0;
			    }
            }  // ## end for step
    
		     // ## replace missing data with floor value, a fix to prevent "ringing" in data gaps 
		     // ## 2-pass approach, forward pass starts at step 1, and looks backward for prefiltered values
		     // ## only replace values that are included by the stride parameter in get_input_data.c to avoid upweighting gaps
		     prev_step = 0;
		     prev_ii = pix_start;
		     for (step = prev_step+1; step < num_steps; step++) 
             {
			     yr = step/MODIS_PERIOD;
			     day = step - yr*MODIS_PERIOD;
			     prev_ii = ii;
                 prev_step = prev_ii-pix_start;
                 ii = pix_start + step;
			     // ## fill evi if the value falls on a stride day and the previous value was snow
			     if ((evi2[ii] == FILL_NBAR) && (kind[prev_ii] > KIND_NSNOW_OBS_MIN) && ((day % stride)==0)) 
                 {
                          evi2[ii] = day_floor[step];
                          kind[ii] = KIND_SNOW_INTERP;
                  } 
			      // ## change the kind flag if the last value was snow
			      if ((evi2[ii] == FILL_NBAR) && (kind[prev_ii] > KIND_NSNOW_OBS_MIN) && ((day % stride)!=0)) 
                  {
                          kind[ii] = KIND_SNOW_INTERP_OFF;
			       }
		        }  // ## end step for first pass
            
		        // ## 2-pass approach, reverse pass starts at the last step - 1, and looks foreward for prefiltered values
		        // ## only replace values that are included by the stride parameter in get_input_data.c to avoid upweighting gaps
		        for (step = (ii-1-pix_start); step >= 0; step--) 
                { 
                    // ii is initialized in the previous loop
                    yr = step/MODIS_PERIOD;
                    day = step - yr*MODIS_PERIOD;
			        prev_ii = ii;
                    prev_step = prev_ii-pix_start;
                    ii = pix_start + step;
			        // ## fill evi if the value falls on a stride day and the previous value was snow
			        if ((evi2[ii] == FILL_NBAR) && (kind[prev_ii] > KIND_NSNOW_OBS_MIN) && ((day % stride)==0)) 
                    {
                          evi2[ii] = day_floor[step];
                          kind[ii] = KIND_SNOW_INTERP;
                    } 
			        // ## change the kind flag if the last value was snow
			        if ((evi2[ii] == FILL_NBAR) && (kind[prev_ii] > KIND_NSNOW_OBS_MIN) && ((day % stride)!=0)) 
                    {
                               kind[ii] = KIND_SNOW_INTERP_OFF;
			        }
		        } // ## end step for second pass

        } // ## end for pix

        // free local arrays
        free(yr_floor);
        free(day_floor);
	    free(day_ceil);
        free(day_mon);

        return;
    }

    /***************************************************************************//**
     * @brief Get specified quantile of input data where ndsi < NDSI_THRESHOLD_MAX
     * @param [in] data 1D array containing data values  
     * @param [in] ndsi 1D array containing NDSI index values  
     * @param [in] numel Number of elements in input arrays
     * @param [in] quantile Scalar, quantile to be returned
     * @param [in] default_quantile Scalar, default value to return if not enough points
     * @returns quantile where ndsi < NDSI_THRESHOLD_MAX
     ******************************************************************************/
    static int16 masked_quantile(int16 *data, int16 *ndsi, size_t numel, 
            float quantile, int16 default_quantile) {

        int16 quantile_value, *valid_data;
        size_t num_valid, ii;

        // allocate temporary arrays        
        assert((valid_data = malloc(numel*sizeof(int16))) != NULL);

        // copy valid (i.e. not snowy and not fill) data to temporary array and count it
        num_valid = 0;
        for (ii = 0; ii < numel; ii++) {
            // ## if ndsi < threshold then it is valid observation
            if ((ndsi[ii] < NDSI_THRESHOLD_MAX) && (data[ii] != FILL_NBAR) && (ndsi[ii] != FILL_NBAR) && (data[ii] > 1400)) {
                valid_data[num_valid] = data[ii]; 
                num_valid++;
            }
        }

        if (num_valid < QUANTILE_MIN_NUM_VALID) {
            // use default if there are two few data
            quantile_value = default_quantile;

        } else {
            // compute the requested quantile, decimal ranks are truncated by casting to int
            qsort(valid_data, num_valid, sizeof(int16), compare); 
            quantile_value = valid_data[ (size_t)(quantile*num_valid) ];
	    if(quantile_value < default_quantile) { quantile_value = default_quantile; }
        }

        // free temporary arrays
        free(valid_data);

        return quantile_value;
    }

    /***************************************************************************//**
     * @brief Load climate data for row chunk
     * @param [in] pp config information
     * @param [out] tmin 2-D Array minimum temperature for each month
     * @param [in] row the current row in the file to be loaded
     ******************************************************************************/
    short* get_clim_dat(config pp, int row, int mon) 
    {
        short *tmin;
        size_t i;
        char full_path[LONG_STR];
        FILE *clim_file;
        long f_loc;

        // ## this is the full name to the climate file - needs to be structured exactly as it is here
        sprintf(full_path,"%s/tmin%d_tiles/tmin%d_%s.bin",pp.clim_dir,mon+1,mon+1,pp.tile_id);
        tmin = malloc(sizeof(short)*pp.num_elem_step);
	
        // ## try to open file
        if( (clim_file = fopen(full_path, "r")) != NULL  )
        {
            // ## floc is the location in the file we will start reading from
            f_loc = row*MODIS_NROW*sizeof(short);
            // ## we seek forward to the location in the file to start reading - f_loc is originally 0
            assert( (fseek(clim_file, f_loc, SEEK_SET)) == 0 );
            // ## we read a row chunk into tmin array
            assert( fread(tmin,sizeof(short),pp.num_elem_step,clim_file) );
            // ## adjust bits because they come out big endian and we need them in little
            for(i=0;i<pp.num_elem_step;i++)
            {
                adjust_bit(&tmin[i]);
            }

            if(pp.verbose == 1)
            {
                printf("Climate file %s was found and read in.\n",full_path);			
                printf("Size of row chunk = %ld, location = %ld, month = %d\n",pp.num_elem_step*sizeof(short),row*MODIS_NROW*sizeof(short),mon);
            }
            // ## only close file if has been read in
            fclose(clim_file);
        }
        // ## if doesnt exist then we fill and print out warning
        else
        {
            printf("Climate file %s not found. Will FILL.\n",full_path);		
            // ## fille for files that don't exist
            for(i=0;i<pp.num_elem_step;i++)
            {
                tmin[i]=FILL_NBAR;
            }	
        }

        return(tmin);
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

    /***************************************************************************//**
     * @brief This function is used to switch big and small endian
     * @param [in/out] a value to be switched
     ******************************************************************************/
    static void adjust_bit(short* data)
    {
        short a,b,c,d=255;
        a=*data;
        b=*data;
        a=a>>8;
        b=b<<8;
        if (a<0)
            a=a&d;
        c=a|b;
        *data=c;	
    }

/***************************************************************************//**
 * @brief Read previous year's splined EVI2 data if exists to get a minimum value
 * @param [in] pp Program configuration options
 * @returns 1D array of fill values
 ******************************************************************************/
static void read_prev_spline(config pp, int row, int16 *snow_vals) {
    
	FILE *file_spline,*file_flag;
	char spline_name[LONG_STR],flag_name[LONG_STR];
	int yr,find_flag;
	size_t ii,jj,chunk_sz,cur_loc,cur_ii; 
	long f_loc;
	int16 *spline_dat,*flag_dat;
    
	yr = pp.start_year; // ## we want the spline of the start year if it exists
    
        // the EVI2 index - should be saved in the out_dir
	sprintf(spline_name, "%s/%s.%s.%d.evi2.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr);
	sprintf(flag_name, "%s/%s.%s.%d.evi2.flag.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr);

	chunk_sz = pp.num_elem_step*MODIS_PERIOD;
	spline_dat = malloc(sizeof(int16)*chunk_sz);
	flag_dat = malloc(sizeof(int16)*chunk_sz);

        // ## try to open file
        if( (file_spline = fopen(spline_name, "r")) != NULL  )
        {
		// ## floc is the location in the file we will start reading from
		f_loc = row*MODIS_NROW*MODIS_PERIOD*sizeof(int16);
	//	## for debugging
	//	printf("%ld, %zu, %d\n",f_loc,chunk_sz,row);
		// ## we seek forward to the location in the file to start reading - f_loc is originally 0
		assert( (fseek(file_spline, f_loc, SEEK_SET)) == 0 );
		// ## we read a row chunk into tmin array
		assert( fread(spline_dat,sizeof(int16),chunk_sz,file_spline) );
		fclose(file_spline);

	 } // ## if doesnt exist then we fill and print out warning
        else
        {
		printf("Previous spline file %s not found. Will FILL.\n",spline_name);		
		// ## fille for files that don't exist
		for(ii=0;ii<chunk_sz;ii++)
		{
 			spline_dat[ii]=FILL_NBAR;
		}	
        }

	if( (file_flag = fopen(flag_name, "r")) != NULL  )
        {
		// ## floc is the location in the file we will start reading from
		f_loc = row*MODIS_NROW*MODIS_PERIOD*sizeof(int16);
		// ## we seek forward to the location in the file to start reading - f_loc is originally 0
		assert( (fseek(file_flag, f_loc, SEEK_SET)) == 0 );
		// ## we read a row chunk into tmin array
		assert( fread(flag_dat,sizeof(int16),chunk_sz,file_flag) );
		fclose(file_flag);
	} else {
		printf("Previous flag file %s not found. Will FILL.\n",flag_name);		
		// ## fille for files that don't exist
		for(ii=0;ii<chunk_sz;ii++)
		{
 			flag_dat[ii]=-1;
		}
        }
	// ## this is simple we just take the last snow fill value as the starting snow-fill value for the next spline
	find_flag = 0;
	for(jj=0;jj<pp.num_elem_step;jj++)
	{
		cur_loc = jj*MODIS_PERIOD;
		cur_ii = MODIS_PERIOD-1;
		for(ii=(MODIS_PERIOD-1);ii>0;ii--)
		{
			if(flag_dat[cur_loc + ii] > KIND_NSNOW_OBS_MIN)
			{
				cur_ii=ii;				
				find_flag=1;				
				break;
			}
		}
		// ## only use it if we have a valid snow min value to use
		if(find_flag==1) {
			snow_vals[jj]=spline_dat[cur_loc + cur_ii];
		}
		else 
		{
			snow_vals[jj]=FILL_NBAR;
		}
	}
	free(spline_dat);
	free(flag_dat);
}

