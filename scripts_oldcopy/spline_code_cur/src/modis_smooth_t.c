/***************************************************************************//**
 * @file
 * @brief MODIS_SMOOTH_T main program
 * @author Damien Sulla-Menashe
 * @author Katia Oleinik
 * @author Keith Ma
 * 
 * Time-series smoothing for stacked MODIS C6 data
 *
 * Style notes: 
 *   - only "public" functions make use of special data formats (config, layer arrays)
 *   - "private" functions are marked as static
 ******************************************************************************/

#include "modis_smooth_t.h"

// local prototypes
static void allocate_chunk_arrays(config, int16***, uint8***, int16***, int16***, int16***);
static void free_chunk_arrays(config, int16**, uint8**, int16**, int16**, int16**);
static void set_size_const(config*);
static int read_pixels(config, long**);
static void print_extract_to_file(config, int16**,int16**,long*,FILE*, int, int);
static int get_line(FILE*,char*);

/***************************************************************************//**
 * @brief Main function for modis_smooth_t
 * @param [in] argc Command line argument cound
 * @param [in] argv Array of command line argument values
 ******************************************************************************/
int main(int argc, char **argv)
{
    
	config param;
    	data_loc nbar_loc, qa_loc;
    	FILE ***smooth_out, ***residual_out, ***flag_out, **quant_out;
    	FILE *extract_out;
    	int row0, row1,num_pix;
    	int16 **data_chunk, **residual_chunk, **index_chunk, **flag_chunk;
    	uint8 **qa_chunk;
    	int extract_flag;
   	    char extract_out_path[LONG_STR];
    	long *pix_list;

    	// ## read program options from defaults, config file, and command line
    	get_options(argc, argv, &param);
    	set_size_const(&param);

    	// ## read in extract file if this is a desired option
	extract_flag = 0;
	if (strcmp(param.extract_path, "NONE") != 0) {
	
		// ## read the pixels for this tile
		num_pix = read_pixels(param, &pix_list);
		// ## only create output if there are pixels for this tile
		if(num_pix>0)
		{
			extract_flag = 1;
			// ## open file for write out
			sprintf(extract_out_path,"%s/%s.%s.%d-%d.extract.txt", param.out_dir, param.out_stem, 
                        param.tile_id, param.start_year,param.start_year+param.num_years-1);
			assert( (extract_out = fopen(extract_out_path, "w")) != NULL);
		}
   	 }

    	// get HDF stacks for NBAR and QA
   	get_input_files(param, &nbar_loc, &qa_loc);

    	// create yearly output files, skipping first and last years
    	smooth_out   = get_output_files(param, 0);
    	residual_out = get_output_files(param, 1); 
    	flag_out = get_output_files(param, 2); 

    	if(param.quantile == 1) {
        	quant_out = get_output_files_quant(param);
    	}

	    // allocate chunk arrays
    	allocate_chunk_arrays(param, &data_chunk, &qa_chunk, &residual_chunk, &flag_chunk, &index_chunk); 

    	// process chunk-by-chunk
    	row1 = -1;
  //      row1 = 1200;
    	while (row1 < MODIS_NROW-1) 
        {
    
        	// define chunk - based on location in while loop
        	row0 = row1+1;
        	row1 = min(row0+param.chunk_num_row-1, MODIS_NROW-1);
        	param.chunk_num_row = row1-row0+1; // only modified for the last chunk, which may be truncated  
        	set_size_const(&param);

        	// read input nbar and qa data
        	wallclock_timer(1);
        	get_chunk(param, nbar_loc, qa_loc, row0, data_chunk, qa_chunk);
        	printf("Read rows %d to %d ... %f seconds\n", row0, row1, wallclock_timer(0));

        	// compute spectral indices EVI2 and NDSI for splining
        	wallclock_timer(1);
        	get_indices(param, data_chunk, qa_chunk);
        	printf("Calculated EVI2 and NDSI indices, rows %d to %d ... %f seconds\n",row0, row1, wallclock_timer(0));

        	// prefilter 
       		wallclock_timer(1);
        	prefilter(param, data_chunk, qa_chunk, flag_chunk, row0);
        	printf("Prefiltered rows %d to %d ... %f seconds\n", row0, row1, wallclock_timer(0));

		    // ## if extract flag is on then we output the pre-splined values for each band in a text file
	 	    if (extract_flag==1) {
			    // ## figure out which pixels we need to read in for this row chunk and write data to file
			    print_extract_to_file(param, data_chunk, flag_chunk, pix_list, extract_out, num_pix, row0);
    		}

        	// smooth chunk, returning smoothed values and residuals
        	wallclock_timer(1);
        	apply_spline(param, data_chunk, qa_chunk, residual_chunk); 
        	printf("Smoothed rows %d to %d ... %f seconds\n", row0, row1, wallclock_timer(0)); 

		    // ## if extract flag is on then we output the splined values for each band in a text file
		    if (extract_flag==1) {
			    // ## figure out which pixels we need to read in for this row chunk post splining and write data to file
			    print_extract_to_file(param, data_chunk, flag_chunk, pix_list ,extract_out ,num_pix, row0);
    		}

        	// write chunk to yearly output files, skipping first and last years
        	wallclock_timer(1);
        	write_chunk(param, param.chunk_num_row, data_chunk, smooth_out, 0);
        	write_chunk(param, param.chunk_num_row, residual_chunk, residual_out, 1);
		    write_chunk(param, param.chunk_num_row, flag_chunk, flag_out, 2);
        	printf("Wrote spline results for rows %d to %d ... %f seconds\n", row0, row1, wallclock_timer(0));

        	// compute and write splined indices and quantiles
        	if(param.quantile == 1) 
            {
            		wallclock_timer(1);
            		get_indices_spl(param, data_chunk, index_chunk);
            		calc_all_quantiles(param, data_chunk, index_chunk, flag_chunk[EVI2], quant_out);
            		printf("Calculated indices and quantiles, rows %d to %d ... %f seconds\n", row0, row1, wallclock_timer(0));
        	}
           
    	}  // ## end while loop

   	// free allocated memory
    	close_output_files(param, 0, smooth_out);
    	close_output_files(param, 1, residual_out);

    	if(param.quantile == 1) {
        	close_output_files_quant(param, quant_out);
    	}
	    // ## close input files and free all the data chunks
        free_chunk_arrays(param, data_chunk, qa_chunk, residual_chunk, flag_chunk, index_chunk);

    	close_loc(&nbar_loc);
    	close_loc(&qa_loc);
    	

     	if (extract_flag==1) {
		    free(pix_list);
		    // ## close file for writing out extract
		    fclose(extract_out);
    	}

    	return EXIT_SUCCESS;
    }  // ## end main

    /***************************************************************************//**
     * @brief allocate memory for data, qa, and residual chunk arrays
     * @param [in] pp program configuration options
     * @param [out] data 2D data array to be allocated
     * @param [out] qa 2D QA array to be allocated
     * @param [out] residual 2D residual array to be allocated
     * @param [out] kind 2D kind flags array to be allocated
     ******************************************************************************/
    static void allocate_chunk_arrays(config pp, int16 ***data, uint8 ***qa, \
            int16 ***residual, int16 ***kind, int16 ***index) {

        int ii;

        assert(((*data)     = malloc(NUM_DATA_LAYERS*sizeof(int16*))) != NULL);  
        assert(((*residual) = malloc(NUM_DATA_LAYERS*sizeof(int16*))) != NULL);  
        assert(((*qa)       = malloc(NUM_DATA_LAYERS*sizeof(uint8*))) != NULL);  
        assert(((*kind)     = malloc(NUM_DATA_LAYERS*sizeof(int16*))) != NULL);  
        assert(((*index) = malloc(NUM_SPL_INDICES*sizeof(int16*))) != NULL); 

        for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {

            if (pp.read[ii] == 1) {
                assert(((*data)[ii] = malloc(pp.num_elem_chunk*sizeof(int16))) != NULL);
                assert(((*qa)[ii]   = malloc(pp.num_elem_chunk*sizeof(uint8))) != NULL);
            }
            if (pp.residual[ii] == 1) {
                assert(((*residual)[ii] = malloc(pp.num_elem_chunk*sizeof(int16))) != NULL);
            }
            if (pp.prefilter[ii] == 1) {
                assert(((*kind)[ii]   = malloc(pp.num_elem_chunk*sizeof(uint16))) != NULL);
            }
        }

        if (pp.quantile == 1) {
            for (ii = 0; ii < NUM_SPL_INDICES; ii++) {
                assert(((*index)[ii]   = malloc(pp.num_elem_chunk*sizeof(int16))) != NULL);
            }
        }

        return;
    }

    /***************************************************************************//**
     * @brief free data, qa, and residual chunk arrays
     * @param [in] pp program configuration options
     * @param [in,out] data 2D data array to be freed
     * @param [in,out] qa 2D QA array to be freed
     * @param [in,out] residual 2D residual array to be freed
     * @param [in,out] kind 2D kind flag array to be freed
     ******************************************************************************/
    void free_chunk_arrays(config pp, int16 **data, uint8 **qa, int16 **residual, \
            int16 **kind, int16 **index) {

        int ii;

        for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {
            if (pp.read[ii] == 1) {
                free(data[ii]);
                free(qa[ii]);
            }
            if (pp.residual[ii] == 1) {
                free(residual[ii]);
            }
            if (pp.prefilter[ii] == 1) {
                free(kind[ii]);
            }
        }
        if (pp.quantile == 1) {
            for (ii = 0; ii < NUM_SPL_INDICES; ii++) {
                free(index[ii]);
            }
        }
        free(data);
        free(qa);
        free(residual);
        free(kind);
        free(index);

        return;
    }


    /***************************************************************************//**
     * @brief Compute array size parameters
     * @param [in, out] pp program configuration options
     ******************************************************************************/
    static void set_size_const(config *pp) {

        pp->num_step       = pp->num_years     * MODIS_PERIOD; 
        pp->num_elem_step  = pp->chunk_num_row * MODIS_NCOL;
        pp->num_elem_chunk = pp->num_step * pp->num_elem_step; 

        return;
    }

    /***************************************************************************//**
     * @brief Get list of pixels if necessary
     * @param [in] pp program configuration options
     * @param [in] pix_list long int array of pixels to extract for that tile
     ******************************************************************************/
    static int read_pixels(config pp, long** pix_list) {
	FILE *extract_in;
	int count,x,y,tile_x,tile_y,id;
	char record_string[SHORT_STR];
	char *token_string;

	// ## open the file and read each row
	assert( (extract_in=fopen(pp.extract_path,"r")) != NULL);
	// ## read through once to just count the number of pixels we will be extracting for this tile
	count=0;
	while(get_line(extract_in,record_string))
	{
		// ## take in the pixels into an array
		token_string=strtok(record_string,",");
		sscanf(token_string,"%d",&tile_x);       //get the tile_x
		token_string=strtok(NULL,",");
		sscanf(token_string,"%d",&tile_y);       // get the tile_y
		token_string=strtok(NULL,",");
		sscanf(token_string,"%d",&x);       //get the pix col
		token_string=strtok(NULL,",");
		sscanf(token_string,"%d",&y);       // get the pix row

		// ## count how many pixels belong to this tile
		if(tile_x==pp.tile_column && tile_y==pp.tile_row && x<MODIS_NCOL && y<MODIS_NROW)
		{
			count++;
		}
	}
	fclose(extract_in);

	// ## open a second time and read in the data we need
	if(count > 0)
	{
		*pix_list = malloc(sizeof(long)*count);

		extract_in=fopen(pp.extract_path,"r");
		count=0;
		while(get_line(extract_in,record_string))
		{
			// ## take in the pixels into an array
			token_string=strtok(record_string,",");
			sscanf(token_string,"%d",&tile_x);       //get the tile_x
			token_string=strtok(NULL,",");
			sscanf(token_string,"%d",&tile_y);       // get the tile_y
			token_string=strtok(NULL,",");
			sscanf(token_string,"%d",&x);       //get the pix col
			token_string=strtok(NULL,",");
			sscanf(token_string,"%d",&y);       // get the pix row
			token_string=strtok(NULL,",");
			sscanf(token_string,"%d",&id);       // get the pix row
		
			if(tile_x==pp.tile_column && tile_y==pp.tile_row)
			{
				(*pix_list)[count] = (long)(x + (y * MODIS_NCOL));
				// ## debug statement
			//	printf("%d,%ld\n",count,(*pix_list)[count]);
				count++;
				
			}
		}
		fclose(extract_in);
	
	}
	return(count);
    }

      /***************************************************************************//**
     * @brief read a valid line from the ASCII file
     * @param [in] file a file with an open file pointer to be read from
     * @param [out] string_input the returned line of text
     ******************************************************************************/
   int get_line(FILE *file,char* string_input)
   {
	int indicator=1;
	char char_input;
	int string_count=0;
	int blank_count=0;
	int initial_sign=0;
	while ((char_input=fgetc(file))!=EOF)
	{
	
		if (initial_sign==0)  /***the first letter of a line****/
		{
			if (char_input=='#')
				indicator=0;
			else
				indicator=1;
			initial_sign=1;
		}
		if (char_input=='\n')  /****the end of a line *****/ 
		{
			if (indicator==1 && string_count>0)  /**** there is valid line ******/
			{
				if (blank_count==string_count)  /***** one line contains only space, ignore it *****/
				{	
					string_count=0;
					blank_count=0;
					initial_sign=0;
					continue;
				}
				string_input[string_count]='\0';
				return 1;
			}
			initial_sign=0;
			continue;
		}
		
		if (indicator==1)
		{
			if (char_input==' ')
				blank_count++;
			string_input[string_count]=char_input;
			string_count++;
		}
	}
	return 0;
   }

   static void print_extract_to_file(config pp, int16** data, int16** flag, long* pix_list, FILE* extract_out, int num_pix, int row0)
   {
	size_t jj, ii, kk, start,end,start_data;
	
	start = row0*MODIS_NCOL;
	end = (row0+pp.chunk_num_row)*MODIS_NCOL;	
	// ## debug statement	
//	printf("Inside print to file: %d,%zu,%zu\n",num_pix,start,end);
	fflush(stdout);
	// main loop, spline each pixel in each layer 
	for (ii = 0; ii < (size_t)num_pix; ii++)
	{
		// ## for the current pixel to output we make sure it is within the row chunk of interest
		if((pix_list[ii] >= (long)start) && (pix_list[ii] < (long)end))
		{  
			// ## debug statement	
		//	printf("Inside loop: %zu,%zu,%zu,%ld\n",ii,start,end,pix_list[ii]);
	  	//	fflush(stdout);
			// index of first time step for pixel pix_list[n] in data array 
			start_data = (pix_list[ii]-start)*pp.num_step; 				
			for (jj = 0; jj < NUM_DATA_LAYERS; jj++) 
			{
				if (pp.smooth[jj] == 0) continue; // skip if layer is not requested
				
				fprintf(extract_out,"%ld,%zu,%zu,",pix_list[ii],ii,jj);
				for (kk = 0; kk < (pp.num_step-1); kk++) {
					fprintf(extract_out,"%d,",data[jj][start_data + kk]);
				}
				fprintf(extract_out,"%d\n",data[jj][start_data + pp.num_step-1]);

			}
			// ## for snow flags
			fprintf(extract_out,"%ld,%zu,255,",pix_list[ii],ii);
			for (kk = 0; kk < (pp.num_step-1); kk++) {
				fprintf(extract_out,"%d,",flag[EVI2][start_data + kk]);
			}
			fprintf(extract_out,"%d\n",flag[EVI2][start_data + pp.num_step-1]);
		} // ## end if cur_loc < tot_pixels

	} // ## end ii loop
	

   } // ## end function

