/***************************************************************************//**
 * @file
 * @brief Functions for finding HDF file stack for a particular MODIS tile
 ******************************************************************************/
 
#include "modis_smooth_t.h"
#include <string.h> // strcmp()
#include <dirent.h> // dirent type, opendir(), etc

// local macros
#define MAX_REDUNDANT_FILES 10 /**<Number of redundant files retained by find_file*() */

// local prototypes
void string_sort(char**, int);
void find_file(char*, char*, int, int, char*, char*);

/***************************************************************************//**
 * @brief Locate MODIS HDF files for the desired time and time period
 * @param [in] pp Program configuration, as populated by set_options()
 * @param [out] nbar Struct containing NBAR input file info
 * @param [out] qa Struct containing QA input file info
 ******************************************************************************/
void get_input_files(config pp, data_loc *nbar, data_loc *qa) {

    int ii, yr, day, loc_count, nbar_file_count, qa_file_count;
    char nbar_file[LONG_STR], qa_file[LONG_STR];

    // initialize outputs 
    nbar->num_files = pp.num_years*MODIS_PERIOD;	
	qa->num_files   = pp.num_years*MODIS_PERIOD;	
	nbar->names = malloc(nbar->num_files*sizeof(char*));
	qa->names   = malloc(qa->num_files*sizeof(char*));
	nbar->is_fill = malloc(nbar->num_files*sizeof(int));
	qa->is_fill   = malloc(qa->num_files*sizeof(int));

	for(ii = 0; ii < nbar->num_files; ii++) {
		nbar->names[ii] = malloc(LONG_STR*sizeof(char));
		qa->names[ii]   = malloc(LONG_STR*sizeof(char));
		nbar->is_fill[ii] = 1; 
		qa->is_fill[ii]   = 1; 
	}
    
    // find files
    loc_count = 0;
    nbar_file_count = 0;
    qa_file_count = 0;
    for(yr = pp.start_year; yr<(pp.start_year+pp.num_years); yr++) { 
		for (day = 1; day<(MODIS_PERIOD+1); day++) {
			// JMG July, 2019; q&d fix so we only need half a year of data past year of interest:
			// sub "FILL" for nbar/qa files when in last year when day is July or past
            if ( ((day-1) % pp.stride) == 0 ) {
				// step is divisible by stride, check if we're past halfway in final year
				if((yr == pp.start_year + pp.num_years) & day > 181) {
					// we're past halway, we'll use FILL for these files
					sprintf(nbar_file, "FILL");
					sprintf(qa_file, "FILL");
				} else {
					// step is divisible by stride and we're not past halfway in final year
					// look for files
					find_file(pp.nbar_dir, pp.tile_id, yr, day, INPUT_STEM_NBAR, nbar_file);
					find_file(pp.qa_dir,   pp.tile_id, yr, day, INPUT_STEM_QA,   qa_file);
				}
            } else {
                // step is NOT divisible by stride, skip
                sprintf(nbar_file, "FILL");
                sprintf(qa_file, "FILL");
            }

            if( strcmp(nbar_file, "FILL") != 0 ) {
                sprintf(nbar->names[loc_count], "%s", nbar_file);
                nbar->is_fill[loc_count] = 0;
                nbar_file_count++;
            }
            if( strcmp(qa_file, "FILL") != 0 ) {
                sprintf(qa->names[loc_count], "%s", qa_file);
                qa->is_fill[loc_count] = 0;
                qa_file_count++;
            }

            if (pp.verbose == 1) {
                printf("get_input_files: NBAR: %s\n", nbar_file);
                printf("get_input_files: QA  : %s\n", qa_file);
            }

            loc_count++;
        }
	}
	
	// report results
	printf("get_input_files: using %d NBAR files for %d days\n", nbar_file_count, nbar->num_files);
	printf("get_input_files: using %d QA files for %d days\n",  qa_file_count,   qa->num_files);
	
	return;
}

/***************************************************************************//**
 * @brief Deallocate data_loc struct
 ******************************************************************************/
void close_loc(data_loc *loc) {
   
    int ii; 

	for(ii = 0; ii < loc->num_files; ii++) {
        free(loc->names[ii]);
    }
    free(loc->names);
    free(loc->is_fill);

    return;
}

/***************************************************************************//**
 * @brief Find MODIS data file given parts of the path 
 * @param [in]  top_dir top-level directory containing input files
 * @param [in]  tile tile ID from config struct
 * @param [in]  year Observation year
 * @param [in]  day Observation day
 * @param [in]  stem File name stem/prefix
 * @param [out] out_str File name, or if no file is found, "FILL"
 ******************************************************************************/
void find_file(char *top_dir, char *tile, int year, int day, char *stem, char* out_str) {

	struct dirent entry, *result;
	DIR *dp;
	int count, ii, status;
    int match_stem, match_date, match_tile;
	char **list;
	char temp_name[LONG_STR], copy_name[LONG_STR];
	char *cmp_str; 
	char tok[] = "."; // token that will break up the file names for counting
    char dir[LONG_STR];
    char date[SHORT_STR];

    // generate file name parts
    sprintf(dir, "%s/%d/%03d", top_dir, year, day);
    sprintf(date, "A%d%03d", year, day);	

	// allocate space for filenames
	list = malloc(MAX_REDUNDANT_FILES*sizeof(char*));
	for(ii = 0 ; ii < MAX_REDUNDANT_FILES; ii++) {
		list[ii] = malloc(LONG_STR*sizeof(char));
	}	

	// access the directory of interest
  	dp = opendir(dir);
  	if (dp == NULL) { // file not found
    	sprintf(out_str, "FILL");
		return;
	}

 	count = 0;
  	while(dp)
	{
		status = readdir_r(dp, &entry, &result);

		if(result!=NULL && status == 0) 
		{
			// assign a temp name to keep the current item
			sprintf(temp_name, "%s", entry.d_name);

			// copy the location of the current name to another temp name
			sprintf(copy_name, "%s", temp_name);

			// Parse the file name and attempt to match all segments 
            		// ...first token is the name stem
            		match_stem = 0;
           		 cmp_str = strtok(copy_name,tok); 
			if( cmp_str != NULL && strcmp(cmp_str, stem) == 0 ) {
				match_stem = 1;
            }
            // ...second token is the date string
            match_date = 0;
            cmp_str = strtok(NULL,tok);  
            if ( cmp_str != NULL && strcmp(cmp_str, date) == 0 ) {
                match_date = 1;
            }
            // ...third token is the tile id string
            match_tile = 0;
            cmp_str = strtok(NULL,tok);  
            if ( cmp_str != NULL && strcmp(cmp_str, tile) == 0 ) {
                match_tile = 1;
            }
				
			// if we find a match, copy to list and go to the next filename
			if( match_stem && match_date && match_tile ) {
				sprintf(list[count], "%s", temp_name);
				count++;
				continue;
			}
			
		} else {
			break;
		}
		
	}
	
	// close the dir
	if(dp) { 
		closedir(dp);
	}

	// extract file name from list	
	if(count==0) { 
		// too few files found, just fill
		sprintf(out_str, "FILL");
	} else if (count == 1) { 
		// we have our single file
		sprintf(out_str, "%s/%s", dir, list[0]);
	} else if (count > 1) {
		// sort list to get most recent
		string_sort(list, count);
		sprintf(out_str, "%s/%s", dir, list[0]);	
	}
		
	// free allocated variables	
	for (ii=0; ii < MAX_REDUNDANT_FILES; ii++) {
		free(list[ii]);
	}
	free(list);

	return;
}

/***************************************************************************//**
 * @brief Sort array of strings in descending order 
 ******************************************************************************/
void string_sort(char **a, int size)
{
	char *temp = NULL;
   	int i = 0;
  	int sorted = FALSE;
   	while(!sorted) {
		sorted = TRUE;
     	for( i = 0; i<size-1; i++) {
       		// sort alphabetically descending = put largest first
			// only switch if the first is smaller than the second
			if(strcmp(a[i], a[i+1]) < 0) {
         		sorted = FALSE;
         		temp = a[i+1];
         		a[i+1] = a[i];
         		a[i]  = temp;
       		}
   		}
	}
	return;
}
