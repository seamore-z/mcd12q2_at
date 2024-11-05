/***************************************************************************//**
 * @file
 * @brief Functions for creating output files and writing data to them
 ******************************************************************************/

#include "modis_smooth_t.h"
#include "mfhdf.h" // HDF types (e.g. int16)

// local macros
#define PIX_LEN 463.312716525                       /**<ENVI header file parameter*/
#define MAP_X -20015109.354                         /**<ENVI header file parameter*/
#define MAP_Y 10007554.677                          /**<ENVI header file parameter*/
#define PROJ "{16, 6371007.181, 0, 0, 0,MODIS Sin}" /**<ENVI header file parameter*/

// local prototypes
FILE* create_file(const char*, size_t, int);
void create_envi_header(const char*, config, int, int);

/***************************************************************************//**
 * @brief Get array of file pointers for output files
 * @param [in] pp Program configuration options
 * @param [in] kind flag indicating file type. Valid options are 0 = smoothed
 *     data layers,  1 = residual data layers, 2 = flags encoding processing info
 * @returns 2D array of file pointers to output files, with dimensions [layer][year]
 *
 * Each output file contains 1 year of observations for 1 layer. The first and
 * last years are not written - they serve as buffers to exclude unpredictable
 * behavior of the smoothing spline at the data edges. Disk space is allocated
 * upon file creation to avoid out-of-space errors that could arise during
 * processing. Outputs are assumed to be int16.
 *
 * This function works for spline, residual, and flag outputs. Note that the
 * flag output is enabled only for prefiltered layers.
 ******************************************************************************/
FILE*** get_output_files(config pp, int kind) {
    
    FILE ***files;
    char head_name[LONG_STR], file_name[LONG_STR], suffix[SHORT_STR];
    int ii, jj, yr, yr0, yr1;
    size_t sz; 
    
    // select suffix and element size
    if (kind == 0) {
        suffix[0] = 0; // empty string
        sz = sizeof(int16);
    } else if (kind == 1) {
        sprintf(suffix, ".residual");
        sz = sizeof(int16);
    } else if (kind == 2) {
        sprintf(suffix, ".flag");
	sz = sizeof(int16);
     //   sz = sizeof(uint8);
    } else {
        printf("get_output_files: Invalid value for kind (%d), exiting.\n", kind);
        exit(EXIT_FAILURE);
    }

    yr0 = pp.start_year+1; // omit initial buffer year
    yr1 = pp.start_year+pp.num_years-2; // omit final buffer year

    files = malloc(NUM_DATA_LAYERS*sizeof(FILE**));
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {

        // skip if layer is not requested
        if (kind == 0 && pp.smooth[ii] == 0) continue; 
        if (kind == 1 && pp.residual[ii] == 0) continue; 
        if (kind == 2 && pp.prefilter[ii] == 0) continue; 

	// ## dsm added this - only want to output EVI smoothed splines and residuals
	if (ii != EVI2) continue;

        // allocate FILE pointers for each year in band
        files[ii] = malloc( (pp.num_years-2)*sizeof(FILE*) );

        // create files and headers
        jj = 0;
        for (yr = yr0; yr <= yr1; yr++) {

            if (ii <= SWIR3) {
                // MODIS bands
                sprintf(head_name, "%s/%s.%s.%d.b%d%s.hdr", pp.out_dir, pp.out_stem, pp.tile_id, yr, ii+1, suffix);
                sprintf(file_name, "%s/%s.%s.%d.b%d%s.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr, ii+1, suffix);
            } else if (ii == EVI2) {
                // the EVI2 index
                sprintf(head_name, "%s/%s.%s.%d.evi2%s.hdr", pp.out_dir, pp.out_stem, pp.tile_id, yr, suffix);
                sprintf(file_name, "%s/%s.%s.%d.evi2%s.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr, suffix);
            } else if (ii == NDSI) {
                // the NDSI index
                sprintf(head_name, "%s/%s.%s.%d.ndsi%s.hdr", pp.out_dir, pp.out_stem, pp.tile_id, yr, suffix);
                sprintf(file_name, "%s/%s.%s.%d.ndsi%s.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr, suffix);
            }

            create_envi_header(head_name, pp, sz, MODIS_PERIOD); 
            files[ii][jj] = create_file(file_name, MODIS_NROW*MODIS_NCOL*MODIS_PERIOD, sz); 
            jj++; // increment year index

            if (pp.verbose == 1) {
                printf("get_output_files: created %s\n", file_name);
            }
        }
    }

    return files; 
}

/***************************************************************************//**
 * @brief Get array of file pointers for output files
 * @param [in] pp Program configuration options
 * @returns 1D array of file pointers to output files, with dimensions [year]
 *
 * Each output file contains 1 year of observations for 1 layer. The first and
 * last years are not written - they serve as buffers to exclude unpredictable
 * behavior of the smoothing spline at the data edges. Disk space is allocated
 * upon file creation to avoid out-of-space errors that could arise during
 * processing. Outputs are assumed to be int16.
 ******************************************************************************/
FILE** get_output_files_quant(config pp) {
    
    FILE **files;
    char head_name[LONG_STR], file_name[LONG_STR];
    int yr, yr0, yr1,jj;
    int num_bands;

    yr0 = pp.start_year+1; // omit initial buffer year
    yr1 = pp.start_year+pp.num_years-2; // omit final buffer year

     // allocate FILE pointers for each year
    files = malloc((pp.num_years-2)*sizeof(FILE*));

     num_bands = (NUM_QUANT_BANDS*NUM_QUANTS) + 1;
     // create files and headers
     jj = 0;
     for (yr = yr0; yr <= yr1; yr++) {
          sprintf(head_name, "%s/%s.%s.%d.quant.hdr", pp.out_dir, pp.out_stem, pp.tile_id, yr);
          sprintf(file_name, "%s/%s.%s.%d.quant.bip", pp.out_dir, pp.out_stem, pp.tile_id, yr);
            
          create_envi_header(head_name, pp, sizeof(int16),num_bands); 
          files[jj] = create_file(file_name, MODIS_NROW*MODIS_NCOL*num_bands, sizeof(int16)); 
          jj++; // increment year index

          if (pp.verbose == 1) {
              printf("get_output_files: created %s\n", file_name);
          }
      }

    return files; 
}

/***************************************************************************//**
 * @brief Create header file used by ENVI to read the output BIP binary                                                                               
 ******************************************************************************/
void create_envi_header(const char *name, config pp, int num_bytes, int num_bands) {

    FILE *fp;
    double ulxmap, ulymap;

    ulxmap = MAP_X+PIX_LEN*MODIS_NCOL*pp.tile_column;
    ulymap = MAP_Y-PIX_LEN*MODIS_NROW*pp.tile_row;

    assert( (fp = fopen(name, "w")) != NULL);

    fprintf(fp, "ENVI\n");
    fprintf(fp, "samples = %d\n", MODIS_NCOL);
    fprintf(fp, "lines = %d\n", MODIS_NROW);
    fprintf(fp, "bands = %d\n", num_bands);
    fprintf(fp, "header offset = 0\n");
    fprintf(fp, "file type = ENVI Standard\n");
    fprintf(fp, "data type = %d\n", num_bytes); // will be 2 for a 16-bit signed integer
    fprintf(fp, "interleave = bip\n"); 
    fprintf(fp, "byte order = 0\n");
    fprintf(fp, "map info = {Sinusoidal, 1, 1, %12.9lf, %12.9lf, %12.9lf, %12.9lf, " 
                "units=Meters}\n", ulxmap, ulymap, PIX_LEN, PIX_LEN);
    fprintf(fp, "projection info = " PROJ "\n");

    assert( fclose(fp) != EOF );

    return;
}

/***************************************************************************//**
 * @brief Create an output file of a specified size
 * @param [in] path Name of file to be created, including path
 * @param [in] num_elem Number of array elements to be written
 * @param [in] elem_bytes Bytes per element to be written
 *
 * Creates a new file, seeks to the desired endpoint and writes a junk
 * character to allocate, then closes and reopens to return a fresh file pointer
 ******************************************************************************/
FILE* create_file(const char *path, size_t num_elem, int elem_bytes) {
    
    FILE *fp;

    assert( (fp = fopen(path, "wb")) != NULL );
    assert( (fseek(fp, num_elem*elem_bytes-1, SEEK_SET)) == 0 );
    assert(  fputc(0, fp) == 0 );
    rewind(fp);

    return fp;
}

/***************************************************************************//**
 * @brief Write output chunk to yearly output files
 * @param [in] pp Program configuration options
 * @param [in] chunk_num_rows Number of rows in this chunk, may be truncated
 * @param [in] data 2D array of data to be written to files
 * @param [in, out] files 2D array of file pointers, which are modified by fwrite calls
 * @param [in] kind flag indicating file type. Valid options are 0 = smoothed
 *     data layers,  1 = residual data layers. 2 = 
 *
 * Write output chunk to yearly output files. As in get_output_files(), the
 * first and last years are skipped.
 ******************************************************************************/
void write_chunk(config pp, size_t chunk_num_rows, int16 **data, FILE ***files, int kind) { 

    int ii, year;
    size_t pixel, idx, sz;

    // select element size
    if        (kind == 0) {
        sz = sizeof(int16);
    } else if (kind == 1) {
        sz = sizeof(int16);
    } else if (kind == 2) {
	sz = sizeof(int16);
     //   sz = sizeof(uint8);
    } else {
        printf("write_chunk: Invalid value for kind, exiting.\n");
        exit(EXIT_FAILURE);
    }

    // loop over layers
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {

        // skip if layer is not requested
        if (kind == 0 && pp.smooth[ii] == 0) continue; 
        if (kind == 1 && pp.residual[ii] == 0) continue; 
        if (kind == 2 && pp.prefilter[ii] == 0) continue; 

	// ## dsm added this - only want to output EVI smoothed splines and residuals
	if (ii != EVI2) continue;

        if (pp.verbose == 1) {
            printf("write_chunk: layer %d\n", ii+1);
        }

        for (pixel = 0; pixel < (chunk_num_rows*MODIS_NCOL); pixel++) {
            for(year = 1; year < (pp.num_years-1); year++) {
                idx = pixel*MODIS_PERIOD*pp.num_years + year*MODIS_PERIOD;
                assert( fwrite((&data[ii][idx]), sz, MODIS_PERIOD, files[ii][year-1]) == MODIS_PERIOD);
            }
        }
    }

    return;
}

/***************************************************************************//**
 * @brief Close output files and deallocate list
 * @param [in] pp Program configuration options
 * @param [in] kind flag indicating file type. Valid options are 0 = smoothed
 *     data layers,  1 = residual data layers
 * @param [in,out] files List of filenames to close and deallocate
 ******************************************************************************/
void close_output_files(config pp, int kind, FILE ***files) {

    int ii, jj;
    
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {

        // skip if layer is not requested
        if (kind == 0 && pp.smooth[ii] == 0) continue; 
        if (kind == 1 && pp.residual[ii] == 0) continue; 
	
	// ## dsm added this - only want to output EVI smoothed splines and residuals
	if (ii != EVI2) continue;

        // close files and deallocate arrays
        for (jj = 0; jj < pp.num_years-2; jj++) {
            assert( fclose(files[ii][jj]) != EOF );
        }
        free(files[ii]);
    }
    free(files);

    return;
}

/***************************************************************************//**
 * @brief Close output files and deallocate list
 * @param [in] pp Program configuration options
 * @param [in,out] files List of filenames to close and deallocate
 ******************************************************************************/
void close_output_files_quant(config pp, FILE **files) {

    int jj;

    // close files and deallocate arrays
    for (jj = 0; jj < pp.num_years-2; jj++) {
          assert( fclose(files[jj]) != EOF );
     }
       
    free(files);

    return;
}
