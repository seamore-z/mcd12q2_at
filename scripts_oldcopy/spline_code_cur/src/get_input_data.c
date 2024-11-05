/***************************************************************************//**
 * @file
 * @brief Functions for reading input data from MODIS HDF stack
 ******************************************************************************/

#include "modis_smooth_t.h"

// local macros
#define QA_MASK 15 /**<Bit mask used to extract 4-bit QA flags from packed data*/
#define MSG_CORRUPT_FILE "Corrupt input file - using fill value and continuing\n" /**<Warning message for corrupt input files*/

// local prototypes
static void read_step_nbar(int, char*, char*, int, size_t, int16*);
#ifdef C5
static void read_step_qa_c5(int, char*, char*, int, size_t, int, uint8*);
#else
static void read_step_qa_c6(int, char*, char*, int, size_t, uint8*);
#endif

/***************************************************************************//**
 * @brief Read a specific chunk from the MODIS HDF stack
 * @param [in] pp Program configuration options
 * @param [in] nbar_loc Location and info for NBAR input files
 * @param [in] qa_loc Location and info for QA input files
 * @param [in] row0 index of initial row to read in for this chunk
 * @param [out] nbar NBAR data for this chunk
 * @param [out] qa QA data for this chunk
 *
 * Reads in "chunks" of NBAR and QA data from the MODIS HDF stack. Chunks
 * include a fixed number of rows, all columns, and all timesteps.
 ******************************************************************************/
void get_chunk(config pp, data_loc nbar_loc, data_loc qa_loc, int row0, 
        int16 **nbar, uint8 **qa) {

    int bb;
    size_t step, elem, ii;
    int16 *nbar_step;
    uint8 *qa_step;

    assert( (nbar_step = malloc(pp.num_elem_step*sizeof(int16))) != NULL);
    assert( (qa_step   = malloc(pp.num_elem_step*sizeof(uint8))) != NULL);

    for (bb = 0; bb < MODIS_NUM_BANDS; bb++) {

        if (pp.read[bb] == 0) continue; // skip if band is not requested
        
        for (step = 0; step < pp.num_step; step++) {

            if (pp.verbose == 1) {
                printf("get_chunk: band %d, file %zu of %zu\n", bb+1, step+1, pp.num_step);
            }

            // read NBAR and QA step from file, or use fill
            read_step_nbar(nbar_loc.is_fill[step], nbar_loc.names[step], pp.nbar_sds_name[bb], row0, pp.chunk_num_row, nbar_step);
            #ifdef C5
            read_step_qa_c5(qa_loc.is_fill[step], qa_loc.names[step], pp.qa_sds_name[bb], row0, pp.chunk_num_row, bb+1, qa_step);
            #else
            read_step_qa_c6(qa_loc.is_fill[step], qa_loc.names[step], pp.qa_sds_name[bb], row0, pp.chunk_num_row, qa_step);
            #endif
            
            // copy step to main arrays
            for (elem = 0; elem < pp.num_elem_step; elem++) {
                ii = step + elem*pp.num_step;
                nbar[bb][ii] = nbar_step[elem];
                qa[bb][ii] = qa_step[elem];
            }
        }
    }

    free(nbar_step);
    free(qa_step);

    return;
}

/***************************************************************************//**
 * @brief Read a block of NBAR data from a single MODIS HDF file 
 * @param [in] is_fill Flag indicating if input file is fill (1) or not (0)
 * @param [in] name Filename of HDF file to be read
 * @param [in] sds_name Standard subdataset name for this band in the HDF file 
 * @param [in] row0 Initial row in chunk to be read
 * @param [in] num_row Number of rows in chunk to be read
 * @param [out] data Data read from file, returned by reference
 *
 * Data will be read from file if is_fill == 0, otherwise the output array is
 * set to a constant fill value. The function handles read failures by printing
 * a warning message and using the fill value.
 ******************************************************************************/
static void read_step_nbar(int is_fill, char *name, char* sds_name, int row0, size_t num_row, int16 *data) {

    int32 sd_id, sds_index, sds_id, start[2], edge[2];
    size_t ii;
    int read_fail;

    read_fail = 0;
 
    if (is_fill == 0) {
        // attempt to read data from file

        // open file and band 
        sd_id = SDstart(name, DFACC_READ);
        sds_index = SDnametoindex(sd_id, sds_name); // get band by name
        sds_id = SDselect(sd_id, sds_index); 

        if ((sd_id == FAIL) || (sds_index == FAIL) || (sds_id == FAIL)) {
            printf("read_step_nbar: " MSG_CORRUPT_FILE);
            printf("read_step_nbar: filename = %s\n", name);
            read_fail = 1;
        }

        if (read_fail == 0) {
            // read data
            start[0] = row0;
            edge[0]  = num_row;
            start[1] = 0;
            edge[1]  = MODIS_NCOL;
            assert(SDreaddata(sds_id, start, NULL, edge, (VOIDP)&(*data)) != FAIL);

            // close file and band
            assert(SDendaccess(sds_id) != FAIL);
            assert(SDend(sd_id) != FAIL);
        }
    }

    if ((is_fill == 1) || (read_fail == 1)) {
        // use fill value
        for (ii = 0; ii < MODIS_NCOL*num_row; ii++) {
            data[ii] = FILL_NBAR;
        }
    }
    
    return;
}

#ifdef C5
/***************************************************************************//**
 * @brief Read a block of QA data from a single MODIS C5 HDF file 
 * @param [in] is_fill Flag indicating if input file is fill (1) or not (0)
 * @param [in] name Filename of HDF file to be read
 * @param [in] sds_name Standard subdataset name for this band in the HDF file 
 * @param [in] row0 Initial row in chunk to be read
 * @param [in] num_row Number of rows in chunk to be read
 * @param [in] band_num ID for the band to be read, expected range is 1-7
 * @param [out] data Data read from file, returned by reference
 *
 * Data will be read from file if is_fill == 0, otherwise the output array is
 * set to a constant fill value. The function handles read failures by printing
 * a warning message and using the fill value.
 *
 * C5 QA data is provided as 4-bit integers for all 7 bands packed into a single
 * 32-bit integer. This function uses the band number to extract the correct
 * bits by shifting and masking.
 ******************************************************************************/
static void read_step_qa_c5(int is_fill, char *name, char *sds_name, int row0, size_t num_row, int band_num, uint8 *data) {

    size_t ii, num_elem;
    int32 sd_id, sds_index, sds_id, start[2], edge[2];
    uint32 *raw_data, tmp;
    int read_fail;

    read_fail = 0;

    if (is_fill == 0) {
        // attempt to read data from file

        // open file and band 
        sd_id = SDstart(name, DFACC_READ);
        sds_index = SDnametoindex(sd_id, sds_name); // get band by name
        sds_id = SDselect(sd_id, sds_index); 

        if ((sd_id == FAIL) || (sds_index == FAIL) || (sds_id == FAIL)) {
            printf("read_step_qa_c6: " MSG_CORRUPT_FILE);
            printf("read_step_qa_c6: filename = %s\n", name);
            read_fail = 1;
        }

        if (read_fail == 0) {
            // allocate space for packed QA data
            num_elem = MODIS_NCOL*num_row; // elements in one timestep of chunk
            assert( (raw_data = malloc(num_elem*sizeof(uint32))) != NULL);

            // read raw data
            start[0] = row0;
            edge[0]  = num_row;
            start[1] = 0;
            edge[1]  = MODIS_NCOL;
            assert(SDreaddata(sds_id, start, NULL, edge, (VOIDP)&(*raw_data)) != FAIL);

            // close file and band
            assert(SDendaccess(sds_id) != FAIL);
            assert(SDend(sd_id) != FAIL);

            // extract flag for the specified band
            band_num -= 1; // band_num is 1-indexed, file is 0-indexed
            for (ii = 0; ii < num_elem; ii++) {
                tmp = raw_data[ii];
                tmp = tmp>>(4*band_num); 
                data[ii] =(uint8)(tmp&QA_MASK);
            }

            // cleanup memory
            free(raw_data);
        }
    }

    if ((is_fill == 1) || (read_fail == 1)) {
        // use fill value
        for (ii = 0; ii < MODIS_NCOL*num_row; ii++) {
            data[ii] = FILL_QA;
        }
    }

    return;
}

#else

/***************************************************************************//**
 * @brief Read a block of QA data from a single MODIS C6 HDF file 
 * @param [in] is_fill Flag indicating if input file is fill (1) or not (0)
 * @param [in] name Filename of HDF file to be read
 * @param [in] sds_name Standard subdataset name for this band in the HDF file 
 * @param [in] row0 Initial row in chunk to be read
 * @param [in] num_row Number of rows in chunk to be read
 * @param [out] data Data read from file, returned by reference
 *
 * Data will be read from file if is_fill == 0, otherwise the output array is
 * set to a constant fill value. The function handles read failures by printing
 * a warning message and using the fill value.
 ******************************************************************************/
static void read_step_qa_c6(int is_fill, char *name, char *sds_name, int row0, size_t num_row, uint8 *data) {

    int32 sd_id, sds_index, sds_id, start[2], edge[2];
    size_t ii, read_fail;

    read_fail = 0;
 
    if (is_fill == 0) {
        // attempt to read data from file

        // open file
        sd_id = SDstart(name, DFACC_READ);
        sds_index = SDnametoindex(sd_id, sds_name); // get band by name
        sds_id = SDselect(sd_id, sds_index); 

        if ((sd_id == FAIL) || (sds_index == FAIL) || (sds_id == FAIL)) {
            printf("read_step_qa_c6: " MSG_CORRUPT_FILE);
            printf("read_step_qa_c6: filename = %s\n", name);
            read_fail = 1;
        }

        if (read_fail == 0) {
            // read data
            start[0] = row0;
            edge[0]  = num_row;
            start[1] = 0;
            edge[1]  = MODIS_NCOL;
            assert(SDreaddata(sds_id, start, NULL, edge, (VOIDP)&(*data)) != FAIL);

            // close file and band
            assert(SDendaccess(sds_id) != FAIL);
            assert(SDend(sd_id) != FAIL);
        }
    }

    if ((is_fill == 1) || (read_fail == 1)) {
        // use fill value
        for (ii = 0; ii < MODIS_NCOL*num_row; ii++) {
            data[ii] = FILL_QA;
        }
    }
    
    return;
}

#endif
