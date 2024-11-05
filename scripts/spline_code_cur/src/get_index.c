/***************************************************************************//**
 * @file
 * @brief Functions for computing indices derived from MODIS bands 
 ******************************************************************************/

#include "modis_smooth_t.h"

//local prototypes
static void get_evi2(size_t, int16*, int16*, int16*, uint8*, uint8*, uint8*);
static void get_ndsi(size_t, int16*, int16*, int16*, uint8*, uint8*, uint8*);
static void get_ndsi_sm(size_t, int16*, int16*, int16*);
static void get_ndwi_sm(size_t, int16*, int16*, int16*);
static void get_ndii1_sm(size_t, int16*, int16*, int16*);
static void get_ndii2_sm(size_t, int16*, int16*, int16*);

/***************************************************************************//**
 * @brief Compute requested indices before spline smooth
 ******************************************************************************/
void get_indices(config pp, int16 **data, uint8 **qa) {

    if(pp.read[EVI2] == 1) { 
        get_evi2(pp.num_elem_chunk, data[NIR], data[RED], data[EVI2], \
                qa[NIR], qa[RED], qa[EVI2]);
    }

    if(pp.read[NDSI] == 1) { 
        get_ndsi(pp.num_elem_chunk, data[GREEN], data[SWIR2], data[NDSI], \
                qa[GREEN], qa[SWIR2], qa[NDSI]);
    }

    return;
}

/***************************************************************************//**
 * @brief Compute requested indices after spline smooth
 ******************************************************************************/
void get_indices_spl(config pp, int16 **data_in, int16 **data_out) {

    get_ndsi_sm(pp.num_elem_chunk, data_in[GREEN], data_in[SWIR2], data_out[0]);
    get_ndwi_sm(pp.num_elem_chunk, data_in[NIR], data_in[SWIR1], data_out[1]);
    get_ndii1_sm(pp.num_elem_chunk, data_in[NIR], data_in[SWIR2], data_out[2]);
    get_ndii2_sm(pp.num_elem_chunk, data_in[NIR], data_in[SWIR3], data_out[3]);

    return;
}



/***************************************************************************//**
 * @brief Compute EVI2 index from MODIS NIR and red bands, skip fill values
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] nir MODIS near infrared band data array
 * @param [in] red MODIS red band data array
 * @param [out] evi2 Computed index value
 * @param [in] nir_qa QA flags for MODIS near infrared band data array
 * @param [in] red_qa QA flags forMODIS red band data array
 * @param [out] evi2_qa QA flags for computed index value
 ******************************************************************************/
static void get_evi2(size_t num_elem, int16 *nir, int16 *red, int16 *evi2, \
              uint8 *nir_qa, uint8 *red_qa, uint8 *evi2_qa) {

    size_t ii;

    for (ii = 0; ii < num_elem; ii++) {
        if ((nir[ii] == FILL_NBAR) || (red[ii] == FILL_NBAR)) {
            evi2[ii] = FILL_NBAR;
            evi2_qa[ii] = FILL_QA;
        } else {
            evi2[ii] = 10000*2.5*(nir[ii] - red[ii])/(nir[ii] + (2.4*red[ii]) + 10000); // red and nir are >= 0, so denom is != 0
            evi2_qa[ii] =  (nir_qa[ii] > red_qa[ii]) ? nir_qa[ii] : red_qa[ii]; // use max (worst) QA flag
        }
    }

    return;
}

/***************************************************************************//**
 * @brief Compute NDSI index from MODIS bands, skip fill values
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] green MODIS green band data array
 * @param [in] swir2 MODIS shortwave infrared 2 band data array
 * @param [out] ndsi Computed index value
 * @param [in] green_qa QA flags for MODIS green band data array
 * @param [in] swir2_qa QA flags for MODIS shortwave infrared 2 band data array
 * @param [out] ndsi_qa QA flags for computed index value
 ******************************************************************************/
static void get_ndsi(size_t num_elem, int16 *green, int16 *swir2, int16 *ndsi, \
        uint8 *green_qa, uint8 *swir2_qa, uint8 *ndsi_qa) {

    size_t ii;
    
    for (ii = 0; ii < num_elem; ii++) {
        if ((green[ii] == FILL_NBAR) || (swir2[ii] == FILL_NBAR) || ((green[ii]+swir2[ii])==0)) {
            ndsi[ii] = FILL_NBAR;
            ndsi_qa[ii] = FILL_QA;
        } else {
	    ndsi[ii] = 10000*(double)(green[ii]-swir2[ii]) / (double)(green[ii]+swir2[ii]); 
            ndsi_qa[ii] =  (green_qa[ii] > swir2_qa[ii]) ? green_qa[ii] : swir2_qa[ii]; // use max (worst) QA flag
        }
    }

    return;
}

/***************************************************************************//**
 * @brief Compute NDSI index from MODIS bands after smoothing
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] green MODIS green band data array
 * @param [in] swir2 MODIS shortwave infrared 2 band data array
 * @param [out] out Computed index value
 ******************************************************************************/
static void get_ndsi_sm(size_t num_elem, int16 *green, int16 *swir2, int16 *out) {

    size_t ii;

    for (ii = 0; ii < num_elem; ii++) {
	if((green[ii]+swir2[ii])!=0) {
	   	out[ii] = 10000*(double)(green[ii]-swir2[ii]) / (double)(green[ii]+swir2[ii]); 
	  } else {
		out[ii] = 0;
	}
    }

    return;
}

/***************************************************************************//**
 * @brief Compute NDWI index from MODIS bands after smoothing
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] nir MODIS nir band data array
 * @param [in] swir1 MODIS shortwave infrared 1 band data array
 * @param [out] out Computed index value
 ******************************************************************************/
static void get_ndwi_sm(size_t num_elem, int16 *nir, int16 *swir1, int16 *out) {

    size_t ii;

    for (ii = 0; ii < num_elem; ii++) {
	  if((nir[ii]+swir1[ii])!=0) {
	   	out[ii] = 10000*(double)(nir[ii]-swir1[ii]) / (double)(nir[ii]+swir1[ii]); 
	  } else {
		out[ii] = 0;
	}
    }

    return;
}

/***************************************************************************//**
 * @brief Compute NDII1 index from MODIS bands after smoothing
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] nir MODIS nir band data array
 * @param [in] swir2 MODIS shortwave infrared 2 band data array
 * @param [out] out Computed index value
 ******************************************************************************/
static void get_ndii1_sm(size_t num_elem, int16 *nir, int16 *swir2, int16 *out) {

    size_t ii;

    for (ii = 0; ii < num_elem; ii++) {
	if((nir[ii]+swir2[ii])!=0) {
	   	out[ii] = 10000*(double)(nir[ii]-swir2[ii]) / (double)(nir[ii]+swir2[ii]);
	  } else {
		out[ii] = 0;
	}
    }

    return;
}

/***************************************************************************//**
 * @brief Compute NDII2 index from MODIS bands after smoothing
 * @param [in] num_elem Number of elements in input and output arrays
 * @param [in] nir MODIS nir band data array
 * @param [in] swir3 MODIS shortwave infrared 3 band data array
 * @param [out] out Computed index value
 ******************************************************************************/
static void get_ndii2_sm(size_t num_elem, int16 *nir, int16 *swir3, int16 *out) {

    size_t ii;

    for (ii = 0; ii < num_elem; ii++) {
	   if((nir[ii]+swir3[ii])!=0) {
	   	out[ii] = 10000*(double)(nir[ii]-swir3[ii]) / (double)(nir[ii]+swir3[ii]);
	  } else {
		out[ii] = 0;
	}
    }

    return;
}

