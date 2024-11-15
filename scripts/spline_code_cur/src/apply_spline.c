/***************************************************************************//**
 * @file
 * @brief Functions to setup and run splining routine
 ******************************************************************************/

#include "modis_smooth_t.h"
#include <math.h>

// local macros
#define SPLINE_MIN_PTS_PER_YEAR 10 /**<Threshold valid points-per-year below which fill value is used */
#define SPLINE_LOW_PTS_PER_YEAR 25 /**<Threshold valid points-per-year below which mean value is used */

#define WEIGHT_QA_0 1.0     /**<Weight for QA value of 0 */
#define WEIGHT_QA_1 0.9     /**<Weight for QA value of 1 */
#define WEIGHT_QA_2 0.6     /**<Weight for QA value of 2 */
#define WEIGHT_QA_3 0.5     /**<Weight for QA value of 3 */
#define WEIGHT_QA_4 0.0     /**<Weight for QA value of 4 */
#define WEIGHT_NONE 0.0     /**<Weight for ignored points */
#define QA_MIN_INCLUDED 3   /**<Minimum QA value with a non-zero weight */

// local prototypes
int get_num_knots(int);
int count_valid_data(int16*, uint8*, int);
double get_weight(int16, uint8);
void use_fill(int16*, uint8*, int, int);
void get_spline_inputs(int16*, uint8*, int, double*, double*, double*);
void get_knots(double*, double*, int, int);
void write_pix(int, int16*, FILE*);

/***************************************************************************//**
 * @brief Main function to setup and run splining routine                                                                               
 * @param [in] pp Program configuration options
 * @param [in, out] data 2D array of data to be splined and returned in place
 * @param [in] qa 2D array of QA flags corresponding to layers, used to construct spline weights 
 * @param [out] residual 2D array of residual values (smoothed-original)
 * @param [in] pixel ids for that row chunk to write out to an output file
 * @param [out] file for writing pre and post splined data for specific pixels
 ******************************************************************************/
void apply_spline(config pp, int16 **data, uint8 **qa, int16 **residual) {

    size_t ii, jj, kk, num_pixels, num_steps, num_knots, num_valid, start;
    int ld4, ispar, isetup;
    double *knot, *coef, *xwy, *hs0, *hs1, *hs2, *hs3, *sg0, *sg1, *sg2, *sg3, \
        *abd, *xs, *ys, *ws, *sz, spar, lspar;

    num_pixels = pp.chunk_num_row*MODIS_NCOL; // number of pixels in chunk
    num_steps = pp.num_years*MODIS_PERIOD; // num time steps in each pixel time series
    // ## number of knots is related to how smooth the output spline will be
    num_knots = get_num_knots(num_steps)+2;

    // constant parameters
    ispar = 1;
    lspar = -1.5;
     // ## very important smoothing parameter - larger number means more smoothing
    spar = 0.4;
    ld4 = 4;

    // main loop, spline each pixel in each layer 
    for (ii = 0; ii < NUM_DATA_LAYERS; ii++) {

        if (pp.smooth[ii] == 0) continue; // skip if layer is not requested

        if (pp.verbose == 1) {
            if (ii <= SWIR3) {
                printf("apply_spline: band %zu\n", ii+1);
            } else if (ii == EVI2) {
                printf("apply_spline: EVI2\n");
            } else if (ii == NDSI) {
                printf("apply_spline: NDSI\n");
            }
        }
	// ## start parallel processing using OpenMP
        #pragma omp parallel default(none) \
        shared(num_pixels, num_steps, num_knots, ispar, lspar, spar, ld4, ii, \
                pp, data, qa, residual) \
        private(knot, coef, xwy, hs0, hs1, hs2, hs3, sg0, sg1, sg2, sg3, abd, \
                xs, ys, ws, sz, start, num_valid, isetup, kk)
        {

            // verbose output from OpenMP
            if (pp.verbose == 1) {
                printf("apply_spline: layer=%zu, num threads=%d, thread id=%d\n", \
                        ii+1, omp_get_num_threads(), omp_get_thread_num());
            }

            // allocate local arrays
            knot = malloc((num_knots+4)  * sizeof(double));
            coef = malloc( num_knots     * sizeof(double)); 
            xwy  = malloc( num_knots     * sizeof(double));
            hs0  = malloc( num_knots     * sizeof(double));
            hs1  = malloc( num_knots     * sizeof(double));
            hs2  = malloc( num_knots     * sizeof(double));
            hs3  = malloc( num_knots     * sizeof(double));
            sg0  = malloc( num_knots     * sizeof(double));
            sg1  = malloc( num_knots     * sizeof(double));
            sg2  = malloc( num_knots     * sizeof(double));
            sg3  = malloc( num_knots     * sizeof(double));
            abd  = malloc( num_knots*ld4 * sizeof(double));
            xs = malloc(num_steps * sizeof(double));
            ys = malloc(num_steps * sizeof(double));
            ws = malloc(num_steps * sizeof(double));
            sz = malloc(num_steps * sizeof(double));

            #pragma omp for schedule(static)
            for (jj = 0; jj < num_pixels; jj++) {

                start = jj*num_steps; // ## index of first step for pixel jj in data array 

                // ## check number of valid data points, if insufficient data, just fill
                num_valid = count_valid_data(&data[ii][start], &qa[ii][start], num_steps); 

                if (num_valid <= (size_t)(SPLINE_MIN_PTS_PER_YEAR*pp.num_years)) {
                    use_fill(&data[ii][start], &qa[ii][start], num_steps, num_valid);
                    continue;
                }

                // ## get spline inputs 
                get_spline_inputs(&data[ii][start], &qa[ii][start], num_steps, xs, ys, ws);
                get_knots(knot, xs, num_knots, num_steps);

                // ## main splining function
                isetup = 0; // ## reset
		// ## the sz array is the output smooth spline
                assert( sbart(xs, ys, ws, num_steps, knot, num_knots, coef, sz, spar,
                            ispar, &lspar, &isetup, xwy, hs0, hs1, hs2, hs3, sg0, sg1,
                            sg2, sg3, abd, ld4) == 0);

		// ## make sure output values are not below 0 or above 10000
		// ## added dsm 6/30/16
                for (kk = 0; kk < num_steps; kk++) {
			if(sz[kk] != FILL_NBAR) {
				// ## in the case of EVI2 check relative to global minimum default
				if (ii == EVI2 && sz[kk] < EVI2_FLOOR_DEFAULT) {
					sz[kk] = EVI2_FLOOR_DEFAULT;
				}
				if(ii!= EVI2 && ii!=NDSI && sz[kk] < 0) {
					sz[kk] = 0;
				}
				if(sz[kk] > 10000) {
					sz[kk] = 10000;
				}
			}
                } // ## kk loop

                // compute residuals where requested as (spline - original_value)
                if (pp.residual[ii] == 1) {
                    for (kk = 0; kk < num_steps; kk++) {
			if(data[ii][start+kk] != FILL_NBAR) {
                        	residual[ii][start+kk] = (int16)sz[kk] - data[ii][start+kk];
			} else {
				residual[ii][start+kk] = FILL_NBAR;
			}
                    } // ## kk loop
                } 

                // copy smooth results to data array, overwriting input data
                for (kk = 0; kk < num_steps; kk++) {
                    data[ii][start+kk] = (int16) sz[kk];
                } // kk

            } // jj

            // free local arrays
            free(knot); 
            free(coef); 
            free(xwy);  
            free(hs0);  
            free(hs1);  
            free(hs2);  
            free(hs3);  
            free(sg0);  
            free(sg1);  
            free(sg2);  
            free(sg3);  
            free(abd);  
            free(xs);
            free(ys);
            free(ws);
            free(sz);

        } // end omp parallel block
    } // ii

    return;
}

/***************************************************************************//**
 * @brief Populate knot locations for splining 
 * @param [out] knot Knot locations
 * @param [in] xs Spline routine input array
 * @param [in] num_knots Number of knots, defines length of knots array
 * @param [in] num_steps Number of time steps in pixel time series
 ******************************************************************************/
void get_knots(double *knot, double *xs, int num_knots, int num_steps) {

    int ii;

    knot[0] = xs[0];
    knot[1] = xs[0];
    knot[2] = xs[0];

    for (ii = 0; ii < num_knots-3; ii++) {
        knot[ii+3] = xs[ (ii*(num_steps-1)/(num_knots-3)) ];
    }

    knot[num_knots+0] = xs[num_steps-1];
    knot[num_knots+1] = xs[num_steps-1];
    knot[num_knots+2] = xs[num_steps-1];
    knot[num_knots+3] = xs[num_steps-1];

    return;
}

/***************************************************************************//**
 * @brief Populate input arrays for smoothing spline routine
 * @param [in] nbar NBAR data for the pixel time series 
 * @param [in] qa QA flags for the pixel time series 
 * @param [in] num_steps Number of points in pixel time series
 * @param [out] xs Spline input array
 * @param [out] ys Spline input array
 * @param [out] ws Spline input array
 ******************************************************************************/
void get_spline_inputs(int16 *nbar, uint8 *qa, int num_steps, double *xs, double *ys, double *ws) {

    int ii, pos_weights;
    double weight_sum, x_range;

    // initialize spline input arrays
    for (ii = 0; ii < num_steps; ii++) {
        xs[ii] = (double) ii;
        ys[ii] = (double) nbar[ii];
        ws[ii] = get_weight(nbar[ii], qa[ii]);
    }

    // normalize spline input arrays
    weight_sum = 0.0;
    pos_weights = 0;
    for (ii = 0; ii < num_steps; ii++) {
        weight_sum += ws[ii];
        if (ws[ii] > 0) pos_weights++;
    }

    x_range = xs[num_steps-1] - xs[0];
    for (ii = 0; ii < num_steps; ii++) {
        ws[ii] = ws[ii]*(double)pos_weights/weight_sum;
        xs[ii] = (xs[ii]-xs[0])/x_range;
        if (ws[ii] <= 0) ys[ii] = ws[ii]*ys[ii]; 
    }

    return;
}

/***************************************************************************//**
 * @brief Fill NBAR pixel time-series with constant weighted mean or no-data value 
 * @param [in, out] nbar NBAR pixel time-series  
 * @param [in] qa Quality flags for pixel time-series, used to get weights
 * @param [in] num_steps Number of steps in the time-series
 * @param [in] num_valid Number of valid points in the time-series
 ******************************************************************************/
void use_fill(int16 *nbar, uint8 *qa, int num_steps, int num_valid) {

    int ii;
    double weighted_sum, total_weight, weighted_mean, weight;

    if (num_valid == 0) {
        // no data, fill with no-data value

        for (ii = 0; ii < num_steps; ii++) {
            nbar[ii] = FILL_NBAR;
        }

    } else {
        // insufficient data, fill with weighted mean
        weighted_sum = 0.0;
        total_weight = 0.0;
        for (ii = 0; ii < num_steps; ii++) {
            weight = get_weight(nbar[ii], qa[ii]);
            weighted_sum += (double)nbar[ii]*weight;
            total_weight += weight;
        }
        weighted_mean = weighted_sum/total_weight;

        for (ii = 0; ii < num_steps; ii++) {
            nbar[ii] = weighted_mean;
        }
    }

    return;
}

/***************************************************************************//**
 * @brief Convert quality flag into weight using values defined in local macros
 * @param [in] nbar NBAR data for this data point
 * @param [in] qa Quality flag for this data point
 * @return Weight as a scalar double in the range [0, 1]
 ******************************************************************************/
double get_weight(int16 nbar, uint8 qa) {

    double w;

    if (nbar == FILL_NBAR) {
        w = WEIGHT_NONE;
    } else {
        switch(qa) {
            case 0:  w = WEIGHT_QA_0; break;
            case 1:  w = WEIGHT_QA_1; break; 
            case 2:  w = WEIGHT_QA_2; break;
            case 3:  w = WEIGHT_QA_3; break;
            case 4:  w = WEIGHT_QA_4; break;
            default: w = WEIGHT_NONE;
        }
    }

    return w;
}

/***************************************************************************//**
 * @brief Count the number of valid data points a pixel time series
 * @param [in] nbar_pixel NBAR data for a pixel time series
 * @param [in] qa_pixel QA data for a pixel time series
 * @param [in] num_steps Number of timesteps in dataset 
 * @return number of valid data points the pixel time series
 *
 * Counts the number of pixels in the NBAR time series that are not set to the
 * fill value and have non-zero weight.
 ******************************************************************************/
int count_valid_data(int16 *nbar_pixel, uint8* qa_pixel, int num_steps) {

     int ii, num_valid;

     num_valid = 0;
     for (ii = 0; ii<num_steps; ii++) {
         if ((nbar_pixel[ii] != FILL_NBAR) && (qa_pixel[ii] <= QA_MIN_INCLUDED)) {
             num_valid++;
         }
     }

     return num_valid;
}

/***************************************************************************//**
 * @brief Calculate the optimal number of knots in the spline
 * @param [in] n Number of points in the dataset to be splined
 * @return Optimal number of knots for the dataset
 *
 * Calculate the optimal number of knots based on the number of data points in
 * the spline. Code is taked from R package sources, no details on the
 * algorithm.
 ******************************************************************************/
int get_num_knots(int n) {

   int knots;
   if ( n < 50 ){
      return (n);

   } else {
      double a1 = log2(50);
      double a2 = log2(100);
      double a3 = log2(140);
      double a4 = log2(200);

      if ( n < 200 ){
         knots = (int) floor( pow( 2, a1 + ( a2 - a1 ) * ( n - 50 )/150.) );
      } else if ( n < 800 ) {
         knots = (int) floor( pow( 2, a2 + ( a3 - a2 ) * ( n - 200 )/600.) );
      } else if ( n < 3200 ) {
         knots = (int) floor( pow( 2, a3 + ( a4 - a3 ) * ( n - 800 )/2400. ) );
      } else {
         knots = (int) floor( pow( 200. + ( n - 3200 ), 0.2 ) );
      }
   }
   return (knots);
}
