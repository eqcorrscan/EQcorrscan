/*
 * =====================================================================================
 *
 *       Filename:  multi_corr.c
 *
 *        Purpose:  Routines for computing cross-correlations in the time-domain
 *
 *        Created:  03/07/17 02:25:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (Calum Chamberlain),
 *   Organization:  EQcorrscan
 *      Copyright:  EQcorrscan developers.
 *        License:  GNU Lesser General Public License, Version 3
 *                  (https://www.gnu.org/copyleft/lesser.html)
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__linux__) || defined(__linux) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    #include <omp.h>
    #ifndef N_THREADS
        #define N_THREADS omp_get_max_threads()
    #endif
#endif

int multi_normxcorr_time(float*, int, int, float*, int, float*, int);

int multi_normxcorr_time(float *templates, int template_len, int n_templates, float *image,
                         int image_len, float *ccc, int num_threads){
	int i, p, k;
	int steps = image_len - template_len + 1;
	double * running_mean = (double *) calloc(steps, sizeof(double));
	double * auto_b = (double *) calloc(steps, sizeof(double));


    // Pre-compute the running mean
	running_mean[0] = 0;
    for (i = 0; i < template_len; ++i){
		running_mean[0] += image[i];
	}
	running_mean[0] = running_mean[0] / template_len;

    for(i = 1; i < steps; ++i){
        running_mean[i] = running_mean[i - 1] + (image[i + template_len - 1] - image[i - 1]) / template_len;
    }

    // Pre-compute the image autocorrelation
    #pragma omp parallel for num_threads(num_threads)
    for(i = 0; i < steps; ++i){
        double _auto_b = 0.0;
        for(p = 0; p < template_len; ++p){
            _auto_b += ((double) image[p + k] - running_mean[k]) * ((double) image[p + k] - running_mean[k]);
        }
        auto_b[i] = _auto_b;
    }

	for (i = 0; i < n_templates; ++i){
	    double auto_a = 0.0;
	    #pragma omp parallel for reduction(+:auto_a) num_threads(num_threads)
    	for(p = 0; p < template_len; ++p){
	    	auto_a += (double) templates[(template_len * i) + p] * (double) templates[(template_len * i) + p];
    	}

        #pragma omp parallel for num_threads(num_threads)
    	for(k = 0; k < steps; ++k){
	        double numerator = 0.0, denom;
    		for(p = 0; p < template_len; ++p){
	    		numerator += (double) templates[(template_len * i) + p] * ((double) image[p + k] - running_mean[k]);
    		}
    		denom = sqrt(auto_a * auto_b[k]);
		    ccc[((image_len - template_len + 1) * i) + k] = (float) (numerator / denom);
    	}
	}
	free(running_mean);
	free(auto_b);
	return 0;
}

