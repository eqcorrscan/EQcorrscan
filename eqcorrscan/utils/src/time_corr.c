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

// Prototypes

int multi_normxcorr_time(float*, long, long, float*, long, float*, int);

// Code

int multi_normxcorr_time(float *templates, long template_len, long n_templates,
                         float *image, long image_len, float *ccc, int num_threads){
	int i, k, p;
	int steps = image_len - template_len + 1;
	double new_samp, old_samp, sum, numerator = 0.0, denom, auto_a;
	double *mean, *auto_b;

    // Pre-compute the running mean
    mean = (double*) malloc(steps * sizeof(double));
    sum = 0.0;
    for (k=0; k < template_len; ++k){
        sum += (double) image[k];
    }
    mean[0] = sum / template_len;

    for(k = 1; k < steps; ++k){
        new_samp = (double) image[k + template_len - 1];
        old_samp = (double) image[k - 1];
        mean[k] = mean[k - 1] + (new_samp - old_samp) / template_len;
    }

    // Pre-compute the image auto-correlation
    auto_b = (double*) malloc(steps * sizeof(double));
    #pragma omp parallel for num_threads(num_threads) private(p)
    for (k = 0; k < steps; ++k){
        double _auto_b = 0.0;
        for(p = 0; p < template_len; ++p){
            _auto_b += ((double) image[p + k] - mean[k]) * ((double) image[p + k] - mean[k]);
        }
        auto_b[k] = _auto_b;
    }

	for (i = 0; i < n_templates; ++i){
	    auto_a = 0.0;
	    for (p = 0; p < template_len; ++p){
	    	auto_a += (double) templates[(template_len * i) + p] * (double) templates[(template_len * i) + p];
    	}
	    #pragma omp parallel for num_threads(num_threads) private(numerator, denom, p)
	    for(k = 0; k < steps; ++k){
	        numerator = 0.0;
		    for(p = 0; p < template_len; ++p){
			    numerator += (double) templates[(template_len * i) + p] * ((double) image[p + k] - mean[k]);
		    }
		    denom = sqrt(auto_a * auto_b[k]);
		    ccc[((image_len - template_len + 1) * i) + k] = (float) (numerator / denom);
	    }
	}
	free(mean);
	free(auto_b);
	return 0;
}
