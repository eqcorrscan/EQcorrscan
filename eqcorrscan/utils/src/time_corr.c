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

int normxcorr_time_threaded(float*, int, float*, int, float*, int);

int normxcorr_time(float*, int, float*, int, float*);

int multi_normxcorr_time(float*, int, int, float*, int, float*);

int multi_normxcorr_time_threaded(float*, int, int, float*, int, float*, int);

int normxcorr_time_threaded(float *template, int template_len, float *image, int image_len, float *ccc, int num_threads){
    // Time domain cross-correlation - requires zero-mean template
	int p, k;
	int steps = image_len - template_len + 1;
	double * mean = (double *) calloc(steps, sizeof(double));
	double numerator = 0.0, denom;
	double auto_a = 0.0, auto_b = 0.0;

    mean[0] = 0;
    for (k=0; k < template_len; ++k){
		mean[0] += image[k];
	}
	mean[0] = mean[0] / template_len;

    for(k = 1; k < steps; ++k){
        mean[k] = mean[k - 1] + (image[k + template_len - 1] - image[k - 1]) / template_len;
    }
	for(p = 0; p < template_len; ++p){
		auto_a += (double) template[p] * (double) template[p];
	}
	#pragma omp parallel for private(numerator, denom, auto_b, p) num_threads(num_threads)
	for(k = 0; k < steps; ++k){
	    numerator = 0.0;
	    auto_b = 0.0;
		for(p = 0; p < template_len; ++p){
			numerator += (double) template[p] * ((double) image[p + k] - mean[k]);
			auto_b += ((double) image[p + k] - mean[k]) * ((double) image[p + k] - mean[k]);
		}
		denom = sqrt(auto_a * auto_b);
		ccc[k] = (float) (numerator / denom);
	}
	free(mean);
	return 0;
}

int normxcorr_time(float *template, int template_len, float *image, int image_len, float *ccc){
    // Time domain cross-correlation - requires zero-mean template
	int p, k;
	int steps = image_len - template_len + 1;
	double * mean = (double *) calloc(steps, sizeof(double));
	double numerator = 0.0, denom;
	double auto_a = 0.0, auto_b = 0.0;

    mean[0] = 0;
    for (k=0; k < template_len; ++k){
		mean[0] += image[k];
	}
	mean[0] = mean[0] / template_len;

    for(k = 1; k < steps; ++k){
        mean[k] = mean[k - 1] + (image[k + template_len - 1] - image[k - 1]) / template_len;
    }
	for(p = 0; p < template_len; ++p){
		auto_a += (double) template[p] * (double) template[p];
	}
	for(k = 0; k < steps; ++k){
	    numerator = 0.0;
	    auto_b = 0.0;
		for(p = 0; p < template_len; ++p){
			numerator += (double) template[p] * ((double) image[p + k] - mean[k]);
			auto_b += ((double) image[p + k] - mean[k]) * ((double) image[p + k] - mean[k]);
		}
		denom = sqrt(auto_a * auto_b);
		ccc[k] = (float) (numerator / denom);
	}
	free(mean);
	return 0;
}

int multi_normxcorr_time(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc){
	int i;
	for (i = 0; i < n_templates; ++i){
		normxcorr_time(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len + 1) * i]);
	}
	return 0;
}

int multi_normxcorr_time_threaded(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc, int num_threads){
	int i;
	for (i = 0; i < n_templates; ++i){
		normxcorr_time_threaded(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len + 1) * i], num_threads);
	}
	return 0;
}
