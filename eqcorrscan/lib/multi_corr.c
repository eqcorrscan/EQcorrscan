/*
 * =====================================================================================
 *
 *       Filename:  multi_corr.c
 *
 *        Purpose:  Routines for computing cross-correlations
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
#include <fftw3.h>
#if defined(__linux__) || defined(__linux) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    #include <omp.h>
    #ifndef N_THREADS
        #define N_THREADS omp_get_max_threads()
    #endif
#endif

// Prototypes
int normxcorr_fftw(float*, int, int, float*, int, float*, int);

int normxcorr_fftw_threaded(float*, int, int, float*, int, float*, int);

int normxcorr_time(float*, int, float*, int, float*);

int multi_normxcorr_fftw(float*, int, int, int, float*, int, float*, int);

int multi_normxcorr_time(float*, int, int, float*, int, float*);

// Functions
int normxcorr_fftw_threaded(float *templates, int template_len, int n_templates,
					        float *image, int image_len, float *ncc, int fft_len){
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
	templates:      Template signals
	template_len:   Length of template
	n_templates:    Number of templates (n0)
	image:          Image signal (to scan through)
	image_len:      Length of image
	ncc:            Output for cross-correlation - should be pointer to memory -
					must be n_templates x image_len - template_len + 1
	fft_len:        Size for fft (n1)
  */
	int N2 = fft_len / 2 + 1;
	int i, t, startind;
	double mean, stdev, old_mean, new_samp, old_samp, c, var=0.0, sum=0.0, acceptedDiff = 0.0000001;
	double * norm_sums = (double *) calloc(n_templates, sizeof(double));
	double * template_ext = (double *) calloc(fft_len * n_templates, sizeof(double));
	double * image_ext = (double *) calloc(fft_len, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * fft_len * n_templates);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	// Initialize threads
	#ifdef N_THREADS
        fftw_init_threads();
	    fftw_plan_with_nthreads(N_THREADS);
	#endif
	// Plan

	fftw_plan pa = fftw_plan_dft_r2c_2d(n_templates, fft_len, template_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_2d(n_templates, fft_len, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (t = 0; t < n_templates; ++t){
		for (i = 0; i < template_len; ++i)
		{
			template_ext[(t * fft_len) + i] = (double) templates[((t + 1) * template_len) - (i + 1)];
			norm_sums[t] += templates[(t * template_len) + i];
		}
	}
	for (i = 0; i < image_len; ++i)
	{
		image_ext[i] = (double) image[i];
	}
	//  Compute ffts of template and image
	#pragma omp parallel sections
	{
	    {fftw_execute(pa); }
	    #pragma omp section
	    {fftw_execute(pb); }
	}
	//  Compute dot product
	for (t = 0; t < n_templates; ++t){
    	for (i = 0; i < N2; ++i)
	    {
		    out[(t * N2) + i][0] = outa[(t * N2) + i][0] * outb[i][0] - outa[(t * N2) + i][1] * outb[i][1];
    		out[(t * N2) + i][1] = outa[(t * N2) + i][0] * outb[i][1] + outa[(t * N2) + i][1] * outb[i][0];
    	}
    }
	//  Compute inverse fft
	fftw_execute(px);
	//  Procedures for normalisation
	// Compute starting mean, will update this
	for (i=0; i < template_len; ++i){
		sum += image[i];
	}
	mean = sum / template_len;

	// Compute starting standard deviation
	for (i=0; i < template_len; ++i){
		var += pow(image[i] - mean, 2) / (template_len);
	}
	stdev = sqrt(var);
    // Used for centering - taking only the valid part of the cross-correlation
	startind = template_len - 1;
	for (t = 0; t < n_templates; ++t){
    	if (var < acceptedDiff){
	    	ncc[t * (image_len - template_len + 1)] = 0;
    	}
	    else {
		    c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean) / stdev;
    		ncc[t * (image_len - template_len + 1)] = (float) c;
	    }
	}
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (image_len - template_len + 1); ++i){
		// Need to cast to double otherwise we end up with annoying floating
		// point errors when the variance is massive - collecting fp errors.
		new_samp = image[i + template_len - 1];
		old_samp = image[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / template_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (template_len);
		stdev = sqrt(var);
		for (t=0; t < n_templates; ++t){
			if (var > acceptedDiff){
				c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean ) / stdev;
				ncc[(t * (image_len - template_len + 1)) + i] = (float) c;
			}
			else{
				ncc[(t * (image_len - template_len + 1)) + i] = 0.0;
			}
		}
	}
	//  Clean up
	fftw_destroy_plan(pa);
	fftw_destroy_plan(pb);
	fftw_destroy_plan(px);

	fftw_free(out);
	fftw_free(outa);
	fftw_free(outb);
	fftw_free(ccc);

	fftw_cleanup();
	fftw_cleanup_threads();

	free(template_ext);
	free(image_ext);

	return 0;
}


int normxcorr_fftw(float *templates, int template_len, int n_templates,
                   float *image, int image_len, float *ncc, int fft_len){
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
	templates:      Template signals
	template_len:   Length of template
	n_templates:    Number of templates (n0)
	image:          Image signal (to scan through)
	image_len:      Length of image
	ncc:            Output for cross-correlation - should be pointer to memory -
					must be n_templates x image_len - template_len + 1
	fft_len:        Size for fft (n1)
  */
	int N2 = fft_len / 2 + 1;
	int i, t, startind;
	double mean, stdev, old_mean, new_samp, old_samp, c, var=0.0, sum=0.0, acceptedDiff = 0.0000001;
	double * norm_sums = (double *) calloc(n_templates, sizeof(double));
	double * template_ext = (double *) calloc(fft_len * n_templates, sizeof(double));
	double * image_ext = (double *) calloc(fft_len, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * fft_len * n_templates);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	// Plan

	fftw_plan pa = fftw_plan_dft_r2c_2d(n_templates, fft_len, template_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_2d(n_templates, fft_len, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (t = 0; t < n_templates; ++t){
		for (i = 0; i < template_len; ++i)
		{
			template_ext[(t * fft_len) + i] = (double) templates[((t + 1) * template_len) - (i + 1)];
			norm_sums[t] += templates[(t * template_len) + i];
		}
	}
	for (i = 0; i < image_len; ++i)
	{
		image_ext[i] = (double) image[i];
	}
	//  Compute ffts of template and image
	fftw_execute(pa);
	fftw_execute(pb);

	//  Compute dot product
	for (t = 0; t < n_templates; ++t){
    	for (i = 0; i < N2; ++i)
	    {
		    out[(t * N2) + i][0] = outa[(t * N2) + i][0] * outb[i][0] - outa[(t * N2) + i][1] * outb[i][1];
    		out[(t * N2) + i][1] = outa[(t * N2) + i][0] * outb[i][1] + outa[(t * N2) + i][1] * outb[i][0];
    	}
    }
	//  Compute inverse fft
	fftw_execute(px);
	//  Procedures for normalisation
	// Compute starting mean, will update this
	for (i=0; i < template_len; ++i){
		sum += image[i];
	}
	mean = sum / template_len;

	// Compute starting standard deviation
	for (i=0; i < template_len; ++i){
		var += pow(image[i] - mean, 2) / (template_len);
	}
	stdev = sqrt(var);
    // Used for centering - taking only the valid part of the cross-correlation
	startind = template_len - 1;
	for (t = 0; t < n_templates; ++t){
    	if (var < acceptedDiff){
	    	ncc[t * (image_len - template_len + 1)] = 0;
    	}
	    else {
		    c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean) / stdev;
    		ncc[t * (image_len - template_len + 1)] = (float) c;
	    }
	}
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (image_len - template_len + 1); ++i){
		// Need to cast to double otherwise we end up with annoying floating
		// point errors when the variance is massive - collecting fp errors.
		new_samp = image[i + template_len - 1];
		old_samp = image[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / template_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (template_len);
		stdev = sqrt(var);
		for (t=0; t < n_templates; ++t){
			if (var > acceptedDiff){
				c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean ) / stdev;
				ncc[(t * (image_len - template_len + 1)) + i] = (float) c;
			}
			else{
				ncc[(t * (image_len - template_len + 1)) + i] = 0.0;
			}
		}
	}
	//  Clean up
	fftw_destroy_plan(pa);
	fftw_destroy_plan(pb);
	fftw_destroy_plan(px);

	fftw_free(out);
	fftw_free(outa);
	fftw_free(outb);
	fftw_free(ccc);

	fftw_cleanup();
	fftw_cleanup_threads();

	free(template_ext); free(image_ext); free(norm_sums);
	return 0;
}


int normxcorr_time(float *template, int template_len, float *image, int image_len, float *ccc){
    // Time domain cross-correlation - requires zero-mean template and image
	int p, k;
	int steps = image_len - template_len + 1;
	double numerator = 0.0, denom, mean = 0.0;
	double auto_a = 0.0, auto_b = 0.0;

    for (k=0; k < template_len; ++k){
		mean += image[k];
	}
	mean = mean / template_len;

	for(p = 0; p < template_len; ++p){
		auto_a += (double) template[p] * (double) template[p];
		numerator += (double) template[p] * ((double) image[p] - mean);
		auto_b += ((double) image[p] - mean) * ((double) image[p] - mean);
	}
	denom = sqrt(auto_a * auto_b);
	ccc[0] = (float) (numerator / denom);
	for(k = 1; k < steps; ++k){
		mean = mean + (image[k + template_len - 1] - image[k - 1]) / template_len;
	    numerator = 0.0;
	    auto_b = 0.0;
		for(p = 0; p < template_len; ++p){
			numerator += (double) template[p] * ((double) image[p + k] - mean);
			auto_b += ((double) image[p + k] - mean) * ((double) image[p + k] - mean);
		}
		denom = sqrt(auto_a * auto_b);
		ccc[k] = (float) (numerator / denom);
	}
	return 0;
}

int multi_normxcorr_fftw(float *templates, int n_templates, int template_len, int n_channels, float *image, int image_len, float *ccc, int fft_len){
    int i, r=0;
    for (i = 0; i < n_channels; ++i){
        r = normxcorr_fftw_threaded(&templates[n_templates * template_len * i], template_len,
                                    n_templates, &image[image_len * i], image_len,
                                    &ccc[(image_len - template_len + 1) * n_templates * i], fft_len);
    }
    return r;
}

int multi_normxcorr_time(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc){
	int i;
	for (i = 0; i < n_templates; ++i){
		normxcorr_time(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len + 1) * i]);
	}
	return 0;
}