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
#include <omp.h>

// Prototypes
int normxcorr_fftw_1d(float *template, int template_len, float *image, int image_len, float *ncc, int fft_len);

int xcorr_fftw_1d(float *template, int template_len, float *image, int image_len, float *ncc, int fft_len);

int normxcorr_fftw_2d(float *templates, int template_len, int n_templates, float *image, int image_len, float *ncc, int fft_len);

int multi_fftw_normxcorr(float *templates, int n_templates, int template_len, int n_channels, float *image, int image_len, float *ccc, int fft_len);

int run_std_mean(int template_len, float *image, int image_len, float *run_std, float *run_mean);

int xcorr (float *template, int template_len, float *image, int image_len, float *ccc);

int multi_corr (float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc);

int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n);

// Functions
int run_std_mean(int template_len, float *image, int image_len, float *run_std, float *run_mean){
	int i;
	double sum = 0.0, mean, stdev, old_mean, var=0.0, new_samp, old_samp;

	for (i=0; i < template_len; ++i){
		sum += (double) image[i];
	}
	mean = sum / template_len;

	// Compute starting standard deviation
	for (i=0; i < template_len; ++i){
		var += pow(image[i] - mean, 2) / (template_len);
	}
	stdev = sqrt(var);

	run_std[0] = (float) stdev;
	run_mean[0] = (float) mean;
	for(i = 1; i < image_len; ++i){
		new_samp = image[i + template_len - 1];
		old_samp = image[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / template_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (template_len);
		stdev = sqrt(var);
		run_mean[i] = (float) mean;
		run_std[i] = (float) stdev;
	}
	return 0;
}

int normxcorr_fftw_2d(float *templates, int template_len, int n_templates,
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
	int i, t, startind, n_threads;
	double mean, stdev, old_mean, new_samp, old_samp, c, var=0.0, sum=0.0, acceptedDiff = 0.0000001;
	double * norm_sums = (double *) calloc(n_templates, sizeof(double));
	double * template_ext = (double *) calloc(fft_len * n_templates, sizeof(double));
	double * image_ext = (double *) calloc(fft_len, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * fft_len * n_templates);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2 * n_templates);
	// Initialize threads
	fftw_init_threads();
	n_threads = omp_get_max_threads();
	fftw_plan_with_nthreads(n_threads);
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


int normxcorr_fftw_1d(float *template, int template_len, float *image, int image_len,
				  float *ncc, int fft_len){
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
	template:       Template signal
	template_len:   Length of template
	image:          Image signal (to scan through)
	image_len:      Length of image
	ncc:            Output for cross-correlation - should be pointer to memory -
					must be image_len - template_len + 1 long
	fft_len:        Size for fft
  */
	int N2 = fft_len / 2 + 1;
	int i, startind;
	double norm_sum = 0.0, sum = 0.0;
	double mean, stdev, old_mean, new_samp, old_samp, c, var=0.0;
	double acceptedDiff = 0.0000001;
	double * template_ext = (double *) calloc(fft_len, sizeof(double));
	double * image_ext = (double *) calloc(fft_len, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * fft_len);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);

	fftw_plan pa = fftw_plan_dft_r2c_1d(fft_len, template_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_1d(fft_len, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (i = 0; i < template_len; ++i)
	{
		template_ext[i] = (double) template[template_len - (i + 1)];
		norm_sum += template[i];
	}
	for (i = 0; i < image_len; ++i)
	{
		image_ext[i] = (double) image[i];
	}
	//  Compute ffts of template and image
	fftw_execute(pa);
	fftw_execute(pb);

	//  Compute dot product
	for (i = 0; i < N2; ++i)
	{
		out[i][0] = outa[i][0] * outb[i][0] - outa[i][1] * outb[i][1];
		out[i][1] = outa[i][0] * outb[i][1] + outa[i][1] * outb[i][0];
	}
	//  Compute inverse fft
	fftw_execute(px);

	startind = template_len - 1;

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

	if (var < acceptedDiff){
		ncc[0] = 0;
	}
	else {
		c = ((ccc[startind] / fft_len) - norm_sum * mean) / stdev;
		ncc[0] = (float) c;
	}
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (image_len - template_len + 1); ++i){
		// Need to cast to double otherwise we end up with annoying floating
		// point errors when the variance is massive.
		new_samp = image[i + template_len - 1];
		old_samp = image[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / template_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (template_len);
		stdev = sqrt(var);
		if (var > acceptedDiff){
			c = ((ccc[i + startind] / fft_len) - norm_sum * mean ) / stdev;
			ncc[i] = (float) c;
		}
		else{
			ncc[i] = 0.0;
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

	free(template_ext);
	free(image_ext);

	return 0;
}


int xcorr_fftw_1d(float *template, int template_len, float *image, int image_len,
				  float *ncc, int fft_len){
  /*
  Purpose: compute frequency domain cross-correlation of real data using fftw
  Author: Calum J. Chamberlain

  Note: This is NOT normalised

  Date: 12/06/2017
  Args:
	template:       Template signal
	template_len:   Length of template
	image:          Image signal (to scan through)
	image_len:      Length of image
	ncc:            Output for cross-correlation - should be pointer to memory -
					must be image_len - template_len + 1 long
	fft_len:        Size for fft
  */
	int N2 = fft_len / 2 + 1;
	int i, startind;
	double * template_ext = (double *) calloc(fft_len, sizeof(double));
	double * image_ext = (double *) calloc(fft_len, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * fft_len);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);

	fftw_plan pa = fftw_plan_dft_r2c_1d(fft_len, template_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_1d(fft_len, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (i = 0; i < template_len; ++i)
	{
		template_ext[i] = template[template_len - (i + 1)];
	}
	for (i = 0; i < image_len; ++i)
	{
		image_ext[i] = image[i];
	}
	//  Compute ffts of template and image
	fftw_execute(pa);
	fftw_execute(pb);

	//  Compute dot product
	for (i = 0; i < N2; ++i)
	{
		out[i][0] = outa[i][0] * outb[i][0] - outa[i][1] * outb[i][1];
		out[i][1] = outa[i][0] * outb[i][1] + outa[i][1] * outb[i][0];
	}
	//  Compute inverse fft
	fftw_execute(px);

	startind = template_len - 1;

	//  Procedures for normalisation
	ncc[0] = ccc[startind] / fft_len;
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (image_len - template_len + 1); ++i){
		ncc[i] = ccc[i + startind] / fft_len;
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

	free(template_ext);
	free(image_ext);

	return 0;
}


int xcorr(float *template, int template_len, float *image, int image_len, float *ccc){
	int p, k;
	int steps = image_len - template_len + 1;
	float numerator, denom;
	float auto_a = 0.0, auto_b = 0.0;

	for(p = 0; p < template_len; ++p){
		auto_a += template[p] * template[p];
	}
	for(k = 0; k < steps; ++k){
		numerator = 0.0;
		auto_b = 0.0;
		for(p = 0; p < template_len; ++p){
			numerator += template[p] * image[p + k];
		}
		for(p = 0; p < template_len; ++p){
			auto_b += image[p + k] * image[p + k];
		}
		denom = sqrtf(auto_a * auto_b);
		ccc[k] = numerator / denom;
	}
	return 0;
}

int multi_fftw_normxcorr(float *templates, int n_templates, int template_len, int n_channels, float *image, int image_len, float *ccc, int fft_len){
    int i, r=0;
    for (i = 0; i < n_channels; ++i){
        r = normxcorr_fftw_2d(&templates[n_templates * template_len * i], template_len,
                              n_templates, &image[image_len * i], image_len,
                              &ccc[(image_len - template_len + 1) * n_templates * i], fft_len);
    }
    return r;
}

int multi_corr(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc){
	int i;
	for (i = 0; i < n_templates; ++i){
		xcorr(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len + 1) * i]);
	}
	return 0;
}


int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n)
{
	int i, k;
	float mean, old_mean, std, sum=0.0, var=0.0;
	float acceptedDiff = 0.00000001;

	// Compute starting mean, will update this
	for (i=0; i < template_len; ++i)
	{
		sum += image[i];
	}
	mean = sum / template_len;

	// Compute starting standard deviation, will update this
	for (i=0; i < template_len; ++i)
	{
		var += powf(image[i] - mean, 2);
	}
	std = sqrtf(var / (template_len));
	if (std < acceptedDiff)
	{
		for (k=0; k<n; ++k)
		{
			ccc[k] = 0.0;
		}
	}
	else
	{
		for (k=0; k<n; ++k)
		{
			ccc[k] = (ccc[k] - norm_sum[k] * mean ) / std;
		}
	}
	if (ccc[0] != ccc[0])
	{
		ccc[0] = 0.0;
	}
	// Loop through updating as we go
	for(i=1; i<ccc_len; ++i)
	{
		old_mean = mean;
		mean = mean + (image[i + template_len - 1] - image[i - 1]) / template_len;
		// Don't know of a usefully accurate and efficient method :(
		var += (image[i + template_len - 1] - image[i - 1]) * (image[i + template_len -1] - mean + image[i - 1] - old_mean) / template_len;
		std = sqrtf(var / template_len);
		if (std < acceptedDiff)
		{
			for (k=0; k<n; ++k)
			{
				ccc[(i * n) + k] = 0.0;
			}
		}
		else
		{
			// Normalize by the std and mean
			for (k=0; k<n; ++k)
			{
				ccc[(i * n) + k] = (ccc[(i * n) + k] - norm_sum[k] * mean ) / std;
			}
		}
		// Convert nans
		if (ccc[i] != ccc[i])
		{
			ccc[i] = 0.0;
		}
	}
	return 0;
}

