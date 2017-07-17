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

// Prototypes
int normxcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N);

int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N);

int run_std_mean(int a_len, float *signalb, int b_len, float *run_std, float *run_mean);

int xcorr (float *signala, int a_len, float *signalb, int b_len, float *ccc);

int multi_corr (float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc);

int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n);

int run_std_mean(int a_len, float *signalb, int b_len, float *run_std, float *run_mean){
    int i;
	double sum = 0.0, mean, stdev, old_mean, var=0.0, new_samp, old_samp;

	for (i=0; i < a_len; ++i){
		sum += (double) signalb[i];
	}
	mean = sum / a_len;

	// Compute starting standard deviation
	for (i=0; i < a_len; ++i){
		var += pow(signalb[i] - mean, 2) / (a_len);
	}
	stdev = sqrt(var);

	run_std[0] = (float) stdev;
	run_mean[0] = (float) mean;
	for(i = 1; i < b_len; ++i){
	    new_samp = signalb[i + a_len - 1];
	    old_samp = signalb[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / a_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (a_len);
		stdev = sqrt(var);
        run_mean[i] = (float) mean;
        run_std[i] = (float) stdev;
	}
	return 0;
}

// Functions
int normxcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
				  float *ncc, int N){
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
	signala:  Template signal
	a_len:    Length of signala
	signalb:  Image signal (to scan through)
	b_len:    Length of signalb
	ncc:      Output for cross-correlation - should be pointer to memory -
			  must be b_len - a_len + 1 long
	N:        Size for fft
  */
	int N2 = N / 2 + 1;
	int i, startind;
	double norm_sum = 0.0, sum = 0.0;
	double mean, stdev, old_mean, new_samp, old_samp, c, var=0.0;
	float acceptedDiff = 0.0000001;
	double * signala_ext = (double *) calloc(N, sizeof(double));
	double * signalb_ext = (double *) calloc(N, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * N);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);

	fftw_plan pa = fftw_plan_dft_r2c_1d(N, signala_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(N, signalb_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_1d(N, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (i = 0; i < a_len; ++i)
	{
		signala_ext[i] = (double) signala[a_len - (i + 1)];
		norm_sum += signala[i];
	}
	for (i = 0; i < b_len; ++i)
	{
		signalb_ext[i] = (double) signalb[i];
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

	startind = a_len - 1;

	//  Procedures for normalisation
	// Compute starting mean, will update this
	for (i=0; i < a_len; ++i){
		sum += signalb[i];
	}
	mean = sum / a_len;

	// Compute starting standard deviation
	for (i=0; i < a_len; ++i){
		var += pow(signalb[i] - mean, 2) / (a_len);
	}
	stdev = sqrt(var);

	if (var < acceptedDiff){
        ncc[0] = 0;
	}
	else {
	    c = ((ccc[startind] / N) - norm_sum * mean) / stdev;
	    ncc[0] = (float) c;
	}
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (b_len - a_len + 1); ++i){
	    // Need to cast to double otherwise we end up with annoying floating
	    // point errors when the variance is massive.
		new_samp = signalb[i + a_len - 1];
	    old_samp = signalb[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / a_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (a_len);
		stdev = sqrt(var);
		if (var > acceptedDiff){
		    c = ((ccc[i + startind] / N) - norm_sum * mean ) / stdev;
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

	free(signala_ext);
	free(signalb_ext);

	return 0;
}


// Functions
int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
				  float *ncc, int N){
  /*
  Purpose: compute frequency domain cross-correlation of real data using fftw
  Author: Calum J. Chamberlain

  Note: This is NOT normalised

  Date: 12/06/2017
  Args:
	signala:  Template signal
	a_len:    Length of signala
	signalb:  Image signal (to scan through)
	b_len:    Length of signalb
	ncc:      Output for cross-correlation - should be pointer to memory -
			  must be b_len - a_len + 1 long
	N:        Size for fft
  */
	int N2 = N / 2 + 1;
	int i, startind;
	double * signala_ext = (double *) calloc(N, sizeof(double));
	double * signalb_ext = (double *) calloc(N, sizeof(double));
	double * ccc = (double *) fftw_malloc(sizeof(double) * N);
	fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N2);

	fftw_plan pa = fftw_plan_dft_r2c_1d(N, signala_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(N, signalb_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_1d(N, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (i = 0; i < a_len; ++i)
	{
		signala_ext[i] = signala[a_len - (i + 1)];
	}
	for (i = 0; i < b_len; ++i)
	{
		signalb_ext[i] = signalb[i];
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

	startind = a_len - 1;

	//  Procedures for normalisation
	ncc[0] = ccc[startind] / N;
	// Center and divide by length to generate scaled convolution
	for(i = 1; i < (b_len - a_len + 1); ++i){
		ncc[i] = ccc[i + startind] / N;
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

	free(signala_ext);
	free(signalb_ext);

	return 0;
}


int xcorr(float *signala, int a_len, float *signalb, int b_len, float *ccc){
	int p, k;
	int steps = b_len - a_len + 1;
	float numerator, denom;
	float auto_a = 0.0, auto_b = 0.0;

	for(p = 0; p < a_len; ++p){
		auto_a += signala[p] * signala[p];
	}
	for(k = 0; k < steps; ++k){
		numerator = 0.0;
		auto_b = 0.0;
		for(p = 0; p < a_len; ++p){
			numerator += signala[p] * signalb[p + k];
		}
		for(p = 0; p < a_len; ++p){
			auto_b += signalb[p + k] * signalb[p + k];
		}
		denom = sqrtf(auto_a * auto_b);
		ccc[k] = numerator / denom;
	}
	return 0;
}


int multi_corr(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc){
	int i;
	#pragma omp parallel for
	for (i = 0; i < n_templates; ++i){
		xcorr(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len) * i]);
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

