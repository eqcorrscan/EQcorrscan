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
int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N);

int xcorr (float *signala, int a_len, float *signalb, int b_len, float *ccc);

int multi_corr (float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc);

int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n);

// Functions

int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
                  float *ncc, int N){
  /*
  Purpose: compute frequency domain cross-correlation of real data using fftw
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
	int i, j, startind;
	float norm_sum = 0.0, sum = 0.0, mean, var = 0.0, stdev;
	float acceptedDiff = 0.0000001;
	float * signala_ext = (float *) calloc(N, sizeof(float));
	float * signalb_ext = (float *) calloc(N, sizeof(float));
	float * ccc = (float *) fftwf_malloc(sizeof(float) * N);
	fftwf_complex * outa = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2);
	fftwf_complex * outb = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2);
	fftwf_complex * out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2);

	fftwf_plan pa = fftwf_plan_dft_r2c_1d(N, signala_ext, outa, FFTW_ESTIMATE);
	fftwf_plan pb = fftwf_plan_dft_r2c_1d(N, signalb_ext, outb, FFTW_ESTIMATE);
	fftwf_plan px = fftwf_plan_dft_c2r_1d(N, out, ccc, FFTW_ESTIMATE);

	// zero padding - and flip template
	for (i = 0; i < a_len; ++i)
	{
	    signala_ext[i] = signala[a_len - (i + 1)];
	    norm_sum += signala[i];
	}
	for (i = 0; i < b_len; ++i)
	{
	    signalb_ext[i] = signalb[i];
	}
    //  Compute ffts of template and image
	fftwf_execute(pa);
	fftwf_execute(pb);

    //  Compute dot product
    for (i = 0; i < N2; ++i)
    {
        out[i][0] = outa[i][0] * outb[i][0] - outa[i][1] * outb[i][1];
        out[i][1] = outa[i][0] * outb[i][1] + outa[i][1] * outb[i][0];
    }
    //  Compute inverse fft
	fftwf_execute(px);

    startind = a_len - 1;

    //  Procedures for normalisation
	// Compute starting mean, will update this
	for (i=0; i < a_len; ++i){
	    sum += signalb[i];
	}
	mean = sum / a_len;

	// Compute starting standard deviation
	for (i=0; i < a_len; ++i){
	    var += powf(signalb[i] - mean, 2);
	}
	stdev = sqrtf(var / (a_len));

	if (var < acceptedDiff){
        ncc[0]=0;
    }
    else {
        ncc[0]=((ccc[startind] / N) - norm_sum * mean) / stdev;
    }
    // Center and divide by length to generate scaled convolution
	for(i = 1; i < (b_len - a_len + 1); ++i){
		mean = mean + (signalb[i + a_len - 1] - signalb[i - 1]) / a_len;
		// Don't know of a usefully accurate and efficient method of running std
		var = 0.0;
		for (j=0; j<a_len; ++j)
		{
			var += powf(signalb[i + j] - mean, 2);
		}
		stdev = sqrtf(var / a_len);
		if (var > acceptedDiff){
		    ncc[i] = ((ccc[i + startind] / N) - norm_sum * mean ) / stdev;
		}
		else{
		    ncc[i] = 0.0;
		}
	}
    //  Clean up
	fftwf_destroy_plan(pa);
	fftwf_destroy_plan(pb);
	fftwf_destroy_plan(px);

	fftwf_free(out);
	fftwf_free(outa);
	fftwf_free(outb);
	fftwf_free(ccc);

	fftwf_cleanup();

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
	int i, j, k;
	float mean, std, sum=0.0, var=0.0;
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
		mean = mean + (image[i + template_len - 1] - image[i - 1]) / template_len;
		// Don't know of a usefully accurate and efficient method :(
		var = 0.0;
		for (j=0; j<template_len; ++j)
		{
			var += powf(image[i + j] - mean, 2);
		}
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

