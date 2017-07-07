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
#include <complex.h>
#include <math.h>
#include <fftw3.h>

// Prototypes
int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N, float *norm_sums);

int normxcorr_fftw_loop(float *templates, int a_len, float *signalb, int b_len, float *ncc, int N, int n_templates);

int xcorr (float *signala, int a_len, float *signalb, int b_len, float *ccc);

int multi_corr (float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc);

int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n);

int center(float *a, int inlen, float *out, int outlen);


// Functions

int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
                  float *ncc, int N, float *norm_sums){
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
    N:        Size for fft (and of result)
  */
	int N2 = N / 2 + 1;
	int i;
	float norm_sum = 0.0;
	float * signala_ext = (float *) calloc(N, sizeof(float));
	float * signalb_ext = (float *) calloc(N, sizeof(float));
	float * result = (float *) calloc(b_len + a_len - 1, sizeof(float));
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
	norm_sums[0] = norm_sum;
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
        out[i] = outa[i] * outb[i];
    }
    //  Compute inverse fft
	fftwf_execute(px);

	for(i = 0; i < (a_len + b_len - 1); ++i){
	    result[i] = ccc[i] / N;
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

    // Still do do: center then normalise
    center(&result[0], a_len + b_len -1, &ncc[0], b_len - a_len + 1);

    free(result);
	return 0;
}


int normxcorr_fftw_loop(float *templates, int a_len, float *signalb, int b_len,
                        float *ncc, int N, int n_templates){
    float * norm_sums = (float *) calloc(n_templates, sizeof(float));
    float * ccc_reshape = (float *) calloc(n_templates * (b_len - a_len + 1), sizeof(float));
    int i, j;
    int ccc_len = b_len - a_len + 1;

    for (i=0; i<n_templates; ++i){
        xcorr_fftw_1d(&templates[a_len * i], a_len, &signalb[0], b_len, &ncc[i * ccc_len], N, &norm_sums[i]);
    }
    // Need to reshape ncc before going to normalise.
    for (i=0; i < ccc_len; ++i){
        for (j=0; j < n_templates; ++j){
            ccc_reshape[(i * n_templates) + j] = ncc[i + (j * ccc_len)];
        }
    }
    multi_normalise(&ccc_reshape[0], b_len - a_len + 1, &signalb[0], &norm_sums[0], a_len, n_templates);
    free(norm_sums);
    for (i=0; i < ccc_len; ++i){
        for (j=0; j < n_templates; ++j){
            ncc[i + (j * ccc_len)] = ccc_reshape[(i * n_templates) + j];
        }
    }
    free(ccc_reshape);
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
	if (std == 0.0)
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
		if (std == 0.0)
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

int center(float *in, int inlen, float *out, int outlen){
    int startind, i;

    if (outlen > inlen){
        // Can't cope with this!
        fprintf(stderr, "Out is bigger than in, aborting\n");
        exit(EXIT_FAILURE);
    }
    startind = (inlen - outlen) / 2;
    for (i=0; i < outlen; ++i){
        out[i] = in[i + startind];
    }
    return 0;
}