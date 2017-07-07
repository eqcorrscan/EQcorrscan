/*
 * =====================================================================================
 *
 *       Filename:  norm.c
 *
 *        Purpose:  Normalise cross-correlations in a memory efficient way.
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

int normalise(float *ccc, int ccc_len, float *image, float norm_sum, int template_len)
{
	int i, j;
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
		ccc[0] = 0.0;
	}
	else
	{
		ccc[0] = (ccc[0] - norm_sum * mean ) / std;
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
			ccc[i] = 0.0;
		}
		else
		{
			// Normalize by the std and mean
			ccc[i] = (ccc[i] - norm_sum * mean ) / std;
		}
	}
	return 0;
}


int multi_normalise_loop(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n)
{
	int i;

	for (i=0; i<n; ++i)
	{
		normalise(&ccc[i * ccc_len], ccc_len, image, norm_sum[i], template_len);
	}
	return 0;
}

//int main(void)
////Used for testing
//{
//	int template_len = 600;
//	int ccc_len = 8640000;
//	int n = 3;
//	float * norm_sum = (float *) calloc(n, sizeof(float));
//	float * ccc = (float *) calloc(ccc_len * n, sizeof(float));
//	float * image = (float *) calloc(ccc_len + template_len, sizeof(float));
//	int i;
//
//	for (i=0; i<ccc_len * n; ++i)
//	{
//		ccc[i] = rand() % 100;
//		ccc[i] -= 50;
//		/* printf("%f\n", ccc[i]); */
//	}
//	for (i=0; i<ccc_len + template_len; ++i)
//	{
//		image[i] = rand() % 100;
//		image[i] -= 50;
//	}
//	for (i=0; i<n; ++i)
//	{
//	    norm_sum[i] = rand() % 1;
//	}
//	multi_normalise(ccc, ccc_len, image, norm_sum, template_len, n);
//	return 0;
//}