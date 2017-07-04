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
