/*
 * =====================================================================================
 *
 *       Filename:  norm.c
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
	float mean, std=0.0, sum=0.0, var=0.0;

	for (i=0; i < template_len; ++i)
	{
		sum += image[i];
	}
	mean = sum / template_len;

	for (i=0; i < template_len; ++i)
	{
		var += pow(image[i] - mean, 2);
	}
	std = sqrt(var / template_len);

	ccc[0] = (ccc[0] - norm_sum * mean ) / std;

	for(i=1; i<ccc_len; ++i)
	{
	    sum -= image[i-1];
	    sum += image[i+template_len-1];
		mean = sum / template_len;

//		I think there should be a better way to do this
        var = 0.0;
		for (j=0; j < template_len; ++j)
		{
		    var += pow(image[i + j] - mean, 2);
		}
		std = sqrt(var / template_len);
//		std = sqrt((var - pow(image[i-1] - mean, 2) + pow(image[i+template_len-1] - mean, 2)) / template_len);
		ccc[i] = (ccc[i] - norm_sum * mean ) / std;
	}
	return 0;
}

int main()
//Used for testing
{
	int template_len = 100;
	int ccc_len = 864000;
	float norm_sum = 5.3;
	float * ccc = (float *) calloc(ccc_len, sizeof(float));
	float * image = (float *) calloc(ccc_len + template_len, sizeof(float));
	int i;

	for (i=0; i<ccc_len; ++i)
	{
		ccc[i] = rand() % 100;
		ccc[i] -= 50;
		/* printf("%f\n", ccc[i]); */
	}
	for (i=0; i<ccc_len + template_len; ++i)
	{
		image[i] = rand() % 100;
		image[i] -= 50;
	}
	normalise(ccc, ccc_len, image, norm_sum, template_len);
	return 0;
}
