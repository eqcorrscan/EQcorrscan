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
#include <string.h>
#include <math.h>
#if (defined(_MSC_VER))
    #if (_MSC_VER < 1800)
        #include <float.h>
        #define isnanf(x) _isnanf(x)
    #endif
    #define inline __inline
#endif
#if (defined(__APPLE__) && !isnanf)
    #define isnanf isnan
#endif
#include <fftw3.h>
#if defined(__linux__) || defined(__linux) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    #include <omp.h>
    #ifndef N_THREADS
        #define N_THREADS omp_get_max_threads()
    #endif
#endif

// Prototypes
int normxcorr_fftw(float*, int, int, float*, int, float*, int, int*, int*);

static inline int set_ncc(int t, int i, int template_len, int image_len, float value, int *used_chans, int *pad_array, float *ncc);

int normxcorr_fftw_main(float*, int, int, float*, int, float*, int, double*, double*, double*,
        fftw_complex*, fftw_complex*, fftw_complex*, fftw_plan, fftw_plan, fftw_plan, int*, int*);

int normxcorr_fftw_threaded(float*, int, int, float*, int, float*, int, int*, int*);

void free_fftw_arrays(int, double**, double**, double**, fftw_complex**, fftw_complex**, fftw_complex**);

int multi_normxcorr_fftw(float*, int, int, int, float*, int, float*, int, int*, int*);

// Functions
int normxcorr_fftw_threaded(float *templates, int template_len, int n_templates,
					        float *image, int image_len, float *ncc, int fft_len,
                            int *used_chans, int *pad_array) {
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
	int i, t, startind, status = 0;
	double mean, stdev, old_mean, new_samp, old_samp, var=0.0, sum=0.0, acceptedDiff = 0.0000001;
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
    if (var >= acceptedDiff) {
        for (t = 0; t < n_templates; ++t){
            double c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean) / stdev;
            status += set_ncc(t, 0, template_len, image_len, (float) c, used_chans, pad_array, ncc);
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
        if (var > acceptedDiff) { // TODO: above is '>=', should they be the same?
            for (t = 0; t < n_templates; ++t){
                double c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean ) / stdev;
                status += set_ncc(t, i, template_len, image_len, (float) c, used_chans, pad_array, ncc);
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

	return status;
}


int normxcorr_fftw(float *templates, int template_len, int n_templates,
                   float *image, int image_len, float *ncc, int fft_len,
                   int *used_chans, int *pad_array){
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
  Notes:
    This is a wrapper around `normxcorr_fftw_main`, allocating required memory and plans
    for that function. We have taken this outside the main function because creating plans
    is not thread-safe and we want to call the main function from within an OpenMP loop.
  */
	int status = 0;
	int N2 = fft_len / 2 + 1;
	// All memory allocated with `fftw_malloc` to ensure 16-byte aligned
	double * template_ext = fftw_alloc_real(fft_len * n_templates);
	double * image_ext = fftw_alloc_real(fft_len);
	double * ccc = fftw_alloc_real(fft_len * n_templates);
	fftw_complex * outa = fftw_alloc_complex(N2 * n_templates);
	fftw_complex * outb = fftw_alloc_complex(N2);
	fftw_complex * out = fftw_alloc_complex(N2 * n_templates);
	// Plan
	fftw_plan pa = fftw_plan_dft_r2c_2d(n_templates, fft_len, template_ext, outa, FFTW_ESTIMATE);
	fftw_plan pb = fftw_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
	fftw_plan px = fftw_plan_dft_c2r_2d(n_templates, fft_len, out, ccc, FFTW_ESTIMATE);

	// Initialise to zero
	memset(template_ext, 0, fft_len * n_templates * sizeof(double));
	memset(image_ext, 0, fft_len * sizeof(double));

	// Call the function to do the work
	status = normxcorr_fftw_main(templates, template_len, n_templates, image, image_len,
			ncc, fft_len, template_ext, image_ext, ccc, outa, outb, out, pa, pb, px,
            used_chans, pad_array);

	// free memory and plans
	fftw_destroy_plan(pa);
	fftw_destroy_plan(pb);
	fftw_destroy_plan(px);

	fftw_free(out);
	fftw_free(outa);
	fftw_free(outb);
	fftw_free(ccc);
	fftw_free(template_ext);
	fftw_free(image_ext);

	fftw_cleanup();
	fftw_cleanup_threads();

	return status;
}


int normxcorr_fftw_main(float *templates, int template_len, int n_templates,
                        float *image, int image_len, float *ncc, int fft_len,
                        double *template_ext, double *image_ext, double *ccc,
                        fftw_complex *outa, fftw_complex *outb, fftw_complex *out,
                        fftw_plan pa, fftw_plan pb, fftw_plan px, int *used_chans,
                        int *pad_array) {
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
                    It is assumed that ncc will be initialised to zero before
                    passing into this function
    fft_len:        Size for fft (n1)
    template_ext:   Input FFTW array for template transform (must be allocated)
    image_ext:      Input FFTW array for image transform (must be allocated)
    ccc:            Output FFTW array for reverse transform (must be allocated)
    outa:           Output FFTW array for template transform (must be allocatd)
    outb:           Output FFTW array for image transform (must be allocated)
    out:            Input array for reverse transform (must be allocated)
    pa:             Forward plan for templates
    pb:             Forward plan for image
    px:             Reverse plan
  */
	int N2 = fft_len / 2 + 1;
	int i, t, startind, status = 0;
	double mean, stdev, old_mean, new_samp, old_samp, var=0.0, sum=0.0, acceptedDiff = 0.0000001;
	double * norm_sums = (double *) calloc(n_templates, sizeof(double));

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
	fftw_execute_dft_r2c(pa, template_ext, outa);
	fftw_execute_dft_r2c(pb, image_ext, outb);

	//  Compute dot product
	for (t = 0; t < n_templates; ++t){
    	for (i = 0; i < N2; ++i)
	    {
		    out[(t * N2) + i][0] = outa[(t * N2) + i][0] * outb[i][0] - outa[(t * N2) + i][1] * outb[i][1];
    		out[(t * N2) + i][1] = outa[(t * N2) + i][0] * outb[i][1] + outa[(t * N2) + i][1] * outb[i][0];
    	}
    }
	//  Compute inverse fft
	fftw_execute_dft_c2r(px, out, ccc);

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
    if (var >= acceptedDiff) {
        for (t = 0; t < n_templates; ++t){
            double c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean) / stdev;
            status += set_ncc(t, 0, template_len, image_len, (float) c, used_chans, pad_array, ncc);
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
        if (var > acceptedDiff) { // TODO: above is '>=', should they be the same?
            for (t = 0; t < n_templates; ++t){
                double c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean ) / stdev;
                status += set_ncc(t, i, template_len, image_len, (float) c, used_chans, pad_array, ncc);
			}
		}
	}
	//  Clean up
	free(norm_sums);

	return status;
}


static inline int set_ncc(int t, int i, int template_len, int image_len, float value, int *used_chans, int *pad_array, float *ncc) {
    int status = 0;

    if (used_chans[t] && (i >= pad_array[t])) {
        size_t ncc_index = t * (image_len - template_len + 1) + i - pad_array[t];

        if (isnanf(value)) {
            // set NaNs to zero
            value = 0.0;
        }
        else if (fabsf(value) > 1.01) {
            // this will raise an exception when we return to Python
            status = 1;
        }
        else if (value > 1.0) {
            value = 1.0;
        }
        else if (value < -1.0) {
            value = -1.0;
        }

        #pragma omp atomic
        ncc[ncc_index] += value;
    }

    return status;
}


void free_fftw_arrays(int size, double **template_ext, double **image_ext, double **ccc,
        fftw_complex **outa, fftw_complex **outb, fftw_complex **out) {
    int i;

    /* free memory */
    for (i = 0; i < size; i++) {
        fftw_free(template_ext[i]);
        fftw_free(image_ext[i]);
        fftw_free(ccc[i]);
        fftw_free(outa[i]);
        fftw_free(outb[i]);
        fftw_free(out[i]);
    }
    free(template_ext);
    free(image_ext);
    free(ccc);
    free(outa);
    free(outb);
    free(out);
}


int multi_normxcorr_fftw(float *templates, int n_templates, int template_len, int n_channels,
        float *image, int image_len, float *ncc, int fft_len, int *used_chans, int *pad_array){
    int i;
    int r = 0;
    size_t N2 = (size_t) fft_len / 2 + 1;
    double **template_ext = NULL;
    double **image_ext = NULL;
    double **ccc = NULL;
    fftw_complex **outa = NULL;
    fftw_complex **outb = NULL;
    fftw_complex **out = NULL;
    fftw_plan pa, pb, px;
    int num_threads = 1;


    #ifdef N_THREADS
    /* set the number of threads - the minimum of the numbers of channels and threads */
    num_threads = (N_THREADS > n_channels) ? n_channels : N_THREADS;
    #endif

    /* allocate memory for all threads here */
    template_ext = (double**) malloc(num_threads * sizeof(double*));
    if (template_ext == NULL) {
        printf("Error allocating template_ext\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    image_ext = (double**) malloc(num_threads * sizeof(double*));
    if (image_ext == NULL) {
        printf("Error allocating image_ext\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    ccc = (double**) malloc(num_threads * sizeof(double*));
    if (ccc == NULL) {
        printf("Error allocating ccc\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    outa = (fftw_complex**) malloc(num_threads * sizeof(fftw_complex*));
    if (outa == NULL) {
        printf("Error allocating outa\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    outb = (fftw_complex**) malloc(num_threads * sizeof(fftw_complex*));
    if (outb == NULL) {
        printf("Error allocating outb\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    out = (fftw_complex**) malloc(num_threads * sizeof(fftw_complex*));
    if (out == NULL) {
        printf("Error allocating out\n");
        free_fftw_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }

    // All memory allocated with `fftw_malloc` to ensure 16-byte aligned.
    for (i = 0; i < num_threads; i++) {
        /* initialise all to NULL so that freeing on error works */
        template_ext[i] = NULL;
        image_ext[i] = NULL;
        ccc[i] = NULL;
        outa[i] = NULL;
        outb[i] = NULL;
        out[i] = NULL;

        /* allocate template_ext arrays */
        template_ext[i] = fftw_alloc_real((size_t) fft_len * n_templates);
        if (template_ext[i] == NULL) {
            printf("Error allocating template_ext[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate image_ext arrays */
        image_ext[i] = fftw_alloc_real(fft_len);
        if (image_ext[i] == NULL) {
            printf("Error allocating image_ext[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate ccc arrays */
        ccc[i] = fftw_alloc_real((size_t) fft_len * n_templates);
        if (ccc[i] == NULL) {
            printf("Error allocating ccc[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate outa arrays */
        outa[i] = fftw_alloc_complex((size_t) N2 * n_templates);
        if (outa[i] == NULL) {
            printf("Error allocating outa[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate outb arrays */
        outb[i] = fftw_alloc_complex((size_t) N2);
        if (outb[i] == NULL) {
            printf("Error allocating outb[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate out arrays */
        out[i] = fftw_alloc_complex((size_t) N2 * n_templates);
        if (out[i] == NULL) {
            printf("Error allocating out[%d]\n", i);
            free_fftw_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }
    }

    //TODO: touch all arrays - NUMA first touch???

    // We create the plans here since they are not thread safe.
    pa = fftw_plan_dft_r2c_2d(n_templates, fft_len, template_ext[0], outa[0], FFTW_ESTIMATE);
    pb = fftw_plan_dft_r2c_1d(fft_len, image_ext[0], outb[0], FFTW_ESTIMATE);
    px = fftw_plan_dft_c2r_2d(n_templates, fft_len, out[0], ccc[0], FFTW_ESTIMATE);

    /* loop over the channels */
    #pragma omp parallel for reduction(+:r) num_threads(num_threads)
    for (i = 0; i < n_channels; ++i){
        int tid = 0; /* each thread has its own workspace */

        #ifdef N_THREADS
        /* get the id of this thread */
        tid = omp_get_thread_num();
        #endif

        /* initialise memory to zero */
        memset(template_ext[tid], 0, (size_t) fft_len * n_templates * sizeof(double));
        memset(image_ext[tid], 0, (size_t) fft_len * sizeof(double));

        /* call the routine */
        r += normxcorr_fftw_main(&templates[(size_t) n_templates * template_len * i], template_len,
                                 n_templates, &image[(size_t) image_len * i], image_len, ncc, fft_len,
                                 template_ext[tid], image_ext[tid], ccc[tid], outa[tid], outb[tid], out[tid],
                                 pa, pb, px, &used_chans[(size_t) i * n_templates],
                                 &pad_array[(size_t) i * n_templates]);
    }

    /* free fftw memory */
    free_fftw_arrays(num_threads, template_ext, image_ext, ccc, outa, outb, out);
    fftw_destroy_plan(pa);
    fftw_destroy_plan(pb);
    fftw_destroy_plan(px);
    fftw_cleanup();

    return r;
}
