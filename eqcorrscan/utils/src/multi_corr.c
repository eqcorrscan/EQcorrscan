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
    #include <float.h>
    #define isnanf(x) _isnan(x)
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
#ifndef OUTER_SAFE
    #if defined(__linux__) || defined(__linux)
        #define OUTER_SAFE 1
    #else
        #define OUTER_SAFE 0
    #endif
#else
    #define OUTER_SAFE 1
#endif
// Define minimum variance to compute correlations - requires some signal
#define ACCEPTED_DIFF 1e-10 //1e-15
// Define difference to warn user on
#define WARN_DIFF 1e-8 //1e-10

// Prototypes
int normxcorr_fftw(float*, long, long, float*, long, float*, long, int*, int*, int*);

static inline int set_ncc(long t, long i, long template_len, long image_len, float value, int *used_chans, int *pad_array, float *ncc);

int normxcorr_fftw_main(float*, long, long, float*, long, float*, long, float*, float*, float*,
        fftwf_complex*, fftwf_complex*, fftwf_complex*, fftwf_plan, fftwf_plan, fftwf_plan, int*, int*, int, int*);

int normxcorr_fftw_threaded(float*, long, long, float*, long, float*, long, int*, int*, int*);

void free_fftwf_arrays(int, float**, float**, float**, fftwf_complex**, fftwf_complex**, fftwf_complex**);

void free_fftw_arrays(int, double**, double**, double**, fftw_complex**, fftw_complex**, fftw_complex**);

int multi_normxcorr_fftw(float*, long, long, long, float*, long, float*, long, int*, int*, int, int, int*);

// Functions
int normxcorr_fftw_threaded(float *templates, long template_len, long n_templates,
                            float *image, long image_len, float *ncc, long fft_len,
                            int *used_chans, int *pad_array, int *variance_warning) {
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
    long N2 = fft_len / 2 + 1;
    long i, t, startind;
    int status = 0;
    int flatline_count = 0;
    double mean, stdev, old_mean, new_samp, old_samp, var=0.0, sum=0.0;
    float * norm_sums = (float *) calloc(n_templates, sizeof(float));
    float * template_ext = (float *) calloc(fft_len * n_templates, sizeof(float));
    float * image_ext = (float *) calloc(fft_len, sizeof(float));
    float * ccc = (float *) fftwf_malloc(sizeof(float) * fft_len * n_templates);
    fftwf_complex * outa = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2 * n_templates);
    fftwf_complex * outb = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2);
    fftwf_complex * out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2 * n_templates);
    // Initialize threads
    #ifdef N_THREADS
        fftwf_init_threads();
        fftwf_plan_with_nthreads(N_THREADS);
    #endif
    // Plan

    fftwf_plan pa = fftwf_plan_dft_r2c_2d(n_templates, fft_len, template_ext, outa, FFTW_ESTIMATE);
    fftwf_plan pb = fftwf_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
    fftwf_plan px = fftwf_plan_dft_c2r_2d(n_templates, fft_len, out, ccc, FFTW_ESTIMATE);

    // zero padding - and flip template
    for (t = 0; t < n_templates; ++t){
        for (i = 0; i < template_len; ++i)
        {
            template_ext[(t * fft_len) + i] = templates[((t + 1) * template_len) - (i + 1)];
            norm_sums[t] += templates[(t * template_len) + i];
        }
    }
    for (i = 0; i < image_len; ++i)
    {
        image_ext[i] = image[i];
    }
    //  Compute ffts of template and image
    #pragma omp parallel sections
    {
        {fftwf_execute(pa); }
        #pragma omp section
        {fftwf_execute(pb); }
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
    fftwf_execute(px);
    //  Procedures for normalisation
    // Compute starting mean, will update this
    for (i=0; i < template_len; ++i){
        sum += (double) image[i];
    }
    mean = sum / template_len;

    // Compute starting standard deviation
    for (i=0; i < template_len; ++i){
        var += pow((double) image[i] - mean, 2) / (template_len);
    }
    stdev = sqrt(var);
    // Used for centering - taking only the valid part of the cross-correlation
    startind = template_len - 1;
    if (var >= ACCEPTED_DIFF) {
        for (t = 0; t < n_templates; ++t){
            float c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean) / stdev;
            status += set_ncc(t, 0, template_len, image_len, (float) c, used_chans, pad_array, ncc);
        }
        if (var <= WARN_DIFF){
            variance_warning[0] = 1;
        }
    }
    // Center and divide by length to generate scaled convolution
    for(i = 1; i < (image_len - template_len + 1); ++i){
        // Need to cast to double otherwise we end up with annoying floating
        // point errors when the variance is massive - collecting fp errors.
        new_samp = (double) image[i + template_len - 1];
        old_samp = (double) image[i - 1];
        old_mean = mean;
        mean = mean + (new_samp - old_samp) / template_len;
        var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (template_len);
        if (new_samp == (double) image[i + template_len - 2]) {
            flatline_count += 1;
         }
         else {
            flatline_count = 0;
        }
        stdev = sqrt(var);
        if (var >= ACCEPTED_DIFF && flatline_count < template_len - 1 && stdev * mean >= ACCEPTED_DIFF) {
            for (t = 0; t < n_templates; ++t){
                float c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean ) / stdev;
                status += set_ncc(t, i, template_len, image_len, (float) c, used_chans, pad_array, ncc);
            }
            if (var <= WARN_DIFF){
                variance_warning[0] += 1;
            }
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
    fftwf_cleanup_threads();

    free(template_ext);
    free(image_ext);

    return status;
}


int normxcorr_fftw(float *templates, long template_len, long n_templates,
                   float *image, long image_len, float *ncc, long fft_len,
                   int *used_chans, int *pad_array, int *variance_warning){
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
    long N2 = fft_len / 2 + 1;
    // All memory allocated with `fftw_malloc` to ensure 16-byte aligned
    float * template_ext = (float*) fftwf_malloc(fft_len * n_templates * sizeof(float));
    float * image_ext = (float*) fftwf_malloc(fft_len * sizeof(float));
    float * ccc = (float*) fftwf_malloc(fft_len * n_templates * sizeof(float));
    fftwf_complex * outa = (fftwf_complex*) fftwf_malloc(N2 * n_templates * sizeof(fftwf_complex));
    fftwf_complex * outb = (fftwf_complex*) fftwf_malloc(N2 * sizeof(fftwf_complex));
    fftwf_complex * out = (fftwf_complex*) fftwf_malloc(N2 * n_templates * sizeof(fftwf_complex));
    // Plan
    fftwf_plan pa = fftwf_plan_dft_r2c_2d(n_templates, fft_len, template_ext, outa, FFTW_ESTIMATE);
    fftwf_plan pb = fftwf_plan_dft_r2c_1d(fft_len, image_ext, outb, FFTW_ESTIMATE);
    fftwf_plan px = fftwf_plan_dft_c2r_2d(n_templates, fft_len, out, ccc, FFTW_ESTIMATE);

    // Initialise to zero
    memset(template_ext, 0, (size_t) fft_len * n_templates * sizeof(float));
    memset(image_ext, 0, (size_t) fft_len * sizeof(float));

    // Call the function to do the work
    // Note: forcing inner threads to 1 for now (could be passed from Python)
    status = normxcorr_fftw_main(templates, template_len, n_templates, image, image_len,
            ncc, fft_len, template_ext, image_ext, ccc, outa, outb, out, pa, pb, px,
            used_chans, pad_array, 1, variance_warning);

    // free memory and plans
    fftwf_destroy_plan(pa);
    fftwf_destroy_plan(pb);
    fftwf_destroy_plan(px);

    fftwf_free(out);
    fftwf_free(outa);
    fftwf_free(outb);
    fftwf_free(ccc);
    fftwf_free(template_ext);
    fftwf_free(image_ext);

    fftwf_cleanup();
    fftwf_cleanup_threads();

    return status;
}


int normxcorr_fftw_main(float *templates, long template_len, long n_templates,
                        float *image, long image_len, float *ncc, long fft_len,
                        float *template_ext, float *image_ext, float *ccc,
                        fftwf_complex *outa, fftwf_complex *outb, fftwf_complex *out,
                        fftwf_plan pa, fftwf_plan pb, fftwf_plan px, int *used_chans,
                        int *pad_array, int num_threads, int *variance_warning) {
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
    long N2 = fft_len / 2 + 1;
    long i, t, startind;
    int status = 0, unused_corr = 0;
    int * flatline_count = (int *) calloc(image_len - template_len + 1, sizeof(int));
    double *mean, *var;
    double new_samp, old_samp, sum=0.0;
    float * norm_sums = (float *) calloc(n_templates, sizeof(float));

    if (norm_sums == NULL) {
        printf("Error allocating norm_sums in normxcorr_fftw_main\n");
        return 1;
    }

    // zero padding - and flip template
    for (t = 0; t < n_templates; ++t){
        for (i = 0; i < template_len; ++i)
        {
            template_ext[(t * fft_len) + i] = templates[((t + 1) * template_len) - (i + 1)];
            norm_sums[t] += templates[(t * template_len) + i];
        }
    }
    for (i = 0; i < image_len; ++i)
    {
        image_ext[i] = image[i];
    }

    //  Compute fft of template
    fftwf_execute_dft_r2c(pa, template_ext, outa);
    
    // Compute fft of image
    fftwf_execute_dft_r2c(pb, image_ext, outb);

    //  Compute dot product
    #pragma omp parallel for num_threads(num_threads) private(i)
    for (t = 0; t < n_templates; ++t){
        for (i = 0; i < N2; ++i)
        {
            out[(t * N2) + i][0] = outa[(t * N2) + i][0] * outb[i][0] - outa[(t * N2) + i][1] * outb[i][1];
            out[(t * N2) + i][1] = outa[(t * N2) + i][0] * outb[i][1] + outa[(t * N2) + i][1] * outb[i][0];
        }
    }

    //  Compute inverse fft
    fftwf_execute_dft_c2r(px, out, ccc);
 
    // Allocate mean and var arrays
    mean = (double*) malloc((image_len - template_len + 1) * sizeof(double));
    if (mean == NULL) {
        printf("Error allocating mean in normxcorr_fftw_main\n");
        free(norm_sums);
        return 1;
    }
    var = (double*) malloc((image_len - template_len + 1) * sizeof(double));
    if (var == NULL) {
        printf("Error allocating var in normxcorr_fftw_main\n");
        free(norm_sums);
        free(mean);
        return 1;
    }
    
    //  Procedures for normalisation
    // Compute starting mean, will update this
    sum = 0.0;
    for (i=0; i < template_len; ++i){
        sum += (double) image[i];
    }
    mean[0] = sum / template_len;

    // Compute starting standard deviation
    sum = 0.0;
    for (i=0; i < template_len; ++i){
        sum += pow((double) image[i] - mean[0], 2) / (template_len);
    }
    var[0] = sum;

    // Used for centering - taking only the valid part of the cross-correlation
    startind = template_len - 1;
    if (var[0] >= ACCEPTED_DIFF) {
        double stdev = sqrt(var[0]);
        for (t = 0; t < n_templates; ++t){
            double c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean[0]);
            c /= stdev;
            status += set_ncc(t, 0, template_len, image_len, (float) c, used_chans, pad_array, ncc);

        }
        if (var[0] <= WARN_DIFF){
            variance_warning[0] = 1;
        }
    } else {
        unused_corr = 1;
    }
 
    // pre-compute the mean and var so we can parallelise the calculation
    for(i = 1; i < (image_len - template_len + 1); ++i){
        // Need to cast to double otherwise we end up with annoying floating
        // point errors when the variance is massive - collecting fp errors.
        new_samp = (double) image[i + template_len - 1];
        old_samp = (double) image[i - 1];
        mean[i] = mean[i - 1] + (new_samp - old_samp) / template_len;
        var[i] = var[i - 1] + (new_samp - old_samp) * (new_samp - mean[i] + old_samp - mean[i - 1]) / (template_len);
        if (new_samp == (double) image[i + template_len - 2]) {
            flatline_count[i] = flatline_count[i - 1] + 1;
        }
        else {
            flatline_count[i] = 0;
        }
    }

    // Center and divide by length to generate scaled convolution
    #pragma omp parallel for reduction(+:status,unused_corr) num_threads(num_threads) private(t)
    for(i = 1; i < (image_len - template_len + 1); ++i){
        if (var[i] >= ACCEPTED_DIFF && flatline_count[i] < template_len - 1) {
            double stdev = sqrt(var[i]);
            double meanstd = fabs(mean[i] * stdev);
            if (meanstd >= ACCEPTED_DIFF){
                for (t = 0; t < n_templates; ++t){
                    double c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean[i]);
                    c /= stdev;
                    status += set_ncc(t, i, template_len, image_len, (float) c, used_chans, pad_array, ncc);
                }
            }
            else {
                unused_corr = 1;
            }
            if (var[i] <= WARN_DIFF){
                variance_warning[0] += 1;
            }
        } else {
            unused_corr = 1;
        }
    }
    if (unused_corr == 1){
        if (status == 0){
            status = 999;
        }
    }

    //  Clean up
    free(norm_sums);
    free(mean);
    free(var);
    free(flatline_count);
    return status;
}


static inline int set_ncc(long t, long i, long template_len, long image_len, float value, int *used_chans, int *pad_array, float *ncc) {

    int status = 0;

    if (used_chans[t] && (i >= pad_array[t])) {
        size_t ncc_index = t * ((size_t) image_len - template_len + 1) + i - pad_array[t];

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
        // prev_ncc = ncc[ncc_index];
        #pragma omp atomic
        ncc[ncc_index] += value;
    }

    return status;
}

void free_fftwf_arrays(int size, float **template_ext, float **image_ext, float **ccc,
        fftwf_complex **outa, fftwf_complex **outb, fftwf_complex **out) {
    int i;

    /* free memory */
    for (i = 0; i < size; i++) {
        fftwf_free(template_ext[i]);
        fftwf_free(image_ext[i]);
        fftwf_free(ccc[i]);
        fftwf_free(outa[i]);
        fftwf_free(outb[i]);
        fftwf_free(out[i]);
    }
    free(template_ext);
    free(image_ext);
    free(ccc);
    free(outa);
    free(outb);
    free(out);
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


int multi_normxcorr_fftw(float *templates, long n_templates, long template_len, long n_channels,
        float *image, long image_len, float *ncc, long fft_len, int *used_chans, int *pad_array,
        int num_threads_outer, int num_threads_inner, int *variance_warning) {
    int i;
    int r=0;
    size_t N2 = (size_t) fft_len / 2 + 1;
    float **template_ext = NULL;
    float **image_ext = NULL;
    float **ccc = NULL;
    int * results = (int *) calloc(n_channels, sizeof(int));
    fftwf_complex **outa = NULL;
    fftwf_complex **outb = NULL;
    fftwf_complex **out = NULL;
    fftwf_plan pa, pb, px;

    #ifdef N_THREADS
    /* num_threads_outer cannot be greater than the number of channels */
    num_threads_outer = (num_threads_outer > n_channels) ? n_channels : num_threads_outer;

    /* Outer loop parallelism seems to cause issues on OSX */
    if (OUTER_SAFE != 1 && num_threads_outer > 1){
        printf("WARNING\tMULTI_NORMXCORR_FFTW\tOuter loop threading disabled for this system\n");
        num_threads_inner *= num_threads_outer;
        printf("WARNING\tMULTI_NORMXCORR_FFTW\tSetting inner threading to %i and outer threading to 1\n", num_threads_inner);
        num_threads_outer = 1;
    }
    if (num_threads_inner > 1) {
        /* initialise FFTW threads */
        fftwf_init_threads();
        fftwf_plan_with_nthreads(num_threads_inner);

        if (num_threads_outer > 1) {
            /* explicitly enable nested OpenMP loops */
            omp_set_nested(1);
        }
    }

    /* warn if the total number of threads is higher than the number of cores */
    if (num_threads_outer * num_threads_inner > N_THREADS) {
        printf("Warning: requesting more threads than available - this could negatively impact performance\n");
    }
    #else
    /* threading/OpenMP is disabled */
    num_threads_outer = 1;
    num_threads_inner = 1;
    #endif

    /* allocate memory for all threads here */
    template_ext = (float**) malloc(num_threads_outer * sizeof(float*));
    if (template_ext == NULL) {
        printf("Error allocating template_ext\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    image_ext = (float**) malloc(num_threads_outer * sizeof(float*));
    if (image_ext == NULL) {
        printf("Error allocating image_ext\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    ccc = (float**) malloc(num_threads_outer * sizeof(float*));
    if (ccc == NULL) {
        printf("Error allocating ccc\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    outa = (fftwf_complex**) malloc(num_threads_outer * sizeof(fftwf_complex*));
    if (outa == NULL) {
        printf("Error allocating outa\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    outb = (fftwf_complex**) malloc(num_threads_outer * sizeof(fftwf_complex*));
    if (outb == NULL) {
        printf("Error allocating outb\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }
    out = (fftwf_complex**) malloc(num_threads_outer * sizeof(fftwf_complex*));
    if (out == NULL) {
        printf("Error allocating out\n");
        free_fftwf_arrays(0, template_ext, image_ext, ccc, outa, outb, out);
        return -1;
    }

    // All memory allocated with `fftw_malloc` to ensure 16-byte aligned.
    for (i = 0; i < num_threads_outer; i++) {
        /* initialise all to NULL so that freeing on error works */
        template_ext[i] = NULL;
        image_ext[i] = NULL;
        ccc[i] = NULL;
        outa[i] = NULL;
        outb[i] = NULL;
        out[i] = NULL;

        /* allocate template_ext arrays */
        template_ext[i] = (float*) fftwf_malloc((size_t) fft_len * n_templates * sizeof(float));
        if (template_ext[i] == NULL) {
            printf("Error allocating template_ext[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate image_ext arrays */
        image_ext[i] = (float*) fftwf_malloc(fft_len * sizeof(float));
        if (image_ext[i] == NULL) {
            printf("Error allocating image_ext[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate ccc arrays */
        ccc[i] = (float*) fftwf_malloc((size_t) fft_len * n_templates * sizeof(float));
        if (ccc[i] == NULL) {
            printf("Error allocating ccc[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate outa arrays */
        outa[i] = (fftwf_complex*) fftwf_malloc((size_t) N2 * n_templates * sizeof(fftwf_complex));
        if (outa[i] == NULL) {
            printf("Error allocating outa[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate outb arrays */
        outb[i] = (fftwf_complex*) fftwf_malloc((size_t) N2 * sizeof(fftwf_complex));
        if (outb[i] == NULL) {
            printf("Error allocating outb[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }

        /* allocate out arrays */
        out[i] = (fftwf_complex*) fftwf_malloc((size_t) N2 * n_templates * sizeof(fftwf_complex));
        if (out[i] == NULL) {
            printf("Error allocating out[%d]\n", i);
            free_fftwf_arrays(i + 1, template_ext, image_ext, ccc, outa, outb, out);
            return -1;
        }
    }

    // We create the plans here since they are not thread safe.
    pa = fftwf_plan_dft_r2c_2d(n_templates, fft_len, template_ext[0], outa[0], FFTW_ESTIMATE);
    pb = fftwf_plan_dft_r2c_1d(fft_len, image_ext[0], outb[0], FFTW_ESTIMATE);
    px = fftwf_plan_dft_c2r_2d(n_templates, fft_len, out[0], ccc[0], FFTW_ESTIMATE);

    /* loop over the channels */
    #pragma omp parallel for num_threads(num_threads_outer)
    for (i = 0; i < n_channels; ++i){
        int tid = 0; /* each thread has its own workspace */

        #ifdef N_THREADS
        /* get the id of this thread */
        tid = omp_get_thread_num();
        #endif
        /* initialise memory to zero */
        memset(template_ext[tid], 0, (size_t) fft_len * n_templates * sizeof(float));
        memset(image_ext[tid], 0, (size_t) fft_len * sizeof(float));

        /* call the routine */
        results[i] = normxcorr_fftw_main(&templates[(size_t) n_templates * template_len * i], template_len,
                                 n_templates, &image[(size_t) image_len * i], image_len, ncc, fft_len,
                                 template_ext[tid], image_ext[tid], ccc[tid], outa[tid], outb[tid], out[tid],
                                 pa, pb, px, &used_chans[(size_t) i * n_templates],
                                 &pad_array[(size_t) i * n_templates], num_threads_inner, &variance_warning[i]);
    }

    // Conduct error handling
    for (i = 0; i < n_channels; ++i){
        if (results[i] != 999 && results[i] != 0){
            // Some error internally, must catch this
            r += results[i];
        } else if (results[i] == 999 && r == 0){
            // First time unused correlation raised and no prior errors
            r = results[i];
        } else if (r == 999 && results[i] == 999){
            // Unused correlations raised multiple times
            r = 999;
        } else if (r == 999 && results[i] != 999){
            // Some error internally.
            r += results[i];
        } else if (r != 0){
            // Any other error
            r += results[i];
        }
    }
    free(results);
    /* free fftw memory */
    free_fftwf_arrays(num_threads_outer, template_ext, image_ext, ccc, outa, outb, out);
    fftwf_destroy_plan(pa);
    fftwf_destroy_plan(pb);
    fftwf_destroy_plan(px);
    if (num_threads_inner > 1) {
        fftwf_cleanup_threads();
    }
    fftwf_cleanup();

    return r;
}
