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

#include <libutils.h>


static inline int set_ncc(
    long t, long i, int chan, int n_chans, long template_len, long image_len,
    float value, int *used_chans, int *pad_array, float *weight_array, float *ncc,
    int stack_option);

// Functions

// Single-channel functions
int normxcorr_fftw_threaded(float *templates, long template_len, long n_templates,
                            float *image, long image_len, float *ncc, long fft_len,
                            int *used_chans, int *pad_array, int *variance_warning,
                            int *missed_corr) {
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

   Note - no weighting is applied
  */
    long N2 = fft_len / 2 + 1;
    long i, t, startind;
    int status = 0;
    int flatline_count = 0, unused_corr = 0;
    double mean, stdev, old_mean, new_samp, old_samp, var=0.0, sum=0.0;
    float * norm_sums = (float *) calloc(n_templates, sizeof(float));
    float * template_ext = (float *) calloc(fft_len * n_templates, sizeof(float));
    float * image_ext = (float *) calloc(fft_len, sizeof(float));
    float * ccc = (float *) fftwf_malloc(sizeof(float) * fft_len * n_templates);
    float * weight_array = (float *) calloc(n_templates, sizeof(float));
    fftwf_complex * outa = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2 * n_templates);
    fftwf_complex * outb = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2);
    fftwf_complex * out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N2 * n_templates);

    // Initialize weights to one
    // memset(weight_array, 1, (size_t) n_templates * sizeof(float));
    for (t = 0; t < n_templates; ++t){
        weight_array[t] = 1;
    }

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
            status += set_ncc(t, 0, 0, 1, template_len, image_len, (float) c, used_chans, pad_array, weight_array, ncc, 0);
        }
        if (var <= WARN_DIFF){
            variance_warning[0] = 1;
        }
    } else {unused_corr += 1;}
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
                status += set_ncc(t, i, 0, 1, template_len, image_len, (float) c, used_chans, pad_array, weight_array, ncc, 0);
            }
            if (var <= WARN_DIFF){
                variance_warning[0] += 1;
            }
        } else {unused_corr += 1;}
    }
    missed_corr[0] = unused_corr;
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
                   int *used_chans, int *pad_array, int *variance_warning,
                   int *missed_corr){
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

    Weights are not applied.
  */
    int status = 0, t;
    long N2 = fft_len / 2 + 1;
    // All memory allocated with `fftw_malloc` to ensure 16-byte aligned
    float * template_ext = (float*) fftwf_malloc(fft_len * n_templates * sizeof(float));
    float * image_ext = (float*) fftwf_malloc(fft_len * sizeof(float));
    float * ccc = (float*) fftwf_malloc(fft_len * n_templates * sizeof(float));
    float * weight_array = (float*) calloc(n_templates, sizeof(float));
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
    // Initialize weights to one
    for (t = 0; t < n_templates; ++t){
        weight_array[t] = 1;
    }
    // memset(weight_array, 1, (size_t) n_templates * sizeof(float));

    // Call the function to do the work
    // Note: forcing inner threads to 1 for now (could be passed from Python)
    status = normxcorr_fftw_main(
        templates, template_len, n_templates, image, image_len, 0, 1, ncc,
        fft_len, template_ext, image_ext, ccc, outa, outb, out, pa, pb, px,
        used_chans, pad_array, weight_array, 1, variance_warning, missed_corr, 0);

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


// Functions to multiple channels
int normxcorr_fftw_main(
    float *templates, long template_len, long n_templates, float *image,
    long image_len, int chan, int n_chans, float *ncc, long fft_len,
    float *template_ext, float *image_ext, float *ccc, fftwf_complex *outa,
    fftwf_complex *outb, fftwf_complex *out, fftwf_plan pa, fftwf_plan pb,
    fftwf_plan px, int *used_chans, int *pad_array, float *weight_array,
    int num_threads, int *variance_warning, int *missed_corr, int stack_option) {
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  for a single-channel
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
    templates:      Template signals
    template_len:   Length of template
    n_templates:    Number of templates (n0)
    image:          Image signal (to scan through)
    image_len:      Length of image
    ncc:            Output for cross-correlation - should be pointer to memory.
                    Shapes and output determined by stack_option:
        1:          Output stack correlograms, ncc must be
                    (n_templates x image_len - template_len + 1) long.
        0:          Output individual channel correlograms, ncc must be
                    (n_templates x image_len - template_len + 1) long and initialised
                    to zero before passing into this function.
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
    used_chans:     Array to fill with number of channels used per template - must
                    be n_templates long
    pad_array:      Array of pads, should be n_templates long
    weight_array:   Array of weights, should be n_templates long
    num_threads:    Number of threads to parallel internal calculations over
    variance_warning: Pointer to array to store warnings for variance issues
    missed_corr:    Pointer to array to store warnings for unused correlations
    stack_option:   Whether to stack correlograms (1) or leave as individual channels (0),
  */
//    double tic, toc, super_tic, super_toc;
    long i, t, chunk, n_chunks, chunk_len, startind = template_len - 1, offset, step_len;
    int status = 0, N2 = fft_len / 2 + 1, unused_corr = 0;
    float * norm_sums = (float *) calloc(n_templates, sizeof(float));
    int * flatline_count = (int *) calloc(image_len - template_len + 1, sizeof(int));
    double * mean = (double*) malloc((image_len - template_len + 1) * sizeof(double));
    double * var = (double*) malloc((image_len - template_len + 1) * sizeof(double));

    if (norm_sums == NULL) {
        printf("ERROR: Error allocating norm_sums in normxcorr_fftw_main\n");
        return -1;
    }

    if (stack_option > 1) {
        printf("ERROR: stack_option %i is not known\n", stack_option);
        return -1;
    }

    // zero padding - and flip template
    for (t = 0; t < n_templates; ++t){
        for (i = 0; i < template_len; ++i)
        {
            template_ext[(t * fft_len) + i] = templates[((t + 1) * template_len) - (i + 1)];
            norm_sums[t] += templates[(t * template_len) + i];
        }
    }

    //  Compute fft of template
//    tic = omp_get_wtime();
    fftwf_execute_dft_r2c(pa, template_ext, outa);
//    toc = omp_get_wtime();
//    printf("Template ffts took \t\t%f s\n", toc - tic);

    if (fft_len >= image_len){
        n_chunks = 1;
        chunk_len = image_len;
        step_len = chunk_len;
    } else {
        chunk_len = fft_len;
        step_len = fft_len - (template_len - 1);
        n_chunks = (image_len - chunk_len) / step_len + ((image_len - chunk_len) % step_len > 0);
        if (n_chunks * step_len < image_len){n_chunks += 1;}
    }

    // Procedures for normalisation
    // TODO: Run this as a parallel section
//    tic = omp_get_wtime();
    running_mean_var(mean, var, flatline_count, image, image_len, template_len);
//    toc = omp_get_wtime();
//    printf("Running mean took \t\t%f s\n", toc - tic);

//    super_tic = omp_get_wtime();
    for (chunk = 0; chunk < n_chunks; ++chunk){
        offset = chunk * step_len;
        if (offset + chunk_len > image_len){
            chunk_len = image_len - offset;}

        memset(image_ext, 0, (size_t) fft_len * sizeof(float));
        for (i = 0; i < chunk_len; ++i){image_ext[i] = image[offset + i];}

        // Forward FFT
//        tic = omp_get_wtime();
        fftwf_execute_dft_r2c(pb, image_ext, outb);
//        toc = omp_get_wtime();
//        printf("Chunk FFT took \t\t%f s\n", toc - tic);

        // Dot product
//        tic = omp_get_wtime();
        #pragma omp parallel for num_threads(num_threads) private(i)
        for (t = 0; t < n_templates; ++t){
            for (i = 0; i < N2; ++i)
            {
                out[(t * N2) + i][0] = outa[(t * N2) + i][0] * outb[i][0] - outa[(t * N2) + i][1] * outb[i][1];
                out[(t * N2) + i][1] = outa[(t * N2) + i][0] * outb[i][1] + outa[(t * N2) + i][1] * outb[i][0];
            }
        }
//        toc = omp_get_wtime();
//        printf("Dot product took \t\t%f s\n", toc - tic);

        //  Compute inverse fft
//        tic = omp_get_wtime();
        fftwf_execute_dft_c2r(px, out, ccc);
//        toc = omp_get_wtime();
//        printf("Inverse FFT took \t\t%f s\n", toc - tic);

        // Centre and normalise

//        tic = omp_get_wtime();
        if (var[offset] >= ACCEPTED_DIFF) {
            double stdev = sqrt(var[offset]);
            for (t = 0; t < n_templates; ++t){
                double c = ((ccc[(t * fft_len) + startind] / (fft_len * n_templates)) - norm_sums[t] * mean[offset]);
                c /= stdev;
                status += set_ncc(t, offset, chan, n_chans, template_len, image_len,
                                  (float) c, used_chans, weight, array, pad_array, ncc, stack_option);
            }
            if (var[offset] <= WARN_DIFF){
                variance_warning[0] = 1;
            }
        } else {
            unused_corr += 1;
        }

        // Center and divide by length to generate scaled convolution
        #pragma omp parallel for reduction(+:status,unused_corr) num_threads(num_threads) private(t)
        for(i = 1; i < (chunk_len - template_len + 1); ++i){
            if (var[offset + i] >= ACCEPTED_DIFF && flatline_count[offset + i] < template_len - 1) {
                double stdev = sqrt(var[offset + i]);
                double meanstd = fabs(mean[offset + i] * stdev);
                if (meanstd >= ACCEPTED_DIFF){
                    for (t = 0; t < n_templates; ++t){
                        double c = ((ccc[(t * fft_len) + i + startind] / (fft_len * n_templates)) - norm_sums[t] * mean[offset + i]);
                        c /= stdev;
                        status += set_ncc(t, i + offset, chan, n_chans, template_len,
                                          image_len, (float) c, used_chans, weight_array,
                                          pad_array, ncc, stack_option);
                    }
                }
                else {
                    unused_corr += 1;
                }
                if (var[offset + i] <= WARN_DIFF){
                    variance_warning[0] += 1;
                }
            } else {
                unused_corr += 1;
            }
        }
        missed_corr[0] += unused_corr;
//        toc = omp_get_wtime();
//        printf("Normalising took \t\t%f s\n", toc - tic);
    }
//    super_toc = omp_get_wtime();
//    printf("Looping over chunks took \t\t%f s\n", super_toc - super_tic);
    free(mean);
    free(var);
    free(flatline_count);
    free(norm_sums);
    return status;
}


int running_mean_var(
    double *mean, double *var, int *flatline_count, float *image, long image_len,
    long template_len)
{
    long i;
    double sum, new_samp, old_samp;

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
    return 0;

}

static inline int set_ncc(
    long t, long i, int chan, int n_chans, long template_len, long image_len,
    float value, int *used_chans, int *pad_array, float *weight_array,
    float *ncc, int stack_option){

    int status = 0;

    if (used_chans[t] && (i >= pad_array[t])) {
        size_t ncc_index = (t * n_chans * ((size_t) image_len - template_len + 1)) +
            (chan * ((size_t) image_len - template_len + 1) + i - pad_array[t]);

        if (isnanf(value)) {
            // set NaNs to zero
            value = 0.0;
        }
        else if (fabsf(value) > 1.01) {
            // this will raise a warning when we return to Python
            printf("WARNING: Correlation value=%f:\tncc_index: %ld\ttemplate: %ld\tchannel: %i\tindex: %ld\nSETTING TO ZERO.\n",
                   value, ncc_index, t, chan, i);
            value = 0.0;
            status = 1;
        }
        else if (value > 1.0) {
            value = 1.0;
        }
        else if (value < -1.0) {
            value = -1.0;
        }
        // Apply weight after checks to allow for weights > 1.0
        value *= weight_array[t];
        if (stack_option == 1){
            #pragma omp atomic
            ncc[ncc_index] += value;
        } else if (stack_option == 0){ncc[ncc_index] = value;}
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
                         float *image, long image_len, float *ncc, long fft_len, int *used_chans,
                         int *pad_array, float *weight_array, int num_threads_inner, int num_threads_outer,
                         int *variance_warning, int *missed_corr, int stack_option)
    {
    int i, chan, n_chans;
    int r = 0;
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

    // Ensure stack-option is within supported range
    if (stack_option > 1){
        printf("Stack option %i unsupported, returning\n", stack_option);
        return -1;
    }

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
        printf("WARNING: Requested %d inner and %d outer = %d total, but %d are available\n", num_threads_inner, num_threads_outer, num_threads_outer * num_threads_inner, N_THREADS);
        printf("WARNING: requesting more threads than available - this could negatively impact performance\n");
    }
    #else
    /* threading/OpenMP is disabled */
    num_threads_outer = 1;
    num_threads_inner = 1;
    #endif
//    printf("Using %d outer threads and %d inner threads\n", num_threads_outer, num_threads_inner);

    /* Check that stack-type is within range (0-1) */
    if (stack_option > 1) {
        printf("ERROR: stack_option %i is not supported\n", stack_option);
        return -1;
    }

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
        // Done internally now. memset(image_ext[tid], 0, (size_t) fft_len * sizeof(float));

        if (stack_option == 1){
            chan = 0;
            n_chans = 1;
        } else {
            chan = i;
            n_chans = n_channels;
        }
        /* call the routine */
        results[i] = normxcorr_fftw_main(&templates[(size_t) n_templates * template_len * i], template_len,
                                 n_templates, &image[(size_t) image_len * i], image_len, chan, n_chans, ncc, fft_len,
                                 template_ext[tid], image_ext[tid], ccc[tid], outa[tid], outb[tid], out[tid],
                                 pa, pb, px, &used_chans[(size_t) i * n_templates],
                                 &pad_array[(size_t) i * n_templates],
                                 &weight_array[(size_t) i * n_templates],
                                 num_threads_inner, &variance_warning[i],
                                 &missed_corr[i], stack_option);
        if (results[i] != 0){
            printf("WARNING: %i out-of-range correlations on channel %i\n", results[i], i);
        }
    }

    // Conduct error handling
    for (i = 0; i < n_channels; ++i){
        r += results[i];
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
