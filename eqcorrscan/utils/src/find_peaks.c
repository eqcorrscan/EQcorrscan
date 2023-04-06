/*
 * =====================================================================================
 *
 *       Filename:  find_peaks.c
 *
 *        Purpose:  Routines for finding peaks in noisy data
 *
 *        Created:  03/07/17 02:25:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Calum Chamberlain
 *   Organization:  EQcorrscan
 *      Copyright:  EQcorrscan developers.
 *        License:  GNU Lesser General Public License, Version 3
 *                  (https://www.gnu.org/copyleft/lesser.html)
 *
 * =====================================================================================
 */
#include <libutils.h>

// Decluster in distance and time
int decluster_dist_time(float *arr, long *indexes, float *distances,
                        long len, float thresh, long trig_int,
                        float dist_thresh, unsigned int *out){
    // Takes a sorted array, with indexes as the time between events, and the
    // distances as a distance matrix sorted in the same way.
    long i, j, step, distance_index;
    int keep;

    if (fabs(arr[0]) < thresh){return 0;}

    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
        if (fabs(arr[i]) < thresh){
            break;
        }
        for (j = 0; j < i; ++j){
            distance_index = (i * len) + j;
            step = labs(indexes[i] - indexes[j]);

            if (trig_int >= step && out[j] == 1 && distances[distance_index] < dist_thresh){
                keep = 0;
                break;
            }
        }
        out[i] = keep;
    }
    return 0;
}

int decluster_dist_time_ll(float *arr, long long *indexes, float *distances,
                           long long len, float thresh, long long trig_int,
                           float dist_thresh, unsigned int *out){
    // Takes a sorted array, with indexes as the time between events, and the
    // distances as a distance matrix sorted in the same way.
    long long i, j, step, distance_index;
    int keep;

    if (fabs(arr[0]) < thresh){return 0;}

    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
        if (fabs(arr[i]) < thresh){
            break;
        }
        for (j = 0; j < i; ++j){
            distance_index = (i * len) + j;
            step = llabs(indexes[i] - indexes[j]);
            if (trig_int >= step && out[j] == 1 && distances[distance_index] < dist_thresh){
                keep = 0;
                break;
            }
        }
        out[i] = keep;
    }
    return 0;
}


// Functions for long longs
int decluster_ll(float *arr, long long *indexes, long long len,
                 float thresh, long long trig_int, unsigned int *out){
    // Takes a sorted array and the indexes
    long long i, j, step;
    int keep;

    if (fabs(arr[0]) < thresh){return 0;}

    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
        if (fabs(arr[i]) < thresh){
            break;
        }
        for (j = 0; j < i; ++j){
            step = llabs(indexes[i] - indexes[j]);
            if (trig_int >= step && out[j] == 1){
                keep = 0;
                break;
            }
        }
        if (keep == 1){
            out[i] = 1;
         }
        else {out[i] = 0;}
    }
    return 0;
}

int multi_decluster_ll(float *arr, long long *indices,
                       long long *lengths, int n, float *thresholds,
                       long long trig_int, unsigned int *out, int threads){
    int i, ret_val = 0;
    long long * start_inds = (long long *) calloc(n, sizeof(long long));
    long long start_ind = 0;

    for (i = 0; i < n; ++i){
        start_inds[i] = start_ind;
        start_ind += lengths[i];
    }

    #ifdef N_THREADS
    if (threads > N_THREADS){
        printf("MULTI-DECLUSTER-LL: Setting threads to %i. OMP found %i threads\n", N_THREADS, omp_get_max_threads());
        threads = N_THREADS;
    }
    #else
    threads = 1;
    #endif

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n; ++i){
        ret_val += decluster_ll(
            &arr[start_inds[i]], &indices[start_inds[i]], lengths[i], thresholds[i],
            trig_int, &out[start_inds[i]]);
    }
    free(start_inds);
    return ret_val;
}


// Functions for longs - should be the same logic as above
int decluster(float *arr, long *indexes, long len,
              float thresh, long trig_int, unsigned int *out){
    // Takes a sorted array and the indexes
    long i, j, step;
    int keep;

    if (fabs(arr[0]) < thresh){return 0;}

    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
        if (fabs(arr[i]) < thresh){
            break;
        }
        for (j = 0; j < i; ++j){
            step = labs(indexes[i] - indexes[j]);
            if (trig_int >= step && out[j] == 1){
                keep = 0;
                break;
            }
        }
        if (keep == 1){
            out[i] = 1;
         }
        else {out[i] = 0;}
    }
    return 0;
}

int multi_decluster(float *arr, long *indices,
                    long *lengths, int n, float *thresholds,
                    long trig_int, unsigned int *out, int threads){
    int i, ret_val = 0;
    long * start_inds = (long *) calloc(n, sizeof(long));
    long start_ind = 0;

    for (i = 0; i < n; ++i){
        start_inds[i] = start_ind;
        start_ind += lengths[i];
    }
    #ifdef N_THREADS
    if (threads > N_THREADS){
//        printf("MULTI-DECLUSTER: Setting threads to %i. OMP found %i threads\n", N_THREADS, omp_get_max_threads());
//        threads = N_THREADS;
        printf("MULTI-DECLUSTER\tMore threads requested than available (%i > %i). Caution\n", threads, N_THREADS);
;    }
    #else
    threads = 1;
    #endif

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n; ++i){
        ret_val += decluster(
            &arr[start_inds[i]], &indices[start_inds[i]], lengths[i], thresholds[i],
            trig_int, &out[start_inds[i]]);
    }
    free(start_inds);
    return ret_val;
}


int find_peaks(float *arr, long len, float thresh, unsigned int *peak_positions){
    // Find peaks in noisy data above some threshold and at-least
    // trig-int samples apart. Sets all other values in array to 0
    float prev_value = 0, value, next_value;
    long i;

    for (i = 0; i < len - 1; ++i){
        value = arr[i];
        next_value = arr[i + 1];
        if (fabs(value) > thresh){
            if ((next_value - value) * (prev_value - value) > 0){
                peak_positions[i] = 1;
            }
        }
        prev_value = value;
    }
    // Do separately for the last value
    i = len - 1;
    value = arr[i];
    next_value = 0;
    if (fabs(value) > thresh){
        if ((next_value - value) * (prev_value - value) > 0){
            peak_positions[i] = 1;
        }
    }
    return 0;
}
