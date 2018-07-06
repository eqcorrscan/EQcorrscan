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


// Functions
// Longs could be unsigned ints...
int decluster(float *arr, long *indexes, int len, float thresh, long trig_int,
               unsigned int *out){
    // Takes a sorted array an the indexes
    int i, j, keep;
    float step;
    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
        if (fabs(arr[i]) < thresh){
            break;
        }
        for (j = 0; j < len; ++j){
            step = indexes[i] - indexes[j];
            if (step < 0){step *= -1;}
            if (i != j && trig_int > step && out[j] == 1){
                keep = 0;
                break;
            }
        }
        if (keep == 1){out[i] = 1;}
        else {out[i] = 0;}
    }
    return 0;
}


int find_peaks(float *arr, long len, float thresh, long trig_int, long number_of_peaks){
    // Find peaks in noisy data above some threshold and at-least
    // trig-int samples apart. Sets all other values in array to 0
    int above_threshold = 0, biggest_in_window_position = 0;
    float biggest_in_window = 0, prev_value = 0, value, next_value;
    int i;
    int * peak_positions = (int *) calloc(len, sizeof(int));
    long number_above_threshold = 0;

    // This tactic won't work if the threshold is low.
    for (i = 0; i < len - 1; ++i){
        value = arr[i];
        next_value = arr[i + 1];
        if (fabs(value) > thresh && i - biggest_in_window_position < trig_int){
            above_threshold = 1;
            number_above_threshold += 1;
            if (fabs(value) > fabs(biggest_in_window)){
                biggest_in_window = value;
                biggest_in_window_position = i;
            }
        } else if (fabs(value) > thresh && i - biggest_in_window_position > trig_int){
            // Start a new window
            above_threshold = 1;
            number_above_threshold += 1;
            if (prev_value < fabs(value) && next_value < fabs(value)){
                biggest_in_window = value;
                biggest_in_window_position = i;
            }
        } else if (above_threshold && fabs(value) < thresh){
            // We have dropped below the threshold
            above_threshold = 0;
            peak_positions[biggest_in_window_position] = 1;
            biggest_in_window = 0;
            biggest_in_window_position = 0;
        }
        prev_value = value;
    }
    // Do separately for the last value
    i = len;
    value = arr[-1];
    next_value = 0;
    if (fabs(value) > thresh && i - biggest_in_window_position < trig_int){
        above_threshold = 1;
        number_above_threshold += 1;
        if (fabs(value) > fabs(biggest_in_window)){
            biggest_in_window = value;
            biggest_in_window_position = i;
        }
    } else if (fabs(value) > thresh && i - biggest_in_window_position > trig_int){
        // Start a new window
        above_threshold = 1;
        number_above_threshold += 1;
        if (prev_value < fabs(value) && next_value < fabs(value)){
            biggest_in_window = value;
            biggest_in_window_position = i;
        }
    } else if (above_threshold && fabs(value) < thresh){
        // We have dropped below the threshold
        above_threshold = 0;
        peak_positions[biggest_in_window_position] = 1;
    }
    if (above_threshold){
        peak_positions[biggest_in_window_position] = 1;
    }
    if (number_above_threshold == len){
        return 1;
    }
    for (i = 0; i < len; ++i){
        if (peak_positions[i] != 1){
            arr[i] = 0;
            number_of_peaks += 1;
        }
    }
    free(peak_positions);
    return 0;
}


int multi_find_peaks(float *arr, long len, int n, float *thresholds,
                     long trig_int, long *number_of_peaks, int threads){
    int i, ret_val=0;

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n; ++i){
        ret_val += find_peaks(&arr[i * len], len, thresholds[i], trig_int, number_of_peaks[i]);
    }
    return ret_val;
}