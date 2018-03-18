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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

 // Prototypes
int find_peaks(float*, float*, int, float, float, unsigned int*);

// Functions
// Longs could be unsigned ints...
int find_peaks(float *arr, float *indexes, int len, float thresh, float trig_int,
               unsigned int *out){
    // Takes a sorted array an the indexes
    int i, j, keep;
    float step;
    // Take first (highest) peak
    out[0] = 1;
    for (i = 1; i < len; ++i){
        keep = 1;
        // Threshold is for absolute values
//        if (-1 * thresh < arr[i] && arr[i] < thresh){
//            break;
//        }
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
