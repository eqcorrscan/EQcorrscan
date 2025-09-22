/*
 * =====================================================================================
 *
 *       Filename:  distance_cluster.c
 *
 *        Purpose:  Routines for calculating distance matrices between points
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
#if defined(__linux__) || defined(__linux) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    #include <omp.h>
    #ifndef N_THREADS
        #define N_THREADS omp_get_max_threads()
    #endif
#endif

#define EARTH_RADIUS 6371.009

float dist_calc(float, float, float, float, float, float);

int distance_array(float, float, float, float*, float*, float*, long, float*, int);

int distance_matrix(float*, float*, float*, long, float*, int);

int remove_unclustered(float*, float*, float*, long, unsigned char*, float, int);

float dist_calc(float lat1, float lon1, float depth1, float lat2, float lon2, float depth2){
//    Function to calculate the distance in km between two points.
//
//    Uses the
//    `haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_
//    to calculate great circle distance at the Earth's surface, then uses
//    trig to include depth.
//
//    :type lat1: float - Latitude of point 1 in radians
//    :type lon1: float - Longitude of point 1 in radians
//    :type depth1: float - Depth of point 1 in km (positive down)
//    :type lat2: float - Latitude of point 2 in radians
//    :type lon2: float - Longitude of point 2 in radians
//    :type depth2: float - Depth of point 2 in km (positive down)
//    :type distance: float - Distance between two points (output, km)


    double central_angle;
    float dlat, dlong, ddepth, distance;

    dlat = lat1 - lat2;
    dlong = lon1 - lon2;
    ddepth = depth1 - depth2;

    central_angle = 2 * asin(sqrt(pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlong / 2), 2)));

    distance = EARTH_RADIUS * central_angle;
    distance = sqrt(pow(distance, 2) + pow(ddepth, 2));
    return distance;
}

int distance_array(float lat1, float lon1, float depth1, float *latitudes,
                   float *longitudes, float *depths, long n_locs,
                   float *dist_array, int n_threads){
    /* Calculate the distances of all locations relative to one location.
    *
    *  :type lat1: Latitude of core location in radians
    *  :type lon1: Longitude of core location in radians
    *  :type depth1: Depth of core location in km (positive down)
    *  :type latitudes: Array of floats of latitudes in radians
    *  :type longitudes: Array of floats of longitudes in radians
    *  :type depths: Array of floats of depths in km (positive down)
    *  :type n_locs: Number of locations
    *  :type dist_array: Array of floats of for output - should be initialized as zeros, and should be n_locs long
    */
    int out = 0;
    long n;

    #pragma omp parallel for num_threads(n_threads)
    for (n = 0; n < n_locs; ++n){
        dist_array[n] = dist_calc(
            lat1, lon1, depth1, latitudes[n], longitudes[n], depths[n]);
    }
    return out;
}


int distance_matrix(float *latitudes, float *longitudes, float *depths, long n_locs,
                    float *dist_mat, int n_threads){
    /* Calculate the distance matrix for a set of locations.
    *
    *  :type latitudes: Array of floats of latitudes in radians
    *  :type longitudes: Array of floats of longitudes in radians
    *  :type depths: Array of floats of depths in km (positive down)
    *  :type n_locs: Number of locations
    *  :type dist_mat: Array of floats of for output - should be initialized as zeros, and should be n_locs * n_locs
    */
    int out = 0;
    long n;

    #pragma omp parallel for num_threads(n_threads)
    for (n = 0; n < n_locs * (n_locs + 1) / 2; ++n){
        long j = n / (n_locs + 1), i = n % (n_locs + 1);
        if (i > j){
            j = n_locs - j - 1;
            i = n_locs - i;
        }
        dist_mat[(i * n_locs) + j] = dist_calc(
            latitudes[i], longitudes[i], depths[i], latitudes[j], longitudes[j], depths[j]);
    }
    return out;
}

int remove_unclustered(float *latitudes, float *longitudes, float *depths, long n_locs,
                       unsigned char *mask, float distance_cutoff, int n_threads){
    /* Check whether locations have any other locations within distance_cutoff and return 0 if not and 1 if true.
    *
    *  :type latitudes: Array of floats of latitudes in radians
    *  :type longitudes: Array of floats of longitudes in radians
    *  :type depths: Array of floats of depths in km (positive down)
    *  :type n_locs: Int: Number of locations
    *  :type mask: Array of uint8 which will be filled as bools - should be initialised as zeros
    *  :type distance_cutoff: float, cutoff distance in km
    *  :type n_threads: int Number of threads to parallel over
    */
    int out = 0;
    long i;

    #pragma omp parallel for num_threads(n_threads)
    for (i = 0; i < n_locs; ++i){
        long j;
        float dist;
        if (mask[i] != 0){continue;}
        for (j = 0; j < n_locs; ++j){
            dist = dist_calc(
                latitudes[j], longitudes[j], depths[j],
                latitudes[i], longitudes[i], depths[i]);
            if (j != i && dist < distance_cutoff){
                mask[i] = 1;
                mask[j] = 1;
                // No more calculation needs to be done.
                break;
            }
        }
    }
    return out;
}
