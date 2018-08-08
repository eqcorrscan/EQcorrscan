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

int distance_matrix(float*, float*, float*, long, float*, int);

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

int distance_matrix(float *latitudes, float *longitudes, float *depths, long n_locs, float* dist_mat, int num_threads){
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

    #pragma omp parallel for num_threads(num_threads)
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