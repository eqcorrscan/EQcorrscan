/*-------------------------------------------------------------------------
# Filename: sliding_xcorr.c
#  Purpose: Sliding window, normalized cross-correlation using openCV
#   Author: Calum John Chamberlain
#  Licence: LGPL, Verison 3
#--------------------------------------------------------------------------*/
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <float.h>

static PyObject *
xcorr(PyObject *self, PyObject *args)
//float *arr1, float *arr2, float *ccs, const int shift_len, const int arr1_len, const int arr2_len)
{
    int i;
    int shift;
    int min_len;
    int arr1_start;
    int arr1_end;
    int arr2_start;
    int arr2_end;
    float *arr1
    float *arr2
    float *ccs
    int shift_len
    int arr1_len
    int arr2_len

    // Extract the arguments and convert to useful C things
    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;


    printf("Array 1 length: %d \n", arr1_len);
    printf("Array 2 length: %d \n", arr2_len);
    if (arr1_len <= arr2_len)
    {
        min_len = arr1_len;
    }
    else
    {
        min_len = arr2_len;
    }
    printf("Minimum_length = %d \n", min_len);
    arr1_start = 0;
    arr2_end = arr2_len;
//    for (i=0;i<min_len;i++)
//    {
//        arr1_end = arr1_len - (shift_len - i);
//        arr2_start = arr2_end - (arr1_end - arr1_start);
//        printf("arr1_start: %d \n", arr1_start);
//        printf("arr1_end: %d \n", arr1_end);
//        printf("arr2_start: %d \n", arr2_start);
//        printf("arr2_end: %d \n", arr2_end);
//    }
    return 0;
}

static PyMethodDef XcorrMethods[] = {
    {"xcorr",  xcorr, METH_VARARGS,
     "Compute the sliding window cross correlation"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initsliding_xcorrlib(void)
{
    (void) Py_InitModule("sliding_xcorrlib", XcorrMethods);
}