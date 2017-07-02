/*------------------------------------------------------------------------------
# Filename: corr_par_cv.c
#  Purpose: Compute cross-correlation of templates with data in parallel
#   Author: Calum Chamberlain
#     Date: 17/2/17
#-----------------------------------------------------------------------------*/
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <stdio.h>

using namespace cv;
using namespace std;

extern "C" int Corr_cv(float* tr1, float* tr2, float* ccc, int l1, int l2, int clen);
  int Corr_cv(float* tr1, float* tr2, float* ccc, int l1, int l2, int clen){
    Mat tr1_mat(1, l1, CV_32F, tr1);
    Mat tr2_mat(1, l2, CV_32F, tr2);
    Mat ccc_mat(1, clen, CV_32F, ccc);

    // Ideally this would accept the FFT of tr2 and compute the FFT
    // of tr1 and convolve them to compute the cross-correlation
    // Currently tr2 is FFT'd a lot!
    matchTemplate(tr1_mat, tr2_mat, ccc_mat, CV_TM_CCOEFF_NORMED);
    // cout << tr1_mat;
    // printf("\n");
    for(int a=0; a < clen; a++){
      // printf("Moving index %i\n", a);
      ccc[a] = ccc_mat.at<float>(0, a);
    }
    return 0;
  }


extern "C" int multiCorr_cv(float* templates, float* tr2, float* cccs, int l1, int l2, int clen, int ntemp);
  int multiCorr_cv(float* templates, float* tr2, float* cccs, int l1, int l2, int clen, int ntemp){
    #pragma omp parallel for
    for(int i=0; i < ntemp; i++){
      Corr_cv(&templates[l1 * i], tr2, &cccs[clen * i], l1, l2, clen);
    }
    return 0;
  }
