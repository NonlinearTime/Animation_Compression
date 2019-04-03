//
// Created by haines on 3/31/19.
//

#ifndef ANIMATION_COMPRESSION_EIGEN_H
#define ANIMATION_COMPRESSION_EIGEN_H

#include <valarray>
#include "Matrix.h"
#include <math.h>
#include <iostream>

using namespace std;

int QR(valarray<double>& main_diag, valarray<double>& minor_diag, Matrix &q, double accuracy, int l);
int Householder(Matrix m, Matrix &q, valarray<double>& main_diag, valarray<double>& minor_diag);



#endif //ANIMATION_COMPRESSION_EIGEN_H
