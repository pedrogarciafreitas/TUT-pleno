// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17 
//
// Module vptd -- View Prediction/Transform/Synthesis

#ifndef vptd_hpp
#define vptd_hpp

#include <stdio.h>

// Predicts the intermediate views.
// Reads:
//  - VIN -- the original input data as LFgamma.
//  - DRV -- the decoded reference views as hatLFgamma.
//
//  Writes:
//  - PRCO -- predictor coefficients
//  - PRMA -- predictor masks
void vptd();

#endif /* vptd_hpp */

