// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module vpts -- View Prediction/Transform/Synthesis

#ifndef vpts_hpp
#define vpts_hpp

#include <stdio.h>

// Finds predictors that predict regions from reference views.
// Reads:
//  - VIN -- the true views as LFgamma[iview,ir,ic,icomp]
//  - DRV -- decoded reference views as hatLFgamma.
// Writes:
//  - PRMA -- Predictor masks.
//  - PRCO -- Predictor coefficients.
void vpts(int Ms);

#endif /* vpts_hpp */
