// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Predictor Transform Coefficients Encoder.

#ifndef ptce_h

// Writes:
//  bs_Pred_ThetaW.txt
// Reads:
//  PRCO -- Predictor coefficients
void ptce_encode_predictor_coefficients();

// Writes:
//  bs_Pred_RegWPacked.txt
// Reads:
//  PRMA -- Predictor masks
void ptce_encode_predictor_mask();


#define ptce_h


#endif /* ptce_h */
