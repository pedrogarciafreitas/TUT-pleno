// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Scene/View analysis and Partitioning module.


#ifndef svap_hpp
#define svap_hpp

#include <stdio.h>

// Reads:
//  the input files (e.g. 003_002.ppm) from directory "input_directory"
// Writes:
//  VIN -- 4D LFgamma with sizes [nviews, nr, nc, ncomp]
void svap(const char* input_directory);

#endif /* svap_hpp */
