// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.
//
//
// Module svpd -- Inverse of Scene View Analysis and Partitioning.

#ifndef svpd_hpp
#define svpd_hpp

#include <stdio.h>

// Reads:
//  VIND -- The decoded LFgamma structure.
// Writes:
//   the images in format "003_003.ppm" into the directory "input_directory".
void svpd(const char* input_directory);

#endif /* svpd_hpp */
