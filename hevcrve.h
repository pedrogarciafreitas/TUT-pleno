// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module HEVC Reference View Encoder.

#ifndef hevcrve_hpp
#define hevcrve_hpp

#include <stdio.h>

void hevcrve(const char* input_directory, double bit_budget,
	const char* ffmpeg_path, const char* x265_encoder_path, double bitrate_factor);

#endif /* hevcrve_hpp */
