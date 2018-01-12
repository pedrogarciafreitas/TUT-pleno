// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module dise -- Disparity encoder.

#ifndef dise_hpp
#define dise_hpp

#include <stdio.h>

// Writes:
// SEGM -- SEGgamma that is a 3D variable with sizes nviews, nr, nc.
//         the segmentation propagated to all views.
// Reads:
// DISP -- The displacements.
// SEGC -- The center view segmentation.
int dise(const char* filename_depth);

// Writes:
//  bs_segm.bin -- the bitstream for center view.
// Reads:
//  SEGC -- 2D nr times nc center view segmentation
void dise_encode_depth_map(const char* filename_depth);

// Writes:
//  bs_Pred_Reg.txt -- the bitrstream for displacements
// Reads:
//  DISP -- The displacements needed to warp the center view to different views.
void dise_encode_displacements();


#endif /* dise_hpp */
