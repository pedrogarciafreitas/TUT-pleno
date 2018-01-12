// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module disd -- Disparity decoder.

#ifndef disd_hpp
#define disd_hpp

#include <stdio.h>

// Writes:
// SEGM -- SEGgamma that is a 3D variable with sizes nviews, nr, nc.
//         the segmentation propagated to all views.
// Reads:
// DISP -- The displacements.
// SEGC -- The center view segmentation.
void disd();

// Reads:
//  bs_Pred_Reg.txt -- the bitrstream for displacements
// Writes:
//  DISP -- The displacements needed to warp the center view to different views.
void disd_decode_displacements();

// Reads:
//    "bs_segm.bin" --  the bitstream for center view segmentation.
// Writes:
//    SEGC -- 2D nr times nc center view segmentation
void disd_decode_depth_map();

#endif /* disd_hpp */
