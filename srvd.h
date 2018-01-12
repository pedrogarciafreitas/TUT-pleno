// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Scalable reference view decoder.

#ifndef srvd_hpp
#define srvd_hpp


// Reads:
//  f_hatLI.jp2 and f_hatLIs.jp2
// Writes:
//  DRV -- Decoded reference views as 4D hatLFgamma
void srvd(char* infile, char* infile2);

#endif /* srvd_hpp */
