// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Wraps the usage of JPEG 2000.


#ifndef jp2_wrapper_h
#define jp2_wrapper_h
#include <stdint.h>

class jp2_wrapper {
public:
    jp2_wrapper();
    ~jp2_wrapper();
    void compress_with_jp2(int* A_ij, const int nr, const int nc, char* filename, float bpp, int precision=16);
    void compress_with_jp2_lossless(int* A_ij, const int nr, const int nc, char* filename, float bpp);
    void compress_with_jp2_lossless8(int* A_ij, const int nr, const int nc, char* filename, float bpp);
    void decompress_with_jp2(int* A_ij,  int* NR,  int* NC, char* filename);
    void decompress_with_jp2_lossless(int* A_ij, const int NR, const int NC, char* filename);
    void decompress_with_jp2_lossless8(int* A_ij, const int NR, const int NC, char* filename);
};
#endif
