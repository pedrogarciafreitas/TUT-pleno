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

#include <cstdio>
#include <cstdlib>
#include <string>
#include "svap.h"
#include "gen_types.hh"


void svap(const char* input_directory) {
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    int ncol = 3;
    
    char* imagename = (char*)malloc(100*sizeof(char));
    
    // create big lenslet image
    int* LI = (int*)malloc(nvr*nr*nvc*nc*ncol*sizeof(int));
    
	#pragma omp parallel for
    for( int i = 0; i < 13; ++i ) {
        for( int j = 0; j < 13; ++j ) {
            sprintf(imagename, "%03d_%03d.ppm", j+1, i+1);
            std::string filename = std::string(input_directory) + "/" + std::string(imagename);
            int* img = (int*) malloc( nr*nc*ncol*sizeof(int));
            FILE* filept = fopen(filename.c_str(), "rb");
            aux_read16ppm(filept, nc, nr, img);
            
            for(int icomp = 0; icomp < ncol; ++icomp ) {
                for(int ir = 0; ir < nr; ++ir ) {
                    for( int jc = 0; jc < nc; ++jc) {
                        int iLI = i+ir*nvr;
                        int jLI = j+jc*nvc;
                        LI[iLI+jLI*nr*nvr+icomp*nr*nvr*nc*nvc] = img[ir+jc*nr+icomp*nr*nc];
                    }
                }
            }
            
            fclose(filept);
            free(img);
            
        }
    }
    
    // go to 4D LFgamma structure
    int* LFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    aux_Construct_LFgamma(LI, 0, LFgamma);
    
    // write to disk
    FILE* f_RV = fopen("VIN", "wb");
    fwrite( &nviews ,sizeof(int),1,f_RV);
    fwrite( &nr ,sizeof(int),1,f_RV);
    fwrite( &nc ,sizeof(int),1,f_RV);
    fwrite( &ncol, sizeof(int), 1, f_RV);
    fwrite( LFgamma ,sizeof(int),nviews*nr*nc*ncol,f_RV);
    fclose(f_RV);
    
    free(LI);
    free(LFgamma);
}
