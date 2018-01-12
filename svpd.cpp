// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.

#include <cstdlib>
#include <string>
#include "svpd.h"
#include "gen_types.hh"
#include <iostream>


void svpd(const char* input_directory) {
    
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    int ncol = 3;
    
    char* imagename = (char*)malloc(100*sizeof(char));
    
    int nrr = 0;
    int ncc = 0;
    FILE* f_RV = fopen("VIND", "rb");
    fread(&nviews, sizeof(int), 1, f_RV);
    fread(&nr, sizeof(int), 1, f_RV);
    fread(&nc, sizeof(int), 1, f_RV);
    fread(&ncol, sizeof(int), 1, f_RV);
    int* LFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    fread(LFgamma, sizeof(int), nviews*nr*nc*ncol, f_RV);
    fclose(f_RV);

    int* LI = (int*)malloc(nvr*nr*nvc*nc*ncol*sizeof(int));
    aux_Inverse_LFgamma_TO_LI(LI, LFgamma, LI);

    for( int i = 0; i < 13; ++i ) {
        for( int j = 0; j < 13; ++j ) {
            sprintf(imagename, "%03d_%03d.ppm", j+1, i+1);
            std::string filename = std::string(input_directory) + "/" + std::string(imagename);
            unsigned short int* img = (unsigned short*) malloc( nr*nc*ncol*sizeof(unsigned short));
            FILE* filept = fopen(filename.c_str(), "wb");
            
            for(int icomp = 0; icomp < ncol; ++icomp ) {
                for(int ir = 0; ir < nr; ++ir ) {
                    for( int jc = 0; jc < nc; ++jc) {
                        int iLI = i+ir*nvr;
                        int jLI = j+jc*nvc;
                        img[ir+jc*nr+icomp*nr*nc] = static_cast<unsigned short int>(LI[iLI+jLI*nr*nvr+icomp*nr*nvr*nc*nvc]);
                    }
                }
            }
            
            
            aux_write16ppm(filept, nc, nr, img);

			std::cout << filename << std::endl;

            fclose(filept);
            free(img);
            
            
            
        }
    }
    
    free(LI);

    
}
