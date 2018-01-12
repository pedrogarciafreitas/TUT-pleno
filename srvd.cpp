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

#include <stdio.h>
#include <cstdlib>
#include "srvd.h"
#include "jp2_wrapper.h"
#include "gen_types.hh"

void srvd(char* infile, char* infile2) {
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int nrr = nvr*nr;
    int ncc = nvc*nc;
    int ncol = 3;
    int nbits = 10;
    
    int* LI_JP = (int*)malloc(nrr*ncc*3*sizeof(int));

    jp2_wrapper jp2;
    jp2.decompress_with_jp2(LI_JP, &nrr, &ncc, infile);
    for(int i = 0; i < nrr*ncc*ncol; ++i ) {
        LI_JP[i] /= 64;
    }

    int nrefi = 5;
    int nrefj = 5;
    int indrefi[] = {0, 3, 6, 9, 12};
    int indrefj[] = {0, 3, 6, 9, 12};
    
    int* LIs = (int*)malloc(nrefi*nr*nrefj*nc*ncol*sizeof(int));
    int sizer = 0;
    int sizec = 0;
    jp2.decompress_with_jp2(LIs, &sizer, &sizec, infile2);
    for(int i = 0; i < nrefi*nr*nrefj*nc*ncol; ++i ) {
        LIs[i] /= 64;
    }
    
    for(int iref = 0; iref < nrefi; ++iref) {
        for(int jref = 0; jref < nrefj; ++jref) {
            for(int i = 0; i < nr; ++i ) {
                for(int j = 0; j < nc; ++j ) {
                    for(int icomp = 0; icomp < ncol; ++icomp ) {
                        int small_ind_i = iref+i*nrefi;
                        int small_ind_j = jref+j*nrefj;
                        
                        int big_ind_i = indrefi[iref]+i*nvr;
                        int big_ind_j = indrefj[jref]+j*nvr;
                        LI_JP[big_ind_i+big_ind_j*nvr*nr+icomp*nvr*nr*nvc*nc] = LIs[small_ind_i+small_ind_j*nrefi*nr+icomp*nrefi*nr*nrefj*nc];
                    }
                }
            }
        }
    }
    
    free(LIs);
    
    
    int* hatLFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    aux_Construct_LFgamma(LI_JP, 0, hatLFgamma);
    
    FILE* f_RV = fopen("DRV", "w");
    fwrite( &nviews ,sizeof(int),1,f_RV);
    fwrite( &nr ,sizeof(int),1,f_RV);
    fwrite( &nc ,sizeof(int),1,f_RV);
    fwrite( &ncol, sizeof(int), 1, f_RV);
    fwrite( hatLFgamma ,sizeof(int),nviews*nr*nc*ncol,f_RV);
    fclose(f_RV);
    
    free(LI_JP);

}
