// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Scalable Reference View Encoder.

#include <cstdio>
#include <cstdlib>
#include "srve.h"
#include "jp2_wrapper.h"
#include "gen_types.hh"

void srve(double bit_budget, double bitrate_factor) {
    
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    int nrr = nr*nvr;
    int ncc = nc*nvc;
    
    int nbits = 10;

    // Read the reference views from disk
    FILE* f_RV = fopen("VIN", "r");
    fread(&nviews, sizeof(int), 1, f_RV);
    fread(&nr, sizeof(int), 1, f_RV);
    fread(&nc, sizeof(int), 1, f_RV);
    int ncol = 0;
    fread(&ncol, sizeof(int), 1, f_RV);
    int* LFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    fread(LFgamma, sizeof(int), nrr*ncc*ncol, f_RV);
    fclose(f_RV);
    
    int* LI = (int*)malloc(nrr*ncc*ncol*sizeof(int));
    
    aux_Inverse_LFgamma_TO_LI(LI, LFgamma, LI);
    
    
    free(LFgamma);
    
    // Replace the corners with another view. The prediction step will bring them close to original.
    for(int i = 0; i < nr; ++i ) {
        for(int j = 0; j < nc; ++j ) {
            for(int icomp = 0; icomp < ncol; ++icomp ) {
                int ind_i = 0+i*nvr;
                int ind_j = 0+j*nvc;
                int ind_ir = 1+i*nvr;
                int ind_jr = 0+j*nvc;
                LI[ind_i+ind_j*nvr*nr+icomp*nvr*nr*nvc*nc] = LI[ind_ir+ind_jr*nvr*nr+icomp*nvr*nr*nvc*nc];
                ind_i = 0+i*nvr;
                ind_j = 12+j*nvc;
                ind_ir = 0+i*nvr;
                ind_jr = 11+j*nvc;
                LI[ind_i+ind_j*nvr*nr+icomp*nvr*nr*nvc*nc] = LI[ind_ir+ind_jr*nvr*nr+icomp*nvr*nr*nvc*nc];
                ind_i = 12+i*nvr;
                ind_j = 0+j*nvc;
                ind_ir = 11+i*nvr;
                ind_jr = 0+j*nvc;
                LI[ind_i+ind_j*nvr*nr+icomp*nvr*nr*nvc*nc] = LI[ind_ir+ind_jr*nvr*nr+icomp*nvr*nr*nvc*nc];
                ind_i = 12+i*nvr;
                ind_j = 12+j*nvc;
                ind_ir = 11+i*nvr;
                ind_jr = 12+j*nvc;
                LI[ind_i+ind_j*nvr*nr+icomp*nvr*nr*nvc*nc] = LI[ind_ir+ind_jr*nvr*nr+icomp*nvr*nr*nvc*nc];
            }
        }
    }
    
    //double bitrate_factor = .2; // How big a portion is allocated to the "baseline" image. The rest is for references.
    double bits_to_background_image = bitrate_factor*bit_budget;
    double bits_to_references = (1-bitrate_factor)*bit_budget;
    // Compress the reference views
    jp2_wrapper jp2;
    if( ncol == 3 && nbits == 10) {
        int* LI16 = (int*)malloc(nrr*ncc*ncol*sizeof(int));
        for(int i = 0; i < nrr*ncc*ncol; ++i ) {
            LI16[i] = LI[i]*64;
        }
        jp2.compress_with_jp2(LI16, nrr, ncc, "f_hatLI.jp2", bits_to_background_image/(double)nrr/(double)ncc);
        free(LI16);
    }
    
    
    // Take a set of reference views separately. This information is known as "SELECTOR MATRIX" in the diagram.
    int nrefi = 5;
    int nrefj = 5;
    int indrefi[] = {0, 3, 6, 9, 12};
    int indrefj[] = {0, 3, 6, 9, 12};
    
    int* LIs = (int*)malloc(nrefi*nr*nrefj*nc*ncol*sizeof(int));
    
    for(int iref = 0; iref < nrefi; ++iref) {
        for(int jref = 0; jref < nrefj; ++jref) {
            for(int i = 0; i < nr; ++i ) {
                for(int j = 0; j < nc; ++j ) {
                    for(int icomp = 0; icomp < ncol; ++icomp ) {
                        int small_ind_i = iref+i*nrefi;
                        int small_ind_j = jref+j*nrefj;
                        
                        int big_ind_i = indrefi[iref]+i*nvr;
                        int big_ind_j = indrefj[jref]+j*nvc;
                        LIs[small_ind_i+small_ind_j*nrefi*nr+icomp*nrefi*nr*nrefj*nc] = LI[big_ind_i+big_ind_j*nvr*nr+icomp*nvr*nr*nvc*nc];
                    }
                }
            }
        }
    }
    
    if( ncol == 3 && nbits == 10) {
        int* LI16 = (int*)malloc(nrefi*nr*nrefj*nc*ncol*sizeof(int));
        for(int i = 0; i < nrefi*nr*nrefj*nc*ncol; ++i ) {
            LI16[i] = LIs[i]*64;
        }
        jp2.compress_with_jp2(LI16, nrefi*nr, nrefj*nc, "f_hatLIs.jp2", bits_to_references/(double)(nrefi*nr)/(double)(nrefj*nc));
        free(LI16);
    }
}
