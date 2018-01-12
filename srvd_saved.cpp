// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Saved Reference View Decoder.

#include <iostream>
#include "srvd_saved.h"
#include "gen_types.hh"

void srvd_saved(const char* filepath) {
    
    // hevc decoding
    
    int status = 0;
    

    
    // ppm to LF
    
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    int ncol = 3;
    
    FILE* f_SM = fopen("selector_matrix.txt", "rb");
    int* SM = (int*)malloc(13*13*sizeof(int));
    fread(SM, sizeof(int), 13*13, f_SM);
    fclose(f_SM);
    
    FILE* f_hmat = fopen("hmat.txt", "rb");
    int* hmat = (int*)malloc(13*13*sizeof(int));
    fread(hmat, sizeof(int), 13*13, f_hmat);
    fclose(f_hmat);
    
    int n_refs = 0;
    int n_nonrefs = 0;
    
    int ordev[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    int ordev_rev[13] = { 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
    
    char* imagename = (char*)malloc(100*sizeof(char));

    
    /*Printing matrix for debug*/
    //	std::cout << std::endl << nsum << "\t" << inn << std::endl;
    for (int i = 0; i < 13; ++i){
        std::cout << std::endl;
        for (int j = 0; j < 13; ++j){
            std::cout << hmat[i + j * 13] << "\t";
        }
    }
    std::cout << std::endl;
    
    int i_ref = 0;
    int i_for_i_ref[225];
    int j_for_i_ref[225];
    
    for( int i = 0; i < 13; ++i ) {
        for (int jj = 0; jj < 13; ++jj) {
            
            int j = 0;
            if (i % 2){
                j = ordev_rev[jj];
            }
            else{
                j = ordev[jj];
            }
            
            if( SM[i+j*13] == 1 ) {
                i_for_i_ref[i_ref] = i;
                j_for_i_ref[i_ref] = j;
                i_ref++;
            }
        }
    }

    
    // create big lenslet image
    int* LI_hevc = (int*)malloc(nvr*nr*nvc*nc*ncol*sizeof(int));
    
    for( int i = 0; i < 13; ++i ) {
        for (int jj = 0; jj < 13; ++jj) {
            int using_same_ref = 0;
            
            int j = 0;
            if (i % 2){
                j = ordev_rev[jj];
            }
            else{
                j = ordev[jj];
            }
            //            if ((!(i % 2) && !(j % 2)) || ((i % 2) && (j % 2))) { // SELECTOR MATRIX
            if( SM[i+j*13] == 1 ) {
                
                int tmpi = hmat[i + j * 13]; /*When we want to use only the reference views in nearest neighbour array*/
                int this_i = i_for_i_ref[tmpi-1];
                int this_j = j_for_i_ref[tmpi-1];
                sprintf(imagename, "%s/%03d_%03d.ppm", filepath, this_j+1, this_i+1);
                n_refs++;
            }
            else{
                int tmpi = hmat[i + j * 13]; /*When we want to use only the reference views in nearest neighbour array*/
                int this_i = i_for_i_ref[tmpi-1];
                int this_j = j_for_i_ref[tmpi-1];
                sprintf(imagename, "%s/%03d_%03d.ppm", filepath, this_j+1, this_i+1);
                using_same_ref = 1;
            }
            
            std::string filename = std::string(imagename);
            std::cout << filename.c_str() << std::endl;
            int* img = (int*) malloc( nr*nc*ncol*sizeof(int));
            FILE* filept = fopen(filename.c_str(), "rb");
            aux_read16ppm(filept, nc, nr, img);
            
            for(int icomp = 0; icomp < ncol; ++icomp ) {
                for(int ir = 0; ir < nr; ++ir ) {
                    for( int jc = 0; jc < nc; ++jc) {
                        int iLI = i+ir*nvr;
                        int jLI = j+jc*nvc;
                        LI_hevc[iLI + jLI*nr*nvr + icomp*nr*nvr*nc*nvc] = img[ir + jc*nr + icomp*nr*nc];
                        if( using_same_ref ) {
                            LI_hevc[iLI + jLI*nr*nvr + icomp*nr*nvr*nc*nvc] += rand() % 2;
                        }
                    }
                }
            }
            
            fclose(filept);
            free(img);
            
        }
    }
    
    // go to 4D LFgamma structure
    int* hatLFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    aux_Construct_LFgamma(LI_hevc, 0, hatLFgamma);
    
    // write to disk
    FILE* f_RV = fopen("DRV", "wb");
    fwrite( &nviews ,sizeof(int),1,f_RV);
    fwrite( &nr ,sizeof(int),1,f_RV);
    fwrite( &nc ,sizeof(int),1,f_RV);
    fwrite( &ncol, sizeof(int), 1, f_RV);
    fwrite(hatLFgamma, sizeof(int), nviews*nr*nc*ncol, f_RV);
    fclose(f_RV);
    
    /* Testing purposes. 
     Save the intermediate hevc encoded views for comparison with 
     actual encoding with disparity compensated sparse prediction. */
//    svpd_hevc_tmp(); 
    
    free(LI_hevc);
    free(hatLFgamma);
    
}



