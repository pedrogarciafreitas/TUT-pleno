// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
//  Module HEVC Reference View Decoder.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "hevcrvd.h"
#include "gen_types.hh"

#include <sstream>
#include <iostream>
#include <fstream>

void svpd_hevc_tmp( void ) {

	const char* input_directory = "hevcrve_decoded";

	int nr = 0;
	int nc = 0;
	int nvr = 0;
	int nvc = 0;
	aux_read_header(&nr, &nc, &nvr, &nvc);
	int nviews = nvr*nvc;
	int ncol = 3;

	char* imagename = (char*)malloc(100 * sizeof(char));

	int nrr = 0;
	int ncc = 0;
	FILE* f_RV = fopen("DRV", "rb");
	fread(&nviews, sizeof(int), 1, f_RV);
	fread(&nr, sizeof(int), 1, f_RV);
	fread(&nc, sizeof(int), 1, f_RV);
	fread(&ncol, sizeof(int), 1, f_RV);
	int* LFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
	fread(LFgamma, sizeof(int), nviews*nr*nc*ncol, f_RV);
	fclose(f_RV);

	int* LI = (int*)malloc(nvr*nr*nvc*nc*ncol*sizeof(int));
	aux_Inverse_LFgamma_TO_LI(LI, LFgamma, LI);

	for (int i = 0; i < 13; ++i) {
		for (int j = 0; j < 13; ++j) {
			sprintf(imagename, "%03d_%03d.ppm", j + 1, i + 1);
			std::string filename = std::string(input_directory) + "/" + std::string(imagename);
			unsigned short int* img = (unsigned short*)malloc(nr*nc*ncol*sizeof(unsigned short));
			FILE* filept = fopen(filename.c_str(), "wb");
			for (int icomp = 0; icomp < ncol; ++icomp) {
				for (int ir = 0; ir < nr; ++ir) {
					for (int jc = 0; jc < nc; ++jc) {
						int iLI = i + ir*nvr;
						int jLI = j + jc*nvc;
						img[ir + jc*nr + icomp*nr*nc] = static_cast<unsigned short int>(LI[iLI + jLI*nr*nvr + icomp*nr*nvr*nc*nvc]);
					}
				}
			}


			aux_write16ppm(filept, nc, nr, img);

			fclose(filept);
			free(img);



		}
	}

	free(LI);


}

void hevcrvd(const char* ffmpeg_path, const char* x265_decoder_path) {
    
    // hevc decoding
    
    int status = 0;

	std::stringstream x265_decoder_b, x265_decoder_ref, x265_decoder_nonref;

	x265_decoder_b << x265_decoder_path;

	x265_decoder_ref << x265_decoder_b.str();
	x265_decoder_ref << " -b ref-views.x265 -o ref-views-dec.yuv";

	x265_decoder_nonref << x265_decoder_b.str();
	x265_decoder_nonref << " -b nonref-views.x265 -o nonref-views-dec.yuv";
    
	status = system( x265_decoder_ref.str().c_str() );
	status = system( x265_decoder_nonref.str().c_str() );

    // yuv to ppm

	std::stringstream ffmpeg_b, ffmpeg_nonref, ffmpeg_ref;

	ffmpeg_b << ffmpeg_path;

	ffmpeg_nonref << ffmpeg_b.str();
	ffmpeg_nonref << " -y -f rawvideo -vcodec rawvideo -s 626x434 -r 30 -pix_fmt yuv422p10le -i ref-views-dec.yuv %03d-refview.ppm";

	ffmpeg_ref << ffmpeg_b.str();
	ffmpeg_ref << " -y -f rawvideo -vcodec rawvideo -s 626x434 -r 30 -pix_fmt yuv422p10le -i nonref-views-dec.yuv %03d-nonrefview.ppm";

    
	status = system(ffmpeg_nonref.str().c_str());
	status = system(ffmpeg_ref.str().c_str());

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

//	/* Create indexing matrix for inserting nearest neighbour from reference views */
//	int* hmat = (int*)malloc(13 * 13 * sizeof(int));
//    for (int i = 0; i < 13; ++i){
//        std::cout << std::endl;
//        for (int j = 0; j < 13; ++j){
//            hmat[i+j*13] = 1;
//        }
//    }
//	int inn = 0, nsum = 0;
//	for (int i = 0; i < 13; ++i){
//		for (int jj = 0; jj < 13; ++jj){
//
//			int j = 0;
//			if (i % 2){
//				j = ordev_rev[jj];
//			}
//			else{
//				j = ordev[jj];
//			}
//
//            if ((!(i % 2) && !(j % 2)) || ((i % 2) && (j % 2))) { // SELECTOR MATRIX
//				for (int rr = -1; rr < 1; ++rr){
//					for (int cc = 0; cc < 1; ++cc){
//						int ir, ic;
//						ir = i + rr;
//						ic = j + cc;
//						if (ir >= 0 && ir < 13 && ic >= 0 && ic < 13){
//							hmat[ir + ic*13] = inn;
//							++nsum;
//						}
//					}
//				}
//				++inn;
//			}
//		}
//	}

	/*Printing matrix for debug*/
//	std::cout << std::endl << nsum << "\t" << inn << std::endl;   
	for (int i = 0; i < 13; ++i){
		std::cout << std::endl;
		for (int j = 0; j < 13; ++j){
			std::cout << hmat[i + j * 13] << "\t";
		}
	}
	std::cout << std::endl;
	
    // create big lenslet image
    int* LI_hevc = (int*)malloc(nvr*nr*nvc*nc*ncol*sizeof(int));
    
	std::cout << "reading files \n";

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

                sprintf(imagename, "%03d-refview.ppm", n_refs+1);
                n_refs++;
            }
            else{
				int tmpi = hmat[i + j * 13]; /*When we want to use only the reference views in nearest neighbour array*/
				sprintf(imagename, "%03d-refview.ppm", tmpi); // should not be tmp+1
                using_same_ref = 1;
//				std::cout << imagename << std::endl;
//                sprintf(imagename, "%03d-nonrefview.ppm", n_nonrefs+1);
//                n_nonrefs++;
            }

            std::string filename = std::string(imagename);
            std::cout << filename << std::endl;
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
    
	std::cout << "done reading files \n";

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
	//svpd_hevc_tmp(); 
    
    free(LI_hevc);
    free(hatLFgamma);
    
}

























