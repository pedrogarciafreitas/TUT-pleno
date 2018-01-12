// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.
//
//
// This is the main file that combines all the modules of JPEG-PLENO encoder and decoder.

#include <iostream>
#include <stdio.h>
#include "gen_types.hh"
#include "svap.h"
#include "hevcrve.h"
#include "hevcrvd.h"
//#include "srve.h"
//#include "srvd.h"
#include "srvd_saved.h"
#include "dise.h"
#include "vpts.h"
#include "ptce.h"
#include "disd.h"
#include "srvd.h"
#include "ptcd.h"
#include "vptd.h"
#include "svpd.h"

#include "omp.h"

int main(int argc, char** argv) {

	int encode = 1;

	const char* filename_depth = argv[4]; // depth map. the indexing of regions should start from 1.
	const char* filepath_orig = argv[5];
	const char* filepath_dec = argv[6];
	const char* filepath_saved_refs = argv[7];
	//    const char filename_depth[] = "1_dm.ppm"; // depth map. the indexing of regions should start from 1.
	//    const char filepath_orig[] = "/Users/helinp/JPEG_PLENO/02_Subaperture_PPM_10_bit_CORRECTED/I01_Bikes";
	//    const char filepath_dec[] = "/Users/helinp/JPEG_PLENO/02_Subaperture_PPM_10_bit_dec/I01_Bikes";

	// define what reference encoder to use
	int use_hevc = 1;
	//    int use_jpeg2000 = 0;
	int use_saved_references = 0;

	// if using hevc, locations for external binaries of x265 encoder and decoder need be defined as well as ffmpeg binary
	const char x265_encoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/x265.exe";
	const char x265_decoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/TAppDecoder.exe";
	const char ffmpeg_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/ffmpeg/bin/ffmpeg.exe";

	omp_set_num_threads(6);

	if (encode) {

		std::cout << "aaaa\n";

		double bitrate = atof(argv[2]);
		int Ms = atof(argv[3]); // prediction order.

		int maxiS = 2;

		// write header that includes the sizes of images
		int nr = 434; // number of rows in one view
		int nc = 626; // number of columns in one view
		int nvr = 13; // number of horizontal views
		int nvc = 13; // number of vertical views
		aux_write_header(nr, nc, nvr, nvc);

		/*Scene view analysis and partitioning.
		 Creates VIN that has a 4D-matrix LFgamma[iview,ir,ic,icomp]*/
		std::cout << "svap\n";
		svap(filepath_orig);
		std::cout << "svap done\n";
		/*Disparity encoder module.
		 Creates SEGgamma similar to LFgamma into file SEGM and displacements into DISP.*/
		maxiS = dise(filename_depth);

		/*Calls cerv for DISP.*/
		dise_encode_displacements();

		/*Calls cerv for the depth map.*/
		dise_encode_depth_map(filename_depth);


		double bppDM = (double)aux_GetFileSize("bs_segm.bin")*8.0 / (double)num_pixels_in_lenslet;
		double bppPred = ((double)nvr*(double)nvc*(double)maxiS)*(19.0*(double)Ms + 3.0) / (double)num_pixels_in_lenslet;
		double bit_budget = (bitrate - bppDM - bppPred - 0.0004)*(double)num_pixels_in_lenslet;

		double bitrate_factor = atof(argv[1]);

		/* Reference view encoder. */
		if (use_hevc) {
			hevcrve(filepath_orig, bit_budget, ffmpeg_path, x265_encoder_path, bitrate_factor);
			hevcrvd(ffmpeg_path, x265_decoder_path); // decode
		}

		/*if( use_jpeg2000 ) {
			srve(bit_budget);
			//Reference view decoder.

			char infile[] = "f_hatLI.jp2";
			char infile2[] = "f_hatLIs.jp2";
			srvd(infile, infile2);
			} */

		if (use_saved_references) {

			srvd_saved(filepath_saved_refs);
		}


		/*View prediction or transform synthesis.
		Creates predictor coefficients and mask
		PRCO nviews*maxiS*63
		PRMA nviews*maxiS*63*/
		vpts(Ms);

		// Encode PRCO
		ptce_encode_predictor_coefficients();
		// Encode PRMA
		ptce_encode_predictor_mask();

		// Decoding.

		vptd();
		svpd(filepath_dec);
	}
	else {

		disd_decode_depth_map();
		disd_decode_displacements();

		disd();

		ptcd_decode_predictor_mask();
		ptcd_decode_predictor_coefficients();

		//char infile[] = "f_hatLI.jp2";
		//char infile2[] = "f_hatLIs.jp2";
		//srvd(infile, infile2);

		vptd();

		svpd(filepath_dec);

	}
	exit(0);
}
