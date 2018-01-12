// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
//  Module HEVC Reference View Encoder.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "hevcrve.h"
#include "gen_types.hh"

#include <iostream>
#include <fstream>
#include <sstream>

void hevcrve(const char* input_directory, double bit_budget,
	const char* ffmpeg_path, const char* x265_encoder_path, double bitrate_factor) {

	int nr = 0;
	int nc = 0;
	int nvr = 0;
	int nvc = 0;
	aux_read_header(&nr, &nc, &nvr, &nvc);
	int nviews = nvr*nvc;
	int ncol = 3;

	char* imagename = (char*)malloc(100 * sizeof(char));

	std::ofstream ref_list_file;
	ref_list_file.open("ref_list.txt");

	std::ofstream nonref_list_file;
	nonref_list_file.open("nonref_list.txt");
    
    FILE* f_SM = fopen("selector_matrix.txt", "rb");
    int* SM = (int*)malloc(13*13*sizeof(int));
    fread(SM, sizeof(int), 13*13, f_SM);
    fclose(f_SM);

	double n_nonref = 0;
	double n_ref = 0;

	int ordev[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	int ordev_rev[13] = { 12,11,10,9,8,7,6,5,4,3,2,1,0 };

	for (int i = 0; i < 13; ++i) {
		for (int jj = 0; jj < 13; ++jj) {
            // for snake
			int j = 0;
			if (i % 2){
				j = ordev_rev[jj];
			}
			else{
				j = ordev[jj];
			}
			sprintf(imagename, "%03d_%03d.ppm", j + 1, i + 1);
			std::string filename = std::string(input_directory) + "/" + std::string(imagename);
			// here writing lists
//			if (!(i % 3) && !(j % 3)) { // SELECTOR MATRIX
            
//                if ((!(i % 2) && !(j % 2)) || ((i % 2) && (j % 2))) { // SELECTOR MATRIX
            if( SM[i+j*13] == 1 ) {
                    
				ref_list_file << "file """ << filename << """\n" << "duration 1\n";
				n_ref++;
			}
			else {
				nonref_list_file << "file """ << filename << """\n" << "duration 1\n";
				n_nonref++;
			}
		}
	}

	ref_list_file.close();
	nonref_list_file.close();

	// now run ffmpeg

	std::stringstream ffmpeg_b, ffmpeg_nonref, ffmpeg_ref;

	ffmpeg_b << ffmpeg_path;

	ffmpeg_ref << ffmpeg_b.str();
	ffmpeg_ref << " -y -r 30 -f concat -safe 0 -i ref_list.txt -s 626x434 -framerate 30 -c:v rawvideo -pix_fmt yuv422p10le ref-views.yuv";

	ffmpeg_nonref << ffmpeg_b.str();
	ffmpeg_nonref << " -y -r 30 -f concat -safe 0 -i nonref_list.txt -s 626x434 -framerate 30 -c:v rawvideo -pix_fmt yuv422p10le nonref-views.yuv";

	int status = 0;

	status = system(ffmpeg_ref.str().c_str());
	status = system(ffmpeg_nonref.str().c_str());

	// hevc encoding

//	double bitrate_factor = .3; // How big a portion is allocated to the "baseline" image. The rest is for references.
	double bits_to_nonrefences = bitrate_factor*bit_budget;
	double bits_to_references = (1 - bitrate_factor)*bit_budget;
    std::cout << "bits: " << bit_budget << " " << bitrate_factor << " " << bits_to_references << std::endl;
    double rate_nonref = floor( (bits_to_nonrefences/1000) / (n_nonref / 30) ); // kilobits / seconds
    double rate_ref = floor( (bits_to_references/1000) / (n_ref / 30) );

	std::stringstream x265_str_b, x265_ref, x265_nonref;


	x265_str_b << x265_encoder_path << " --input-depth 10 --input-csp i422 --fps 30 --input-res 626x434 --output-depth 10 --profile main422-10";

	if(1){ // heavy compression, slow ~0.5fps
		/*x265_str_b << " --rd 5 --weightb";
		x265_str_b << " --rd 5 --me full --max-merge 5 --subme 7 --weightb";*/
		//x265_str_b << " --rd 6 --ref 8 --allow-non-conformance --me full --max-merge 5 --subme 7 --weightb";
        x265_str_b << " --preset 8 --tune psnr";
	}

	std::cout << x265_str_b.str();

	x265_ref << x265_str_b.str();
	x265_ref << " --stats stats-ref.log --input ref-views.yuv --output ref-views.x265 --bitrate " << (int)rate_ref;

	x265_nonref << x265_str_b.str();
	x265_nonref << " --stats stats-nonref.log --input nonref-views.yuv --output nonref-views.x265 --bitrate " << (int)rate_nonref;

    std::cout << "rate: " << rate_ref << std::endl;
    std::cout << x265_ref.str() << std::endl;
    
	//std::cout << x265_nonref.str();

	for (int ijj = 0; ijj < 3; ijj++)
	{
		std::stringstream x265_tmp_ref;
		x265_tmp_ref << x265_ref.str();
		x265_tmp_ref << " --pass " << ijj + 1;

		status = system( x265_tmp_ref.str().c_str() );

		std::stringstream x265_tmp_nonref;
		x265_tmp_nonref << x265_nonref.str();
		x265_tmp_nonref << " --pass " << ijj + 1;

		status = system(x265_tmp_nonref.str().c_str());
	}

}

























