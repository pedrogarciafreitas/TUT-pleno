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

void hevcrve(const char* input_directory, int nbits,
	const char* ffmpeg_path, const char* x265_encoder_path, const char* output_path,
	const char* reflist_path, const char* selectormatrix_path, 
	const char* yuv_path,
	const int nr, const int nc, const int nvr, const int nvc,
	const int width, const int height, const int xx, const int yy) {


	char* imagename = (char*)malloc(100 * sizeof(char));

	std::ofstream ref_list_file;
	ref_list_file.open( reflist_path );

	FILE* f_SM = fopen( selectormatrix_path , "rb");
	int* SM = (int*)malloc( nvr * nvc * sizeof(int));
	fread(SM, sizeof(int), nvr * nvc, f_SM);
	fclose(f_SM);

	int n_ref = 0;

	int curindex = 1;
	bool added_frame = 1;

	while ( added_frame ){

		added_frame = 0;

		for (int ii = 0; ii < nvr; ++ii) {
			for (int jj = 0; jj < nvc; ++jj) {
				if (SM[ii + jj * nvc] == curindex) {

					sprintf(imagename, "%03d_%03d.ppm", jj, ii);
					std::string filename = std::string(input_directory) + "/" + std::string(imagename);

					ref_list_file << "file """ << filename << """\n" << "duration 1\n";
					n_ref++;
					added_frame = 1;

				}
			}
		}
	}

	/*
	//int ordev[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	//int ordev_rev[13] = { 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };

	int *ordev = new int[nvc];
	int *ordev_rev = new int[nvc];

	for (int ij = 0; ij < nvc; ij++)
		ordev[ij] = ij;


	for (int ij = nvc-1; ij > 0; ij--)
		ordev_rev[ij] = ij;

	for (int i = 0; i < nvr; ++i) {
		for (int jj = 0; jj < nvc; ++jj) {
			// for snake we flip left-right every other row
			int j = 0;
			if (i % 2){
				j = ordev_rev[jj];
			}
			else{
				j = ordev[jj];
			}
			sprintf(imagename, "%03d_%03d.ppm", j, i );
			std::string filename = std::string(input_directory) + "/" + std::string(imagename);

			if (SM[i + j * nvc] == 1) {

				ref_list_file << "file """ << filename << """\n" << "duration 1\n";
				n_ref++;

			}
		}
	}
	*/

	ref_list_file.close();

	std::stringstream ffmpeg_b, ffmpeg_ref;

	ffmpeg_b << ffmpeg_path;

	ffmpeg_ref << ffmpeg_b.str();
	ffmpeg_ref << " -y -r 30 -f concat -safe 0 -i " << reflist_path << " -s " << nc << "x" << nr;
	
	if (width > 0)
	{
		ffmpeg_ref << " filter:v ""crop=""" << width << ":" << height << ":" << ":" << xx << ":" << yy << """";
	}

	ffmpeg_ref << " -framerate 30 -c:v rawvideo -pix_fmt yuv422p10le " << yuv_path;

	int status = 0;

	status = system(ffmpeg_ref.str().c_str());

	// hevc encoding

	double bits_to_references = nbits;
	std::cout << "bits: " << bits_to_references << std::endl;
	double rate_ref = floor((bits_to_references / 1000) / ((double)n_ref / 30));

	std::stringstream x265_str_b, x265_ref;

	x265_str_b << x265_encoder_path << " --input-depth 10 --input-csp i422 --fps 30 --input-res " << nc << "x" << nr << " --output-depth 10 --profile main422-10";

	if (1)
	{ // heavy compression, slow ~0.5fps
		/*x265_str_b << " --rd 5 --weightb";
		x265_str_b << " --rd 5 --me full --max-merge 5 --subme 7 --weightb";*/
		//x265_str_b << " --rd 6 --ref 8 --allow-non-conformance --me full --max-merge 5 --subme 7 --weightb";
		x265_str_b << " --preset 8 --tune psnr";
	}

	std::cout << x265_str_b.str();

	x265_ref << x265_str_b.str();
	x265_ref << " --stats stats-ref.log --input " << yuv_path << " --output " << output_path << " --bitrate " << (int)rate_ref;

	std::cout << "rate: " << rate_ref << std::endl;
	std::cout << x265_ref.str() << std::endl;

	int NPASS = 1;

	for (int ijj = 0; ijj < NPASS; ijj++)
	{
		std::stringstream x265_tmp_ref;
		x265_tmp_ref << x265_ref.str();
		x265_tmp_ref << " --pass " << ijj + 1;
		status = system(x265_tmp_ref.str().c_str());
	}

	return;
}

int main(int argc, char** argv) {

	const int nbits = atoi(argv[1]);

	const char* filepath_orig = argv[2];

	const char* x265_encoder_path = argv[3];
	const char* x265_decoder_path = argv[4];
	const char* ffmpeg_path = argv[5];
	const char* output_path = argv[6];

	const char* reflist_path = argv[7];

	const char* selectormatrix_path = argv[8];

	const char* yuv_path = argv[9];

	const int nr = atoi(argv[10]);
	const int nc = atoi(argv[11]);
	const int nvr = atoi(argv[12]);
	const int nvc = atoi(argv[13]);

	const int width = atoi(argv[14]);
	const int height = atoi(argv[15]);
	const int xx = atoi(argv[16]);
	const int yy = atoi(argv[17]);

	hevcrve(
		filepath_orig,
		nbits,
		ffmpeg_path,
		x265_encoder_path,
		output_path,
		reflist_path,
		selectormatrix_path,
		yuv_path,
		nr,
		nc,
		nvr,
		nvc,
		width,
		height,
		xx,
		yy);

	exit(0);
}
