#include <stdlib.h>
#include <stdio.h>
#include "jp2_wrapper.h"
#include <iostream>
#include <math.h>

int main() {

	//const int NR = 434;
	//const int NC = 541;
	const int NR = 5368;
	const int NC = 7728;
	const int im_size = NR*NC*3;

	//double im[im_size];
	uint16_t* im16 = (uint16_t*)malloc(im_size*sizeof(uint16_t));
	FILE* fid = fopen("li.bin", "rb");
	fread(im16, sizeof(uint16_t), im_size, fid);
	fclose(fid);

	//std::cout << im[400] << std::endl;
	for( int i = 0; i < im_size; ++i ) {
		im16[i] /= 64;
		//im16[i] = static_cast<uint16_t>(1023*im[i]);
	}

	std::cout << im16[400] << std::endl;

	uint16_t* im16_dec = (uint16_t*)malloc(im_size*sizeof(uint16_t));

	compress_with_jp2(im16, NR, NC, "compressed.jp2");
	
	decompress_with_jp2(im16_dec, NR, NC, "compressed.jp2");

	double MSE = 0;
	for(int i = 0; i < im_size; ++i ) {
		//std::cout << im16[i] << " " << im16_dec[i] << std::endl;
		MSE += (double)(((int)im16[i]-(int)im16_dec[i])*((int)im16[i]-(int)im16_dec[i]));
	}
	MSE /= im_size;
	double PSNR = 10*log10(1023*1023/MSE);
	std::cout << "PSRN: " << PSNR << std::endl;

	//// lossless check
	//for( int i = 0; i < im_size; ++i ) {
	//	if(im16[i] == im16_dec[i]) {
	//		std::cout << "works" << std::endl;
	//	}
	//	else{
	//		std::cout << "not" << std::endl;
	//	}
	//}


	free(im16);
	free(im16_dec);
}
