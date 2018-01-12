// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Wraps the usage of JPEG 2000.

#include <stdio.h>
#include <stdlib.h>
#include "openjpeg.h"
#include "jp2_wrapper.h"


// JPEG 2000 stuff
/**
sample error debug callback expecting no client object
*/
static void error_callback(const char *msg, void *client_data) {
	(void)client_data;
	//fprintf(stdout, "[ERROR] %s", msg);
}
/**
sample warning debug callback expecting no client object
*/
static void warning_callback(const char *msg, void *client_data) {
	(void)client_data;
	//fprintf(stdout, "[WARNING] %s", msg);
}
/**
sample debug callback expecting no client object
*/
static void info_callback(const char *msg, void *client_data) {
	(void)client_data;
	//fprintf(stdout, "[INFO] %s", msg);
}

jp2_wrapper::jp2_wrapper() {}
jp2_wrapper::~jp2_wrapper() {}

void jp2_wrapper::compress_with_jp2(int* A_ij, const int NR, const int NC, char* filename, float bpp, int precision) {

  // JP2 stuff
  opj_cparameters_t l_param;
  opj_set_default_encoder_parameters(&l_param);
  l_param.tcp_numlayers = 1;
  /* tile definitions parameters */
  /* position of the tile grid aligned with the image */
  l_param.cp_tx0 = 0;
  l_param.cp_ty0 = 0;
  /* tile size, we are using tile based encoding */
  l_param.tile_size_on = OPJ_TRUE;
  l_param.cp_tdx = (NC);
  l_param.cp_tdy = (NR);
  /** progression order to use*/
  /** OPJ_LRCP, OPJ_RLCP, OPJ_RPCL, PCRL, CPRL */
  l_param.prog_order = OPJ_LRCP;

  if( bpp > 98 ) {
	  l_param.cp_fixed_quality = 1;
  }
  else {
	  //l_param.cp_fixed_quality = 1;
	  //l_param.tcp_distoratio[0] = 30; //PSNR
	  OPJ_FLOAT32 l_mct [] =
	  {
		  0.2568,    0.5041,    0.0979,
		  -0.1482,   -0.2910,    0.4392,
		  0.4392,   -0.3678,   -0.0714
	  };
	  OPJ_INT32 l_offsets [] =
	  {
		  4112,       32896,       32896
	  };
	  l_param.cp_disto_alloc=1;
	  float rate = 3*16/bpp;
	  l_param.tcp_rates[0] = rate; // rate jossain yksikossa
	  //l_param.max_cs_size = 24000; 
	  l_param.irreversible = 1; 


	  /** number of resolutions */
	  l_param.numresolution = 7; //wavelet 

	  opj_set_MCT(&l_param,l_mct,l_offsets,3);
  } 



  opj_image_cmptparm_t l_params [3];
  opj_image_cmptparm_t * l_current_param_ptr;
  l_current_param_ptr = l_params;
  for (int i=0;i<3;++i) {
     l_current_param_ptr->dx = 1;
     l_current_param_ptr->dy = 1;

     l_current_param_ptr->h = (OPJ_UINT32)(NR);
     l_current_param_ptr->w = (OPJ_UINT32)(NC);

        l_current_param_ptr->sgnd = 0;
     l_current_param_ptr->prec = (OPJ_UINT32)precision;

     l_current_param_ptr->x0 = 0;
     l_current_param_ptr->y0 = 0;

     ++l_current_param_ptr;
  }

  opj_codec_t * l_codec = opj_create_compress(OPJ_CODEC_JP2);

  opj_image_t * l_image = opj_image_tile_create(3,l_params,OPJ_CLRSPC_SRGB);
  l_image->x0 = 0;
  l_image->y0 = 0;
  l_image->x1 = (OPJ_UINT32)(NC);
  l_image->y1 = (OPJ_UINT32)(NR);
  l_image->color_space = OPJ_CLRSPC_SRGB;

  opj_setup_encoder(l_codec,&l_param,l_image);

  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_FALSE);

  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_start_compress(l_codec,l_image,l_stream);

  int l_data_size;
  OPJ_BYTE* l_data;
 // if(precision==8) {
 //   l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
 //   l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
 //   int i = 0;
 //   for( int icomp = 0; icomp < 3; ++icomp) {
 //     for(int iir = 0; iir < A_ij->size().height; ++iir) {
 //       for(int iic = 0; iic < A_ij->size().width ; ++iic) {
 //         l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
 //         ++i;
 //       }
 //     }
 //   }
 // }
  if(precision==16) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3 *(OPJ_UINT32)2; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    for( int icomp = 0; icomp < 3; ++icomp) {
    for(int iir = 0; iir < NR; ++iir) {
      for(int iic = 0; iic < NC; ++iic) {
        uint16_t val = (uint16_t)A_ij[iir+iic*NR+icomp*(NR*NC)];
        l_data[i+1] = (OPJ_BYTE)(val >> 8);
        l_data[i] = (OPJ_BYTE)(val & 0x00ff);
        i += 2;
      }
    }
    }
  }
  if(precision==8) {
	  l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
	  l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
	  int i = 0;
	  for( int icomp = 0; icomp < 3; ++icomp) {
		  for(int iir = 0; iir < NR; ++iir) {
			  for(int iic = 0; iic < NC; ++iic) {
				  uint8_t val = (uint8_t)A_ij[iir+iic*NR+icomp*(NR*NC)];
				  l_data[i] = (OPJ_BYTE)(val & 0x00ff);
				  i += 1;
			  }
		  }
	  }
  }
  int i=0;
  opj_write_tile(l_codec,i,l_data,l_data_size,l_stream);

  opj_end_compress(l_codec,l_stream);
  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(l_image);
  free(l_data);
}



void jp2_wrapper::decompress_with_jp2(int* A_ij, int* NR, int* NC, char* filename) {
   opj_image_t * l_image;
  opj_dparameters_t l_param;
  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_TRUE);
  /* Set the default decoding parameters */
  opj_set_default_decoder_parameters(&l_param);

  /* */
  //l_param.decod_format = JP2_CFMT; 

  /** you may here add custom decoding parameters */
  /* do not use layer decoding limitations */
  l_param.cp_layer = 0;

  /* do not use resolutions reductions */
  l_param.cp_reduce = 0;
  opj_codec_t* l_codec = opj_create_decompress(OPJ_CODEC_JP2);
  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_setup_decoder(l_codec, &l_param);
  opj_read_header(l_stream, l_codec, &l_image);
  OPJ_UINT32 da_x0=0;
  OPJ_UINT32 da_y0=0;
  OPJ_UINT32 da_x1=l_image->x1;
  OPJ_UINT32 da_y1=l_image->y1;
  opj_set_decode_area(l_codec, l_image, da_x0, da_y0,da_x1, da_y1);
  OPJ_UINT32 l_tile_index, l_data_size, l_nb_comps;
  OPJ_INT32 l_current_tile_x0, l_current_tile_y0, l_current_tile_x1, l_current_tile_y1, l_go_on;
  opj_read_tile_header( l_codec,
        l_stream,
        &l_tile_index,
        &l_data_size,
        &l_current_tile_x0,
        &l_current_tile_y0,
        &l_current_tile_x1,
        &l_current_tile_y1,
        &l_nb_comps,
        &l_go_on);

    *NR = da_y1;
    *NC = da_x1;
  OPJ_BYTE* l_data = (OPJ_BYTE *) malloc(l_data_size);
  opj_decode_tile_data(l_codec,l_tile_index,l_data,l_data_size,l_stream);
        /** now should inspect image to know the reduction factor and then how to behave with data */
  int i = 0;
if( l_data_size/(*NR)/(*NC)/3 == 2 ) {
  //cout << (int)l_data[0] << " ";
  for( int icomp = 0; icomp < 3; ++icomp) {
     for(int iir = 0; iir < (*NR); ++iir) {
     //for(int iir = A->size().height-1; iir > 0; --iir) {
        for(int iic = 0; iic < (*NC); ++iic) {
           //l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
	uint8_t val1 = l_data[i];
	uint8_t val2 = l_data[i+1];
	uint16_t val = uint16_t(val2<<8)+val1;
	A_ij[iir+iic*(*NR)+icomp*(*NR)*(*NC)] = val;
	i+=2;
        }
     }
  }
} else {
  for( int icomp = 0; icomp < 3; ++icomp) {
     for(int iir = 0; iir < (*NR); ++iir) {
        for(int iic = 0; iic < (*NC); ++iic) {
	uint8_t val = l_data[i];
	A_ij[iir+iic*(*NR)+icomp*(*NR)*(*NC)] = val;
	i+=1;
        }
     }
  }
}

  free(l_data);
}

void jp2_wrapper::compress_with_jp2_lossless8(int* A_ij, const int NR, const int NC, char* filename, float bpp) {
  int precision = 8;
  // JP2 stuff
  opj_cparameters_t l_param;
  opj_set_default_encoder_parameters(&l_param);
  l_param.tcp_numlayers = 1;
  l_param.cp_fixed_quality = 1;
  //l_param.tcp_distoratio[0] = 30; //PSNR
  //OPJ_FLOAT32 l_mct [] =
  //{
  //        0.2568,    0.5041,    0.0979,
  //        -0.1482,   -0.2910,    0.4392,
  //        0.4392,   -0.3678,   -0.0714
  //};
  //OPJ_INT32 l_offsets [] =
  //{
  //        4112,       32896,       32896
  //};
  //l_param.cp_disto_alloc=1;
  //float rate = 3*16/bpp;
  //l_param.tcp_rates[0] = rate; // rate jossain yksikossa
  //l_param.max_cs_size = 24000; 
  //l_param.irreversible = 1; 

  /* tile definitions parameters */
  /* position of the tile grid aligned with the image */
  l_param.cp_tx0 = 0;
  l_param.cp_ty0 = 0;
  /* tile size, we are using tile based encoding */
  l_param.tile_size_on = OPJ_TRUE;
  l_param.cp_tdx = (NC);
  l_param.cp_tdy = (NR);

  /** number of resolutions */
  l_param.numresolution = (int)bpp; //wavelet 

  /** progression order to use*/
  /** OPJ_LRCP, OPJ_RLCP, OPJ_RPCL, PCRL, CPRL */
  l_param.prog_order = OPJ_LRCP;

  //opj_set_MCT(&l_param,l_mct,l_offsets,3);



  opj_image_cmptparm_t l_params [1];
  opj_image_cmptparm_t * l_current_param_ptr;
  l_current_param_ptr = l_params;
  for (int i=0;i<1;++i) {
     l_current_param_ptr->dx = 1;
     l_current_param_ptr->dy = 1;

     l_current_param_ptr->h = (OPJ_UINT32)(NR);
     l_current_param_ptr->w = (OPJ_UINT32)(NC);

        l_current_param_ptr->sgnd = 0;
     l_current_param_ptr->prec = (OPJ_UINT32)precision;

     l_current_param_ptr->x0 = 0;
     l_current_param_ptr->y0 = 0;

     ++l_current_param_ptr;
  }

  opj_codec_t * l_codec = opj_create_compress(OPJ_CODEC_JP2);

  opj_image_t * l_image = opj_image_tile_create(1,l_params,OPJ_CLRSPC_GRAY);
  l_image->x0 = 0;
  l_image->y0 = 0;
  l_image->x1 = (OPJ_UINT32)(NC);
  l_image->y1 = (OPJ_UINT32)(NR);
  l_image->color_space = OPJ_CLRSPC_GRAY;

  opj_setup_encoder(l_codec,&l_param,l_image);

  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_FALSE);

  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_start_compress(l_codec,l_image,l_stream);

  int l_data_size;
  OPJ_BYTE* l_data;
 // if(precision==8) {
 //   l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
 //   l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
 //   int i = 0;
 //   for( int icomp = 0; icomp < 3; ++icomp) {
 //     for(int iir = 0; iir < A_ij->size().height; ++iir) {
 //       for(int iic = 0; iic < A_ij->size().width ; ++iic) {
 //         l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
 //         ++i;
 //       }
 //     }
 //   }
 // }
  if(precision==8) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR); 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    //for( int icomp = 0; icomp < 3; ++icomp) {
    for(int iir = 0; iir < NR; ++iir) {
      for(int iic = 0; iic < NC; ++iic) {
        uint8_t val = (uint8_t)A_ij[iir+iic*NR];
        l_data[i] = (OPJ_BYTE)(val) ;
        i++;
      }
    }
    //}
  }
  int i=0;
  opj_write_tile(l_codec,i,l_data,l_data_size,l_stream);

  opj_end_compress(l_codec,l_stream);
  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(l_image);
  free(l_data);
}

void jp2_wrapper::compress_with_jp2_lossless(int* A_ij, const int NR, const int NC, char* filename, float bpp) {
  int precision = 16;
  // JP2 stuff
  opj_cparameters_t l_param;
  opj_set_default_encoder_parameters(&l_param);
  l_param.tcp_numlayers = 1;
  l_param.cp_fixed_quality = 1;
  //l_param.tcp_distoratio[0] = 30; //PSNR
  //OPJ_FLOAT32 l_mct [] =
  //{
  //        0.2568,    0.5041,    0.0979,
  //        -0.1482,   -0.2910,    0.4392,
  //        0.4392,   -0.3678,   -0.0714
  //};
  //OPJ_INT32 l_offsets [] =
  //{
  //        4112,       32896,       32896
  //};
  //l_param.cp_disto_alloc=1;
  //float rate = 3*16/bpp;
  //l_param.tcp_rates[0] = rate; // rate jossain yksikossa
  //l_param.max_cs_size = 24000; 
  //l_param.irreversible = 1; 

  /* tile definitions parameters */
  /* position of the tile grid aligned with the image */
  l_param.cp_tx0 = 0;
  l_param.cp_ty0 = 0;
  /* tile size, we are using tile based encoding */
  l_param.tile_size_on = OPJ_TRUE;
  l_param.cp_tdx = (NC);
  l_param.cp_tdy = (NR);

  /** number of resolutions */
  l_param.numresolution = (int)bpp; //wavelet 

  /** progression order to use*/
  /** OPJ_LRCP, OPJ_RLCP, OPJ_RPCL, PCRL, CPRL */
  l_param.prog_order = OPJ_LRCP;

  //opj_set_MCT(&l_param,l_mct,l_offsets,3);



  opj_image_cmptparm_t l_params [1];
  opj_image_cmptparm_t * l_current_param_ptr;
  l_current_param_ptr = l_params;
  for (int i=0;i<1;++i) {
     l_current_param_ptr->dx = 1;
     l_current_param_ptr->dy = 1;

     l_current_param_ptr->h = (OPJ_UINT32)(NR);
     l_current_param_ptr->w = (OPJ_UINT32)(NC);

        l_current_param_ptr->sgnd = 0;
     l_current_param_ptr->prec = (OPJ_UINT32)precision;

     l_current_param_ptr->x0 = 0;
     l_current_param_ptr->y0 = 0;

     ++l_current_param_ptr;
  }

  opj_codec_t * l_codec = opj_create_compress(OPJ_CODEC_JP2);

  opj_image_t * l_image = opj_image_tile_create(1,l_params,OPJ_CLRSPC_GRAY);
  l_image->x0 = 0;
  l_image->y0 = 0;
  l_image->x1 = (OPJ_UINT32)(NC);
  l_image->y1 = (OPJ_UINT32)(NR);
  l_image->color_space = OPJ_CLRSPC_GRAY;

  opj_setup_encoder(l_codec,&l_param,l_image);

  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_FALSE);

  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_start_compress(l_codec,l_image,l_stream);

  int l_data_size;
  OPJ_BYTE* l_data;
 // if(precision==8) {
 //   l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
 //   l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
 //   int i = 0;
 //   for( int icomp = 0; icomp < 3; ++icomp) {
 //     for(int iir = 0; iir < A_ij->size().height; ++iir) {
 //       for(int iic = 0; iic < A_ij->size().width ; ++iic) {
 //         l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
 //         ++i;
 //       }
 //     }
 //   }
 // }
  if(precision==16) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)*(OPJ_UINT32)2; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    //for( int icomp = 0; icomp < 3; ++icomp) {
    for(int iir = 0; iir < NR; ++iir) {
      for(int iic = 0; iic < NC; ++iic) {
        uint16_t val = (uint16_t)A_ij[iir+iic*NR];
        l_data[i+1] = (OPJ_BYTE)(val >> 8);
        l_data[i] = (OPJ_BYTE)(val & 0x00ff);
        i += 2;
      }
    }
    //}
  }
  int i=0;
  opj_write_tile(l_codec,i,l_data,l_data_size,l_stream);

  opj_end_compress(l_codec,l_stream);
  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(l_image);
  free(l_data);
}



void jp2_wrapper::decompress_with_jp2_lossless(int* A_ij, const int NR, const int NC, char* filename) {
   opj_image_t * l_image;
  opj_dparameters_t l_param;
  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_TRUE);
  /* Set the default decoding parameters */
  opj_set_default_decoder_parameters(&l_param);

  /* */
  //l_param.decod_format = JP2_CFMT; 

  /** you may here add custom decoding parameters */
  /* do not use layer decoding limitations */
  l_param.cp_layer = 0;

  /* do not use resolutions reductions */
  l_param.cp_reduce = 0;
  opj_codec_t* l_codec = opj_create_decompress(OPJ_CODEC_JP2);
  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_setup_decoder(l_codec, &l_param);
  opj_read_header(l_stream, l_codec, &l_image);
  OPJ_UINT32 da_x0=0;
  OPJ_UINT32 da_y0=0;
  OPJ_UINT32 da_x1=l_image->x1;
  OPJ_UINT32 da_y1=l_image->y1;
  opj_set_decode_area(l_codec, l_image, da_x0, da_y0,da_x1, da_y1);
  OPJ_UINT32 l_tile_index, l_data_size, l_nb_comps;
  OPJ_INT32 l_current_tile_x0, l_current_tile_y0, l_current_tile_x1, l_current_tile_y1, l_go_on;
  opj_read_tile_header( l_codec,
        l_stream,
        &l_tile_index,
        &l_data_size,
        &l_current_tile_x0,
        &l_current_tile_y0,
        &l_current_tile_x1,
        &l_current_tile_y1,
        &l_nb_comps,
        &l_go_on);

  OPJ_BYTE* l_data = (OPJ_BYTE *) malloc(l_data_size);
  opj_decode_tile_data(l_codec,l_tile_index,l_data,l_data_size,l_stream);
        /** now should inspect image to know the reduction factor and then how to behave with data */
  int i = 0;
  //cout << (int)l_data[0] << " ";
     for(int iir = 0; iir < NR; ++iir) {
     //for(int iir = A->size().height-1; iir > 0; --iir) {
        for(int iic = 0; iic < NC; ++iic) {
           //l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
	uint8_t val1 = l_data[i];
	uint8_t val2 = l_data[i+1];
	uint16_t val = uint16_t(val2<<8)+val1;
	A_ij[iir+iic*NR] = val; 
	i+=2;
        }
  }

  free(l_data);
}

void jp2_wrapper::decompress_with_jp2_lossless8(int* A_ij, const int NR, const int NC, char* filename) {
   opj_image_t * l_image;
  opj_dparameters_t l_param;
  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename, OPJ_TRUE);
  /* Set the default decoding parameters */
  opj_set_default_decoder_parameters(&l_param);

  /* */
  //l_param.decod_format = JP2_CFMT; 

  /** you may here add custom decoding parameters */
  /* do not use layer decoding limitations */
  l_param.cp_layer = 0;

  /* do not use resolutions reductions */
  l_param.cp_reduce = 0;
  opj_codec_t* l_codec = opj_create_decompress(OPJ_CODEC_JP2);
  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_setup_decoder(l_codec, &l_param);
  opj_read_header(l_stream, l_codec, &l_image);
  OPJ_UINT32 da_x0=0;
  OPJ_UINT32 da_y0=0;
  OPJ_UINT32 da_x1=l_image->x1;
  OPJ_UINT32 da_y1=l_image->y1;
  opj_set_decode_area(l_codec, l_image, da_x0, da_y0,da_x1, da_y1);
  OPJ_UINT32 l_tile_index, l_data_size, l_nb_comps;
  OPJ_INT32 l_current_tile_x0, l_current_tile_y0, l_current_tile_x1, l_current_tile_y1, l_go_on;
  opj_read_tile_header( l_codec,
        l_stream,
        &l_tile_index,
        &l_data_size,
        &l_current_tile_x0,
        &l_current_tile_y0,
        &l_current_tile_x1,
        &l_current_tile_y1,
        &l_nb_comps,
        &l_go_on);

  OPJ_BYTE* l_data = (OPJ_BYTE *) malloc(l_data_size);
  opj_decode_tile_data(l_codec,l_tile_index,l_data,l_data_size,l_stream);
        /** now should inspect image to know the reduction factor and then how to behave with data */
  int i = 0;
  //cout << (int)l_data[0] << " ";
     for(int iir = 0; iir < NR; ++iir) {
     //for(int iir = A->size().height-1; iir > 0; --iir) {
	     for(int iic = 0; iic < NC; ++iic) {
		     //l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
		     uint8_t val1 = l_data[i];
		     A_ij[iir+iic*NR] = val1; 
		     i++;
        }
  }

  free(l_data);
}
