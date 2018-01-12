// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Predictor Transform Coefficients Encoder.

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "ptce.h"
#include "cerv/cerv.h"
#include "golomb_coder.hh"
#include "gen_types.hh"

void ptce_encode_predictor_coefficients() {
    FILE* f_main_header = fopen("HDR", "rb");
    int nr = 0;
    fread(&nr, sizeof(int), 1, f_main_header);
    int nc = 0;
    fread(&nc, sizeof(int), 1, f_main_header);
    int nvr = 0;
    fread(&nvr, sizeof(int), 1, f_main_header);
    int nvc = 0;
    fread(&nvc, sizeof(int), 1, f_main_header);
    fclose(f_main_header);

    
    int nviews = nvr*nvc;
    
    FILE* f_IntPred_ThetaW = fopen("PRCO", "rb");
    fread(&nviews, sizeof(int), 1, f_IntPred_ThetaW);
    int maxiS = 0;
    fread(&maxiS, sizeof(int), 1, f_IntPred_ThetaW);
    int Ms = 0;
    fread(&Ms, sizeof(int), 1, f_IntPred_ThetaW);
    int* IntPred_ThetaW = (int*)malloc(nviews*maxiS*63*sizeof(int));
    fread(IntPred_ThetaW, sizeof(int), nviews*maxiS*63, f_IntPred_ThetaW);
    fclose(f_IntPred_ThetaW);
    
    GolombCoder golomb_coder("bs_Pred_ThetaW.txt", 0);
    int lin_ind = 0;
    for(int i = 0; i < nviews; ++i ) {
        std::vector<int> symbols;
        for(int iR = 0; iR < maxiS; ++iR ) {
            for( int iM = 0; iM < Ms; ++iM ) {
                int val = IntPred_ThetaW[i+nviews*(iR+maxiS*iM)];// IntPred_ThetaW[iview, iR, i];
                if( abs(val) > 10000 ) {
                    val = 0;
                }
                symbols.push_back(val);
                //std::cout << i << " " << iR << " " << iM << " " << val << std::endl;
            }
        }
        golomb_coder.encode_symbols(symbols, 10);
    }

}

void ptce_encode_predictor_mask() {
    int i, j, ik, i1;
    long packed, len_vec, len_rep;
    char *repack;
    FILE *fid = NULL;
    
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    
    FILE* f_Pred_RegW = fopen("PRMA", "rb");
    fread(&nviews, sizeof(int), 1, f_Pred_RegW);
    int maxiS = 0;
    fread(&maxiS, sizeof(int), 1, f_Pred_RegW);
    int Ms = 0;
    fread(&Ms, sizeof(int), 1, f_Pred_RegW);
    int* Pred_RegW = (int*)malloc(nviews*maxiS*63*sizeof(int));
    fread(Pred_RegW, sizeof(int), nviews*maxiS*63, f_Pred_RegW);
    fclose(f_Pred_RegW);


    len_vec = nviews*(long)Ms*maxiS;
    ik = len_vec -(int)(((int)len_vec)/4)*4;
    len_rep = ik + (((int)len_vec)/4)*3;
    repack = (char*)malloc( len_rep*sizeof(char) );
    for (i=0; i<len_rep; i++)
        repack[i] = 0;
     
    
    i=0; j = 0;

    while( (i+3) < len_vec )
    {
        packed = (Pred_RegW[i+3] << 18) | (Pred_RegW[i+2] << 12) |  (Pred_RegW[i+1] << 6) | Pred_RegW[i];
        i = i+4;
        // repack on bytes
        repack[j] = (char)(packed & (long)0xff);
        repack[j+1] = (char)( (packed & (long)0xff00) >> 8 );
        repack[j+2] = (char)( (packed & (long)0xff0000) >> 16 );
        j = j+3;
        
        //printf("i\t %d \t j\t %d \n",i,j);
    }
    
    //printf("a2\n");
    
    // treat the last part, if any
    for( i1 = (i); i1<len_vec; i1++)
    {
        repack[j+(i1-i)] = (char)Pred_RegW[i1];
    }
    //printf("%d %d \n",len_rep-1,j+len_vec-1-i);
    if(len_rep-1 != j+len_vec-1-i)
        printf("Error of counting chars");
    
    fid = fopen ( "bs_Pred_RegWPacked.txt" , "wb" );
    if (fid == NULL)
    { fputs("File error2",stderr);  exit(1); }
    fwrite(repack,sizeof(char),len_rep ,fid);
    fclose(fid);
    
    if( repack != NULL )
        free(repack);
}

