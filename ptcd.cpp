// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module Predictor Transform Coefficients Decoder.

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "ptcd.h"
#include "gen_types.hh"
#include "golomb_coder.hh"

void ptcd_decode_predictor_coefficients() {
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int Ms = 0;
    int maxiS = 0;
    aux_read_pred_header(&Ms, &maxiS);
    
    // Read file Pred_ThetaW1
    int* Pred_ThetaW1 = alocaVector((int)nviews*63*maxiS);
    
    GolombCoder golomb_coder("bs_Pred_ThetaW.txt", 1);
    int lin_ind = 0;
    for(int i = 0; i < nviews; ++i ) {
        std::vector<int> symbols;
        golomb_coder.decode_symbols(symbols, 10);
        lin_ind = 0;
        for(int iR = 0; iR < maxiS; ++iR ) {
            for( int iM = 0; iM < Ms; ++iM ) {
                Pred_ThetaW1[i+nviews*(iR+maxiS*iM)] = symbols[lin_ind];// Pred_RegW[iview, iR, i] = PredRegr[i];
                lin_ind++;
            }
        }
    }
    
    FILE* f_IntPred_ThetaW = fopen("PRCO", "wb");
    fwrite(&nviews, sizeof(int), 1, f_IntPred_ThetaW);
    fwrite(&maxiS, sizeof(int), 1, f_IntPred_ThetaW);
    fwrite(&Ms, sizeof(int), 1, f_IntPred_ThetaW);
    fwrite(Pred_ThetaW1, sizeof(int), nviews*maxiS*63, f_IntPred_ThetaW);
    fclose(f_IntPred_ThetaW);
    
    
//    for ( i=0; i<nviews*63*maxiS; i++) {
//        Pred_ThetaW[i]  = (double) ((double)Pred_ThetaW1[i]/pow((double)2,12));
//    }
}

void ptcd_decode_predictor_mask() {
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int Ms = 0;
    int maxiS = 0;
    aux_read_pred_header(&Ms, &maxiS);
    
    int i, j, ik, i1;
    unsigned long unpacked, len_vec, len_rep;
    unsigned char *repack;
    FILE *fid = NULL;
    
    len_vec = nviews*(long)Ms*maxiS;
    ik = len_vec -(int)(((int)len_vec)/4)*4;
    len_rep = ik + (((int)len_vec)/4)*3;
    repack = (unsigned char*)malloc( len_rep*sizeof(char) );
    
    fid = fopen ( "bs_Pred_RegWPacked.txt" , "rb" );
    if (fid == NULL)
    { fputs("File error",stderr);  exit(1); }
    //Neigh9_165 = alocaVector((int)nviews*9);
    fread(repack,sizeof(unsigned char),len_rep,fid);
    fclose(fid);
    
    //for (i=0; i<len_rep; i++)
    //                         repack[i] = 0;
    
    //printf("a1\t %d \n",len_rep);
    
    int* Pred_RegW = (int*)malloc(nviews*maxiS*63*sizeof(int));

    
    i=0; j = 0;
    while( (j+2) < len_rep )
    {
        unpacked = (unsigned long) ( ((unsigned long)repack[j+2]) << 16 |  ((unsigned long)repack[j+1]) << 8 | (unsigned long)repack[j] );
        
        //printf("%i \n", unpacked);
        
        j = j+3;
        // repack on bytes
        Pred_RegW[i] = (int)(unpacked & (unsigned long)0x3f);
        Pred_RegW[i+1] = (int)( (unpacked>>6 & (unsigned long)0x3f) );
        Pred_RegW[i+2] = (int)( (unpacked>>12 & (unsigned long)0x3f) );
        Pred_RegW[i+3] = (int)( (unpacked>>18 & (unsigned long)0x3f) );
        
        //printf("%d %d %d \n",(unsigned char)repack[j],(unsigned char)repack[j+1],(unsigned char)repack[j+2]);
        //printf("%d %d %d %d \n", Pred_RegW[i], Pred_RegW[i+1], Pred_RegW[i+2], Pred_RegW[i+3] );
        //
        //printf("%d ",unpacked & (unsigned long)0x3f);
        //printf("%d ",(unpacked>>6 & (unsigned long)0x3f));
        //printf("%d ",(unpacked>>12 & (unsigned long)0x3f));
        //printf("%d ",(unpacked>>18 & (unsigned long)0x3f));
        //printf("\n i\t %d \t j\t %d \n",i,j);
        
        i = i+4;
        
        //exit(1);
    }
    
    //printf("a2\n");
    
    // treat the last part, if any
    for( i1 = (i); i1<len_vec; i1++)
    {
        Pred_RegW[i1] = (int)repack[j+(i1-i)];
        //repack[j+(i1-i)] = (char)Pred_RegW[i1];
    }
    //printf("%d %d \n",len_rep-1,j+len_vec-1-i);
    if(len_rep-1 != j+len_vec-1-i)
        printf("Error of counting chars");
    
    //fid = fopen ( "Pred_RegWPacked.txt" , "wb" );
    //if (fid == NULL)
    //{ fputs("File error2",stderr);  exit(1); }
    //fwrite(repack,sizeof(char),len_rep ,fid);
    //fclose(fid);
    
    FILE* f_Pred_RegW = fopen ( "PRMA" , "wb" );
    fwrite( &nviews, sizeof(int), 1, f_Pred_RegW);
    fwrite( &maxiS ,sizeof(int),1,f_Pred_RegW);
    fwrite( &Ms ,sizeof(int),1,f_Pred_RegW);
    fwrite(Pred_RegW,sizeof(int),nviews*63*maxiS,f_Pred_RegW);
    fclose(f_Pred_RegW);
    
    if( repack != NULL )
        free(repack);
    
}
