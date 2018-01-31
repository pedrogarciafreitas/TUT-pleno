// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module disd -- Disparity decoder.

#include <cstdio>
#include <cstdlib>
#include "disd.h"
#include "gen_types.hh"
#include "cerv/cerv.h"

static void CC_CreateStructForSegmentPart2(int SEGgamma[], int SEGMFINAL[], int Pred_Reg[], int maxiS)
{
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int MIr, MIc, MIr1, MIc1, in, iview, i, iR, iT ;
    int  icount1, icount2, max_reg,i1,i2;
    FILE *f_SEGgamma;
    
    // Set SEGgamma(0,:,:) = SEGMFINAL
    for (MIc=0; MIc<nc; MIc++)
    {
        for (MIr=0; MIr<nr; MIr++)
        {
            // hint: linear index LFgamma(ig,MIr,MIc,icomp) :  ((icomp*nc+MIc)*nr+MIr)*nviews+ig
            // linear index for SEGgamma(ig,MIr,MIc) : (MIc*nr+MIr)*nviews+ig
            // linear index for SEGMFINAL(MIr,MIc) : (MIc*nr+MIr)
            in = MIc*nr+MIr;
            SEGgamma[in*nviews] = SEGMFINAL[in];
        }
    }
    //printf(" Passed 1  \n" );
    
    // Take all side views in turn to find their segmentation
    for (iview=1; iview<nviews; iview++)
    {
        // printf(" iview [%d]\n",iview);
        // Take all regions, from farthest to closest
        for ( iR=(maxiS-1); iR>=0; iR--)
        {
            in = iR*nviews+iview ; // Pred_Reg[iview,iR]
            iT = Pred_Reg[in]; // imaxi
            // printf(" Passed 1 iR [%d]\n",iR);
            for (MIr=0; MIr<(nr); MIr++)  // (MIr=0; MIr<nr; MIr++)
            {
                for (MIc=0; MIc<(nc); MIc++)	// (MIc=0; MIc<nc;MIc++)
                {
                    if(1)//MIr%2 == 1)
                    {
                        MIr1 = MIr + hexag_even_R_const[iT];
                        MIc1 = MIc + hexag_even_C_const[iT];
                    }
                    else
                    {
                        MIr1 = MIr + hexag_even_R_const[iT];
                        MIc1 = MIc + hexag_odd_C_const[iT];
                    }
                    if( MIr1>0 && MIr1<nr && MIc1>0 && MIc1<nc )
                        // we can access SEGMFINAL[MIr1,MIc1]
                        if( SEGMFINAL[MIr1+MIc1*nr] == (iR+1) )
                        { // here we should warp reg iR, but only if it is not occupied with a closer region
                            // since we scan from farthest to closest, there is no worry
                            
                            // linear index for SEGgamma(ig,MIr,MIc) : (MIc*nr+MIr)*nviews+ig
                            // linear index for SEGMFINAL(MIr,MIc) : (MIc*nr+MIr)
                            in = MIc*nr+MIr;
                            SEGgamma[in*nviews+iview] = iR+1; // mark it to the final structure
                        }
                }
            }
        }
    }
    // A missing value exists when SEGgamma is 0
    icount1 = 0; icount2 = 0;
    for (iview=1; iview<nviews; iview++)
    {
        for (MIr=0; MIr<(nr); MIr++)  // (MIr=0; MIr<nr; MIr++)
        {
            for (MIc=0; MIc<(nc); MIc++)	// (MIc=0; MIc<nc;MIc++)
            {
                in = MIc*nr+MIr;
                if( SEGgamma[in*nviews+iview] == 0 )
                {
                    icount1 = icount1+1;
                    // Missing pixel; Check if in the previous view it was set to a region
                    if( SEGgamma[in*nviews+iview-1] != 0 )
                    {
                        max_reg = SEGgamma[in*nviews+iview-1];
                        for(i1=0; i1<(2); i1++)
                            for(i2=0; i2<(2); i2++)
                            {
                                MIr1 = MIr + i1;
                                MIc1 = MIc + i2;
                                if( MIr1>0 && MIr1<nr && MIc1>0 && MIc1<nc )
                                    if( SEGgamma[(MIc1*nr+MIr1)*nviews+iview-1] > max_reg )
                                        max_reg = SEGgamma[(MIc1*nr+MIr1)*nviews+iview-1];
                            }
                        SEGgamma[in*nviews+iview] = max_reg;
                        icount2 = icount2 +1;
                    }
                }
            }
        }
        
        
        printf("   iview, icount1,icount2  [%d][%d][%d]\n", iview, icount1,icount2);
        //if( iR==1 && iview==8 )
        //{
        //
        //	printf(" PHIdiag[imaxi],PSI[imaxi],ilenR [%lf][%lf][%d] \n",PHIdiag[imaxi],PSI[imaxi],ilenR);
        //}
    } //for (iview=1; iview<nviews; iview++)
    
    // write SEGgamma on file
    //f_SEGgamma = fopen ( "f_SEGgamma.txt" , "wb" );
    //if (f_SEGgamma == NULL)
    //{ fputs("File error",stderr);  exit(1); }
    
    //fwrite(SEGgamma,sizeof(int),(long)nr*(long)nc*nviews,f_SEGgamma);
    //fclose(f_SEGgamma);
    
}


void disd() {
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    // read center view segmentation
    FILE* f_SEG = fopen("SEGC", "rb");
    fread(&nr, sizeof(int), 1, f_SEG);
    fread(&nc, sizeof(int), 1, f_SEG);
    int* SEGMFINAL = (int*) malloc(nr*nc*sizeof(int));
    fread(SEGMFINAL, sizeof(int), nr*nc, f_SEG);
    fclose(f_SEG);
    
    // Find number of regions
    int maxiS = 0;
    for (int i=0; i <nr*nc; i++) {
        if( SEGMFINAL[i] > maxiS ) {
            maxiS = SEGMFINAL[i];
        }
    }
    
    // Read displacements
    FILE* f_DISP = fopen("DISP", "rb");
    fread(&nviews, sizeof(int), 1, f_DISP);
    fread(&maxiS, sizeof(int), 1, f_DISP);
    int* Pred_Reg = (int*)malloc(nviews*maxiS*sizeof(int));
    fread(Pred_Reg, sizeof(int), nr*nc, f_DISP);
    fclose(f_DISP);


    // get SEGgamma
    int* SEGgamma = alocaVector(nr*nc*nviews);
    CC_CreateStructForSegmentPart2(SEGgamma, SEGMFINAL, Pred_Reg, maxiS);
    
    FILE* f_SEGM = fopen("SEGM", "wb");
    fwrite( &nviews ,sizeof(int),1,f_SEGM);
    fwrite( &nr ,sizeof(int),1,f_SEGM);
    fwrite( &nc ,sizeof(int),1,f_SEGM);
    fwrite( SEGgamma ,sizeof(int), nr*nc*nviews,f_SEGM);
    fclose(f_SEGM);
    
    free(SEGMFINAL);
    free(SEGgamma);
    free(Pred_Reg);
}

void disd_decode_displacements() {
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    // read center view segmentation
    FILE* f_SEG = fopen("SEGC", "rb");
    fread(&nr, sizeof(int), 1, f_SEG);
    fread(&nc, sizeof(int), 1, f_SEG);
    int* SEGMFINAL = (int*) malloc(nr*nc*sizeof(int));
    fread(SEGMFINAL, sizeof(int), nr*nc, f_SEG);
    fclose(f_SEG);
    
    // Find number of regions
    int maxiS = 0;
    for (int i=0; i <nr*nc; i++) {
        if( SEGMFINAL[i] > maxiS ) {
            maxiS = SEGMFINAL[i];
        }
    }
    
    int* Pred_Reg = (int*)malloc(nviews*maxiS*sizeof(int));
    // Read displacements from Pred_Reg
    int** Pred_Reg_im = (int**)malloc(nviews*sizeof(int*));
    for(int i = 0; i < nviews; ++i ) {
        Pred_Reg_im[i] = (int*)malloc(maxiS*sizeof(int));
    }
    cerv_decode(Pred_Reg_im, nviews, maxiS, "bs_Pred_Reg.txt");
    for(int i = 0; i < nviews; ++i ) {
        for(int j = 0; j < maxiS; ++j ) {
            Pred_Reg[i+j*nviews] = Pred_Reg_im[i][j];
        }
    }
    for(int i = 0; i < nviews; ++i ) {
        free(Pred_Reg_im[i]);
    }
    free(Pred_Reg_im);
    
    FILE* f_D = fopen("DISP", "wb");
    fwrite( &nviews ,sizeof(int),1,f_D);
    fwrite( &maxiS ,sizeof(int),1,f_D);
    fwrite( Pred_Reg ,sizeof(int), maxiS*nviews,f_D);
    fclose(f_D);

    free(SEGMFINAL);
    free(Pred_Reg);
}

void disd_decode_depth_map() {
    
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    
    int* SEGMFINAL = alocaVector((int)nr*nc);
    //fread(SEGMFINAL,sizeof(int),(int)nr*nc,f_SEGMFINAL);
    int** SEGM2D = (int**)malloc(nr*sizeof(int*));
    for(int i=0; i < nr; ++i) { SEGM2D[i] = (int*)malloc(nc*sizeof(int));}
    for(int i=0; i < nr; ++i) {
        for(int j = 0; j < nc; ++j) {
            SEGM2D[i][j] = 0;
        }
    }
    cerv_decode(SEGM2D, nr, nc, "bs_segm.bin");
    for(int i=0; i < nr; ++i) {
        for(int j=0; j < nc; ++j) {
            SEGMFINAL[i+j*nr] = SEGM2D[i][j];
        }
    }
    for(int i=0; i < nr; ++i) { free(SEGM2D[i]); } free(SEGM2D);
    FILE* f_SEGMFINAL = fopen ( "SEGC" , "wb" );
    fwrite(&nr, sizeof(int), 1, f_SEGMFINAL);
    fwrite(&nc, sizeof(int), 1, f_SEGMFINAL);
    fwrite(SEGMFINAL,sizeof(int),nr*nc,f_SEGMFINAL);
    fclose(f_SEGMFINAL);
    
    free(SEGMFINAL);

}
