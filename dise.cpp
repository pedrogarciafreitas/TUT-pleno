// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module dise -- Disparity encoder.

#include <time.h>
#include "dise.h"
#include "gen_types.hh"
#include "cerv/cerv.h"

static void CC_CreateStructForSegmentPart1(int SEGgamma[], int LFgamma[], int SEGMFINAL[], int Pred_Reg[], int maxiS)
{
    
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int MIr, MIc, MIr1, MIc1, in, iview, i, iR, iT, icomp, maxi, yd;
    int iT1, in1, phi2, imaxi, ilenR;
    double cmaxi, cc;
    FILE *f_Pred_Reg;
    
    double *phi, *PHI, *PSI, *PHIdiag, *ydi;
    ydi = alocaDoubleVector((int)3);
    phi = alocaDoubleVector((int)241*3);
    PHI= alocaDoubleVector((int)241*241);
    PSI= alocaDoubleVector((int)241);
    PHIdiag= alocaDoubleVector((int)241);
    
    
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
    // find the maximum value in SEGMFINAL; there are so many potential regions
    
    for ( iR=(maxiS-1); iR>=0; iR--)
        Pred_Reg[nviews*iR] = 120;  // Pred_Reg[0,iR] no displacements
    //printf(" Passed 1 maxiS [%d]\n",maxiS);
    // Take all side views in turn to find their segmentation
    for (iview=1; iview<nviews; iview++)
    {
        // printf(" iview [%d]\n",iview);
        // Take all regions, from farthest to closest
        for ( iR=(maxiS-1); iR>=0; iR--)
        {	 ilenR = 0;
            for (iT=0; iT<241;iT++)
            {
                PSI[iT] = 0.;
                PHIdiag[iT] = 0.;
            }
            // printf(" Passed 1 iR [%d]\n",iR);
            for (MIr=9; MIr<(nr-9); MIr++)  // (MIr=0; MIr<nr; MIr++)
            {
                for (MIc=9; MIc<(nc-9);MIc++)	// (MIc=0; MIc<nc;MIc++)
                {
                    if( SEGMFINAL[MIr+MIc*nr] == (iR+1) )
                    { // In central view at (MIr,MIc);
                        // imsv(MIr,MIc,icomp)=LFgamma(iview,MIr,MIc,icomp)
                        ilenR = ilenR +1;
                        for (icomp=0; icomp<3;icomp++)
                        {
                            in = ((icomp*nc+MIc)*nr+MIr)*nviews+iview;
                            ydi[icomp] = (double)LFgamma[in]/(double)1023;
                        }
                        //collect the vector phi from im0(MIr,MIc,icomp)
                        // im0(MIr,MIc,icomp)=LFgamma(0,MIr,MIc,icomp)
                        for (iT=0; iT<241;iT++)
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
                            for (icomp=0; icomp<3;icomp++)
                            {
                                in = ((icomp*nc+MIc1)*nr+MIr1)*nviews;					// LFgamma[0,MIr1,MIc1,icomp]
                                phi[iT+icomp*241] = (double)LFgamma[in]/(double)1023;    // phi[iT,icomp]
                            }
                        }
                        for (iT=0; iT<241;iT++)
                        {
                            for (icomp=0; icomp<3;icomp++)
                            {
                                PSI[iT] = PSI[iT] + (double)phi[iT+icomp*241]*(double)ydi[icomp];
                                PHIdiag[iT] = PHIdiag[iT] + (double)phi[iT+icomp*241]*(double)phi[iT+icomp*241];
                            } // for (icomp=0; icomp<3;icomp++)
                            //for (iT1=0; iT1<241;iT1++)
                            //{
                            //	phi2 = phi[iT+icomp*241] *phi[iT1+icomp*241];
                            //	in = iT*241+iT1;  // PHI[iT,iT1]
                            //	PHI[in] = PHI[in] + phi2;
                            //	in1 = iT1*241+iT;  // PHI[iT,iT1]
                            //	PHI[in1] = PHI[in1] + phi2;
                            //}
                        }
                    } // if( SEGMFINAL[MIr+MIc*nr] == (iR+1) )
                }
            } //for (MIr=0; MIr<nr; MIr++)
            
            // printf(" Passed 1 iR [%d]\n",iR);
            // Find the best basis
            cmaxi = 0; imaxi = 0;
            for (iT=0; iT<241;iT++)
            {
                cc =  (PSI[iT]*PSI[iT])/PHIdiag[iT];
                if( cc>cmaxi)
                {
                    cmaxi = cc;
                    imaxi = iT;
                }
            }
            //if( iR==1 && iview==8 )
            //{
            //	printf(" cmaxi imaxi , iR , iview [%lf][%d][%d][%d]\n",cmaxi,imaxi,iR, iview);
            //	printf(" PHIdiag[imaxi],PSI[imaxi],ilenR [%lf][%lf][%d] \n",PHIdiag[imaxi],PSI[imaxi],ilenR);
            //}
            
            in = iR*nviews+iview ; // Pred_Reg[iview,iR]
            Pred_Reg[in] = imaxi;
        }// for ( iR=maxi; iR>0; iR--)
    } //for (iview=1; iview<nviews; iview++)
    //f_Pred_Reg = fopen ( "f_Pred_Reg.txt" , "wb" );
    //if (f_Pred_Reg == NULL) 
    //{ fputs("File error",stderr);  exit(1); }
    //fwrite(Pred_Reg, sizeof(int), maxiS*nviews, f_Pred_Reg);
    //fclose(f_Pred_Reg);
    
//    int** Pred_Reg_im = (int**)malloc(nviews*sizeof(int*));
//    for(i = 0; i < nviews; ++i ) {
//        Pred_Reg_im[i] = (int*)malloc(maxiS*sizeof(int));
//    }
//    for(i = 0; i < nviews; ++i ) {
//        for(int j = 0; j < maxiS; ++j ) {
//            Pred_Reg_im[i][j] = Pred_Reg[i+j*nviews];
//        }
//    }
//    cerv_encode(Pred_Reg_im, nviews, maxiS, "Pred_Reg.txt");
//    for(i = 0; i < nviews; ++i ) {
//        free(Pred_Reg_im[i]);
//    }
//    free(Pred_Reg_im);
    
    
    //free memory for double *phi, *PHI, *PSI, *PHIdiag, *ydi;
    if(phi != NULL)
        free(phi);
    if(PHI != NULL)
        free(PHI);
    if(PSI != NULL)
        free(PSI);
    if(PHIdiag != NULL)
        free(PHIdiag);
    if(ydi != NULL)
        free(ydi);
    
}

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

int dise(const char* filename_depth) {
    
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    int nrr = nr*nvr;
    int ncc = nc*nvc;
    int ncols = 3;
    
    // read the color views
    FILE* f_RV = fopen("VIN", "rb");
    fread(&nviews, sizeof(int), 1, f_RV);
    fread(&nr, sizeof(int), 1, f_RV);
    fread(&nc, sizeof(int), 1, f_RV);
    int ncol = 0;
    fread(&ncol, sizeof(int), 1, f_RV);
    int* LFgamma = (int*)malloc(nviews*nr*nc*ncol*sizeof(int));
    fread(LFgamma, sizeof(int), nviews*nr*nc*ncol, f_RV);
    fclose(f_RV);
    
    // read depth
    int* SEGMFINAL = alocaVector(nr*nc);
    FILE* f_SEGMFINAL = fopen ( filename_depth, "rb" );
    aux_read8ppm(f_SEGMFINAL, nr, nc, SEGMFINAL);
    fclose(f_SEGMFINAL);
    
    // Find number of regions
    int maxiS = 0;
    for (int i=0; i <nr*nc; i++) {
        if( SEGMFINAL[i] > maxiS ) {
            maxiS = SEGMFINAL[i];
        }
    }
    //printf("Main maxiS [%d]\n",maxiS);
    
    
    // create SEGgamma
    int* SEGgamma = (int*)malloc(nr*nc*nviews*sizeof(int));
    
    int* Pred_Reg = alocaVector(maxiS*nviews);
    // Find displacements in Pred_Reg and save them on disk
    time_t begin = clock();
    CC_CreateStructForSegmentPart1(SEGgamma, LFgamma, SEGMFINAL, Pred_Reg, maxiS);
    time_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time for CreateStructForSegmentPart1: %f\n", elapsed_secs);
    
    begin = clock();
    CC_CreateStructForSegmentPart2(SEGgamma, SEGMFINAL, Pred_Reg, maxiS);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time for CreateStructForSegmentPart2: %f\n", elapsed_secs);
    
    FILE* f_D = fopen("DISP", "wb");
    fwrite( &nviews ,sizeof(int),1,f_D);
    fwrite( &maxiS ,sizeof(int),1,f_D);
    fwrite( Pred_Reg ,sizeof(int), maxiS*nviews,f_D);
    fclose(f_D);
    
    FILE* f_SEGM = fopen("SEGM", "wb");
    fwrite( &nviews ,sizeof(int),1,f_SEGM);
    fwrite( &nr ,sizeof(int),1,f_SEGM);
    fwrite( &nc ,sizeof(int),1,f_SEGM);
    fwrite( SEGgamma ,sizeof(int), nr*nc*nviews,f_SEGM);
    fclose(f_SEGM);
    
    f_SEGM = fopen("LFGAMMA", "wb");
    fwrite( &nviews ,sizeof(int),1,f_SEGM);
    fwrite( &nr ,sizeof(int),1,f_SEGM);
    fwrite( &nc ,sizeof(int),1,f_SEGM);
    fwrite( &ncols, sizeof(int), 1, f_SEGM);
    fwrite( LFgamma ,sizeof(int), nr*nc*nviews*ncols,f_SEGM);
    fclose(f_SEGM);
    
    free(LFgamma);
    free(SEGgamma);
    free(Pred_Reg);
    
    return maxiS;

}


void dise_encode_displacements() {
    //    int nr = 0;
    //    int nc = 0;
    //    int nvr = 0;
    //    int nvc = 0;
    //    aux_read_header(&nr, &nc, &nvr, &nvc);
    //    int nviews = nvr*nvc;
    
    FILE* f_Pred_Reg = fopen("DISP", "rb");
    int nviews = 0;
    fread(&nviews, sizeof(int), 1, f_Pred_Reg);
    int maxiS = 0;
    fread(&maxiS, sizeof(int), 1, f_Pred_Reg);
    int* Pred_Reg = (int*)malloc(nviews*maxiS*sizeof(int));
    fread(Pred_Reg, sizeof(int), nviews*maxiS, f_Pred_Reg);
    fclose(f_Pred_Reg);
    
    int** Pred_Reg_im = (int**)malloc(nviews*sizeof(int*));
    for(int i = 0; i < nviews; ++i ) {
        Pred_Reg_im[i] = (int*)malloc(maxiS*sizeof(int));
    }
    for(int i = 0; i < nviews; ++i ) {
        for(int j = 0; j < maxiS; ++j ) {
            Pred_Reg_im[i][j] = Pred_Reg[i+j*nviews];
        }
    }
    cerv_encode(Pred_Reg_im, nviews, maxiS, "bs_Pred_Reg.txt");
    for(int i = 0; i < nviews; ++i ) {
        free(Pred_Reg_im[i]);
    }
    free(Pred_Reg_im);
}

void dise_encode_depth_map(const char* filename_depth) {
    
    // read header
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    
    
    // read depth
    int* SEGMFINAL = alocaVector(nr*nc);
    FILE* f_SEGMFINAL = fopen ( filename_depth, "rb" );
    aux_read8ppm(f_SEGMFINAL, nr, nc, SEGMFINAL);
    fclose(f_SEGMFINAL);
    
    int** SEGM2D = (int**)malloc(nr*sizeof(int*));
    for(int i=0; i < nr; ++i) { SEGM2D[i] = (int*)malloc(nc*sizeof(int));}
    for(int i=0; i < nr; ++i) {
        for(int j = 0; j < nc; ++j) {
            SEGM2D[i][j] = SEGMFINAL[i+j*nr];
        }
    }
    clock_t begin = clock();
    cerv_encode(SEGM2D, nr, nc, "bs_segm.bin");
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time for cerv: %f\n", elapsed_secs);
    for(int i=0; i < nr; ++i) { free(SEGM2D[i]); } free(SEGM2D);
}



