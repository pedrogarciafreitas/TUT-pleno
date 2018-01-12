// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module vptd -- View Prediction/Transform/Synthesis

#include <cstdlib>
#include <math.h>
#include "vptd.h"
#include "gen_types.hh"

static void CC_WarpAllEncoder4Part2(int SEGgamma[], int Neigh9_165[], int hatLFgamma[], int hatLFgamma2[],
                             int Pred_RegW[], double Pred_ThetaW[], int maxiS, int Ms)
{
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int MIr, MIc, MIr1, MIc1, in, in1, iview, i, iR, iT, iT1, iT2, iT3 ;
    int icomp, Mtrue;
    int  icount1, icount2, max_reg,i1,i2,j,found_it,idone;
    int iNN, nn1, ITT;
    
    
    int hexag_even_R1[7] = { -1, -1, 0, 0, 0,  1, 1 };
    int hexag_even_C1[7] = {  0, 1, -1, 0, 1,  0, 1 };
    int hexag_odd_C1[7] = {  0, 1, -1, 0, 1,  0, 1 };
    //    int hexag_odd_C1[7]  = { -1, 0, -1, 0, 1, -1, 0 };
    
    int *PredRegr0,  mTheta, MT, *Pred_ThetaW1;
    double *PredTheta0;
    double  *ycrt, phi2, yd2;
    FILE *f_hatLFgamma, *f_SEGgamma, *f_Neigh9_165, *f_Pred_RegW, *f_Pred_ThetaW;
    
    ycrt = alocaDoubleVector((int)3);
    PredTheta0= alocaDoubleVector((int)63);
    PredRegr0= alocaVector((int)63);
    
    
    for( i=0; i<nviews*nr*nc*3; i++)
        hatLFgamma2[i] = hatLFgamma[i];
    
    // Take all side views in turn to find their warping
    for (iview=0; iview<nviews; iview++)
    {
        
        
        // Take this side view to find its warping
        
        //printf(" Enc4 iview [%d]\n",iview);
        // Take all regions, from farthest to closest
        for ( iR=(maxiS-1); iR>=0; iR--)
        {
            //printf("did it 1\n");
            
            for( i=0; i<Ms;i++ )
            {
                PredRegr0[i] = Pred_RegW[iview+nviews*(iR+maxiS*i)];// Pred_RegW[iview, iR, i] = PredRegr[i];
                PredTheta0[i] = Pred_ThetaW[iview+nviews*(iR+maxiS*i)];
            }
            for (MIr=2; MIr<(nr-2); MIr++)  // (MIr=0; MIr<nr; MIr++)
            {
                for (MIc=2; MIc<(nc-2);MIc++)	// (MIc=0; MIc<nc;MIc++)
                {
                    if( SEGgamma[(MIr+MIc*nr)*nviews+iview] == (iR+1) ) // SEGgamma[iview,MIr,MIc)
                    {
                        //printf("did it 2\n");
                        //ilenR = ilenR +1; // We are in iview, at region (iR+1)
                        // Target is side view at (MIr,MIc);
                        // imsv(MIr,MIc,icomp)=LFgamma(iview,MIr,MIc,icomp)
                        for (icomp=0; icomp<3;icomp++)
                        {
                            ycrt[icomp] = 0.;
                        }
                        //collect the vector phi from neighbor views
                        ITT = 0; // index in phi
                        for (iNN=0; iNN<9;iNN++)
                        {
                            //printf("did it 3\n");
                            // ATTENTION!!! Neigh9_165 is given in (1..165). Take one out in C.
                            nn1 = Neigh9_165[iNN*nviews+iview];   // Neigh9_165[iview,iNN]
                            // printf(" Enc4 nn1 iview [%d][%d] \n", nn1,iview);
                            if( nn1 > 0)
                            { // Use neighbor nn1 in constructing phi
                                // im0(MIr,MIc,icomp)=LFgamma(0,MIr,MIc,icomp)
                                nn1 = nn1-1; // This is the view index in C !!!!!!!!!!
                                // printf(" Enc4 nn1 iview [%d][%d] \n", nn1,iview);
                                for (iT=0; iT<7;iT++)
                                {
                                    found_it = 0;
                                    for( idone=0; idone<Ms; idone++)
                                    {
                                        if(ITT == PredRegr0[idone])
                                        {
                                            found_it =1;
                                            break;
                                        }
                                    }
                                    if( found_it == 1 )
                                    {
                                        // multiply the regressor by PredTheta0[idone]
                                        if(1)//(MIr%2 == 1)
                                        {
                                            MIr1 = MIr + hexag_even_R1[iT];
                                            MIc1 = MIc + hexag_even_C1[iT];
                                        }
                                        else
                                        {
                                            MIr1 = MIr + hexag_even_R1[iT];
                                            MIc1 = MIc + hexag_odd_C1[iT];
                                        }
                                        for (icomp=0; icomp<3;icomp++)
                                        {
                                            in = ((icomp*nc+MIc1)*nr+MIr1)*nviews+nn1;		 // hatLFgamma[nn1,MIr1,MIc1,icomp]
                                            // phi[icomp*63+ITT] = (double)hatLFgamma[in]/(double)256 ;    // phi[ITT,icomp]
                                            ycrt[icomp] = ycrt[icomp] + PredTheta0[idone]*(double)hatLFgamma[in];
                                        }
                                    }
                                    ITT = ITT + 1;
                                } // for (iT=0; iT<7;iT++)
                            }
                        } // for (iNN=0; iNN<9;iNN++)
                        // printf("j [%d]",j);
                        //printf("did it 4\n");
                        for (icomp=0; icomp<3;icomp++)
                        {
                            in = ((icomp*nc+MIc)*nr+MIr)*nviews+iview;
                            hatLFgamma2[in]= (int)ycrt[icomp]; // LFgamma(iview,MIr,MIc,icomp)
                        }
                        //printf("did it 5\n");
                    } // if( SEGMFINAL[MIr+MIc*nr] == (iR+1) )
                } //for (MIc=2; MIc<(nc-2);MIc++)
            } //for (MIr=2; MIr<nr; MIr++)
        } // for ( iR=(maxiS-1); iR>=0; iR--)
        
        // NOW SAVE THE BETTER VERSION AT THE LOCATION OF THIS VIEW
        
        
        for (MIr=2; MIr<nr; MIr++)
            for (MIc=2; MIc<(nc-2);MIc++)
                for (icomp=0; icomp<3;icomp++)
                {
                    in = ((icomp*nc+MIc)*nr+MIr)*nviews+iview;
                    hatLFgamma[in] = hatLFgamma2[in];  // hatLFgamma[iview,MIr,MIc,icomp] iview+nviews*(MIr+nr*(MIc+nc*icomp))
                }
        
        // END OF PART COPIED IDENTICALLY FROM DECODING (FROM PART2)
    } // for (iview=1; iview<nviews; iview++)
    
    // take out the estimated better image
    for( i=0; i<nviews*nr*nc*3; i++)
        hatLFgamma2[i] = hatLFgamma[i];
    //printf("did it [%d]\n",sizeof(long));
    
    //free memory for double *phi, *PHI, *PSI, *PHIdiag, *ydi;
    
    if(PredRegr0 != NULL)
        free(PredRegr0);
    if(PredTheta0 != NULL)
        free(PredTheta0);
    //if(Pred_ThetaW1 != NULL)
    //    free(Pred_ThetaW1);
    if(ycrt != NULL)
        free(ycrt);
    
    //if(hexag_even_R1 != NULL)
    //	free(hexag_even_R1);
    //if(hexag_even_C1 != NULL)
    //	free(hexag_even_C1);
    //if(hexag_odd_C1 != NULL)
    //	free(hexag_odd_C1);
    
    //printf("did it too[%d]\n",sizeof(long));
    
}


void vptd() {
    int nr = 0;
    int nc = 0;
    int nvr = 0;
    int nvc = 0;
    aux_read_header(&nr, &nc, &nvr, &nvc);
    int nviews = nvr*nvc;
    
    int Ms = 0;
    int maxiS = 0;
    aux_read_pred_header(&Ms, &maxiS);
    
    // read segmentation for all views
    FILE* f_SEGM = fopen("SEGM", "rb");
    fread(&nviews, sizeof(int), 1, f_SEGM);
    fread(&nr, sizeof(int), 1, f_SEGM);
    fread(&nc, sizeof(int), 1, f_SEGM);
    int* SEGgamma = (int*) malloc(nr*nc*nviews*sizeof(int));
    fread(SEGgamma, sizeof(int), nr*nc*nviews, f_SEGM);
    fclose(f_SEGM);
    
    // Read predictor coefficients
    FILE* f_IntPred_ThetaW = fopen("PRCO", "rb");
    fread(&nviews, sizeof(int), 1, f_IntPred_ThetaW);
    fread(&maxiS, sizeof(int), 1, f_IntPred_ThetaW);
    fread(&Ms, sizeof(int), 1, f_IntPred_ThetaW);
    int* IntPred_ThetaW = (int*)malloc(nviews*maxiS*63*sizeof(int));
    fread(IntPred_ThetaW, sizeof(int), nviews*maxiS*63, f_IntPred_ThetaW);
    fclose(f_IntPred_ThetaW);
    
    double* Pred_ThetaW = (double*)malloc(nviews*maxiS*63*sizeof(double));
    for (int i=0; i<nviews*63*maxiS; i++) {
        Pred_ThetaW[i]  = (double) ((double)IntPred_ThetaW[i]/pow((double)2,12));
    }
    
    // read predictor masks
    FILE* f_Pred_RegW = fopen("PRMA", "rb");
    fread(&nviews, sizeof(int), 1, f_Pred_RegW);
    fread(&maxiS, sizeof(int), 1, f_Pred_RegW);
    fread(&Ms, sizeof(int), 1, f_Pred_RegW);
    int* Pred_RegW = (int*)malloc(nviews*maxiS*63*sizeof(int));
    fread(Pred_RegW, sizeof(int), nviews*maxiS*63, f_Pred_RegW);
    fclose(f_Pred_RegW);
    

    int nrr = 0;
    int ncc = 0;
    int ncols = 0;
    // get estimate for LFgamma
    FILE* f_DRV = fopen("DRV", "rb");
    fread(&nviews, sizeof(int), 1, f_DRV);
    fread(&nr, sizeof(int), 1, f_DRV);
    fread(&nc, sizeof(int), 1, f_DRV);
    fread(&ncols, sizeof(int), 1, f_DRV);
    int* hatLFgamma = (int*) malloc(nr*nc*nviews*ncols*sizeof(int));
    fread(hatLFgamma, sizeof(int), nr*nc*nviews*ncols, f_DRV);
    fclose(f_DRV);
    
    int* Neigh9_165 = alocaVector((int)nviews*9);
    for(int i = 0; i < nviews*9; ++i ) {
        Neigh9_165[i] = Neigh9_165p_const[i];
    }
    
    int* hatLFgamma2 = (int*) malloc(nr*nc*nviews*ncols*sizeof(int));
    
    CC_WarpAllEncoder4Part2(SEGgamma, Neigh9_165, hatLFgamma, hatLFgamma2, Pred_RegW, Pred_ThetaW, maxiS, Ms);
    
    unsigned short* reconsLI = (unsigned short*) malloc(nrr*ncc*ncols*sizeof( unsigned short));

//    aux_Inverse_LFgamma_TO_LI(reconsLI, hatLFgamma2, hatLI);
//    
//    int* LI = (int*) malloc(nrr*ncc*ncols*sizeof(int));
//    for(int i = 0; i < nrr*ncc*ncols; ++i ) {
//        LI[i] = static_cast<int>(reconsLI[i]);
//    }
    int ncol = 3;
    
    
    FILE* f_RV = fopen("VIND", "wb");
    fwrite( &nviews ,sizeof(int),1,f_RV);
    fwrite( &nr ,sizeof(int),1,f_RV);
    fwrite( &nc ,sizeof(int),1,f_RV);
    fwrite( &ncol, sizeof(int), 1, f_RV);
    fwrite( hatLFgamma2,sizeof(int),nviews*nr*nc*ncol,f_RV);
    fclose(f_RV);
    
    free(hatLFgamma2);
    free(SEGgamma);
    free(reconsLI);
    free(Neigh9_165);
}
