// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Module vpts -- View Prediction/Transform/Synthesis

#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "vpts.h"
#include "gen_types.hh"

int FastOLS(double PHI[], double PSI[], double yd2, int PredRegr0[], double PredTheta0[], int Ms, int MT, int MPHI)
{
    int mTheta, M, iM, iM1;
    double *B, *C,  sigerr, *Ag, *g;
    double C1, valm1, temp, crit, sabsval;
    int p, j_p, i, j, k, itemp;
    
    // Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
    // Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
    M = MT+1;
    B = alocaDoubleVector(M*M);// ((M+2)*(M+2));
    C = alocaDoubleVector(M*M);// ((M+2)*(M+2));
    Ag = alocaDoubleVector(M*M);// ((M+2)*(M+2));
    g = alocaDoubleVector(M);
    
    
    // Inputs: PHI is MTxMT, PSI is MTx1;
    // Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
    // Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM
    
    B[MT+MT*M] = yd2; //B[MT,MT] = yd2; // we start from B[0,0]
    for( iM=0; iM<MT; iM++)
    {
        PredRegr0[iM] = iM;
        B[iM+MT*M] = PSI[iM];// B[iM,MT] = PSI[iM];
        B[MT+iM*M] = PSI[iM];//B[MT,iM] = PSI[iM];
        for( iM1=0; iM1<MT; iM1++)
        {
            B[iM+iM1*M]=PHI[iM+iM1*MPHI];//B[iM,iM1]=PHI[iM,iM1];
        }
    }
    for( iM=0; iM<M; iM++)
        for( iM1=0; iM1<M; iM1++)
            C[iM+iM1*M] = 0;//C[iM,iM1] = 0
    for( iM=0; iM<MT; iM++)
        C[iM+iM*M] = 1;// C[iM,iM] = 1;
    crit = B[MT+MT*M];//crit = B[MT,MT];
    if(crit < 0.0000001)
    {
        printf("crir, yd2 [%f] [%f] ",crit, yd2 );
        i = 0;
        return i;
    }
    
    
    for( p=0; p<Ms; p++ )
    {
        valm1 = 0; j_p = 0; // pick the max value in next loop
        for( j=p; j<MT; j++)
        {
            //if(B[j+j*M] > 0.00000000000000001)
            sigerr  = B[j+MT*M]*B[j+MT*M]/B[j+j*M];//sigerr  = B[j,MT]*B[j,MT]/B[j,j];
            //else
            //	sigerr = 0;
            if( sigerr > valm1 )
            {
                valm1 = sigerr ;
                j_p = j;
            }
        } // j_p is the index of maximum
        crit = crit-valm1;
        itemp = PredRegr0[j_p]; PredRegr0[j_p] = PredRegr0[p]; PredRegr0[p] = itemp;
        for( j=p; j<M; j++)
        {
            //% interchange B(p:end,j_p) with B(p:end,p)
            temp = B[j+j_p*M]; //temp = B[j,j_p];
            B[j+j_p*M] = B[j+p*M]; //B[j,j_p] = B[j,p];
            B[j+p*M] = temp;// B[j,p] = temp;
        }
        for( j=p; j<M; j++)
        {
            //% interchange B(j_p,p:end) with B(p,p:end)
            temp = B[j_p+j*M]; // temp = B[j_p,j];
            B[j_p+j*M] = B[p+j*M]; //B[j_p,j] = B[p,j];
            B[p+j*M] = temp;//B[p,j] = temp;
        }
        
        //% Fast
        for( j=0; j<=p-1; j++)
        {
            // % interchange C(1:p-1,j_p) with C(1:p-1,p)
            temp = C[j+j_p*M]; // temp = C[j,j_p];
            C[j+j_p*M] = C[j+p*M]; //C[j,j_p] = C[j,p];
            C[j+p*M] = temp; //C[j,p] = temp;
        }
        for( j=(p+1); j<M; j++)
        {
            //if(B[p+p*M] > 0.00000000000000000000001)
            C[p+j*M] = B[p+j*M]/B[p+p*M];//C[p,j] = B[p,j]/B[p,p];
            //else
            //C[p+j*M] = 0;
        }
        
        for( j=(p+1); j<MT; j++)
            for( k=j; k<=MT; k++)
            {
                B[j+k*M] = B[j+k*M]-C[p+j*M]*C[p+k*M]*B[p+p*M];//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
            }
        for( j=(p+1); j<MT; j++)
            for( k=j; k<=MT; k++)
            {
                B[k+j*M] = B[j+k*M];//B[k,j] = B[j,k];
            }
        //for j = (p+1):M
        //    for k = j:(M+1)
        //        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
        //        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
        //        B(k,j) = B(j,k);
        //    end
        //end
//        for( iM=0; iM<Ms; iM++)
//        {
//        for( iM1=0; iM1<Ms; iM1++)
//        	printf("C[%f] ",C[iM+iM1*M]);
//        printf(" \n" );
//        }
        // scanf("%d",&i);
    } //% for( p=0; p<M; p++ )
    
    
    // final triangular backsolving
    for( i=0; i<Ms; i++ )
    {
        g[i] = C[i+MT*M];//g[i] = C[i,MT];
        for( j=0; j<Ms; j++ )
            Ag[i+j*M] = C[i+j*M];//Ag[i,j] = C[i,j];
    }
    PredTheta0[Ms-1] = g[Ms-1];
    for( i=Ms-2; i>=0; i-- )
    {
        PredTheta0[i] = g[i];
        for (j=i+1; j<Ms; j++)
            PredTheta0[i] = PredTheta0[i]- Ag[i+j*M]*PredTheta0[j];//PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
    }
     //printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
    // printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );
    
	if (PredTheta0[0] != PredTheta0[0])
	{// if is nan
		printf("PredTheta0[0]  is NaN\n");
		PredTheta0[0] = 1.0;
		for (i = 1; i < Ms; i++)
		{
			PredTheta0[i] = 0.0;
		}
	}

    sabsval = 0;
    for( i=0; i<Ms; i++ )
    {

		if (PredTheta0[i] != PredTheta0[i])
		{// if is nan
			PredTheta0[i] = 0.0;
			printf("PredTheta0[%i]  is NaN\n",i);
		}

        if(PredTheta0[i] > 0)
            sabsval = sabsval + PredTheta0[i];
        else
            sabsval = sabsval - PredTheta0[i];
    }
    //printf("%f\n", sabsval);
    if( sabsval > 2*Ms ) // if average coefficients are too high forget about intrpolation
    {
        PredTheta0[0] = C[0+MT*M];//g[0] = C[0,MT];
        
		//PredTheta0[0] = 1;

		if (PredTheta0[0] != PredTheta0[0])
		{// if is nan
			printf("C[0+MT*M]  is NaN\n", i);
			PredTheta0[0] = 1.0;
		}

		// fix ???
		//if (abs(PredTheta0[0]) > 100)
			//PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

        for( i=1; i<Ms; i++ ) {
            PredTheta0[i]  = 0.0;
        }
        
        i = 1;
    }
    else {
        i = Ms;
    }
    
    free(B); 
    free(C); 
    free(Ag);
    free(g);
    
    //printf("pred Ms [%d]",Ms);
    return i;
    
}



void CC_WarpAllEncoder4Part1(int SEGgamma[], int Neigh9_165[], int hatLFgamma[], int LFgamma[],
                             int Pred_RegW[], double Pred_ThetaW[], int maxiS, int Ms)
{
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
    
    int MIr, MIc, MIr1, MIc1, in, in1, iview, i, iR, iT, iT1, iT2, iT3 ;
    int icomp, Mtrue;
    int  icount1, icount2, max_reg,i1,i2,j;
    int iNN, nn1, ITT, MPHI,  found_it, idone;
    int *hatLFgamma2 = NULL;
    
    
    int hexag_even_R1[7] = { -1, -1, 0, 0, 0,  1, 1 };
    int hexag_even_C1[7] = {  0, 1, -1, 0, 1,  0, 1 };
     int hexag_odd_C1[7] = {  0, 1, -1, 0, 1,  0, 1 };
//    int hexag_odd_C1[7]  = { -1, 0, -1, 0, 1, -1, 0 };
    
    int *PredRegr0,  mTheta, MT;
    double *PredTheta0;
    double *phi, *PHI, *PSI, *PHIdiag, *ydi, phi2, yd2, *ycrt;
    FILE *f_hatLFgamma, *f_LFgamma, *f_SEGgamma, *f_Neigh9_165;
    
    ydi = alocaDoubleVector((int)3);
    phi = alocaDoubleVector((int)63*3);
    PHI= alocaDoubleVector((int)63*63);
    PSI= alocaDoubleVector((int)63);
    PHIdiag= alocaDoubleVector((int)63);
    PredTheta0= alocaDoubleVector((int)63);
    PredRegr0= alocaVector((int)63);
    hatLFgamma2 = alocaVector((int)nviews*nr*nc*3);
    ycrt = alocaDoubleVector((int)3);
    
    
    MPHI = 63;
    ///////////////////////////////////////////////////////////////////////////////
    
    
    for( i=0; i<nviews*nr*nc*3; i++)
        hatLFgamma2[i] = hatLFgamma[i];
    
    //printf(" Passed 2  \n" );
    // Take all side views in turn to find their warping
    for (iview=0; iview<nviews; iview++)
    {
        printf(" Enc4 iview [%d]\n",iview);
        // Take all regions, from farthest to closest
        for ( iR=(maxiS-1); iR>=0; iR--)
        {	 // initialize System of Equation matrices at region iR
            yd2 =0; j=0;
            for (iT=0; iT<63;iT++)
            {
                PSI[iT] = 0.;
                for (iT1=0; iT1<63;iT1++)
                {
                    PHI[iT*63+iT1] = 0.;
                    PHI[iT1*63+iT] = 0.;
                }
            }
            // printf(" Passed 1 iR [%d]\n",iR);
            ITT = 0; // index in the vector phi associated to the pixel (MIr,MIc)
            
            for (MIr=2; MIr<(nr-2); MIr++)  // (MIr=0; MIr<nr; MIr++)
            {
                for (MIc=2; MIc<(nc-2);MIc++)	// (MIc=0; MIc<nc;MIc++)
                {
                    if( SEGgamma[(MIr+MIc*nr)*nviews+iview] == (iR+1) ) // SEGgamma[iview,MIr,MIc)
                    {
                        //ilenR = ilenR +1; // We are in iview, at region (iR+1)
                        // Target is side view at (MIr,MIc);
                        // imsv(MIr,MIc,icomp)=LFgamma(iview,MIr,MIc,icomp)
                        j = j+1;
                        for (icomp=0; icomp<3;icomp++)
                        {
                            in = ((icomp*nc+MIc)*nr+MIr)*nviews+iview;
                            ydi[icomp] = (double)LFgamma[in]/(double)1024 ; // LFgamma(iview,MIr,MIc,icomp)
                            yd2 = yd2 + ydi[icomp]*ydi[icomp];
                            // printf(" Enc4 ydi[icomp] [%f]\n",ydi[icomp]);
                        }
                        //collect the vector phi from neighbor views
                        ITT = 0; // index in the vector phi associated to the pixel (MIr,MIc)
                        for (iNN=0; iNN<9;iNN++)
                        {
                            // ATTENTION!!! Neigh9_165 is given in (1..165). Take one out in C.
                            nn1 = Neigh9_165[iNN*nviews+iview];   // Neigh9_165[iview,iNN]
                            // printf(" Enc4 nn1 iview [%d][%d] \n", nn1,iview);
                            if( nn1 > 0)
                            { // Use neighbor nn1 in constructing phi
                                // phi(ITT) = hatLFgamma(nn1,MIr+delta1,MIc+delta2,icomp)
                                nn1 = nn1-1; // This is the view index in C !!!!!!!!!!
                                // printf(" Enc4 nn1 iview [%d][%d] \n", nn1,iview);
                                for (iT=0; iT<7;iT++)
                                {
									if (MIr % 2 == 1) /* should this be if (1) ?*/
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
                                        phi[icomp*63+ITT] = (double)hatLFgamma[in]/(double)1024 ;    // phi[ITT,icomp]
                                    }
                                    ITT = ITT + 1; // we added an entry to vector phi
                                } // for (iT=0; iT<7;iT++)
                            } // if( nn1 > 0)
                        } // for (iNN=0; iNN<9;iNN++)
                        //
                        if( ITT > 0)
                            // length of phi is ITT
                            // printf(" Enc4 iview [%f][%d][%d]\n", yd2, iview, ITT);
                            for (iT=0; iT<ITT;iT++)
                            {
                                for (icomp=0; icomp<3;icomp++)
                                {
                                    PSI[iT] = PSI[iT] + (double)phi[iT+icomp*63]*(double)ydi[icomp]; // phi(IT,icomp)
                                    // PHIdiag[iT] = PHIdiag[iT] + (double)phi[iT+icomp*63]*(double)phi[iT+icomp*63];
                                    for (iT1=0; iT1<ITT;iT1++)
                                    {
                                        phi2 = phi[iT+icomp*63] *phi[iT1+icomp*63]; // phi(IT,icomp)
                                        in = iT*MPHI+iT1;  // PHI[iT,iT1]
                                        PHI[in] = PHI[in] + phi2;
                                    }
                                }// for (icomp=0; icomp<3;icomp++)
                            } // for (iT=0; iT<ITT;iT++)
                        
                    } // if( SEGMFINAL[MIr+MIc*nr] == (iR+1) )
                } //for (MIc=2; MIc<(nc-2);MIc++)
            } //for (MIr=2; MIr<nr; MIr++)
            
            //printf("j, ITT [%d [%d]\n",j, ITT);

            
            if(1)
            {
                // printf("PredTheta0 [%f]",yd2);
                //	for(i=0; i<7; i++)
                //	  printf("PredTheta0 [%f][%f][%f][%f][%f][%f]\n",PHI[i+0*63],PHI[i+1*63],PHI[i+2*63],PHI[i+3*63],PHI[i+4*63],PHI[i+5*63]);
                // MT = 63;
                // Ms = 7;
                Mtrue = FastOLS(PHI, PSI, yd2, PredRegr0, PredTheta0, Ms, ITT, MPHI);

                if( (Mtrue>0) && (Mtrue<=Ms))
                {
                    for( i=0; i<Ms;i++ )
                    {
                        Pred_RegW[iview+nviews*(iR+maxiS*i)] = PredRegr0[i];// Pred_RegW[iview, iR, i] = PredRegr[i];
                        Pred_ThetaW[iview+nviews*(iR+maxiS*i)] = PredTheta0[i];
                    }
                }
            }
        } // for ( iR=(maxiS-1); iR>=0; iR--)
        
        
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
                                        if(MIr%2 == 1)
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
                                            // phi[icomp*63+ITT] = (double)hatLFgamma[in]/(double)1024 ;    // phi[ITT,icomp]
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
							if (hatLFgamma2[in]>1023)
								hatLFgamma2[in] = 1023;
							if (hatLFgamma2[in] < 0)
								hatLFgamma2[in] = 0;
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
        
    } // for (iview=1; iview<nviews; iview++)
    
    //printf("did it [%d]\n",sizeof(long));
    
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
    if(PredRegr0 != NULL)
        free(PredRegr0);
    if(PredTheta0 != NULL)
        free(PredTheta0);
    if(hatLFgamma2 != NULL)
        free(hatLFgamma2);
    if(ycrt != NULL)
        free(ycrt);
    
    
}


void vpts(int Ms) {
    
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
    
    // read segmentation for all views
    FILE* f_SEGM = fopen("SEGM", "rb");
    fread(&nviews, sizeof(int), 1, f_SEGM);
    fread(&nr, sizeof(int), 1, f_SEGM);
    fread(&nc, sizeof(int), 1, f_SEGM);
    int* SEGgamma = (int*) malloc(nr*nc*nviews*sizeof(int));
    fread(SEGgamma, sizeof(int), nr*nc*nviews, f_SEGM);
    fclose(f_SEGM);
    
    int maxiS = 0;
    int counter = 0;
    for (int i=0; i <nr; i++) {
        for( int j = 0; j < nc; ++j ) {
            if( SEGgamma[0+i*nviews+j*nr*nviews] > maxiS ) {
                maxiS = SEGgamma[0+i*nviews+j*nr*nviews];
                

            }
            if( SEGgamma[0+i*nviews+j*nr*nviews] == 10 ) {
                counter++;
            }
        }
    }
    
    

    
    // get LFgamma
    FILE* f_VIN = fopen("VIN", "rb");
    fread(&nviews, sizeof(int), 1, f_VIN);
    fread(&nr, sizeof(int), 1, f_VIN);
    fread(&nc, sizeof(int), 1, f_VIN);
    fread(&ncols, sizeof(int), 1, f_VIN);
    int* LFgamma = (int*) malloc(nr*nc*nviews*ncols*sizeof(int));
    fread(LFgamma, sizeof(int), nviews*nr*nc*ncols, f_VIN);
    fclose(f_VIN);
    
    // get estimate for LFgamma
    FILE* f_DRV = fopen("DRV", "rb");
    fread(&nviews, sizeof(int), 1, f_DRV);
    fread(&nr, sizeof(int), 1, f_DRV);
    fread(&nc, sizeof(int), 1, f_DRV);
    fread(&ncols, sizeof(int), 1, f_DRV);
    int* hatLFgamma = (int*) malloc(nr*nc*nviews*ncols*sizeof(int));
    fread(hatLFgamma, sizeof(int), nr*nc*nviews*ncols, f_DRV);
    fclose(f_DRV);
    
    int* Pred_RegW = alocaVector((int)nviews*maxiS*63); // Pred_RegW(iview,iR,Ms)
    memset(Pred_RegW, 0, nviews*maxiS*63);
    double* Pred_ThetaW = alocaDoubleVector((int)nviews*maxiS*63);
    memset(Pred_ThetaW, 0, nviews*maxiS*63);
    int* IntPred_ThetaW = alocaVector((int)nviews*maxiS*63);
    memset(IntPred_ThetaW, 0, nviews*maxiS*63);
    
    int* Neigh9_165 = alocaVector((int)nviews*9);
    for(int i = 0; i < nviews*9; ++i ) {
        Neigh9_165[i] = Neigh9_165p_const[i];
    }

	time_t begin = clock();
	CC_WarpAllEncoder4Part1(SEGgamma, Neigh9_165, hatLFgamma, LFgamma, Pred_RegW, Pred_ThetaW, maxiS, Ms);
	time_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed time for WarpAllEncoder4Part1: %f\n", elapsed_secs);
  
    FILE* f_Pred_RegW = fopen ( "PRMA" , "wb" );
    fwrite( &nviews, sizeof(int), 1, f_Pred_RegW);
    fwrite( &maxiS ,sizeof(int),1,f_Pred_RegW);
    fwrite( &Ms ,sizeof(int),1,f_Pred_RegW);
    fwrite(Pred_RegW,sizeof(int),nviews*63*maxiS,f_Pred_RegW);
    fclose(f_Pred_RegW);
    
    
    for (int i=0; i<nviews*63*maxiS; i++) {
        IntPred_ThetaW[i] = (int) (Pred_ThetaW[i]* pow((double)2,12));
		if (abs(IntPred_ThetaW[i])>100)
			printf("[%i] %i\n", i, IntPred_ThetaW[i]);
    }
    FILE* f_Pred_ThetaW = fopen ( "PRCO" , "wb" );
    fwrite( &nviews, sizeof(int), 1, f_Pred_RegW);
    fwrite( &maxiS ,sizeof(int),1,f_Pred_ThetaW);
    fwrite( &Ms ,sizeof(int),1,f_Pred_ThetaW);
    fwrite( IntPred_ThetaW ,sizeof(int),(long)nviews*63*maxiS,f_Pred_ThetaW);
    fclose(f_Pred_ThetaW);
    
    
    f_SEGM = fopen("LFGAMMA2", "wb");
    fwrite( &nviews ,sizeof(int),1,f_SEGM);
    fwrite( &nr ,sizeof(int),1,f_SEGM);
    fwrite( &nc ,sizeof(int),1,f_SEGM);
    fwrite( &ncols, sizeof(int), 1, f_SEGM);
    fwrite( hatLFgamma ,sizeof(int), nr*nc*nviews*ncols,f_SEGM);
    fclose(f_SEGM);
    
    aux_write_pred_header(Ms, maxiS);
    
    free(hatLFgamma);
    free(LFgamma);
    free(Pred_RegW);
    free(Pred_ThetaW);

}
