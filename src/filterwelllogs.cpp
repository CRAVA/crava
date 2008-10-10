#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include "lib/lib_matr.h"
#include "lib/random.h"
#include "lib/kriging1d.h"
#include "src/welldata.h"
#include "src/filterwelllogs.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/model.h"
#include "nrlib/iotools/logkit.hpp"


FilterWellLogs::FilterWellLogs(FFTGrid *postCovAlpha, FFTGrid *postCovBeta, FFTGrid *postCovRho, FFTGrid *postCrCovAlphaBeta,
                               FFTGrid *postCrCovAlphaRho, FFTGrid *postCrCovBetaRho, int nzp, int nz, WellData ** wells, int nWells, 
                               float lowCut, float highCut, int relative, Simbox *simbox, const Simbox& osimbox, RandomGen *random, 
                               float *corrprior, float ** sigma0, bool faciesprob)
: origSimbox_(osimbox)
{

  faciesprob_ = faciesprob;
  fftw_real *postcova, *postcovb, *postcovr, *postcrab, *postcrar, *postcrbr;

    int rnzp = 2*(nzp/2+1);
    postcova = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    postcovb = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    postcovr = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    postcrab = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    postcrar = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    postcrbr = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    
    if(postCovAlpha->getIsTransformed()==true)
      postCovAlpha->invFFTInPlace();
    
    if(postCovBeta->getIsTransformed()==true)
      postCovBeta->invFFTInPlace();
    if(postCovRho->getIsTransformed()==true)
      postCovRho->invFFTInPlace();

    if(postCrCovAlphaBeta->getIsTransformed()==true)
      postCrCovAlphaBeta->invFFTInPlace();
    if(postCrCovAlphaRho->getIsTransformed()==true)
      postCrCovAlphaRho->invFFTInPlace();
    if(postCrCovBetaRho->getIsTransformed()==true)
      postCrCovBetaRho->invFFTInPlace();

    int i;
    for(i=0;i<nzp;i++)
    {
      int refk;
      if( i < nzp/2+1)
        refk = i;
      else
        refk = nzp - i;
      if(refk < nz)
      {
        postcova[i] = postCovAlpha->getRealValue(0,0,refk);
        postcovb[i] = postCovBeta->getRealValue(0,0,refk);
        postcovr[i] = postCovRho->getRealValue(0,0,refk);
        postcrab[i] = postCrCovAlphaBeta->getRealValue(0,0,refk);
        postcrar[i] = postCrCovAlphaRho->getRealValue(0,0,refk);
        postcrbr[i] = postCrCovBetaRho->getRealValue(0,0,refk);
      }
      else
      {
        postcova[i] = 0.0;
        postcovb[i] = 0.0;
        postcovr[i] = 0.0;
        postcrab[i] = 0.0;
        postcrar[i] = 0.0;
        postcrbr[i] = 0.0;
      }
    }
    
    
    corrprior_ = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
    for(i=0;i<nzp;i++)
      corrprior_[i] =corrprior[i];

    WellData** ppWellData = static_cast<WellData**>(wells);
    doFiltering(ppWellData,nWells,
                           postcova,postcovb,postcovr,
                           postcrab,postcrar,postcrbr, 
                           lowCut, highCut, relative, nz, nzp, simbox, random, sigma0);

    postCovAlpha->fftInPlace();
    postCovBeta->fftInPlace();
    postCovRho->fftInPlace();
    postCrCovAlphaBeta->fftInPlace();
    postCrCovAlphaRho->fftInPlace();
    postCrCovBetaRho->fftInPlace();

    fftw_free(postcova);
    fftw_free(postcovb);
    fftw_free(postcovr);
    fftw_free(postcrab);
    fftw_free(postcrar);
    fftw_free(postcrbr);

}

FilterWellLogs::~FilterWellLogs()
{
 fftw_free(corrprior_);
delete [] alphafiltered_;
  delete [] betafiltered_;
  delete [] rhofiltered_;
  delete [] alphablock_;
  delete [] betablock_;
  delete [] rhoblock_;
  delete [] facieslog_;

}

void FilterWellLogs::doFiltering(WellData **wells, int nWells,
                                fftw_real *postcova, fftw_real *postcovb, fftw_real *postcovr,
                                fftw_real *postcrab, fftw_real *postcrar, fftw_real *postcrbr, 
                                float lowCut, float highCut, int relative, int nz, int nzp, 
                                Simbox *simbox, RandomGen *random, float ** sigma0)
{
  LogKit::LogFormatted(LogKit::LOW,"\nFiltering well logs\n");
  LogKit::LogFormatted(LogKit::MEDIUM,"  Wells available:\n");
  for (int w = 0 ; w < nWells ; w++)
    if(wells[w]->getUseForFaciesProbabilities())
      LogKit::LogFormatted(LogKit::MEDIUM,"    %s\n",wells[w]->getWellname());
  
  float domega  = static_cast<float> (1000.0/(nzp*simbox->getdz()));  //dz in milliseconds
  int w, w1,i,j;
  ndata_ = nWells*nz;
  alphafiltered_ = new float[ndata_];
  betafiltered_ = new float[ndata_];
  rhofiltered_ = new float[ndata_];
  alphablock_ = new float[ndata_];
  betablock_ = new float[ndata_];
  rhoblock_ = new float[ndata_];
  facieslog_ = new int[ndata_];

  for(i=0;i<ndata_;i++)
  {
    alphafiltered_[i] = RMISSING;
    betafiltered_[i] = RMISSING;
    rhofiltered_[i] = RMISSING;
    alphablock_[i] = RMISSING;
    betablock_[i] = RMISSING;
    rhoblock_[i] = RMISSING;
    facieslog_[i] = IMISSING;
  }
  
 /* sigmaPost_[0][0]=postcova[0];
  sigmaPost_[1][1]=postcovb[0];
  sigmaPost_[2][2]=postcovr[0];
  sigmaPost_[0][1]=postcrab[0];
  sigmaPost_[1][0]=postcrab[0];
  sigmaPost_[0][2]=postcrar[0];
  sigmaPost_[2][0]=postcrar[0];
  sigmaPost_[1][2]=postcrbr[0];
  sigmaPost_[2][1]=postcrbr[0];*/


  // Do fourier transform of covariances
  rfftwnd_plan p1 = rfftwnd_create_plan(1, &nzp, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  int rnzp = 2*(nzp/2+1);
  fftw_complex *postcova_cAmp  = reinterpret_cast<fftw_complex*>(postcova);
  fftw_complex *postcovb_cAmp  = reinterpret_cast<fftw_complex*>(postcovb);
  fftw_complex *postcovr_cAmp  = reinterpret_cast<fftw_complex*>(postcovr);
  fftw_complex *postcovar_cAmp = reinterpret_cast<fftw_complex*>(postcrar);
  fftw_complex *postcovab_cAmp = reinterpret_cast<fftw_complex*>(postcrab);
  fftw_complex *postcovbr_cAmp = reinterpret_cast<fftw_complex*>(postcrbr);
  fftw_complex *corrprior_cAmp = reinterpret_cast<fftw_complex*>(corrprior_);

  rfftwnd_one_real_to_complex(p1,corrprior_, corrprior_cAmp);
  rfftwnd_one_real_to_complex(p1,postcova, postcova_cAmp);
  rfftwnd_one_real_to_complex(p1,postcovb, postcovb_cAmp);
  rfftwnd_one_real_to_complex(p1,postcovr, postcovr_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrab, postcovab_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrar, postcovar_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrbr, postcovbr_cAmp);

  fftw_real    * alpha_rAmp = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
  fftw_real    * beta_rAmp  = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
  fftw_real    * rho_rAmp   = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
 
  fftw_complex * alpha_cAmp = reinterpret_cast<fftw_complex*>(alpha_rAmp);
  fftw_complex * beta_cAmp  = reinterpret_cast<fftw_complex*>(beta_rAmp);
  fftw_complex * rho_cAmp   = reinterpret_cast<fftw_complex*>(rho_rAmp);

  int insideOrigSimbox;
  double x,y,z;

  float * vtAlpha = new float[nz];            // vt = vertical trend
  float * vtBeta  = new float[nz];
  float * vtRho   = new float[nz];

  float * vtAlphaBg = new float[nz];          // vt = vertical trend
  float * vtBetaBg  = new float[nz];
  float * vtRhoBg   = new float[nz];
  int   * vtFacies  = new int[nz];

  char * fileName     = new char[MAX_STRING];
  char * fullFileName;

  sprintf(fileName,"BW_filteredlogs.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file = fopen(fullFileName,"w");
  delete [] fullFileName;

  sprintf(fileName,"BW_originallogs.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file2 = fopen(fullFileName,"w");
  delete [] fullFileName;

  sprintf(fileName,"BW_background.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file3 = fopen(fullFileName,"w");
  delete [] fullFileName;
  delete [] fileName;

  for (w1 = 0 ; w1 < nWells ; w1++)
  {
    if(wells[w1]->getUseForFaciesProbabilities())
    { 
      BlockedLogs * bw = new BlockedLogs(wells[w1], simbox, random) ;
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bw->getVerticalTrend(bw->getAlphaBackgroundResolution(), vtAlphaBg); // might be missing data at end because of simbox 
      bw->getVerticalTrend(bw->getBetaBackgroundResolution(), vtBetaBg);
      bw->getVerticalTrend(bw->getRhoBackgroundResolution(), vtRhoBg);
      
      extrapolate(vtAlphaBg,nz);
      extrapolate(vtBetaBg, nz);
      extrapolate(vtRhoBg, nz);
      
      bw->getVerticalTrend(bw->getAlpha(),vtAlpha);
      bw->getVerticalTrend(bw->getBeta(),vtBeta);
      bw->getVerticalTrend(bw->getRho(),vtRho);
      if(faciesprob_==true)
        bw->getVerticalTrend(bw->getFacies(),vtFacies,random);

      const int * ipos = bw->getIpos();
      const int * jpos = bw->getJpos();
      const int * kpos = bw->getKpos();
      
      for(i=0;i<nz;i++)
      {
        simbox->getCoord(ipos[0],jpos[0],kpos[0]+i,x,y,z);
        insideOrigSimbox = origSimbox_.getIndex(x,y,z);
        if(vtAlpha[i]!=RMISSING && vtBeta[i]!=RMISSING && vtRho[i]!=RMISSING && insideOrigSimbox!=IMISSING && faciesprob_==true)
          facieslog_[i+w1*nz] = vtFacies[i];
        else
          facieslog_[i+w1*nz] = IMISSING;
      }
      
      float dz = static_cast<float> (simbox->getdz()*simbox->getRelThick(ipos[0],jpos[0]));
      Kriging1D::krigVector(vtAlpha, vtAlphaBg, nz, dz);
      Kriging1D::krigVector(vtBeta, vtBetaBg, nz, dz);
      Kriging1D::krigVector(vtRho, vtRhoBg, nz, dz);
      simbox->getCoord(ipos[0],jpos[0],kpos[0],x,y,z);

      for(i=0;i<nz;i++)
      {
        alpha_rAmp[i] = vtAlpha[i] - vtAlphaBg[i];
        beta_rAmp[i] = vtBeta[i] - vtBetaBg[i];
        rho_rAmp[i] = vtRho[i] - vtRhoBg[i];
        
        fprintf(file2,"%d %f %f %f %f %d\n",w1, z+i*dz,vtAlpha[i], vtBeta[i], vtRho[i],facieslog_[i+w1*nz]);
        fprintf(file3,"%d %f %f %f %f %d\n",w1, z+i*dz,vtAlphaBg[i], vtBetaBg[i], vtRhoBg[i],facieslog_[i+w1*nz] );
      }
      
      j = nz-1;
      float diffa = (alpha_rAmp[0]-alpha_rAmp[j])/(nzp-j-1);
      float diffb = (beta_rAmp[0]-beta_rAmp[j])/(nzp-j-1);
      float diffr = (rho_rAmp[0]-rho_rAmp[j])/(nzp-j-1);
      //  for(i=sBW[w1]->GetNumberOfObs();i<nzp_;i++) // do padding
      for(i=nz;i<nzp;i++)
      {
        alpha_rAmp[i] = alpha_rAmp[j]+(i-j)*diffa;
        beta_rAmp[i] = beta_rAmp[j]+(i-j)*diffb;
        rho_rAmp[i] = rho_rAmp[j]+(i-j)*diffr;
      }
      // do fourier transform of alphadiff,betadiff,rhodiff
      
      rfftwnd_one_real_to_complex(p1, alpha_rAmp, alpha_cAmp);
      rfftwnd_one_real_to_complex(p1, beta_rAmp, beta_cAmp);
      rfftwnd_one_real_to_complex(p1, rho_rAmp, rho_cAmp);
      
      //loop over frequencies, create 3*3 matrices, and the filter
      fftw_complex **sigmaK = new fftw_complex*[3];
      fftw_complex **sigmaE = new fftw_complex*[3];
      fftw_complex *paramvec = new fftw_complex[3];
      fftw_complex *help = new fftw_complex[3];
      int ok;
      double **F = new double*[3];
      
      for(i=0;i<3;i++){
        sigmaK[i] = new fftw_complex[3];
        sigmaE[i] = new fftw_complex[3];
        F[i] =new double[3];
      }

      for(w=0;w<nzp/2+1;w++)
      {
        //  if(corrprior_cAmp[w].re<delta)
        //    corrprior_cAmp[w].re = 0.0;
        sigmaK[0][0].re = corrprior_cAmp[w].re*sigma0[0][0];
        //sigmaK[0][0].im = corrprior_cAmp[w].im*sigma0_[0][0];
        sigmaK[0][0].im = 0.0;
        sigmaK[1][1].re = corrprior_cAmp[w].re*sigma0[1][1];
        //sigmaK[1][1].im = corrprior_cAmp[w].im*sigma0_[1][1];
        sigmaK[1][1].im = 0.0;
        sigmaK[2][2].re = corrprior_cAmp[w].re*sigma0[2][2];
        //sigmaK[2][2].im = corrprior_cAmp[w].im*sigma0_[2][2];
        sigmaK[2][2].im = 0.0;
        sigmaK[1][0].re = corrprior_cAmp[w].re*sigma0[1][0];
        //sigmaK[1][0].im = corrprior_cAmp[w].im*sigma0_[1][0];
        sigmaK[1][0].im = 0.0;
        sigmaK[0][1].re = corrprior_cAmp[w].re*sigma0[1][0];
        //sigmaK[0][1].im = -corrprior_cAmp[w].im*sigma0_[1][0];
        sigmaK[0][1].im = 0.0;
        sigmaK[2][0].re = corrprior_cAmp[w].re*sigma0[2][0];
        //sigmaK[2][0].im = corrprior_cAmp[w].im*sigma0_[2][0];
        sigmaK[2][0].im = 0.0;
        sigmaK[0][2].re = corrprior_cAmp[w].re*sigma0[2][0];
        //sigmaK[0][2].im = -corrprior_cAmp[w].im*sigma0_[2][0];
        sigmaK[0][2].im = 0.0;
        sigmaK[1][2].re = corrprior_cAmp[w].re*sigma0[1][2];
        //sigmaK[1][2].im = corrprior_cAmp[w].im*sigma0_[1][2];
        sigmaK[1][2].im = 0.0;
        sigmaK[2][1].re = corrprior_cAmp[w].re*sigma0[1][2];
        //sigmaK[2][1].im = -corrprior_cAmp[w].im*sigma0_[1][2];
        sigmaK[2][1].im = 0.0;

        sigmaE[0][0].re = sigmaK[0][0].re - postcova_cAmp[w].re;
        sigmaE[0][0].im = sigmaK[0][0].im - postcova_cAmp[w].im;
        sigmaE[1][1].re = sigmaK[1][1].re - postcovb_cAmp[w].re;
        sigmaE[1][1].im = sigmaK[1][1].im - postcovb_cAmp[w].im;
        sigmaE[2][2].re = sigmaK[2][2].re - postcovr_cAmp[w].re;
        sigmaE[2][2].im = sigmaK[2][2].im - postcovr_cAmp[w].im;
        sigmaE[0][1].re = sigmaK[0][1].re - postcovab_cAmp[w].re;
        sigmaE[0][1].im = sigmaK[0][1].im - postcovab_cAmp[w].im;
        sigmaE[1][0].re = sigmaE[0][1].re;
        sigmaE[1][0].im = -sigmaE[0][1].im;
        sigmaE[0][2].re = sigmaK[0][2].re - postcovar_cAmp[w].re;
        sigmaE[0][2].im = sigmaK[0][2].im - postcovar_cAmp[w].im;
        sigmaE[2][0].re = sigmaE[0][2].re;
        sigmaE[2][0].im = -sigmaE[0][2].im;
        sigmaE[1][2].re = sigmaK[1][2].re - postcovbr_cAmp[w].re;
        sigmaE[1][2].im = sigmaK[1][2].im - postcovbr_cAmp[w].im;
        sigmaE[2][1].re = sigmaE[1][2].re;
        sigmaE[2][1].im = -sigmaE[1][2].im;
        
       
        
        //Do cholesky of sigmaK
        ok = lib_matrCholCpx(3,sigmaK);
        if(ok==0 && (w*domega > lowCut && w*domega < highCut) )
        {
          calcFilter(sigmaK, sigmaE, F);
          //apply filter on alpha, beta, rho
          paramvec[0] = alpha_cAmp[w];
          paramvec[1] = beta_cAmp[w];
          paramvec[2] = rho_cAmp[w];
          lib_matrProdMatRVecCpx(F, paramvec, 3,3, help);
          alpha_cAmp[w] = help[0];
          beta_cAmp[w] = help[1];
          rho_cAmp[w] = help[2];
        }
        else
        {
          alpha_cAmp[w].re = 0.0;
          beta_cAmp[w].re = 0.0;
          rho_cAmp[w].re = 0.0;
          alpha_cAmp[w].im = 0.0;
          beta_cAmp[w].im = 0.0;
          rho_cAmp[w].im = 0.0;
        }
      }   
      
      // do inverse fourier transform, and add background model
      rfftwnd_plan p2 = rfftwnd_create_plan(1, &nzp, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_complex_to_real(p2, alpha_cAmp, alpha_rAmp);
      rfftwnd_one_complex_to_real(p2, beta_cAmp, beta_rAmp);
      rfftwnd_one_complex_to_real(p2, rho_cAmp, rho_rAmp);
      // Add background model to filtered data if not relative method is used.
      
      simbox->getCoord(ipos[0],jpos[0],kpos[0],x,y,z);
      
      for(i=0;i<nz;i++)
      {
        alphafiltered_[i+w1*nz] = alpha_rAmp[i]/nzp+(1-relative)*vtAlphaBg[i];
        betafiltered_[i+w1*nz] = beta_rAmp[i]/nzp+(1-relative)*vtBetaBg[i];
        rhofiltered_[i+w1*nz] = rho_rAmp[i]/nzp+(1-relative)*vtRhoBg[i];
        alphablock_[i+w1*nz] = vtAlpha[i];
        betablock_[i+w1*nz] = vtBeta[i];
        rhoblock_[i+w1*nz] = vtRho[i];
        fprintf(file,"%d %f %f %f %f %d \n",w1, z+dz*i,
                alphafiltered_[i+w1*nz],
                betafiltered_[i+w1*nz],
                rhofiltered_[i+w1*nz],
                facieslog_[i+w1*nz]);   
      }
      fftwnd_destroy_plan(p2);


      for(i=0;i<3;i++)
      {
        delete [] sigmaK[i];
        delete [] sigmaE[i];
        delete [] F[i];
      }
      delete [] sigmaK;
      delete [] sigmaE;
      delete [] F;
      delete [] help;
      delete [] paramvec;
      delete bw;
    } // end if
  } //end for wells
  // Calculate min and max
  fclose(file);
  fclose(file2);
  fclose(file3);
  fftwnd_destroy_plan(p1);

  delete [] vtAlpha;            // vt = vertical trend
  delete [] vtBeta;
  delete [] vtRho;
  delete [] vtFacies;

  delete [] vtAlphaBg;          // vt = vertical trend
  delete [] vtBetaBg;
  delete [] vtRhoBg;

  fftw_free(alpha_rAmp);
  fftw_free(beta_rAmp);
  fftw_free(rho_rAmp);
}

void
FilterWellLogs::extrapolate(float * log,
                        int     nz) 
{
  //
  // Extrapolate log[] in both ends if needed
  //  
  int i=0;
  while (i<nz && log[i]==RMISSING)
    i++;
  if (i < nz - 1) 
  { 
    int first_nonmissing = i;
    i = nz - 1;
    while (i>0 && log[i]==RMISSING)
      i--;
    int last_nonmissing = i;
    
    for(int i=0 ; i < first_nonmissing ; i++) { 
      log[i] = log[first_nonmissing]; 
    }
    for(int i=last_nonmissing + 1 ; i < nz ; i++) { 
      log[i] = log[last_nonmissing]; 
    }
   
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"WARNING: All ... trend values are missing.\n");
  }
}

void
FilterWellLogs::calcFilter(fftw_complex **sigmaK, fftw_complex **sigmaE, double **F)
{
  int i,j;
  int     * ok1 = new int[1];
  double  * ev  = new double[3];
  double ** M   = new double*[3];
  double ** V   = new double*[3];
  double ** EV  = new double*[3];
  double ** sigmaKreal = new double*[3];
  fftw_complex ** help = new fftw_complex*[3];

  for(i=0;i<3;i++){
    help[i] = new fftw_complex[3];
    M[i] = new double[3];
    V[i] = new double[3];
    EV[i] =new double[3];
    sigmaKreal[i] = new double[3];
  }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmaKreal[i][j] = sigmaK[i][j].re;

  lib_matrLXeqBMatCpx(3, sigmaK, sigmaE, 3);
  lib_matrAdjoint(sigmaE,3,3,help);
  lib_matrLXeqBMatCpx(3, sigmaK, help, 3);
  lib_matrAdjoint(help,3,3,sigmaE); //sigmaE = M
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      M[i][j] = sigmaE[i][j].re;

  lib_matr_eigen(M,3,V,ev,ok1);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      if(i==j && ev[i]>0.0)
        EV[i][j] = sqrt(ev[i]);
      else
        EV[i][j] = 0.0;
    }
  lib_matr_prod(V,EV,3,3,3,M);  // M=V*EV
  lib_matrTranspose(V,3,3,EV);  // EV = V^T
  lib_matr_prod(M,EV,3,3,3,V);  // V = A^K=VD^0.5V^T
  //F = L*V*L^-1
    
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      EV[i][j] = sigmaKreal[i][j];
  EV[0][2] = 0.0;
  EV[0][1] = 0.0;
  EV[1][2] = 0.0;
  lib_matr_prod(EV,V,3,3,3,M); //M=LA
  lib_matrTranspose(M,3,3,V);
  lib_matrLtXeqBR(3, sigmaKreal, V, 3);  //L^T*F^T = V
  lib_matrTranspose(V,3,3,F);     
  
  delete [] ok1;
  delete [] ev;
  for (int i=0 ; i < 3 ; i++)
  {
    delete [] M[i];
    delete [] V[i];
    delete [] EV[i];
    delete [] sigmaKreal[i];
    delete [] help[i];

  }
  delete [] M;
  delete [] V;
  delete [] EV;
  delete [] sigmaKreal;
  delete [] help;
}

