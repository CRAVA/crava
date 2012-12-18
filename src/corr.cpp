#include <stdio.h>
#include <math.h>

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/seismicparametersholder.h"
#include "src/spatialwellfilter.h"
#include "src/modelsettings.h"
#include "src/definitions.h"
#include "src/fftfilegrid.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/tasklist.h"
#include "src/corr.h"
#include "src/io.h"

Corr::Corr(float  ** pointVar0,
           float  ** priorVar0,
           float   * priorCorrT,
           int       n,
           float     dt,
           Surface * priorCorrXY)
  : priorCorrTFiltered_(NULL),
    postVar0_(NULL),
    postCovAlpha00_(NULL),
    postCovBeta00_(NULL),
    postCovRho00_(NULL),
    postCrCovAlphaBeta00_(NULL),
    postCrCovAlphaRho00_(NULL),
    postCrCovBetaRho00_(NULL),
    postCovAlpha_(NULL),
    postCovBeta_(NULL),
    postCovRho_(NULL),
    postCrCovAlphaBeta_(NULL),
    postCrCovAlphaRho_(NULL),
    postCrCovBetaRho_(NULL)
{
  pointVar0_         = pointVar0;
  priorVar0_         = priorVar0;
  priorCorrT_        = priorCorrT;
  priorCorrXY_       = priorCorrXY;
  n_                 = n;
  dt_                = dt;
  common_correlation_ = true;

}

Corr::Corr(int                       n,
           float                     dt,
           SeismicParametersHolder & seismicParameters)
  : errCorr_(NULL)
{
  postCovAlpha_       = seismicParameters.GetCovAlpha();
  postCovBeta_        = seismicParameters.GetCovBeta();
  postCovRho_         = seismicParameters.GetCovRho();
  postCrCovAlphaBeta_ = seismicParameters.GetCrCovAlphaBeta();
  postCrCovAlphaRho_  = seismicParameters.GetCrCovAlphaRho();
  postCrCovBetaRho_   = seismicParameters.GetCrCovBetaRho();
  dt_                 = dt;
  common_correlation_  = false;

  if((n % 2) == 0)
    n_ = n/2+1;
  else
    n_ = n/2;

 createPriorVar0();
}


Corr::~Corr(void)
{
  if(common_correlation_ == true){
    if(pointVar0_!=NULL){
      for(int i=0; i<3; i++)
        delete [] pointVar0_[i];
      delete [] pointVar0_;
    }
    delete priorCorrXY_;
    delete [] priorCorrT_;

    if(priorCorrTFiltered_!=NULL)
      delete [] priorCorrTFiltered_;

  }
  else
    delete errCorr_;

  if(priorVar0_!=NULL){
    for(int i=0; i<3; i++)
      delete [] priorVar0_[i];
    delete [] priorVar0_;
  }

  if(postVar0_!=NULL){
    for(int i=0 ; i<3 ; i++)
      delete[] postVar0_[i];
    delete [] postVar0_;
  }

  if(postCovAlpha00_!=NULL)
    delete [] postCovAlpha00_ ;
  if(postCovBeta00_!=NULL)
    delete [] postCovBeta00_ ;
  if(postCovRho00_!=NULL)
    delete [] postCovRho00_ ;
  if(postCrCovAlphaBeta00_!=NULL)
    delete [] postCrCovAlphaBeta00_ ;
  if(postCrCovAlphaRho00_!=NULL)
    delete [] postCrCovAlphaRho00_ ;
  if(postCrCovBetaRho00_!=NULL)
    delete [] postCrCovBetaRho00_;

  //Covariance grids are deleted in seismicParametersHolder
}

//--------------------------------------------------------------------
void
Corr::createPostGrids(int nx,  int ny,  int nz,
                      int nxp, int nyp, int nzp,
                      bool fileGrid)
{
  if(common_correlation_ == true){
    postCovAlpha_       = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder
    postCovBeta_        = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder
    postCovRho_         = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder
    postCrCovAlphaBeta_ = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder
    postCrCovAlphaRho_  = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder
    postCrCovBetaRho_   = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);   //Deleted in seismicParametersHolder

    postCovAlpha_       ->setType(FFTGrid::COVARIANCE);
    postCovBeta_        ->setType(FFTGrid::COVARIANCE);
    postCovRho_         ->setType(FFTGrid::COVARIANCE);
    postCrCovAlphaBeta_ ->setType(FFTGrid::COVARIANCE);
    postCrCovAlphaRho_  ->setType(FFTGrid::COVARIANCE);
    postCrCovBetaRho_   ->setType(FFTGrid::COVARIANCE);

    postCovAlpha_       ->createRealGrid();     // First used in time domain
    postCovBeta_        ->createRealGrid();     // First used in time domain
    postCovRho_         ->createComplexGrid();  // First used in Fourier domain
    postCrCovAlphaBeta_ ->createComplexGrid();  // First used in Fourier domain
    postCrCovAlphaRho_  ->createComplexGrid();  // First used in Fourier domain
    postCrCovBetaRho_   ->createComplexGrid();  // First used in Fourier domain
  }
  else{
    errCorr_ = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
    errCorr_ ->setType(FFTGrid::COVARIANCE);
    errCorr_ ->createRealGrid();                // First used in time domain
  }
}
//--------------------------------------------------------------------
void
Corr::computeCircCorrT(fftw_real* CircCorrT, int nzp)
{
  int k,refk;

  assert(priorCorrT_[0] != 0);

  for(k = 0 ; k < 2*(nzp/2+1) ; k++ )
  {
    if(k < nzp)
    {
      if( k < nzp/2+1)
        refk = k;
      else
        refk = nzp - k;
      if(refk < n_)
        CircCorrT[k] = priorCorrT_[refk];
      else
        CircCorrT[k] = 0.0;
    }
    else
    {
      CircCorrT[k] = RMISSING;
    }
  }
}
//-------------------------------------------------------------------
fftw_real*
Corr::extractParamCorr(FFTGrid * covAlpha, int nzp)
{
  assert(covAlpha->getIsTransformed() == false);

  covAlpha->setAccessMode(FFTGrid::RANDOMACCESS);

  fftw_real * circCorrT = reinterpret_cast<fftw_real*>(fftw_malloc(2*(nzp/2+1)*sizeof(fftw_real)));
  int         refk;
  float       constant = covAlpha->getRealValue(0,0,0);

  for(int k = 0 ; k < 2*(nzp/2+1) ; k++ ){
    if(k < nzp){
      if( k < nzp/2+1)
        refk = k;
      else
        refk = nzp - k;
      if(refk < n_)
        circCorrT[k] = covAlpha->getRealValue(0,0,refk)/constant;
      else
        circCorrT[k] = 0.0;
    }
    else
      circCorrT[k]=RMISSING;
  }
  covAlpha->endAccess();

  return circCorrT;//fftw_free(circCorrT);
}
//------------------------------------------------------------------
void
Corr::makeCircCorrTPosDef(fftw_real* CircCorrT, int minIntFq, int nzp)
{
  int k;
  fftw_complex* fftCircCorrT;
  fftCircCorrT = FFTGrid::fft1DzInPlace(CircCorrT, nzp);
  for(k = 0; k < nzp/2+1; k++)
  {
    if(k <= minIntFq)
      fftCircCorrT[k].re = 0.0 ;
    else
      fftCircCorrT[k].re = float(sqrt(fftCircCorrT[k].re * fftCircCorrT[k].re +
                                      fftCircCorrT[k].im * fftCircCorrT[k].im ));
    fftCircCorrT[k].im = 0.0;
  }

  CircCorrT   = FFTGrid::invFFT1DzInPlace(fftCircCorrT, nzp);
  //
  // NBNB-PAL: If the number of layers is too small CircCorrT[0] = 0. How
  //           do we avoid this, or how do we flag the problem?
  //
  float scale;
  if (CircCorrT[0] > 1.e-5f) // NBNB-PAL: Temporary solution for above mentioned problem
    scale = float( 1.0/CircCorrT[0] );
  else
  {
    LogKit::LogFormatted(LogKit::Low,"\nERROR: The circular temporal correlation (CircCorrT) is undefined. You\n");
    LogKit::LogFormatted(LogKit::Low,"       probably need to increase the number of layers...\n\nAborting\n");
    exit(1);
  }
  for(k = 0; k < nzp; k++)
  {
    CircCorrT[k] *= scale;
  }
}
//-------------------------------------------------------------------
FFTGrid*
Corr::createFFTGrid(int nx,  int ny,  int nz,
                    int nxp, int nyp, int nzp,
                    bool fileGrid)
{
  FFTGrid * fftGrid;
  if(fileGrid)
    fftGrid = new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fftGrid = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);
  return(fftGrid);
}
//--------------------------------------------------------------------
void
Corr::setPriorVar0(float ** priorVar0)
{
  if(priorVar0_ != NULL) {
    for(int i=0 ; i<3 ; i++)
      delete [] priorVar0_[i];
    delete [] priorVar0_;
  }
  priorVar0_ = priorVar0;
}

//--------------------------------------------------------------------
void
Corr::setPriorCorrTFiltered(float * corrT, int nz, int nzp)
{
  // This is the cyclic and filtered version of CorrT which
  // has one or more zeros in the middle

  priorCorrTFiltered_ = new float[nzp];

  int refk;
  for(int i = 0 ; i < nzp ; i++ )
  {
    if( i < nzp/2+1)
      refk = i;
    else
      refk = nzp - i;

    if(refk < nz && corrT != NULL)
      priorCorrTFiltered_[i] = corrT[refk];
    else
      priorCorrTFiltered_[i] = 0.0;
  }
}

//--------------------------------------------------------------------
void
Corr::invFFT(void)
{
  LogKit::LogFormatted(LogKit::High,"\nBacktransforming correlation grids from FFT domain to time domain...");
  if (postCovAlpha_->getIsTransformed())
    postCovAlpha_->invFFTInPlace();
  if (postCovBeta_->getIsTransformed())
    postCovBeta_->invFFTInPlace();
  if (postCovRho_->getIsTransformed())
    postCovRho_->invFFTInPlace();
  if (postCrCovAlphaBeta_->getIsTransformed())
    postCrCovAlphaBeta_->invFFTInPlace();
  if (postCrCovAlphaRho_->getIsTransformed())
    postCrCovAlphaRho_->invFFTInPlace();
  if (postCrCovBetaRho_->getIsTransformed())
    postCrCovBetaRho_->invFFTInPlace();
  if (common_correlation_==false && errCorr_->getIsTransformed())
    errCorr_->invFFTInPlace();
  LogKit::LogFormatted(LogKit::High,"...done\n");
}

//--------------------------------------------------------------------
void
Corr::FFT(void)
{
  LogKit::LogFormatted(LogKit::High,"Transforming correlation grids from time domain to FFT domain...");
  if (!postCovAlpha_->getIsTransformed())
    postCovAlpha_->fftInPlace();
  if (!postCovBeta_->getIsTransformed())
    postCovBeta_->fftInPlace();
  if (!postCovRho_->getIsTransformed())
    postCovRho_->fftInPlace();
  if (!postCrCovAlphaBeta_->getIsTransformed())
    postCrCovAlphaBeta_->fftInPlace();
  if (!postCrCovAlphaRho_->getIsTransformed())
    postCrCovAlphaRho_->fftInPlace();
  if (!postCrCovBetaRho_->getIsTransformed())
    postCrCovBetaRho_->fftInPlace();
  if (common_correlation_==false && !errCorr_->getIsTransformed())
    errCorr_->fftInPlace();
  LogKit::LogFormatted(LogKit::High,"...done\n");
}

//--------------------------------------------------------------------
void Corr::printPriorVariances(void) const
{
  LogKit::LogFormatted(LogKit::Low,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Variances           ln Vp     ln Vs    ln Rho         \n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"Inversion grid:   %.1e   %.1e   %.1e (used by program)\n",priorVar0_[0][0],priorVar0_[1][1],priorVar0_[2][2]);
  if(pointVar0_ != NULL)
    LogKit::LogFormatted(LogKit::Low,"Well logs     :   %.1e   %.1e   %.1e                  \n",pointVar0_[0][0],pointVar0_[1][1],pointVar0_[2][2]);

  float corr01 = priorVar0_[0][1]/(sqrt(priorVar0_[0][0]*priorVar0_[1][1]));
  float corr02 = priorVar0_[0][2]/(sqrt(priorVar0_[0][0]*priorVar0_[2][2]));
  float corr12 = priorVar0_[1][2]/(sqrt(priorVar0_[1][1]*priorVar0_[2][2]));
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::Low,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::Low,"ln Rho |                     %5.2f \n",1.0f);

  if (std::abs(corr01) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Vs correlation is wrong (%.2f).\n",corr01);
    TaskList::addTask("Check your prior correlations. Corr(Vp,Vs) is out of bounds.");
  }
  if (std::abs(corr02) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Rho correlation is wrong (%.2f).\n",corr02);
    TaskList::addTask("Check your prior correlations. Corr(Vp,Rho) is out of bounds.");
  }
  if (std::abs(corr12) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vs-Rho correlation is wrong (%.2f).\n",corr12);
    TaskList::addTask("Check your prior correlations. Corr(Vs,Rho) is out of bounds.");
  }
}

//--------------------------------------------------------------------
void
Corr::printPostVariances(void) const
{
  LogKit::WriteHeader("Posterior Covariance");

  LogKit::LogFormatted(LogKit::Low,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"               ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-----------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"Variances:   %.1e   %.1e   %.1e    \n",postVar0_[0][0],postVar0_[1][1],postVar0_[2][2]);
  LogKit::LogFormatted(LogKit::Low,"\n");
  float corr01 = postVar0_[0][1]/(sqrt(postVar0_[0][0]*postVar0_[1][1]));
  float corr02 = postVar0_[0][2]/(sqrt(postVar0_[0][0]*postVar0_[2][2]));
  float corr12 = postVar0_[1][2]/(sqrt(postVar0_[1][1]*postVar0_[2][2]));
  LogKit::LogFormatted(LogKit::Low,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::Low,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::Low,"ln Rho |                     %5.2f \n",1.0f);

  if (std::abs(corr01) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Vs correlation is wrong (%.2f).\n",corr01);
    TaskList::addTask("Check your posterior correlations. Corr(Vp,Vs) is out of bounds.");
  }
  if (std::abs(corr02) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Rho correlation is wrong (%.2f).\n",corr02);
    TaskList::addTask("Check your posterior correlations. Corr(Vp,Rho) is out of bounds.");
  }
  if (std::abs(corr12) > 1.0) {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vs-Rho correlation is wrong (%.2f).\n",corr12);
    TaskList::addTask("Check your posterior correlations. Corr(Vs,Rho) is out of bounds.");
  }
}
//--------------------------------------------------------------------
void
Corr::getNextParameterCovariance(fftw_complex **& parVar)
{
  if(common_correlation_ == true) {
    fftw_complex ijkParLam;
    fftw_complex ijkTmp = postCovAlpha_->getNextComplex();

    ijkParLam.re = float ( sqrt(ijkTmp.re * ijkTmp.re));
    ijkParLam.im = 0.0;

    for(int l = 0; l < 3; l++ ){
      for(int m = 0; m < 3; m++ )
      {
        parVar[l][m].re = priorVar0_[l][m] * ijkParLam.re;
        parVar[l][m].im = 0.0;
        // if(l!=m)
        //   parVar[l][m].re *= 0.75;  //NBNB OK DEBUG TEST
      }
    }
  }
  else{
    fftw_complex ii;
    fftw_complex jj;
    fftw_complex kk;
    fftw_complex ij;
    fftw_complex ik;
    fftw_complex jk;

    fftw_complex iiTmp = postCovAlpha_      ->getNextComplex();
    fftw_complex jjTmp = postCovBeta_       ->getNextComplex();
    fftw_complex kkTmp = postCovRho_        ->getNextComplex();
    fftw_complex ijTmp = postCrCovAlphaBeta_->getNextComplex();
    fftw_complex ikTmp = postCrCovAlphaRho_ ->getNextComplex();
    fftw_complex jkTmp = postCrCovBetaRho_  ->getNextComplex();

    ii.re = float( sqrt(iiTmp.re * iiTmp.re));
    ii.im = 0.0;
    jj.re = float( sqrt(jjTmp.re * jjTmp.re));
    jj.im = 0.0;
    kk.re = float( sqrt(kkTmp.re * kkTmp.re));
    kk.im = 0.0;
    ij.re = float( sqrt(ijTmp.re * ijTmp.re));
    ij.im = 0.0;
    ik.re = float( sqrt(ikTmp.re * ikTmp.re));
    ik.im = 0.0;
    jk.re = float( sqrt(jkTmp.re * jkTmp.re));
    jk.im = 0.0;

    parVar[0][0].re = ii.re;
    parVar[0][0].im = ii.im;

    parVar[1][1].re = jj.re;
    parVar[1][1].im = jj.im;

    parVar[2][2].re = kk.re;
    parVar[2][2].im = kk.im;

    parVar[0][1].re = ij.re;
    parVar[0][1].im = ij.im;
    parVar[1][0].re = ij.re;
    parVar[1][0].im = -ij.im;

    parVar[0][2].re = ik.re;
    parVar[0][2].im = ik.im;
    parVar[2][0].re = ik.re;
    parVar[2][0].im = ik.im;

    parVar[1][2].re = jk.re;
    parVar[1][2].im = jk.im;
    parVar[2][1].re = jk.re;
    parVar[2][1].im = -jk.im;
  }
}
//--------------------------------------------------------------------
void
Corr::getNextErrorVariance(fftw_complex **& errVar,
                           fftw_complex * errMult1,
                           fftw_complex * errMult2,
                           fftw_complex * errMult3,
                           int ntheta,
                           float wnc,
                           double ** errThetaCov,
                           bool invert_frequency)
{
  fftw_complex ijkErrLam;
  fftw_complex ijkTmp;
  if(common_correlation_ == true)
    ijkTmp = postCovBeta_->getNextComplex();
  else
    ijkTmp = errCorr_->getNextComplex();

  ijkErrLam.re        = float( sqrt(ijkTmp.re * ijkTmp.re));
  ijkErrLam.im        = 0.0;


  if(invert_frequency) // inverting only relevant frequencies
  {
    for(int l = 0; l < ntheta; l++ ) {
      for(int m = 0; m < ntheta; m++ )
      {        // Note we multiply kWNorm[l] and comp.conj(kWNorm[m]) hence the + and not a minus as in pure multiplication
        errVar[l][m].re  = float( 0.5*(1.0-wnc)*errThetaCov[l][m] * ijkErrLam.re * ( errMult1[l].re *  errMult1[m].re +  errMult1[l].im *  errMult1[m].im));
        errVar[l][m].re += float( 0.5*(1.0-wnc)*errThetaCov[l][m] * ijkErrLam.re * ( errMult2[l].re *  errMult2[m].re +  errMult2[l].im *  errMult2[m].im));
        if(l==m) {
          errVar[l][m].re += float( wnc*errThetaCov[l][m] * errMult3[l].re  * errMult3[l].re);
          errVar[l][m].im   = 0.0;
        }
        else {
          errVar[l][m].im  = float( 0.5*(1.0-wnc)*errThetaCov[l][m] * ijkErrLam.re * (-errMult1[l].re * errMult1[m].im + errMult1[l].im * errMult1[m].re));
          errVar[l][m].im += float( 0.5*(1.0-wnc)*errThetaCov[l][m] * ijkErrLam.re * (-errMult2[l].re * errMult2[m].im + errMult2[l].im * errMult2[m].re));
        }
      }
    }
  }
}

void Corr::initializeCorrelationsSyntWells(SpatialWellFilter                    * spatwellfilter,
                                           std::vector<SyntWellData *>            wells,
                                           int                                    nWells,
                                           int                                    lowIntCut,
                                           float                                  corrGradI,
                                           float                                  corrGradJ)
{
  //fftw_real * corrT = NULL;

  if(common_correlation_ == true){

    postCovAlpha_->fillInParamCorr(this, lowIntCut, corrGradI, corrGradJ);
    //setPriorCorrTFiltered(corrT, nz, nzp); // Can have zeros in the middle

    if(spatwellfilter != NULL){
        postCovAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);

        for(int i=0; i<nWells; i++)
          spatwellfilter->setPriorSpatialCorrSyntWell(postCovAlpha_, wells[i], i);

        postCovAlpha_->endAccess();
    }
  }
  //delete corrT;
}

//--------------------------------------------------------------------
fftw_real*
Corr::initializeCorrelations(SpatialWellFilter   * spatwellfilter, std::vector<WellData *> wells, float corrGradI, float corrGradJ, int lowIntCut, int nWells, int nz, int nzp)
{

  fftw_real * corrT = NULL; // =  fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real));

  if(common_correlation_ == true){
    corrT = postCovAlpha_->fillInParamCorr(this, lowIntCut, corrGradI, corrGradJ);
    setPriorCorrTFiltered(corrT, nz, nzp); // Can have zeros in the middle

    if(spatwellfilter != NULL){
      postCovAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);

      for(int i=0; i<nWells; i++)
        spatwellfilter->setPriorSpatialCorr(postCovAlpha_, wells[i], i);

      postCovAlpha_->endAccess();
    }

    postCovBeta_->fillInErrCorr(this, corrGradI, corrGradJ);
  }
  else{

    corrT = extractParamCorr(postCovAlpha_, nzp);
    setPriorCorrTFiltered(corrT, nz, nzp);

    errCorr_->fillInErrCorrFromCovAlpha(postCovAlpha_, corrGradI, corrGradJ);
  }

  return(corrT);
}
//--------------------------------------------------------------------
void
Corr::initializeAccess(void)
{
  if(common_correlation_ == true){

    postCovAlpha_      ->fftInPlace();
    postCovBeta_       ->fftInPlace();

    postCovAlpha_      ->setAccessMode(FFTGrid::READANDWRITE);  //endAccess() in Crava::computePostMeanResidAndFFTCov()
    postCovBeta_       ->setAccessMode(FFTGrid::READANDWRITE);
    postCovRho_        ->setAccessMode(FFTGrid::WRITE);
    postCrCovAlphaBeta_->setAccessMode(FFTGrid::WRITE);
    postCrCovAlphaRho_ ->setAccessMode(FFTGrid::WRITE);
    postCrCovBetaRho_  ->setAccessMode(FFTGrid::WRITE);
  }
  else{

    FFT();

    errCorr_           ->setAccessMode(FFTGrid::READ);          //endAccess() in Corr::terminateAccess
    postCovAlpha_      ->setAccessMode(FFTGrid::READANDWRITE);  //endAccess() in Crava::computePostMeanResidAndFFTCov()
    postCovBeta_       ->setAccessMode(FFTGrid::READANDWRITE);
    postCovRho_        ->setAccessMode(FFTGrid::READANDWRITE);
    postCrCovAlphaBeta_->setAccessMode(FFTGrid::READANDWRITE);
    postCrCovAlphaRho_ ->setAccessMode(FFTGrid::READANDWRITE);
    postCrCovBetaRho_  ->setAccessMode(FFTGrid::READANDWRITE);
  }
}
//--------------------------------------------------------------------
void
Corr::terminateAccess(void)
{
  if(common_correlation_ == false)
    errCorr_->endAccess();
}
//--------------------------------------------------------------------
void
Corr::createPostVariances(void)
{
  postVar0_ = new float*[3];
  for(int i = 0 ; i < 3 ; i++)
    postVar0_[i] = new float[3];

  postVar0_[0][0] = getOrigin(postCovAlpha_);
  postVar0_[1][1] = getOrigin(postCovBeta_);
  postVar0_[2][2] = getOrigin(postCovRho_);
  postVar0_[0][1] = getOrigin(postCrCovAlphaBeta_);
  postVar0_[1][0] = postVar0_[0][1];
  postVar0_[2][0] = getOrigin(postCrCovAlphaRho_);
  postVar0_[0][2] = postVar0_[2][0];
  postVar0_[2][1] = getOrigin(postCrCovBetaRho_);
  postVar0_[1][2] = postVar0_[2][1];

  postCovAlpha00_       = createPostCov00(postCovAlpha_);
  postCovBeta00_        = createPostCov00(postCovBeta_);
  postCovRho00_         = createPostCov00(postCovRho_);
  postCrCovAlphaBeta00_ = createPostCov00(postCrCovAlphaBeta_);
  postCrCovAlphaRho00_  = createPostCov00(postCrCovAlphaRho_);
  postCrCovBetaRho00_   = createPostCov00(postCrCovBetaRho_);
}

//--------------------------------------------------------------------
float
Corr::getOrigin(FFTGrid * grid) const
{
  grid->setAccessMode(FFTGrid::RANDOMACCESS);
  float value = grid->getRealValue(0,0,0);
  grid->endAccess();
  return value;
}

//--------------------------------------------------------------------
float *
Corr::createPostCov00(FFTGrid * postCov)
{
  int nz = postCov->getNz();
  float * postCov00 = new float[nz];
  postCov->setAccessMode(FFTGrid::RANDOMACCESS);
  for(int k=0 ; k < nz ; k++)
  {
    postCov00[k] = postCov->getRealValue(0,0,k);
  }
  postCov->endAccess();
  return (postCov00);
}
//--------------------------------------------------------------------
void
Corr::createPriorVar0(void)
{
  postCovAlpha_      ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovBeta_       ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovRho_        ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaBeta_->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovBetaRho_  ->setAccessMode(FFTGrid::RANDOMACCESS);

  priorVar0_ = new float*[3];
  for(int i=0 ; i<3 ; i++)
    priorVar0_[i] = new float[3];

  priorVar0_[0][0] = postCovAlpha_      ->getRealValue(0,0,0);
  priorVar0_[1][1] = postCovBeta_       ->getRealValue(0,0,0);
  priorVar0_[2][2] = postCovRho_        ->getRealValue(0,0,0);
  priorVar0_[0][1] = postCrCovAlphaBeta_->getRealValue(0,0,0);
  priorVar0_[1][0] = postCrCovAlphaBeta_->getRealValue(0,0,0);
  priorVar0_[0][2] = postCrCovAlphaRho_ ->getRealValue(0,0,0);
  priorVar0_[2][0] = postCrCovAlphaRho_ ->getRealValue(0,0,0);
  priorVar0_[1][2] = postCrCovBetaRho_  ->getRealValue(0,0,0);
  priorVar0_[2][1] = postCrCovBetaRho_  ->getRealValue(0,0,0);

  postCovAlpha_      ->endAccess();
  postCovBeta_       ->endAccess();
  postCovRho_        ->endAccess();
  postCrCovAlphaBeta_->endAccess();
  postCrCovAlphaRho_ ->endAccess();
  postCrCovBetaRho_  ->endAccess();
}
//--------------------------------------------------------------------
void
Corr::writeFilePriorCorrT(float* corrT, int nzp) const
{
  // This is the cyclic and filtered version of CorrT
  std::string baseName = IO::PrefixPrior() + IO::FileTemporalCorr() + IO::SuffixGeneralData();
  std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed
       << std::right
       << std::setprecision(6)
       << dt_ << "\n";
  for(int i=0 ; i<nzp ; i++) {
    file << std::setw(9) << corrT[i] << "\n";
  }
  file.close();
}

//--------------------------------------------------------------------
void Corr::writeFilePriorVariances(const ModelSettings * modelSettings) const
{
  std::string baseName1 = IO::PrefixPrior() + IO::FileParameterCov() + IO::SuffixCrava();
  std::string baseName2 = IO::PrefixPrior() + IO::FileTemporalCorr() + IO::SuffixCrava();
  std::string baseName3 = IO::PrefixPrior() + IO::FileLateralCorr();
  std::string fileName1 = IO::makeFullFileName(IO::PathToCorrelations(), baseName1);
  std::string fileName2 = IO::makeFullFileName(IO::PathToCorrelations(), baseName2);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName1);
  file << std::fixed
       << std::right
       << std::setprecision(10);
  for(int i=0 ; i<3 ; i++) {
    for(int j=0 ; j<3 ; j++) {
      file << std::setw(13) << priorVar0_[i][j] << " ";
    }
    file << "\n";
  }
  file.close();

  NRLib::OpenWrite(file, fileName2);
  file << std::fixed
       << std::right
       << std::setprecision(8)
       << dt_ << "\n";
  for(int i=0 ; i<n_ ; i++) {
    file << std::setw(11) << priorCorrT_[i] << "\n";
  }
  file.close();

  IO::writeSurfaceToFile(*priorCorrXY_, baseName3, IO::PathToCorrelations(), modelSettings->getOutputGridFormat());
}

//--------------------------------------------------------------------
void
Corr::writeFilePostVariances(void) const
{
  std::string baseName = IO::PrefixPosterior() + IO::FileParameterCov() + IO::SuffixGeneralData();
  std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), baseName);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed;
  file << std::right;
  file << std::setprecision(6);
  for(int i=0 ; i<3 ; i++) {
    for(int j=0 ; j<3 ; j++) {
      file << std::setw(10) << postVar0_[i][j] << " ";
    }
    file << "\n";
  }
  file.close();

  int nz = postCovAlpha_->getNz();
  std::string baseName1 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vp" +IO::SuffixGeneralData();
  std::string baseName2 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vs" +IO::SuffixGeneralData();
  std::string baseName3 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Rho"+IO::SuffixGeneralData();
  writeFilePostCorrT(postCovAlpha00_, nz, IO::PathToCorrelations(), baseName1);
  writeFilePostCorrT(postCovBeta00_ , nz, IO::PathToCorrelations(), baseName2);
  writeFilePostCorrT(postCovRho00_  , nz, IO::PathToCorrelations(), baseName3);
}

//--------------------------------------------------------------------
void
Corr::writeFilePostCorrT(float             * postCov,
                         int                 nz,
                         const std::string & subDir,
                         const std::string & baseName) const
{
  std::string fileName = IO::makeFullFileName(subDir,baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed;
  file << std::setprecision(6);
  file << std::right;
   float c0 = 1.0f/postCov[0];
  for(int k=0 ; k < nz ; k++)
  {
    file << std::setw(9) << postCov[k]*c0 << "\n";
  }
  file.close();
}

//--------------------------------------------------------------------
void
Corr::writeFilePostCovGrids(Simbox const * simbox) const
{
  std::string fileName;
  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp";
  postCovAlpha_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovAlpha_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for Vp");
  postCovAlpha_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs";
  postCovBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovBeta_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for Vs");
  postCovBeta_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho";
  postCovRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for density");
  postCovRho_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs";
  postCrCovAlphaBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaBeta_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vp,Vs)");
  postCrCovAlphaBeta_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho";
  postCrCovAlphaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vp,density)");
  postCrCovAlphaRho_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho";
  postCrCovBetaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovBetaRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vs,density)");
  postCrCovBetaRho_ ->endAccess();
}
//--------------------------------------------------------------------
