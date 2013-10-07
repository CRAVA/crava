#include <stdio.h>

#include "nrlib/flens/nrlib_flens.hpp"

#include "src/seismicparametersholder.h"
#include "src/fftfilegrid.h"
#include "src/fftgrid.h"
#include "src/modelgeneral.h"
#include "src/tasklist.h"
#include "src/modelsettings.h"

SeismicParametersHolder::SeismicParametersHolder(void)
{
  muAlpha_        = NULL;
  muBeta_         = NULL;
  muRho_          = NULL;
  covAlpha_       = NULL;
  covBeta_        = NULL;
  covRho_         = NULL;
  crCovAlphaBeta_ = NULL;
  crCovAlphaRho_  = NULL;
  crCovBetaRho_   = NULL;

  priorVar0_.resize(3,3);
}

//--------------------------------------------------------------------
SeismicParametersHolder::~SeismicParametersHolder(void)
{
  if(covAlpha_!=NULL)
    delete covAlpha_;

  if(covBeta_!=NULL)
    delete covBeta_;

  if(covRho_!=NULL)
    delete covRho_;

  if(crCovAlphaBeta_!=NULL)
    delete crCovAlphaBeta_ ;

  if(crCovAlphaRho_!=NULL)
    delete crCovAlphaRho_ ;

  if(crCovBetaRho_!=NULL)
    delete crCovBetaRho_;

  if(muAlpha_!=NULL)
    delete muAlpha_;

  if(muBeta_!=NULL)
    delete muBeta_;

  if(muRho_!=NULL)
    delete muRho_;

}
//--------------------------------------------------------------------

void
SeismicParametersHolder::setBackgroundParameters(FFTGrid  * muAlpha,
                                                 FFTGrid  * muBeta,
                                                 FFTGrid  * muRho)
{
  muAlpha_   = muAlpha;
  muBeta_    = muBeta;
  muRho_     = muRho;
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::copyBackgroundParameters(FFTGrid  * muAlpha,
                                                  FFTGrid  * muBeta,
                                                  FFTGrid  * muRho)
{
  if(muAlpha_!=NULL)
    delete muAlpha_;

  if(muBeta_!=NULL)
    delete muBeta_;

  if(muRho_!=NULL)
    delete muRho_;

  muAlpha_= new FFTGrid(muAlpha);
  muBeta_ = new FFTGrid(muBeta);
  muRho_  = new FFTGrid(muRho);
}



void
SeismicParametersHolder::setCorrelationParameters(float                    ** priorVar0,
                                                  const std::vector<float>  & priorCorrT,
                                                  Surface                   * priorCorrXY,
                                                  const int                 & minIntFq,
                                                  const float               & corrGradI,
                                                  const float               & corrGradJ,
                                                  const int                 & nx,
                                                  const int                 & ny,
                                                  const int                 & nz,
                                                  const int                 & nxPad,
                                                  const int                 & nyPad,
                                                  const int                 & nzPad)
{
  priorVar0_.resize(3,3);

  priorVar0_(0,0) = static_cast<double>(priorVar0[0][0]);
  priorVar0_(1,0) = static_cast<double>(priorVar0[1][0]);
  priorVar0_(2,0) = static_cast<double>(priorVar0[2][0]);
  priorVar0_(0,1) = static_cast<double>(priorVar0[0][1]);
  priorVar0_(1,1) = static_cast<double>(priorVar0[1][1]);
  priorVar0_(2,1) = static_cast<double>(priorVar0[2][1]);
  priorVar0_(0,2) = static_cast<double>(priorVar0[0][2]);
  priorVar0_(1,2) = static_cast<double>(priorVar0[1][2]);
  priorVar0_(2,2) = static_cast<double>(priorVar0[2][2]);


  // check if covariance is well conditioned and robustify
  NRLib::Vector eVals(3);
  NRLib::Matrix eVec(3,3);
  NRLib::Matrix tmp(3,3);
  tmp=priorVar0_;
  NRLib::ComputeEigenVectors(tmp,eVals,eVec);

  double maxVal = eVals(0);

  for(int i=1;i<3;i++)
    if(maxVal<eVals(i))
      maxVal=eVals(i);

  for(int k=0;k<3;k++){
    if(eVals(k)<maxVal*0.001){
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
          priorVar0_(i,j)+= eVec(k,i)*eVec(k,j)*(maxVal*0.0011-eVals(k));
    }
  }
  //tmp=priorVar0_;
  //NRLib::ComputeEigenVectors(tmp,eVals,eVec);

  createCorrGrids(nx, ny, nz, nxPad, nyPad, nzPad, false);

  initializeCorrelations(priorCorrXY,
                         priorCorrT,
                         corrGradI,
                         corrGradJ,
                         minIntFq,
                         nzPad);
}
//--------------------------------------------------------------------

void
SeismicParametersHolder::createCorrGrids(int  nx,
                                         int  ny,
                                         int  nz,
                                         int  nxp,
                                         int  nyp,
                                         int  nzp,
                                         bool fileGrid)
{
  covAlpha_       = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
  covBeta_        = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
  covRho_         = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
  crCovAlphaBeta_ = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
  crCovAlphaRho_  = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);
  crCovBetaRho_   = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);

  covAlpha_       ->setType(FFTGrid::COVARIANCE);
  covBeta_        ->setType(FFTGrid::COVARIANCE);
  covRho_         ->setType(FFTGrid::COVARIANCE);
  crCovAlphaBeta_ ->setType(FFTGrid::COVARIANCE);
  crCovAlphaRho_  ->setType(FFTGrid::COVARIANCE);
  crCovBetaRho_   ->setType(FFTGrid::COVARIANCE);

  covAlpha_       ->createRealGrid();
  covBeta_        ->createRealGrid();
  covRho_         ->createRealGrid();
  crCovAlphaBeta_ ->createRealGrid();
  crCovAlphaRho_  ->createRealGrid();
  crCovBetaRho_   ->createRealGrid();
}
//-------------------------------------------------------------------
FFTGrid *
SeismicParametersHolder::createFFTGrid(int nx,  int ny,  int nz,
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
//-------------------------------------------------------------------
void
SeismicParametersHolder::initializeCorrelations(const Surface            * priorCorrXY,
                                                const std::vector<float> & priorCorrT,
                                                const float              & corrGradI,
                                                const float              & corrGradJ,
                                                const int                & lowIntCut,
                                                const int                & nzp)
{
  fftw_real * circCorrT = computeCircCorrT(priorCorrT, lowIntCut, nzp);

  covAlpha_      ->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);
  covBeta_       ->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);
  covRho_        ->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);
  crCovAlphaBeta_->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);
  crCovAlphaRho_ ->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);
  crCovBetaRho_  ->fillInParamCorr(priorCorrXY, circCorrT, corrGradI, corrGradJ);

  covAlpha_      ->multiplyByScalar(static_cast<float>(priorVar0_(0,0)));
  covBeta_       ->multiplyByScalar(static_cast<float>(priorVar0_(1,1)));
  covRho_        ->multiplyByScalar(static_cast<float>(priorVar0_(2,2)));
  crCovAlphaBeta_->multiplyByScalar(static_cast<float>(priorVar0_(0,1)));
  crCovAlphaRho_ ->multiplyByScalar(static_cast<float>(priorVar0_(0,2)));
  crCovBetaRho_  ->multiplyByScalar(static_cast<float>(priorVar0_(1,2)));

  fftw_free(circCorrT);
}

//--------------------------------------------------------------------
NRLib::Matrix
SeismicParametersHolder::getPriorVar0(void) const
{
  return priorVar0_;
}

//--------------------------------------------------------------------

void
SeismicParametersHolder::allocateGrids(const int nx, const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad)
{
  createCorrGrids(nx, ny, nz, nxPad, nyPad, nzPad, false);

  muAlpha_ = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  muBeta_  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  muRho_   = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);

  muAlpha_->setType(FFTGrid::PARAMETER);
  muBeta_ ->setType(FFTGrid::PARAMETER);
  muRho_  ->setType(FFTGrid::PARAMETER);

  muAlpha_->createGrid();
  muBeta_ ->createGrid();
  muRho_  ->createGrid();

  muAlpha_->fillInConstant(0.0);
  muBeta_ ->fillInConstant(0.0);
  muRho_  ->fillInConstant(0.0);

  covAlpha_->fillInConstant(0.0);
  covBeta_ ->fillInConstant(0.0);
  covRho_  ->fillInConstant(0.0);

  crCovAlphaBeta_->fillInConstant(0.0);
  crCovAlphaRho_ ->fillInConstant(0.0);
  crCovBetaRho_  ->fillInConstant(0.0);
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::invFFTAllGrids()
{
  LogKit::LogFormatted(LogKit::High,"\nBacktransforming background grids from FFT domain to time domain...");

  if(muAlpha_->getIsTransformed())
    muAlpha_->invFFTInPlace();

  if(muBeta_->getIsTransformed())
    muBeta_->invFFTInPlace();

  if(muRho_->getIsTransformed())
    muRho_->invFFTInPlace();
  LogKit::LogFormatted(LogKit::High,"...done\n");

  invFFTCovGrids();
}

//--------------------------------------------------------------------
void
SeismicParametersHolder::FFTAllGrids()
{
  LogKit::LogFormatted(LogKit::High,"\nTransforming background grids from time domain to FFT domain ...");

  if(!muAlpha_->getIsTransformed())
    muAlpha_->fftInPlace();

  if(!muBeta_->getIsTransformed())
    muBeta_->fftInPlace();

  if(!muRho_->getIsTransformed())
    muRho_->fftInPlace();
  LogKit::LogFormatted(LogKit::High,"...done\n");

  FFTCovGrids();


}
//-----------------------------------------------------------------------------------------

void
SeismicParametersHolder::invFFTCovGrids()
{
  LogKit::LogFormatted(LogKit::High,"\nBacktransforming correlation grids from FFT domain to time domain...");

  if (covAlpha_->getIsTransformed())
    covAlpha_->invFFTInPlace();

  if (covBeta_->getIsTransformed())
    covBeta_->invFFTInPlace();

  if (covRho_->getIsTransformed())
    covRho_->invFFTInPlace();

  if (crCovAlphaBeta_->getIsTransformed())
    crCovAlphaBeta_->invFFTInPlace();

  if (crCovAlphaRho_->getIsTransformed())
    crCovAlphaRho_->invFFTInPlace();

  if (crCovBetaRho_->getIsTransformed())
    crCovBetaRho_->invFFTInPlace();

  LogKit::LogFormatted(LogKit::High,"...done\n");
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::FFTCovGrids()
{
  LogKit::LogFormatted(LogKit::High,"Transforming correlation grids from time domain to FFT domain...");

  if (!covAlpha_->getIsTransformed())
    covAlpha_->fftInPlace();

  if (!covBeta_->getIsTransformed())
    covBeta_->fftInPlace();

  if (!covRho_->getIsTransformed())
    covRho_->fftInPlace();

  if (!crCovAlphaBeta_->getIsTransformed())
    crCovAlphaBeta_->fftInPlace();

  if (!crCovAlphaRho_->getIsTransformed())
    crCovAlphaRho_->fftInPlace();

  if (!crCovBetaRho_->getIsTransformed())
    crCovBetaRho_->fftInPlace();

  LogKit::LogFormatted(LogKit::High,"...done\n");
}
//--------------------------------------------------------------------

void
SeismicParametersHolder::getNextParameterCovariance(fftw_complex **& parVar) const
{

  fftw_complex ii;
  fftw_complex jj;
  fftw_complex kk;
  fftw_complex ij;
  fftw_complex ik;
  fftw_complex jk;

  fftw_complex iiTmp = covAlpha_      ->getNextComplex();
  fftw_complex jjTmp = covBeta_       ->getNextComplex();
  fftw_complex kkTmp = covRho_        ->getNextComplex();
  fftw_complex ijTmp = crCovAlphaBeta_->getNextComplex();
  fftw_complex ikTmp = crCovAlphaRho_ ->getNextComplex();
  fftw_complex jkTmp = crCovBetaRho_  ->getNextComplex();

  if(priorVar0_(0,0) != 0)
    iiTmp.re = iiTmp.re / static_cast<float>(priorVar0_(0,0));

  if(priorVar0_(1,1) != 0)
    jjTmp.re = jjTmp.re / static_cast<float>(priorVar0_(1,1));

  if(priorVar0_(2,2) != 0)
    kkTmp.re = kkTmp.re / static_cast<float>(priorVar0_(2,2));

  if(priorVar0_(0,1) != 0)
    ijTmp.re = ijTmp.re / static_cast<float>(priorVar0_(0,1));

  if(priorVar0_(0,2) != 0)
    ikTmp.re = ikTmp.re / static_cast<float>(priorVar0_(0,2));

  if(priorVar0_(1,2) != 0)
    jkTmp.re = jkTmp.re / static_cast<float>(priorVar0_(1,2));

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

  parVar[0][0].re = ii.re * static_cast<float>(priorVar0_(0,0));
  parVar[0][0].im = ii.im;

  parVar[1][1].re = jj.re * static_cast<float>(priorVar0_(1,1));
  parVar[1][1].im = jj.im;

  parVar[2][2].re = kk.re * static_cast<float>(priorVar0_(2,2));
  parVar[2][2].im = kk.im;

  parVar[0][1].re = ij.re * static_cast<float>(priorVar0_(0,1));
  parVar[0][1].im = ij.im;
  parVar[1][0].re = ij.re * static_cast<float>(priorVar0_(0,1));
  parVar[1][0].im = -ij.im;

  parVar[0][2].re = ik.re * static_cast<float>(priorVar0_(0,2));
  parVar[0][2].im = ik.im;
  parVar[2][0].re = ik.re * static_cast<float>(priorVar0_(0,2));
  parVar[2][0].im = -ik.im;

  parVar[1][2].re = jk.re * static_cast<float>(priorVar0_(1,2));
  parVar[1][2].im = jk.im;
  parVar[2][1].re = jk.re * static_cast<float>(priorVar0_(1,2));
  parVar[2][1].im = -jk.im;
}

//--------------------------------------------------------------------
fftw_real *
SeismicParametersHolder::computeCircCorrT(const std::vector<float> & priorCorrT,
                                          const int                & minIntFq,
                                          const int                & nzp) const
{
  assert(priorCorrT[0] != 0);

  int n = static_cast<int>(priorCorrT.size());

  fftw_real * circCorrT = reinterpret_cast<fftw_real*>(fftw_malloc(2*(nzp/2+1)*sizeof(fftw_real)));

  for(int k = 0; k < 2*(nzp/2+1); k++ ) {
    if(k < nzp) {
      int refk;

      if(k < nzp/2+1)
        refk = k;
      else
        refk = nzp - k;

      if(refk < n)
        circCorrT[k] = priorCorrT[refk];
      else
        circCorrT[k] = 0.0;
    }
    else
      circCorrT[k] = RMISSING;
  }

  makeCircCorrTPosDef(circCorrT, minIntFq, nzp);

  return circCorrT;

}
//--------------------------------------------------------------------

void
SeismicParametersHolder::makeCircCorrTPosDef(fftw_real * circCorrT,
                                             const int & minIntFq,
                                             const int & nzp) const
{
  fftw_complex * fftCircCorrT;
  fftCircCorrT = FFTGrid::fft1DzInPlace(circCorrT, nzp);

  for(int k=0; k<nzp/2+1; k++) {
    if(k <= minIntFq)
      fftCircCorrT[k].re = 0.0 ;
    else
      fftCircCorrT[k].re = float(sqrt(fftCircCorrT[k].re * fftCircCorrT[k].re +
                                      fftCircCorrT[k].im * fftCircCorrT[k].im ));
    fftCircCorrT[k].im = 0.0;
  }

  circCorrT   = FFTGrid::invFFT1DzInPlace(fftCircCorrT, nzp);
  //
  // NBNB-PAL: If the number of layers is too small CircCorrT[0] = 0. How
  //           do we avoid this, or how do we flag the problem?
  //
  float scale;
  if (circCorrT[0] > 1.e-5f) // NBNB-PAL: Temporary solution for above mentioned problem
    scale = float( 1.0f/circCorrT[0] );
  else  {
    LogKit::LogFormatted(LogKit::Low,"\nERROR: The circular temporal correlation (CircCorrT) is undefined. You\n");
    LogKit::LogFormatted(LogKit::Low,"       probably need to increase the number of layers...\n\nAborting\n");
    exit(1);
  }

  for(int k=0; k<nzp; k++)
    circCorrT[k] *= scale;
}


//--------------------------------------------------------------------
float *
SeismicParametersHolder::getPriorCorrTFiltered(int nz, int nzp) const
{
  // This is the cyclic and filtered version of CorrT which
  // has one or more zeros in the middle

  fftw_real * circCorrT = extractParamCorrFromCovAlpha(nzp);

  float * priorCorrTFiltered = new float[nzp];

  for(int i=0; i<nzp; i++ ) {
    int refk;

    if( i < nzp/2+1 )
      refk = i;
    else
      refk = nzp - i;

    if(refk < nz && circCorrT != NULL)
      priorCorrTFiltered[i] = circCorrT[refk];
    else
      priorCorrTFiltered[i] = 0.0;
  }

  fftw_free(circCorrT);

  return priorCorrTFiltered;
}

//--------------------------------------------------------------------
void
SeismicParametersHolder::writeFilePriorCorrT(fftw_real   * priorCorrT,
                                             const int   & nzp,
                                             const float & dt) const
{
  // This is the cyclic and filtered version of CorrT
  std::string baseName = IO::PrefixPrior() + IO::FileTemporalCorr() + IO::SuffixGeneralData();
  std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed
       << std::right
       << std::setprecision(6)
       << dt << "\n";
  for(int i=0 ; i<nzp; i++) {
    file << std::setw(9) << priorCorrT[i] << "\n";
  }
  file.close();
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::writeFilePostCorrT(const std::vector<float> & postCov,
                                            const std::string        & subDir,
                                            const std::string        & baseName) const
{
  std::string fileName = IO::makeFullFileName(subDir,baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed;
  file << std::setprecision(6);
  file << std::right;
  float c0 = 1.0f/postCov[0];

  for(int k=0 ; k < static_cast<int>(postCov.size()) ; k++)
    file << std::setw(9) << postCov[k]*c0 << "\n";

  file.close();
}

//--------------------------------------------------------------------
void
SeismicParametersHolder::writeFilePriorVariances(const ModelSettings      * modelSettings,
                                                 const std::vector<float> & priorCorrT,
                                                 const Surface            * priorCorrXY,
                                                 const float              & dt) const
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
      file << std::setw(13) << priorVar0_(i,j) << " ";
    }
    file << "\n";
  }
  file.close();

  NRLib::OpenWrite(file, fileName2);
  file << std::fixed
       << std::right
       << std::setprecision(8)
       << dt << "\n";
  for(int i=0 ; i<static_cast<int>(priorCorrT.size()); i++) {
    file << std::setw(11) << priorCorrT[i] << "\n";
  }
  file.close();

  IO::writeSurfaceToFile(*priorCorrXY, baseName3, IO::PathToCorrelations(), modelSettings->getOutputGridFormat());
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::printPriorVariances(void) const
{
  LogKit::LogFormatted(LogKit::Low,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Variances           ln Vp     ln Vs    ln Rho         \n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"Inversion grid:   %.1e   %.1e   %.1e (used by program)\n",priorVar0_(0,0),priorVar0_(1,1),priorVar0_(2,2));

  float corr01 = static_cast<float>(priorVar0_(0,1)/(sqrt(priorVar0_(0,0)*priorVar0_(1,1))));
  float corr02 = static_cast<float>(priorVar0_(0,2)/(sqrt(priorVar0_(0,0)*priorVar0_(2,2))));
  float corr12 = static_cast<float>(priorVar0_(1,2)/(sqrt(priorVar0_(1,1)*priorVar0_(2,2))));
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::Low,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::Low,"ln Rho |                     %5.2f \n",1.0f);
  LogKit::LogFormatted(LogKit::Low,"\n");

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
SeismicParametersHolder::printPostVariances(const NRLib::Matrix & postVar0) const
{
  LogKit::WriteHeader("Posterior Covariance");

  LogKit::LogFormatted(LogKit::Low,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"               ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-----------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"Variances:   %.1e   %.1e   %.1e    \n",postVar0(0,0),postVar0(1,1),postVar0(2,2));
  LogKit::LogFormatted(LogKit::Low,"\n");
  float corr01 = static_cast<float>(postVar0(0,1)/(sqrt(postVar0(0,0)*postVar0(1,1))));
  float corr02 = static_cast<float>(postVar0(0,2)/(sqrt(postVar0(0,0)*postVar0(2,2))));
  float corr12 = static_cast<float>(postVar0(1,2)/(sqrt(postVar0(1,1)*postVar0(2,2))));
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
SeismicParametersHolder::writeFilePostVariances(const NRLib::Matrix      & postVar0,
                                                const std::vector<float> & postCovAlpha00,
                                                const std::vector<float> & postCovBeta00,
                                                const std::vector<float> & postCovRho00) const
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
      file << std::setw(10) << postVar0(i,j) << " ";
    }
    file << "\n";
  }
  file.close();

  std::string baseName1 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vp" +IO::SuffixGeneralData();
  std::string baseName2 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vs" +IO::SuffixGeneralData();
  std::string baseName3 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Rho"+IO::SuffixGeneralData();
  writeFilePostCorrT(postCovAlpha00, IO::PathToCorrelations(), baseName1);
  writeFilePostCorrT(postCovBeta00,  IO::PathToCorrelations(), baseName2);
  writeFilePostCorrT(postCovRho00,   IO::PathToCorrelations(), baseName3);
}

//--------------------------------------------------------------------
void
SeismicParametersHolder::writeFilePostCovGrids(Simbox const * simbox) const
{
  std::string fileName;
  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp";
  covAlpha_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  covAlpha_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for Vp");
  covAlpha_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs";
  covBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  covBeta_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for Vs");
  covBeta_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho";
  covRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  covRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior covariance for density");
  covRho_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs";
  crCovAlphaBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  crCovAlphaBeta_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vp,Vs)");
  crCovAlphaBeta_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho";
  crCovAlphaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  crCovAlphaRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vp,density)");
  crCovAlphaRho_ ->endAccess();

  fileName = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho";
  crCovBetaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  crCovBetaRho_ ->writeFile(fileName, IO::PathToCorrelations(), simbox, "Posterior cross-covariance for (Vs,density)");
  crCovBetaRho_ ->endAccess();
}

//-------------------------------------------------------------------
fftw_real *
SeismicParametersHolder::extractParamCorrFromCovAlpha(int nzp) const
{
  assert(covAlpha_->getIsTransformed() == false);

  covAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);

  fftw_real * circCorrT = reinterpret_cast<fftw_real*>(fftw_malloc(2*(nzp/2+1)*sizeof(fftw_real)));
  //int         refk;
  float       constant = covAlpha_->getRealValue(0,0,0);

  for(int k = 0 ; k < 2*(nzp/2+1) ; k++ ){
    if(k < nzp)
      circCorrT[k] = covAlpha_->getRealValue(0,0,k,true)/constant;
    else
      circCorrT[k] = RMISSING;
  }

  covAlpha_->endAccess();

  return circCorrT;//fftw_free(circCorrT);
}
//--------------------------------------------------------------------
void
SeismicParametersHolder::updatePriorVar()
{
  priorVar0_(0,0) = getOrigin(covAlpha_);
  priorVar0_(1,1) = getOrigin(covBeta_);
  priorVar0_(2,2) = getOrigin(covRho_);
  priorVar0_(0,1) = getOrigin(crCovAlphaBeta_);
  priorVar0_(1,0) = priorVar0_(0,1);
  priorVar0_(2,0) = getOrigin(crCovAlphaRho_);
  priorVar0_(0,2) = priorVar0_(2,0);
  priorVar0_(2,1) = getOrigin(crCovBetaRho_);
  priorVar0_(1,2) = priorVar0_(2,1);
}
//--------------------------------------------------------------------
float
SeismicParametersHolder::getOrigin(FFTGrid * grid) const
{
  grid->setAccessMode(FFTGrid::RANDOMACCESS);
  float value = grid->getRealValue(0,0,0);
  grid->endAccess();
  return value;
}

//--------------------------------------------------------------------
std::vector<float>
SeismicParametersHolder::createPostCov00(FFTGrid * postCov) const
{
  int nz = postCov->getNz();
  std::vector<float> postCov00(nz);

  postCov->setAccessMode(FFTGrid::RANDOMACCESS);
  for(int k=0 ; k < nz ; k++)
    postCov00[k] = postCov->getRealValue(0,0,k);

  postCov->endAccess();
  return postCov00;
}
//--------------------------------------------------------------------
