#include <stdio.h>

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "lib/utils.h"

#include "src/definitions.h"
#include "src/fftfilegrid.h"
#include "src/fftgrid.h"
#include "src/model.h"
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
  pointVar0_   = pointVar0;
  priorVar0_   = priorVar0;
  priorCorrT_  = priorCorrT;
  priorCorrXY_ = priorCorrXY;
  n_           = n;
  dt_          = dt;
}

Corr::~Corr(void)
{
  if(pointVar0_ != NULL) {
    for(int i=0;i<3;i++)
      delete [] pointVar0_[i];
    delete [] pointVar0_;
  }

  for(int i=0;i<3;i++)
    delete [] priorVar0_[i];
  delete [] priorVar0_;

  delete priorCorrXY_;

  delete [] priorCorrT_;

  if(priorCorrTFiltered_!=NULL)      
    delete [] priorCorrTFiltered_ ;

  if (postVar0_ != NULL) {
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

  if(postCovAlpha_!=NULL)      
    delete postCovAlpha_ ;
  if(postCovBeta_!=NULL)       
    delete postCovBeta_ ;
  if(postCovRho_!=NULL)        
    delete postCovRho_ ;
  if(postCrCovAlphaBeta_!=NULL)
    delete postCrCovAlphaBeta_ ;
  if(postCrCovAlphaRho_!=NULL) 
    delete postCrCovAlphaRho_ ;
  if(postCrCovBetaRho_!=NULL)  
    delete postCrCovBetaRho_;
}

//--------------------------------------------------------------------
void
Corr::createPostGrids(int nx,  int ny,  int nz,
                      int nxp, int nyp, int nzp,
                      bool fileGrid)
{
  postCovAlpha_       = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);                     
  postCovBeta_        = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);    
  postCovRho_         = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);    
  postCrCovAlphaBeta_ = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid); 
  postCrCovAlphaRho_  = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);                     
  postCrCovBetaRho_   = createFFTGrid(nx,ny,nz,nxp,nyp,nzp,fileGrid);    

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

//--------------------------------------------------------------------
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
float * 
Corr::getPriorCorrT(int &n, float &dt) const
{
  n = n_;
  dt = dt_;
  return priorCorrT_;
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
  LogKit::LogFormatted(LogKit::HIGH,"Backtransforming correlation grids from FFT domain to time domain...");
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
  LogKit::LogFormatted(LogKit::HIGH,"...done\n");
}

//--------------------------------------------------------------------
void
Corr::FFT(void)
{
  LogKit::LogFormatted(LogKit::HIGH,"Transforming correlation grids from time domain to FFT domain...");
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
  LogKit::LogFormatted(LogKit::HIGH,"...done\n");
}

//--------------------------------------------------------------------
void Corr::printPriorVariances(void) const
{
  LogKit::LogFormatted(LogKit::LOW,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"Variances           ln Vp     ln Vs    ln Rho         \n");
  LogKit::LogFormatted(LogKit::LOW,"---------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"Inversion grid:   %.1e   %.1e   %.1e (used by program)\n",priorVar0_[0][0],priorVar0_[1][1],priorVar0_[2][2]);
  if(pointVar0_ != NULL)
    LogKit::LogFormatted(LogKit::LOW,"Well logs     :   %.1e   %.1e   %.1e                  \n",pointVar0_[0][0],pointVar0_[1][1],pointVar0_[2][2]);

  float corr01 = priorVar0_[0][1]/(sqrt(priorVar0_[0][0]*priorVar0_[1][1]));
  float corr02 = priorVar0_[0][2]/(sqrt(priorVar0_[0][0]*priorVar0_[2][2]));
  float corr12 = priorVar0_[1][2]/(sqrt(priorVar0_[1][1]*priorVar0_[2][2]));
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::LOW,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::LOW,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::LOW,"ln Rho |                     %5.2f \n",1.0f);
}

//--------------------------------------------------------------------
void
Corr::printPostVariances(void) const
{
  Utils::writeHeader("Posterior Covariance");
  LogKit::LogFormatted(LogKit::LOW,"\nVariances and correlations for parameter residuals:\n");
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"               ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::LOW,"-----------------------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"Variances:   %.1e   %.1e   %.1e    \n",postVar0_[0][0],postVar0_[1][1],postVar0_[2][2]); 
  LogKit::LogFormatted(LogKit::LOW,"\n");
  float corr01 = postVar0_[0][1]/(sqrt(postVar0_[0][0]*postVar0_[1][1]));
  float corr02 = postVar0_[0][2]/(sqrt(postVar0_[0][0]*postVar0_[2][2]));
  float corr12 = postVar0_[1][2]/(sqrt(postVar0_[1][1]*postVar0_[2][2]));
  LogKit::LogFormatted(LogKit::LOW,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::LOW,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::LOW,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::LOW,"ln Rho |                     %5.2f \n",1.0f);
}
//--------------------------------------------------------------------
void
Corr::getPostVariances(void)
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

  postCovAlpha00_       = getPostCov00(postCovAlpha_);
  postCovBeta00_        = getPostCov00(postCovBeta_);
  postCovRho00_         = getPostCov00(postCovRho_);
  postCrCovAlphaBeta00_ = getPostCov00(postCrCovAlphaBeta_);
  postCrCovAlphaRho00_  = getPostCov00(postCrCovAlphaRho_);
  postCrCovBetaRho00_   = getPostCov00(postCrCovBetaRho_);
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
Corr::getPostCov00(FFTGrid * postCov)
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
Corr::writeFilePriorCorrT(float* corrT, int nzp) const
{
  // This is the cyclic and filtered version of CorrT
  std::string baseName = IO::FilePriorTemporalCorr() + IO::SuffixGeneralData();
  std::string fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), baseName);
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
void Corr::writeFilePriorVariances(void) const
{
  std::string baseName1 = IO::FilePriorParameterCov()           + IO::SuffixGeneralData();
  std::string baseName2 = IO::FilePriorTemporalCorrUnfiltered() + IO::SuffixGeneralData();
  std::string baseName3 = IO::FilePriorLateralCorr()            + IO::SuffixGeneralData();
  std::string fileName1 = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), baseName1);
  std::string fileName2 = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), baseName2);
  std::string fileName3 = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), baseName3);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName1);
  file << std::fixed
       << std::right
       << std::setprecision(6);
  for(int i=0 ; i<3 ; i++) {
    for(int j=0 ; j<3 ; j++) {
      file << std::setw(10) << priorVar0_[i][j] << " ";
    }
    file << "\n";
  }
  file.close();

  NRLib::OpenWrite(file, fileName2);
  file << std::fixed 
       << std::right  
       << std::setprecision(6)
       << dt_ << "\n";
  for(int i=0 ; i<n_ ; i++) {
    file << std::setw(9) << priorCorrT_[i] << "\n";
  }
  file.close();

  NRLib::WriteIrapClassicAsciiSurf(*priorCorrXY_, fileName3);
}

//--------------------------------------------------------------------
void
Corr::writeFilePostVariances(void) const
{
  std::string baseName = IO::FilePosteriorParameterCov() + IO::SuffixGeneralData();
  std::string fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), baseName);

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
  std::string baseName1 = IO::PrefixPosteriorTemporalCorr()+"Vp" +IO::SuffixGeneralData();
  std::string baseName2 = IO::PrefixPosteriorTemporalCorr()+"Vs" +IO::SuffixGeneralData();
  std::string baseName3 = IO::PrefixPosteriorTemporalCorr()+"Rho"+IO::SuffixGeneralData();
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
  std::string fileName = ModelSettings::makeFullFileName2(subDir,baseName);
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
Corr::writeFilePostCovGrids(Simbox * simbox) const
{
  std::string fileName;
  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCov()+"Vp");
  postCovAlpha_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovAlpha_ ->writeFile(IO::PrefixPosteriorCov()+"Vp",simbox, "Posterior covariance for Vp");
  postCovAlpha_ ->endAccess();

  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCov()+"Vs");
  postCovBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovBeta_ ->writeFile(IO::PrefixPosteriorCov()+"Vs",simbox, "Posterior covariance for Vs");
  postCovBeta_ ->endAccess();

  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCov()+"Rho");
  postCovRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovRho_ ->writeFile(IO::PrefixPosteriorCov()+"Rho",simbox, "Posterior covariance for density");
  postCovRho_ ->endAccess();

  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCrCov()+"VpVs");
  postCrCovAlphaBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaBeta_ ->writeFile(IO::PrefixPosteriorCrCov()+"VpVs",simbox, "Posterior cross-covariance for (Vp,Vs)");
  postCrCovAlphaBeta_ ->endAccess();

  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCrCov()+"VpRho");
  postCrCovAlphaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaRho_ ->writeFile(IO::PrefixPosteriorCrCov()+"VpRho",simbox, "Posterior cross-covariance for (Vp,density)");
  postCrCovAlphaRho_ ->endAccess();

  fileName = ModelSettings::makeFullFileName2(IO::PathToCorrelations(), IO::PrefixPosteriorCrCov()+"VsRho");
  postCrCovBetaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovBetaRho_ ->writeFile(IO::PrefixPosteriorCrCov()+"VsRho",simbox, "Posterior cross-covariance for (Vs,density)");
  postCrCovBetaRho_ ->endAccess();
}
