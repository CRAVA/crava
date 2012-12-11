#include <stdio.h>

#include "src/seismicparametersholder.h"
#include "src/fftgrid.h"
#include "src/corr.h"
#include "src/modelgeneral.h"

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
}

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

void
SeismicParametersHolder::setSeismicParameters(FFTGrid  * muAlpha,
                                              FFTGrid  * muBeta,
                                              FFTGrid  * muRho,
                                              Corr     * correlations)
{
  muAlpha_        = muAlpha;
  muBeta_         = muBeta ;
  muRho_          = muRho;
  covAlpha_       = correlations->getPostCovAlpha();
  covBeta_        = correlations->getPostCovBeta();
  covRho_         = correlations->getPostCovRho();
  crCovAlphaBeta_ = correlations->getPostCrCovAlphaBeta();
  crCovAlphaRho_  = correlations->getPostCrCovAlphaRho();
  crCovBetaRho_   = correlations->getPostCrCovBetaRho();
}

void
SeismicParametersHolder::allocateGrids(const int nx, const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad)
{
  muAlpha_ = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  muAlpha_->fillInConstant(0.0);
  muBeta_  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  muBeta_->fillInConstant(0.0);
  muRho_   = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  muRho_->fillInConstant(0.0);
  covAlpha_ = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  covAlpha_->fillInConstant(0.0);
  covBeta_  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  covBeta_->fillInConstant(0.0);
  covRho_   = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  covRho_->fillInConstant(0.0);
  crCovAlphaBeta_ = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovAlphaBeta_->fillInConstant(0.0);
  crCovAlphaRho_  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovAlphaRho_->fillInConstant(0.0);
  crCovBetaRho_   = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovBetaRho_->fillInConstant(0.0);
}

