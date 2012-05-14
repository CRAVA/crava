#include <stdio.h>

#include "src/seismicparametersholder.h"
#include "src/fftgrid.h"
#include "src/corr.h"

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

