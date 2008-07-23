#include "src/bwellpt.h"
#include "lib/global_def.h"
#include "lib/random.h"
#include "nrlib/iotools/logkit.hpp"

#include <stdio.h>
#include <assert.h>
#include <math.h>


CBWellPt::CBWellPt(int i, int j, int k):
i_(i), j_(j), k_(k)
{ 
  alpha_ = beta_ = rho_ = alphaOrig_ = betaOrig_ = rhoOrig_ = 0.0; 
  noValidObsAlpha_ = noValidObsBeta_ = noValidObsRho_ = 0;
}

CBWellPt::~CBWellPt(void)
{
}

void CBWellPt::AddLog(float alpha, float beta, float rho) 
{
  if (alpha != RMISSING) {
    alpha_ += float(log(alpha)); alphaOrig_ = alpha_;
    noValidObsAlpha_++;
  }
  if (beta != RMISSING) {
    beta_ += float(log(beta)); betaOrig_ = beta_;
    noValidObsBeta_++;
  }

  if (rho != RMISSING) {
    rho_ += float(log(rho)); rhoOrig_ = rho_;
    noValidObsRho_++;
  }
}

void CBWellPt::SubtractOnly(Gamma type, float gamma) {
  if (gamma == RMISSING)
    return;

  switch (type) {
  case ALPHA :	
    alpha_ = alphaOrig_ - gamma;		
    break;
  case BETA :
    beta_ = betaOrig_ - gamma;		
    break;
  case RHO :
    rho_ = rhoOrig_ - gamma;		
    break;
  default :
    assert(1 == 0);
    break;
  }  
}

void CBWellPt::Divide() 
{
  if (noValidObsAlpha_ > 0) {
    alpha_ /= noValidObsAlpha_; alphaOrig_ = alpha_;
  }
  if (noValidObsBeta_ > 0) {
    beta_ /= noValidObsBeta_; betaOrig_ = beta_;
  }

  if (noValidObsRho_ > 0) {
    rho_ /= noValidObsRho_; rhoOrig_ = rho_;
  }  
}
