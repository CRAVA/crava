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
 faciescount_ = 0;
 facies_ = IMISSING;
 nFacies_ = 0;
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

void CBWellPt::AddLog(Gamma gamma, float val) 
{
  if (gamma == ALPHA && val!= RMISSING) {  
    alpha_ += float(log(val)); alphaOrig_ = alpha_;
    noValidObsAlpha_++;   
  }
  else if (gamma == BETA && val != RMISSING) {
    beta_ += float(log(val)); betaOrig_ = beta_;
    noValidObsBeta_++;   
  }
  else {
    if (val != RMISSING) {
      rho_ += float(log(val)); rhoOrig_ = rho_;
      noValidObsRho_++;
    }
  }
}

void CBWellPt::SubtractOnly(float alpha, float beta, float rho) 
{
  if (alpha != RMISSING){
    alpha_ = alphaOrig_ - alpha;		
  }
  if (beta != RMISSING) {
    beta_ = betaOrig_ - beta;		
  }
  if (rho != RMISSING) {
    rho_ = rhoOrig_ - rho;		
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

void CBWellPt::SetFaciesParam(int nfac, int *fnr)
{
  faciesnr_ = fnr;
  nFacies_ = nfac;
  faciescount_ = new int[nfac];
  int i;
  for(i=0;i<nfac;i++)
    faciescount_[i] = 0;
}

void CBWellPt::AddFaciesLog(int val)
{
  int i;
  for(i=0;i<nFacies_;i++)
    if(faciesnr_[i]==val)
      faciescount_[i]++;
}
void CBWellPt::FindLargestFaciesCount(double unif01)
{
  int i;
  int maxcount = 0;
  for(i=0;i<nFacies_;i++)
  {
    if(maxcount<faciescount_[i])
    {
      maxcount = faciescount_[i];
      facies_ = faciesnr_[i];
    }
  }
  if(maxcount==0) // missing data
	  facies_ = IMISSING;
  else
  {
  // check if more than one has the maxcount, then do a random drawing
  int n = 0;
  int j = 0;
  //RandomGen random = model->getRandomGen();
  for(i=0;i<nFacies_;i++)
  {
    if(faciescount_[i]==maxcount)
      n++;
  }
  if(n>1)
  {
    //j = (int) random.unif01()*n;
    j = (int) unif01*n;
    i = -1;
    int k;
    for(k=0;k<j;k++)
    {
      i++;
      while(faciescount_[i]!=maxcount)
        i++;	
    }
    facies_ = faciesnr_[i];
  }
  }
}
bool CBWellPt::FaciesDefined() const
{
  if(faciescount_ == 0)
 // if(facies_!=IMISSING)
    return 0;
  else
    return 1;
}
bool CBWellPt::FaciesNotMissing() const
{
  if(facies_!=IMISSING)
    return 1;
  else 
    return 0;
}
