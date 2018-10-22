#include "rplib/syntwelldata.h"

#include "src/definitions.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>


SyntWellData::SyntWellData(double                     const  & trend1,
                           double                     const  & trend2,
                           int                                 i,
                           int                                 j,
                           std::vector<float>         const  & alpha,
                           std::vector<float>         const  & beta,
                           std::vector<float>         const  & rho,
                           std::vector<int>           const  & faciesLog,
                           std::vector<std::string>   const  & faciesNames):
nBins_(static_cast<int>(alpha.size())),
nFacies_(static_cast<int>(faciesNames.size())),
trend1Value_(trend1),
trend2Value_(trend2),
trend1Index_(i),
trend2Index_(j){

  assert(alpha.size()==beta.size() && beta.size()==rho.size());

  faciesNames_.resize(nFacies_);
  for(int ii=0; ii<nFacies_; ii++){
    faciesNames_[ii] = faciesNames[ii];
  }

  ipos_ = new int[nBins_];
  jpos_ = new int[nBins_];
  kpos_ = new int[nBins_];
  alpha_.resize(nBins_,0);
  beta_.resize(nBins_,0);
  rho_.resize(nBins_,0);
  faciesLog_.resize(nBins_,-1);

  for(int ii=0; ii<nBins_; ii++){
    alpha_[ii] = alpha[ii];
    beta_[ii] = beta[ii];
    rho_[ii] = rho[ii];
    faciesLog_[ii] = faciesLog[ii];
    ipos_[ii] = 0;
    jpos_[ii] = 0;
    kpos_[ii] = ii; //only vertical
  }

}

SyntWellData::~SyntWellData(){
  delete [] ipos_;
  delete [] jpos_;
  delete [] kpos_;
}

