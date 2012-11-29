#include "rplib/syntwelldata.h"

#include "src/definitions.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>


SyntWellData::SyntWellData(double                     const  & trend1,
                           double                     const  & trend2,
                           std::vector<float>         const  & alpha,
                           std::vector<float>         const  & beta,
                           std::vector<float>         const  & rho,
                           std::vector<int>           const  & faciesLog,
                           std::vector<std::string>   const  & faciesNames):
nBins_(static_cast<int>(alpha.size())),
nFacies_(static_cast<int>(faciesNames.size())),
trend1Value_(trend1),
trend2Value_(trend2){

  assert(alpha.size()==beta.size() && beta.size()==rho.size());

  faciesNames_.resize(nFacies_);
  for(int i=0; i<nFacies_; i++){
    faciesNames_[i] = faciesNames[i];
  }

  ipos_ = new int[nBins_];
  jpos_ = new int[nBins_];
  kpos_ = new int[nBins_];
  alpha_.resize(nBins_,0);
  beta_.resize(nBins_,0);
  rho_.resize(nBins_,0);
  faciesLog_.resize(nBins_,-1);

  for(int i=0; i<nBins_; i++){
    alpha_[i] = alpha[i];
    beta_[i] = beta[i];
    rho_[i] = rho[i];
    faciesLog_[i] = faciesLog[i];
    ipos_[i] = 0;
    jpos_[i] = 0;
    kpos_[i] = i; //only vertical
  }

}

SyntWellData::~SyntWellData(){
  delete [] ipos_;
  delete [] jpos_;
  delete [] kpos_;
}

