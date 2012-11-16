#include "rplib/syntwelldata.h"
#include <string.h>
#include <stdlib.h>

SyntWellData::SyntWellData(double                     const  & trend1,
                           double                     const  & trend2,
                           std::vector<float>         const  & alpha,
                           std::vector<float>         const  & beta,
                           std::vector<float>         const  & rho,
                           std::vector<int>           const  & faciesLog,
                           std::vector<std::string>   const  & faciesNames):
nBins_(static_cast<int>(alpha.size())),
nFacies_(static_cast<int>(faciesLog.size())),
trend1Value_(trend1),
trend2Value_(trend2),
ipos_(NULL),
jpos_(NULL),
kpos_(NULL){

  assert(alpha.size()==beta.size() && beta.size()==rho.size());

  for(int i=0; i<nFacies_; i++){
    faciesLog_[i] = faciesLog[i];
    faciesNames_[i] = faciesNames[i];
  }

  ipos_ = new int[nBins_];
  jpos_ = new int[nBins_];
  kpos_ = new int[nBins_];

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

