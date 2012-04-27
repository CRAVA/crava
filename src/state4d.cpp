#include "src/state4d.h"

State4D::State4D(){

  // Initilizing the 27 grids. All grids initialize to NULL.
  mu_static_.resize(3, NULL);
  mu_dynamic_.resize(3, NULL);
  sigma_static_static_.resize(6, NULL);
  sigma_dynamic_dynamic_.resize(6, NULL);
  sigma_static_dynamic_.resize(9, NULL);
}

State4D::~State4D(){
}

void State4D::SetDynamicMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho){
  mu_dynamic_[0] = vp;
  mu_dynamic_[1] = vs;
  mu_dynamic_[2] = rho;
}

void State4D::SetStaticMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho){
  mu_static_[0] = vp;
  mu_static_[1] = vs;
  mu_static_[2] = rho;
}

void State4D::SetStaticSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho){
  sigma_static_static_[0] = vpvp;
  sigma_static_static_[1] = vpvs;
  sigma_static_static_[2] = vprho;
  sigma_static_static_[3] = vsvs;
  sigma_static_static_[4] = vsrho;
  sigma_static_static_[5] = rhorho;
}

void State4D::SetDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho){
  sigma_dynamic_dynamic_[0] = vpvp;
  sigma_dynamic_dynamic_[1] = vpvs;
  sigma_dynamic_dynamic_[2] = vprho;
  sigma_dynamic_dynamic_[3] = vsvs;
  sigma_dynamic_dynamic_[4] = vsrho;
  sigma_dynamic_dynamic_[5] = rhorho;
}
void State4D::SetStaticDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvp, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhovp, FFTGrid *rhovs, FFTGrid *rhorho){
  sigma_static_dynamic_[0] = vpvp;
  sigma_static_dynamic_[1] = vpvs;
  sigma_static_dynamic_[2] = vprho;
  sigma_static_dynamic_[3] = vsvs;
  sigma_static_dynamic_[4] = vsrho;
  sigma_static_dynamic_[5] = rhorho;
  sigma_static_dynamic_[6] = vsvp;  // OBS note order of elements! Not same as in parameter list.
  sigma_static_dynamic_[7] = rhovp;
  sigma_static_dynamic_[8] = rhovs;
}
