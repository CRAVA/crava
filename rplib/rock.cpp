#include "rplib/rock.h"


Rock::Rock(const Fluid * base_fluid, double rmissing) : base_fluid_(base_fluid), rmissing_(rmissing)      
{
  vp_                     = rmissing_;   
  vs_                     = rmissing_;
  rho_                    = rmissing_;
  poro_                   = rmissing_;
  k_mineral_              = rmissing_;

  mean_cov_is_calculated_ = false;
  mean_seis_params_.resize(3, rmissing_);
  cov_seis_params_.resize(6, rmissing_);
  mean_poro_              = rmissing_;
  var_poro_               = rmissing_;
}

Rock::~Rock()
{
}

void 
Rock::CalculateMeanAndCov()
{
  int n_samples = 1000;

  double vp       = 0.0;
  double vs       = 0.0;
  double rho      = 0.0;
  double vp2      = 0.0;
  double vs2      = 0.0;
  double rho2     = 0.0;
  double vpvs     = 0.0;
  double vprho    = 0.0;
  double vsrho    = 0.0;
  double poro     = 0.0;  // To be updated only if Rock resamples poro_.
  double poro2    = 0.0;  // To be updated only if Rock resamples poro_.
  
  for (int i = 0; i < n_samples; ++i){
    ReSample();
    vp    += vp_;
    vs    += vs_;
    rho   += rho_;
    vp2   += vp_ * vp_;
    vs2   += vs_ * vs_;
    rho2  += rho_ * rho_;
    vpvs  += vp_ * vs_;
    vprho += vp_ * rho_;
    vsrho += vs_ * rho_;

    if (poro_ != rmissing_){  // Rock without porosity does not resample poro_.
      poro  += poro_;
      poro2 += poro_ * poro_;
    }
  }

  vp    /= double(n_samples);
  vs    /= double(n_samples);
  rho   /= double(n_samples);
  vp2   /= double(n_samples);
  vs2   /= double(n_samples);
  rho2  /= double(n_samples);
  vpvs  /= double(n_samples);
  vprho /= double(n_samples);
  vsrho /= double(n_samples);
  if (poro_ != rmissing_){  // Rock without porosity does not resample poro_.
    poro  /= double(n_samples);
    poro2 /= double(n_samples);
  }

  mean_seis_params_.resize(3, 0.0);
  cov_seis_params_.resize(6, 0.0);

  mean_seis_params_[0] = vp;
  mean_seis_params_[1] = vs;
  mean_seis_params_[2] = rho;
  cov_seis_params_[0]  = vp2   - vp * vp;
  cov_seis_params_[1]  = vs2   - vs * vs;
  cov_seis_params_[2]  = rho2  - rho * rho;
  cov_seis_params_[3]  = vpvs  - vp * vs;
  cov_seis_params_[4]  = vprho - vp * rho;
  cov_seis_params_[5]  = vsrho - vs * rho;

  if (poro_ != rmissing_){ // Rock without porosity does not resample poro_.
    mean_poro_ = poro;
    var_poro_  = poro2 - poro * poro;
  }

  mean_cov_is_calculated_ = true;

}

void 
Rock::CalculateMeanSeismicParams(double & mean_vp,
                                 double & mean_vs,
                                 double & mean_rho)  
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  mean_vp  = mean_seis_params_[0];
  mean_vs  = mean_seis_params_[1];
  mean_rho = mean_seis_params_[2];
}
void 
Rock::CalculateVarVp(double & var_vp)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_vp = cov_seis_params_[0];
}

void 
Rock::CalculateVarVs(double & var_vs)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_vs = cov_seis_params_[1];
}

void 
Rock::CalculateVarRho(double & var_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_rho = cov_seis_params_[2];
}

void 
Rock::CalculateCrCovVpVs(double & cr_cov_vp_vs)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vp_vs = cov_seis_params_[3];
}

void 
Rock::CalculateCrCovVpRho(double & cr_cov_vp_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vp_rho = cov_seis_params_[4];
}

void 
Rock::CalculateCrCovVsRho(double & cr_cov_vs_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vs_rho = cov_seis_params_[5];
}


void 
Rock::CalculateMeanPoro(double & mean_poro)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  mean_poro = mean_poro_;
}

void 
Rock::CalculateVarPoro(double & var_poro)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_poro = var_poro_;
}
