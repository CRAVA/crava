#include "rplib/rock.h"


Rock::Rock(const Fluid * base_fluid) : base_fluid_(base_fluid)      
{
  rmissing_               = -999999.999;   // Byttes til Cravas versjon av RMISSING
  vp_                     = rmissing_;   
  vs_                     = rmissing_;
  rho_                    = rmissing_;
  poro_                   = rmissing_;
  k_mineral_              = rmissing_;
  mean_cov_is_calculated_ = false;
  mean_seis_params_.resize(0,rmissing_);
  cov_seis_params_.resize(0,rmissing_);
  mean_poro_              = rmissing_;
  var_poro_               = rmissing_;
}

Rock::~Rock()
{
}

void 
Rock::CalculateMeanAndCov()
{
  int n_samples = 1000;                   // Er det ok å hardkode denne, og hva er da en bra nok verdi?

  ReSample();
  double vp       = vp_;
  double vs       = vs_;
  double rho      = rho_;
  double vp2      = vp_  * vp_;
  double vs2      = vs_  * vs_;
  double rho2     = rho_ * rho_;
  double vpvs     = vp_  * vs_;
  double vprho    = vp_  * rho_;
  double vsrho    = vs_  * rho_;

  double poro  = 0.0;
  double poro2 = 0.0;
  if (poro_ != rmissing_){ // Rock without porosity does not resample poro_.
    poro  = poro_;
    poro2 = poro_ * poro_;
  }

  for (int i = 1; i < n_samples; ++i){
    ReSample();
    vp    = (double(i) * vp    + vp_)  / (double(i+1));
    vs    = (double(i) * vs    + vs_)  / (double(i+1));
    rho   = (double(i) * rho   + rho_) / (double(i+1));
    vp2   = (double(i) * vp2   + vp_ * vp_)   / (double(i+1));
    vs2   = (double(i) * vs2   + vs_ * vs_)   / (double(i+1));
    rho2  = (double(i) * rho2  + rho_ * rho_) / (double(i+1));
    vpvs  = (double(i) * vpvs  + vp_ * vs_)   / (double(i+1));
    vprho = (double(i) * vprho + vp_ * rho_)  / (double(i+1));
    vsrho = (double(i) * vsrho + vs_ * rho_)  / (double(i+1));

    if (poro_ != rmissing_){  // Rock without porosity does not resample poro_.
      poro  = (double(i) * poro + poro_) / (double(i+1));
      poro2 = (double(i) * poro2 + poro_ * poro_) / (double(i+1));
    }
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

}

void 
Rock::GetMeanSeismicParams(double & mean_vp,
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
Rock::GetVarVp(double & var_vp)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_vp = cov_seis_params_[0];
}

void 
Rock::GetVarVs(double & var_vs)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_vs = cov_seis_params_[1];
}

void 
Rock::GetVarRho(double & var_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_rho = cov_seis_params_[2];
}

void 
Rock::GetCrCovVpVs(double & cr_cov_vp_vs)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vp_vs = cov_seis_params_[3];
}

void 
Rock::GetCrCovVpRho(double & cr_cov_vp_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vp_rho = cov_seis_params_[4];
}

void 
Rock::GetCrCovVsRho(double & cr_cov_vs_rho)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  cr_cov_vs_rho = cov_seis_params_[5];
}


void 
Rock::GetMeanPoro(double & mean_poro)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  mean_poro = mean_poro_;
}

void 
Rock::GetVarPoro(double & var_poro)
{
  if (!mean_cov_is_calculated_)
    CalculateMeanAndCov();
  var_poro = var_poro_;
}
