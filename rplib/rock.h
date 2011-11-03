#ifndef ROCK_H
#define ROCK_H

#include "rplib/fluid.h"
#include "rplib/mineral.h"


// Abstract rock class.
class Rock {
public:

  Rock(const Fluid * base_fluid);

  virtual ~Rock();

  const Fluid * GetBaseFluid()                                  const { return(base_fluid_); }
  void GetSeismicParams(double & vp, double & vs, double & rho) const
  {
    vp = vp_;
    vs = vs_;
    rho = rho_;
  }
  void GetPoro(double & poro)                                   const {poro = poro_;} // Crava's RMISSING if rock without porosity.
  void GetKMineral(double & k_mineral)                          const {k_mineral = k_mineral_;} 
  void GetMeanSeismicParams(double & mean_vp,
                            double & mean_vs,
                            double & mean_rho);
  void GetVarVp(double & var_vp);
  void GetVarVs(double & var_vs);
  void GetVarRho(double & var_rho);
  void GetCrCovVpVs(double & cr_cov_vp_vs);
  void GetCrCovVpRho(double & cr_cov_vp_rho);
  void GetCrCovVsRho(double & cr_cov_vs_rho);
  void GetMeanPoro(double & mean_poro); // Crava's RMISSING if rock without porosity.
  void GetVarPoro(double & var_poro);   // Crava's RMISSING if rock without porosity.

  virtual void ReSample() = 0;  // Resamples the class parameters, using base_fluid_. Poro not resampled if non-porous rock.

protected:

  void CalculateMeanAndCov();
	
  const Fluid * base_fluid_;

  double rmissing_;
  double vp_;        // Updated by ReSample()
  double vs_;        // Updated by ReSample()
  double rho_;       // Updated by ReSample()
  double poro_;      // Updated by ReSample() if applicable.
  double k_mineral_; // Bulk modulus for rock if poro_ were 0. Fixed value, not resampled.
  std::vector<double> mean_seis_params_;
  std::vector<double> cov_seis_params_;
  double mean_poro_;
  double var_poro_;
  bool mean_cov_is_calculated_;

};

#endif
