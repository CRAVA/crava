#ifndef ROCK_H
#define ROCK_H

#include "rplib/fluid.h"
#include "rplib/mineral.h"


// Abstract rock class.
class Rock {
public:

  Rock(const Fluid * base_fluid, double rmissing);

  virtual ~Rock();

  const Fluid * GetBaseFluid()                                  const { return(base_fluid_); }
  void GetSeismicParams(double & vp, double & vs, double & rho) const
  {
    vp  = vp_;
    vs  = vs_;
    rho = rho_;
  }
  void GetPoro(double & poro)                                   const {poro = poro_;} // rmissing_ if rock without porosity.
  void GetKMineral(double & k_mineral)                          const {k_mineral = k_mineral_;} 
  void CalculateMeanSeismicParams(double & mean_vp,
                                  double & mean_vs,
                                  double & mean_rho);
  void CalculateVarVp(double & var_vp);
  void CalculateVarVs(double & var_vs);
  void CalculateVarRho(double & var_rho);
  void CalculateCrCovVpVs(double & cr_cov_vp_vs);
  void CalculateCrCovVpRho(double & cr_cov_vp_rho);
  void CalculateCrCovVsRho(double & cr_cov_vs_rho);
  void CalculateMeanPoro(double & mean_poro); // rmissing_ if rock without porosity.
  void CalculateVarPoro(double & var_poro);   // rmissing_ if rock without porosity.

  bool IsRMissing(double x)                                     const {return x == rmissing_;}

  // Diagenesis of rock properties may be implemented in derived classes, 
  // in the form of rules that time develop the sampled values vp_, vs_, rho_, and poro_ (if applicable).
  // The function is not allowed to do changes to the rock object.
  // Return values are seismic parameters and porosity (may be rmissing_) for the times specified in the input.
  // Base class implementation is equivalent to no diagenesis.
  virtual void DoDiagenesis(const std::vector< int > & time_lapse, // Measured in days. 
                            std::vector< double >    & vp, 
                            std::vector< double >    & vs, 
                            std::vector< double >    & rho,
                            std::vector< double >    & poro) const
  {
    vp.resize(time_lapse.size(),vp_);
    vs.resize(time_lapse.size(),vs_);
    rho.resize(time_lapse.size(),rho_);
    poro.resize(time_lapse.size(),poro_);
  }

  virtual void ReSample() = 0;  // Resamples the class parameters, using base_fluid_. Poro not resampled if non-porous rock.

protected:

  void CalculateMeanAndCov();
	
  const Fluid * base_fluid_;

  double rmissing_;
  double vp_;                    // Updated by ReSample().
  double vs_;                    // Updated by ReSample().
  double rho_;                   // Updated by ReSample().
  double poro_;                  // Updated by ReSample() if applicable.
  double k_mineral_;             // Bulk modulus for rock if poro_ were 0. Fixed value, not resampled.
  std::vector<double> mean_seis_params_;
  std::vector<double> cov_seis_params_;
  double mean_poro_;
  double var_poro_;
  bool mean_cov_is_calculated_;  // Initially set to false. 
                                 // The first time one of the public Calculate-functions is used,
                                 // CalculateMeanAndCov() calculates all mean and cov values, and
                                 // mean_cov_is_calculated_ is set to true.

};

#endif
