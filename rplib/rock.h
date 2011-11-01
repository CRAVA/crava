#ifndef ROCK_H
#define ROCK_H

#include "rplib/fluid.h"
#include "rplib/mineral.h"


// Abstract rock class.
class Rock {
public:

  Rock(const Fluid * base_fluid) : base_fluid_(base_fluid)      
  {
    double rmissing_tmp = -999999.999;   // Byttes til Cravas versjon av RMISSING
    vp_                 = rmissing_tmp;   
    vs_                 = rmissing_tmp;
    rho_                = rmissing_tmp;
    poro_               = rmissing_tmp;
  }

  virtual ~Rock(){}

  const Fluid * GetBaseFluid()                                  const { return(base_fluid_); }
  void GetSeismicParams(double & vp, double & vs, double & rho) const
  {
    vp = vp_;
    vs = vs_;
    rho = rho_;
  }
  void GetPoro(double & poro)                                   const {poro = poro_;} 

  virtual void ReSample() = 0;  // Resamples the class parameters, using base_fluid_.

protected:
	
  const Fluid * base_fluid_;

  double vp_;   // Updated by ReSample()
  double vs_;   // Updated by ReSample()
  double rho_;  // Updated by ReSample()
  double poro_; // Updated by ReSample() if applicable.

};

#endif
