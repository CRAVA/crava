#ifndef STATE4D_H
#define STATE4D_H

#include <vector>
#include <map>
#include <nrlib/flens/nrlib_flens.hpp>

class SeismicParametersHolder;
class FFTGrid;
class TimeEvolution;
class TimeLine;
class DistributionsRock;
class Simbox;

// This class holds FFTGrids for the 4D inversion. This includes grids for the static and dynamic variables \mu and \sigma, and the covariance grids between static and dynamic covariances.
// The grids are:
//   mu_static (3 grids)
//   mu_dynamic (3 grids)
//   sigma_static (6 grids)
//   sigma_dynamic (6 grids)
//   sigma_static_dynamic (9 grids, no symmetry)

// NB: Naming convensions open for dissucion.

class State4D
{

public:
  State4D();

  ~State4D();

  void setStaticMu (FFTGrid *vp, FFTGrid *vs, FFTGrid *rho);
  void setDynamicMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho);

  void updateStaticMu(FFTGrid * vp, FFTGrid * vs, FFTGrid * rho);
  void updateDynamicMu(FFTGrid * vp, FFTGrid * vs, FFTGrid * rho);

  void setStaticSigma (FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho);
  void setDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho);

  void setStaticDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvp, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhovp, FFTGrid *rhovs, FFTGrid *rhorho);  // OBS note order of parameters
  void setRelativeGridBase(int nx, int ny, int nz, int nxPad, int nyPad, int nzPad);
  NRLib::Matrix GetFullCov();
  NRLib::Vector GetFullMean000();

  void   merge(SeismicParametersHolder & current_state );
  void   updateWithSingleParameter(FFTGrid * Epost, FFTGrid * CovPost, int parameterNumber);
  void   updateAllignment(FFTGrid* new_vp_relative_to_base);
  void   mergeCov(std::vector<FFTGrid * > & sigma);
  void   split(SeismicParametersHolder & current_state );
  void   evolve(int time_step, const TimeEvolution & timeEvolution );
  std::vector<FFTGrid*> doRockPhysicsInversion(TimeLine                               & time_line,
                                               const std::vector<DistributionsRock *>   rock_distributions,
                                               TimeEvolution                          & timeEvolution);


  bool   isActive() const {return(mu_static_.size() > 0);}

  // Suggestions for get-functions. Not yet used. Naming etc to be decided later when actual in use.
  FFTGrid * getMuVpStatic(void)  const { return mu_static_[0]; }
  FFTGrid * getMuVsStatic(void)  const { return mu_static_[1]; }
  FFTGrid * getMuRhoStatic(void) const { return mu_static_[2]; }

  FFTGrid * getMuVpDynamic(void)  const { return mu_dynamic_[0]; }
  FFTGrid * getMuVsDynamic(void)  const { return mu_dynamic_[1]; }
  FFTGrid * getMuRhoDynamic(void) const { return mu_dynamic_[2]; }

  FFTGrid * getCovVpVpStaticStatic(void)   const { return sigma_static_static_[0]; }
  FFTGrid * getCovVpVsStaticStatic(void)   const { return sigma_static_static_[1]; }
  FFTGrid * getCovVpRhoStaticStatic(void)  const { return sigma_static_static_[2]; }
  FFTGrid * getCovVsVsStaticStatic(void)   const  { return sigma_static_static_[3]; }
  FFTGrid * getCovVsRhoStaticStatic(void)  const { return sigma_static_static_[4]; }
  FFTGrid * getCovRhoRhoStaticStatic(void) const { return sigma_static_static_[5]; }

  FFTGrid * getCovVpVpDynamicDynamic(void)   const { return sigma_dynamic_dynamic_[0]; }
  FFTGrid * getCovVpVsDynamicDynamic(void)   const { return sigma_dynamic_dynamic_[1]; }
  FFTGrid * getCovVpRhoDynamicDynamic(void)  const { return sigma_dynamic_dynamic_[2]; }
  FFTGrid * getCovVsVsDynamicDynamic(void)   const { return sigma_dynamic_dynamic_[3]; }
  FFTGrid * getCovVsRhoDynamicDynamic(void)  const { return sigma_dynamic_dynamic_[4]; }
  FFTGrid * getCovRhoRhoDynamicDynamic(void) const { return sigma_dynamic_dynamic_[5]; }

  // Note the order in private data stucture
  FFTGrid * getCovVpVpStaticDynamic(void)  const { return sigma_static_dynamic_[0];  }
  FFTGrid * getCovVpVsStaticDynamic(void)  const { return sigma_static_dynamic_[1];  }
  FFTGrid * getCovVpRhoStaticDynamic(void) const { return sigma_static_dynamic_[2];  }
  FFTGrid * getCovVsVpStaticDynamic(void)  const { return sigma_static_dynamic_[3];  }
  FFTGrid * getCovVsVsStaticDynamic(void)  const { return sigma_static_dynamic_[4];  }
  FFTGrid * getCovVsRhoStaticDynamic(void) const { return sigma_static_dynamic_[5];  }
  FFTGrid * getCovRhoVpStaticDynamic(void) const { return sigma_static_dynamic_[6];  }
  FFTGrid * getCovRhoVsStaticDynamic(void) const { return sigma_static_dynamic_[7];  }
  FFTGrid * getCovRhoRhoStaticDynamic(void)const { return sigma_static_dynamic_[8];  }
  FFTGrid * getRelativeVelocity(void)      const { return velocity_relative_to_base_;}


  void      FFT();
  void      iFFT();
  void      iFFTMean();
  void      iFFTCov();

private:
  bool allGridsAreTransformed();
  FFTGrid *              velocity_relative_to_base_;  //  V_current/V_initial
  std::vector<FFTGrid *> mu_static_;            // [0] = vp, [1] = vs, [2] = rho
  std::vector<FFTGrid *> mu_dynamic_;           // [0] = vp, [1] = vs, [2] = rho
  std::vector<FFTGrid *> sigma_static_static_;  // [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho (all static)
  std::vector<FFTGrid *> sigma_dynamic_dynamic_;// [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho (all dynamix)
  std::vector<FFTGrid *> sigma_static_dynamic_; // [0] = vpStat_vpDyn, [1] = vpStat_vsDyn, [2] = vpStat_rhoDyn ,[3] = vsStat_vpDyn,
                                                // [4] = vsStat_vsDyn, [5] = vsStat_rhoDyn, [6]= rhoStat_vpDyn, [7] = rhoStat_vsDyn, [8] = rhoStat_rhoDyn
  std::vector<std::vector<double> >    makeSeismicParamsFromrockSample(std::vector<std::vector<std::vector<double> > > rS);
  std::vector<double>                  getRockPropertiesFromRockSample(std::vector<std::vector<std::vector<double> > > rS,int varNumber);
};

#endif

