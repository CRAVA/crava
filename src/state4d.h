#ifndef STATE4D_H
#define STATE4D_H

#include <vector>

class FFTGrid;

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

  void SetStaticMu (FFTGrid *vp, FFTGrid *vs, FFTGrid *rho);
  void SetDynamicMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho);

  void SetStaticSigma (FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho);
  void SetDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho);

  void SetStaticDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvp, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhovp, FFTGrid *rhovs, FFTGrid *rhorho);  // OBS note order of parameters

  /*
  // Suggestions for get-functions. Not yet used. Naming etc to be decided later when actual in use.

  FFTGrid * getVpStatic(void)  { return mu_static_[0]; }
  FFTGrid * getVsStatic(void)  { return mu_static_[1]; }
  FFTGrid * getRhoStatic(void) { return mu_static_[2]; }

  FFTGrid * getVpDynamic(void)  { return mu_dynamic_[0]; }
  FFTGrid * getVsDynamic(void)  { return mu_dynamic_[1]; }
  FFTGrid * getRhoDynamic(void) { return mu_dynamic_[2]; }

  FFTGrid * getCovVpVpStaticStatic(void)   { return sigma_static_static_[0]; }
  FFTGrid * getCovVpVsStaticStatic(void)   { return sigma_static_static_[1]; }
  FFTGrid * getCovVpRhoStaticStatic(void)  { return sigma_static_static_[2]; }
  FFTGrid * getCovVsVsStaticStatic(void)   { return sigma_static_static_[3]; }
  FFTGrid * getCovVsRhoStaticStatic(void)  { return sigma_static_static_[4]; }
  FFTGrid * getCovRhoRhoStaticStatic(void) { return sigma_static_static_[5]; }

  FFTGrid * getCovVpVpDynamicDynamic(void)   { return sigma_dynamic_dynamic_[0]; }
  FFTGrid * getCovVpVsDynamicDynamic(void)   { return sigma_dynamic_dynamic_[1]; }
  FFTGrid * getCovVpRhoDynamicDynamic(void)  { return sigma_dynamic_dynamic_[2]; }
  FFTGrid * getCovVsVsDynamicDynamic(void)   { return sigma_dynamic_dynamic_[3]; }
  FFTGrid * getCovVsRhoDynamicDynamic(void)  { return sigma_dynamic_dynamic_[4]; }
  FFTGrid * getCovRhoRhoDynamicDynamic(void) { return sigma_dynamic_dynamic_[5]; }

  // Note the order in private data stucture
  FFTGrid * getCovVpVpStaticDynamic(void)   { return sigma_static_dynamic_[0]; }
  FFTGrid * getCovVpVsStaticDynamic(void)   { return sigma_static_dynamic_[1]; }
  FFTGrid * getCovVpRhoStaticDynamic(void)  { return sigma_static_dynamic_[2]; }
  FFTGrid * getCovVsVsStaticDynamic(void)   { return sigma_static_dynamic_[3]; }
  FFTGrid * getCovVsRhoStaticDynamic(void)  { return sigma_static_dynamic_[4]; }
  FFTGrid * getCovRhoRhoStaticDynamic(void) { return sigma_static_dynamic_[5]; }
  FFTGrid * getCovVsVpStaticDynamic(void)   { return sigma_static_dynamic_[6]; }
  FFTGrid * getCovRhoVpStaticDynamic(void)  { return sigma_static_dynamic_[7]; }
  FFTGrid * getCovRhoVsStaticDynamic(void)  { return sigma_static_dynamic_[8]; }
  */

private:
  std::vector<FFTGrid *> mu_static_;            // [0] = vp, [1] = vs, [2] = rho
  std::vector<FFTGrid *> mu_dynamic_;           // [0] = vp, [1] = vs, [2] = rho
  std::vector<FFTGrid *> sigma_static_static_;  // [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho
  std::vector<FFTGrid *> sigma_dynamic_dynamic_;
  std::vector<FFTGrid *> sigma_static_dynamic_; // [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho, [6]=vs_vp, [7] = rho_vp, [8] = rho_vs

};

#endif
