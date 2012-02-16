#ifndef RPLIB_ROCK_PHYSICS_STORAGE_HPP
#define RPLIB_ROCK_PHYSICS_STORAGE_HPP

#include "rplib/distributionsrockt0.h"
#include "nrlib/trend/trendstorage.hpp"

class RockPhysicsStorage {
public:
  RockPhysicsStorage();

  virtual ~RockPhysicsStorage();

  virtual DistributionsRockT0 * GenerateRockPhysics(const std::string & /*path*/,
                                                    std::string       & /*errTxt*/) const = 0;
};

class GaussianRockPhysicsStorage : public RockPhysicsStorage {
public:
  GaussianRockPhysicsStorage(NRLib::TrendStorage * mean_vp,
                             NRLib::TrendStorage * mean_vs,
                             NRLib::TrendStorage * mean_density,
                             NRLib::TrendStorage * variance_vp,
                             NRLib::TrendStorage * variance_vs,
                             NRLib::TrendStorage * variance_density,
                             NRLib::TrendStorage * correlation_vp_vs,
                             NRLib::TrendStorage * correlation_vp_density,
                             NRLib::TrendStorage * correlation_vs_density);

  virtual ~GaussianRockPhysicsStorage();

  virtual DistributionsRockT0 * GenerateRockPhysics(const std::string & path,
                                                    std::string       & errTxt) const;

private:

  NRLib::TrendStorage * mean_vp_;
  NRLib::TrendStorage * mean_vs_;
  NRLib::TrendStorage * mean_density_;
  NRLib::TrendStorage * variance_vp_;
  NRLib::TrendStorage * variance_vs_;
  NRLib::TrendStorage * variance_density_;
  NRLib::TrendStorage * correlation_vp_vs_;
  NRLib::TrendStorage * correlation_vp_density_;
  NRLib::TrendStorage * correlation_vs_density_;
};

#endif
