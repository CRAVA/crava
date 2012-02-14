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
  GaussianRockPhysicsStorage(TrendStorage * mean_vp,
                             TrendStorage * mean_vs,
                             TrendStorage * mean_density,
                             TrendStorage * variance_vp,
                             TrendStorage * variance_vs,
                             TrendStorage * variance_density,
                             TrendStorage * correlation_vp_vs,
                             TrendStorage * correlation_vp_density,
                             TrendStorage * correlation_vs_density);

  virtual ~GaussianRockPhysicsStorage();

  virtual DistributionsRockT0 * GenerateRockPhysics(const std::string & path,
                                                    std::string       & errTxt) const;

private:

  TrendStorage * mean_vp_;
  TrendStorage * mean_vs_;
  TrendStorage * mean_density_;
  TrendStorage * variance_vp_;
  TrendStorage * variance_vs_;
  TrendStorage * variance_density_;
  TrendStorage * correlation_vp_vs_;
  TrendStorage * correlation_vp_density_;
  TrendStorage * correlation_vs_density_;
};

#endif
