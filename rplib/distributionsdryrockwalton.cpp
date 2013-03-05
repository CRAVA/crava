#include "rplib/distributionsdryrockwalton.h"

#include "rplib/distributionssolid.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/dryrockwalton.h"
#include "rplib/solid.h"
#include "rplib/demmodelling.h"

#include "nrlib/random/distribution.hpp"

#include <cassert>
#include <typeinfo>
#include <typeinfo>

DistributionsDryRockWalton::DistributionsDryRockWalton(DistributionsSolid                         * distr_solid,
                                                       DistributionWithTrend                      * distr_friction_weight,
                                                       DistributionWithTrend                      * distr_pressure,
                                                       DistributionWithTrend                      * distr_porosity,
                                                       DistributionWithTrend                      * distr_coord_number,
                                                       std::vector<double>                        & alpha)
: DistributionsDryRock()
{
  distr_solid_      = distr_solid->Clone();

  if(distr_friction_weight->GetIsShared() == false)
    distr_friction_weight_ = distr_friction_weight->Clone();
  else
    distr_friction_weight_ = distr_friction_weight;

  if(distr_pressure->GetIsShared() == false)
    distr_pressure_ = distr_pressure->Clone();
  else
    distr_pressure_ = distr_pressure;

  if(distr_porosity->GetIsShared() == false)
    distr_porosity_ = distr_porosity->Clone();
  else
    distr_porosity_ = distr_porosity;

  if(distr_coord_number->GetIsShared() == false)
    distr_coord_number_ = distr_coord_number->Clone();
  else
    distr_coord_number_ = distr_coord_number;

  alpha_  = alpha;   // Order in alpha: friction_weight, pressure, porosity, coord_number

}

DistributionsDryRockWalton::DistributionsDryRockWalton(const DistributionsDryRockWalton & dist)
: DistributionsDryRock(dist)
{
  distr_solid_ = dist.distr_solid_->Clone();

  if(dist.distr_friction_weight_->GetIsShared() == false)
    distr_friction_weight_ = dist.distr_friction_weight_->Clone();
  else
    distr_friction_weight_ = dist.distr_friction_weight_;

  if(dist.distr_pressure_->GetIsShared() == false)
    distr_pressure_ = dist.distr_pressure_->Clone();
  else
    distr_pressure_ = dist.distr_pressure_;

  if(dist.distr_porosity_->GetIsShared() == false)
    distr_porosity_ = dist.distr_porosity_->Clone();
  else
    distr_porosity_ = dist.distr_porosity_;

  if(dist.distr_coord_number_->GetIsShared() == false)
    distr_coord_number_ = dist.distr_coord_number_->Clone();
  else
    distr_coord_number_ = dist.distr_coord_number_;

  alpha_ = dist.alpha_;   // Order in alpha: friction_weight, pressure, porosity, coord_number

}

DistributionsDryRockWalton::~DistributionsDryRockWalton()
{
  delete distr_solid_;

  if(distr_friction_weight_->GetIsShared() == false)
    delete distr_friction_weight_;

  if(distr_pressure_->GetIsShared() == false)
    delete distr_pressure_;

  if(distr_porosity_->GetIsShared() == false)
    delete distr_porosity_;

  if(distr_coord_number_->GetIsShared() == false)
    delete distr_coord_number_;
}

DistributionsDryRock *
DistributionsDryRockWalton::Clone() const
{
  return new DistributionsDryRockWalton(*this);
}

DryRock *
DistributionsDryRockWalton::GenerateSample(const std::vector<double> & trend_params)
{
  Solid * solid     = distr_solid_->GenerateSample(trend_params);

  std::vector<double> u(4);
  for(size_t i = 0; i < u.size(); i++)
    u[i] = NRLib::Random::Unif01();

  DryRock * new_dryrock = GetSample(u, trend_params, solid);

  // Deep copy taken by constructor of DryRockWalton, hence delete solid here
  delete solid;
  return new_dryrock;
}

bool
DistributionsDryRockWalton::HasDistribution() const
{
  if (distr_solid_->HasDistribution() == true || distr_friction_weight_->GetIsDistribution() == true ||
      distr_pressure_->GetIsDistribution() == true || distr_porosity_->GetIsDistribution() == true||
      distr_coord_number_->GetIsDistribution() == true)
        return true;

  return false;
}

std::vector<bool>
DistributionsDryRockWalton::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> friction_weight_trend          = distr_friction_weight_->GetUseTrendCube();
  std::vector<bool> pressure_trend                 = distr_pressure_->GetUseTrendCube();
  std::vector<bool> porosity_trend                 = distr_porosity_ ->GetUseTrendCube();
  std::vector<bool> coord_number_trend             = distr_coord_number_->GetUseTrendCube();

  std::vector<bool> solid_trend                    = distr_solid_->HasTrend();

  for(int i=0; i<2; i++) {
    if(friction_weight_trend[i] == true || pressure_trend[i] == true || porosity_trend[i] == true ||
       coord_number_trend[i]    == true || solid_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

DryRock *
DistributionsDryRockWalton::UpdateSample(double                      corr_param,
                                         bool                        param_is_time,
                                         const std::vector<double> & trend,
                                         const DryRock             * sample)
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

  assert(typeid(*sample) == typeid(DryRockWalton));
  const DryRockWalton * core_sample = dynamic_cast<const DryRockWalton *>(sample);

  Solid * updated_solid = distr_solid_->UpdateSample(corr_param,
                                                     param_is_time,
                                                     trend,
                                                     core_sample->GetSolid());

  DryRock * updated_sample = GetSample(u, trend, updated_solid);

  delete updated_solid;

  return updated_sample;
}

DryRock *
DistributionsDryRockWalton::GetSample(const std::vector<double>    & u,
                                      const std::vector<double>    & trend_params,
                                      const Solid                  * solid)
{

  // Order in u: friction_weight, pressure, porosity, coord_number

  double friction_weight = distr_friction_weight_->GetQuantileValue(u[0], trend_params[0], trend_params[1]);
  double pressure        = distr_pressure_->GetQuantileValue(u[1], trend_params[0], trend_params[1]);
  double porosity        = distr_porosity_->GetQuantileValue(u[2], trend_params[0], trend_params[1]);
  double coord_number    = distr_coord_number_->GetQuantileValue(u[3], trend_params[0], trend_params[1]);

  return new DryRockWalton(solid, friction_weight, pressure, porosity, coord_number, u);
}
