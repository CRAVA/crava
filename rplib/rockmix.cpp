#include "rplib/rockmix.h"

#include "nrlib/exception/exception.hpp"

#include <numeric>



RockMixOfRock::RockMixOfRock(const std::vector<Rock*>      & rock,
                             const std::vector<double>     & volume_fraction,
                             DEMTools::MixMethod             mix_method)
: Rock()
{
  // Deep copy of rock:
  rock_.resize(rock.size());
  for (size_t i = 0; i < rock.size(); ++i) {
    rock_[i] = rock[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;

  vp_   = vs_ = rho_ = porosity_ = 0;

  if (rock.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> vp(rock.size()), vs(rock.size()), rho(rock.size()), porosity(rock.size());
    for (size_t i = 0; i < rock.size(); i++) {
      rock[i]->GetSeismicParams(vp[i], vs[i], rho[i]);
      porosity[i] = rock[i]->GetPorosity();
    }

    switch (mix_method_) {
      case DEMTools::Hill :
        vp_     = DEMTools::CalcEffectiveElasticModuliUsingHill(vp, volume_fraction);
        vs_     = DEMTools::CalcEffectiveElasticModuliUsingHill(vs, volume_fraction);
        break;
      case DEMTools::Reuss :
        vp_     = DEMTools::CalcEffectiveElasticModuliUsingReuss(vp, volume_fraction);     // homogeneous
        vs_     = DEMTools::CalcEffectiveElasticModuliUsingReuss(vs, volume_fraction);
        break;
      case DEMTools::Voigt :
        vp_     = DEMTools::CalcEffectiveElasticModuliUsingVoigt(vp, volume_fraction);
        vs_     = DEMTools::CalcEffectiveElasticModuliUsingVoigt(vs, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid rock mixing algorithm specified.");
    }
    rho_        = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
    porosity_   = DEMTools::CalcEffectivePorosity(porosity, volume_fraction);
  }
}


RockMixOfRock::~RockMixOfRock()
{
  for (size_t i = 0; i < rock_.size(); ++i)
    delete rock_[i];
}

Rock *
RockMixOfRock::Clone() const
{
  // Provide base class variables.
  RockMixOfRock * s = new RockMixOfRock(*this);

  // Provide variables specific to RockMixOfRock.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;
  s->porosity_        = this->porosity_;

  size_t n_rocks = this->rock_.size();
  s->rock_.resize(n_rocks);
  for (size_t i = 0; i < n_rocks; ++i)
    s->rock_[i] = this->rock_[i]->Clone();

  return s;
}

RockMixOfRock& RockMixOfRock::operator=(const RockMixOfRock& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;
    porosity_        = rhs.porosity_;

    size_t n_rocks_old = rock_.size();
    for (size_t i = 0; i < n_rocks_old; ++i)
      delete rock_[i];

    size_t n_rocks = rhs.rock_.size();
    rock_.resize(n_rocks);
    for (size_t i = 0; i < n_rocks; ++i)
      rock_[i] = rhs.rock_[i]->Clone();
  }
  return *this;
}

Rock *
RockMixOfRock::Evolve(const std::vector<int>          & delta_time,
                      const std::vector< Rock * >     & rock) const
{
  size_t n_rocks = rock_.size();
  std::vector<Rock*> rock_new(n_rocks);
  for (size_t i = 0; i < n_rocks; ++i)
    rock_new[i] = rock_[i]->Evolve(delta_time, rock);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Rock * rock_mixed_new = new RockMixOfRock(rock_new, volume_fraction, mix_method_);

  // Deep copy taken by constructor of RockMixOfRock, hence delete rock_new here:
  for (size_t i = 0; i < n_rocks; ++i)
    delete rock_new[i];

  return rock_mixed_new;
}

double
RockMixOfRock::GetPorosity() const
{
  return porosity_;
}

void
RockMixOfRock::SetPorosity(double porosity)
{
  porosity_ = porosity;
}

