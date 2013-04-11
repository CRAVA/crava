#include "rplib/dryrockwalton.h"

#include "rplib/solid.h"

#include "nrlib/math/constants.hpp"

#include <cassert>
#include <cmath>


DryRockWalton::DryRockWalton(const Solid                       * solid,
                             const double                      & friction_weight,
                             const double                      & pressure,
                             const double                      & porosity,
                             const double                      & coord_number,
                             const std::vector<double>         & u)
: DryRock()
{
  u_ = u; // u contains independent samples used in quantiles of friction_weight, pressure, porosity, coord_number

  solid_ = solid->Clone();

  friction_weight_ = friction_weight;
  pressure_        = pressure;
  total_porosity_  = porosity;

  if (coord_number < 0) {
    coord_number_         = CalcCoordNumber();
    calc_coord_number_    = true;
  }
  else {
    coord_number_         = coord_number;
    calc_coord_number_    = false;
  }

  ComputeElasticParams();

}

DryRockWalton::DryRockWalton() : DryRock()
{
  rho_ = k_ = mu_ = 0;
}

DryRockWalton::~DryRockWalton()
{
  delete solid_;

}

DryRockWalton& DryRockWalton::operator=(const DryRockWalton& rhs)
{
  if (this != &rhs) {
    DryRock::operator=(rhs);

    friction_weight_    = rhs.friction_weight_;
    pressure_           = rhs.pressure_;
    coord_number_       = rhs.coord_number_;
    calc_coord_number_  = rhs.calc_coord_number_;

    delete solid_;
    solid_ = rhs.solid_->Clone();

  }
  return *this;
}

DryRock *
DryRockWalton::Clone() const {
  // Provide base class variables.
  DryRockWalton * r = new DryRockWalton(*this);

  // Provide variables specific to DryRockWalton.
  r->solid_             = this->solid_->Clone();          // Deep copy.

  r->friction_weight_   = this->friction_weight_;
  r->pressure_          = this->pressure_;
  r->coord_number_      = this->coord_number_;
  r->calc_coord_number_ = this->calc_coord_number_;

  return r;
}

void
DryRockWalton::SetTotalPorosity(double porosity) {
  DryRock::SetTotalPorosity(porosity);

  if (calc_coord_number_ == true)
    coord_number_ = CalcCoordNumber();
  ComputeElasticParams();
}

void
DryRockWalton::ComputeElasticParams() {
  static const double pi      = NRLib::Pi;
  static const double pi4_inv = 1.0/(4*pi);

  double mu_solid, k_solid, rho_solid;
  solid_->GetElasticParams(k_solid, mu_solid, rho_solid);
  mineral_moduli_k_ = k_solid;

  const double mu_solid_inv = 1.0/mu_solid;
  const double lambda       = k_solid - (2.0*mu_solid)/3.0;

  const double modulus_inv = 1.0/(mu_solid + lambda);

  double a = pi4_inv*(mu_solid_inv - modulus_inv);
  double b = pi4_inv*(mu_solid_inv + modulus_inv);

  double x = std::pow((3*coord_number_*coord_number_*(1 - total_porosity_)*(1 - total_porosity_)*pressure_)/(pi*pi*pi*pi*b*b), 1.0/3.0);

  double k_rough  = x/6.0;
  double mu_rough = (3.0/5.0)*k_rough*(5.0*b + a)/(2.0*b + a);

  double mu_smooth = x/10.0;
  //double k_smooth  = 5.0*mu_smooth/3.0; //== k_rough

  k_    = k_rough;
  mu_   = friction_weight_*mu_rough + (1.0 - friction_weight_)*mu_smooth;

  rho_  = (1-total_porosity_)*rho_solid;
}

double
DryRockWalton::CalcCoordNumber() const {
  static const double coord [] = {14.007, 12.336, 10.843, 9.5078, 8.3147, 7.2517, 6.3108, 5.4878, 4.7826, 4.1988, 3.7440};
  //int    size    = 11;
  double dp      = 0.05;
  double p_min   = 0.2;
  double p_max   = 0.7;

  if (total_porosity_ <= p_min)
    return p_min;
  else if (total_porosity_ >= p_max)
    return p_max;

  double x  = (total_porosity_ - p_min)/dp;
  int    i1 = static_cast<int>(std::floor(x));
  double t  = x - i1;

  double ret = t*coord[i1 + 1] + (1.0 - t)*coord[i1];

  return ret;

}

