#include "rplib/fluid.h"

Fluid::Fluid(std::string name, double k, double rho) 
: name_(name), k_(k), rho_(rho)
{
  elastics_.push_back(&k_);
  elastics_.push_back(&rho_);
}

Fluid::Fluid(const Fluid & rhs)
  : name_(rhs.name_), 
    k_(rhs.k_),
    rho_(rhs.rho_)
{
  elastics_.push_back(&k_);
  elastics_.push_back(&rho_);
}

Fluid::~Fluid() {}

Fluid& Fluid::operator=(const Fluid& rhs) {
  if (this != &rhs) {
    name_ = rhs.name_;
    k_ = rhs.k_;
    rho_ = rhs.rho_;

    elastics_.push_back(&k_);
    elastics_.push_back(&rho_); 
  }
  return *this;
}



