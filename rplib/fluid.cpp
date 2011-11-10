#include "rplib/fluid.h"

Fluid::Fluid(std::string name, double k, double rho) 
: name_(name)
{
  elastics_.resize(0);
  elastics_.push_back(k);
  elastics_.push_back(rho);
}

Fluid::Fluid(const Fluid & rhs)
  : name_(rhs.name_)
{
  elastics_.resize(0);
  elastics_.push_back(rhs.elastics_[0]);
  elastics_.push_back(rhs.elastics_[1]);
}

Fluid::~Fluid() {}

Fluid& Fluid::operator=(const Fluid& rhs) {
  if (this != &rhs) {
    name_ = rhs.name_;
    elastics_.resize(0);
    elastics_.push_back(rhs.elastics_[0]);
    elastics_.push_back(rhs.elastics_[1]);
  }
  return *this;
}



