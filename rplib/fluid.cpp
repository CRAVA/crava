#include "rplib/fluid.h"

Fluid::Fluid(const std::string& name, double k, double rho) 
: name_(name)
{
  elastics_.push_back(k);
  elastics_.push_back(rho);
}

Fluid::Fluid(const Fluid & rhs)
  : name_(rhs.name_), elastics_(rhs.elastics_)
{
 
}

Fluid::~Fluid() {}

Fluid& Fluid::operator=(const Fluid& rhs) {
  if (this != &rhs) {
    name_ = rhs.name_;
    elastics_ = rhs.elastics_;
  }
  return *this;
}



