#include "rplib/mineral.h"

Mineral::Mineral(const std::string& name, double k, double g, double rho) 
: name_(name)
{
  elastics_.push_back(k);
  elastics_.push_back(g);
  elastics_.push_back(rho);
}

Mineral::Mineral(const Mineral & rhs)
  : name_(rhs.name_), elastics_(rhs.elastics_)
{
  
}

Mineral::~Mineral() {}

Mineral& Mineral::operator=(const Mineral& rhs) {
  if (this != &rhs) {
    name_ = rhs.name_;
    elastics_ = rhs.elastics_;
  }
  return *this;
}





