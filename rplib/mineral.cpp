#include "rplib/mineral.h"

Mineral::Mineral(std::string name, double k, double g, double rho) 
: name_(name)
{
  elastics_.resize(0);
  elastics_.push_back(k);
  elastics_.push_back(g);
  elastics_.push_back(rho);
}

Mineral::Mineral(const Mineral & rhs)
  : name_(rhs.name_)
{
  elastics_.resize(0);
  elastics_.push_back(rhs.elastics_[0]);
  elastics_.push_back(rhs.elastics_[1]);
  elastics_.push_back(rhs.elastics_[2]);
}

Mineral::~Mineral() {}

Mineral& Mineral::operator=(const Mineral& rhs) {
  if (this != &rhs) {
    name_ = rhs.name_;
    elastics_.resize(0);
    elastics_.push_back(rhs.elastics_[0]);
    elastics_.push_back(rhs.elastics_[1]);
    elastics_.push_back(rhs.elastics_[2]);
  }
  return *this;
}





