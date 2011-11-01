#include "rplib/mineral.h"

Mineral::Mineral(std::string name, double k, double g, double rho) 
: name_(name), k_(k), g_(g), rho_(rho)
{
  elastics_.push_back(&k_);
  elastics_.push_back(&g_);
  elastics_.push_back(&rho_);
}

Mineral::Mineral(const Mineral & rhs)
  : name_(rhs.name_), 
    k_(rhs.k_),
    g_(rhs.g_),
    rho_(rhs.rho_)
{
  elastics_.push_back(&k_);
  elastics_.push_back(&g_);
  elastics_.push_back(&rho_);
}

Mineral::~Mineral() {}




