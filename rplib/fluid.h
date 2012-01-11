#ifndef FLUID_H
#define FLUID_H

#include <string>
#include <vector>


class Fluid {
public:

  Fluid(const std::string& name, double k, double rho);
  Fluid(const Fluid & rhs);
  ~Fluid();
  Fluid& operator=(const Fluid& rhs);

  const std::string         & GetName()        const {return name_;}
  double                      GetBulkModulus() const {return elastics_[0];}
  double                      GetDensity()     const {return elastics_[1];}
  const std::vector<double> & GetElastics()    const {return elastics_;}


private:
  std::string name_;// Fluid identifier.

  // Elastic properties assumed to be deterministic.
  // Replace vector type double with NRLib::Distribution if changing to stochastic.
  // If stochastic, the Get-functions must change return values,
  // if we want to give access to the whole distribution.
  // The return value should then be a const reference.
  // Example: NRLib::Distribution & const GetBulkModulus() const {return elastics_[0];}
  // Sampled values from the distributions would typically also be requested.
  // Example: double GetBulkModulusSample() const;

  std::vector<double> elastics_; // elastics_[0] = bulk modulus, elastics_[1] = density.

};

#endif
