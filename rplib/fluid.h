#ifndef FLUID_H
#define FLUID_H

#include <string>
#include <vector>


class Fluid {
public:

  Fluid(std::string name, double k, double rho);
  Fluid(const Fluid & rhs);
  ~Fluid();
  Fluid& operator=(const Fluid& rhs);

  std::string                GetName()        const {return name_;}
  double                     GetBulkModulus() const {return k_;}
  double                     GetDensity()     const {return rho_;}
  const std::vector<double*> GetElastics()    const {return elastics_;}


private:
	std::string name_;		// Fluid identifier.

	// Elastic properties assumed to be deterministic. 
	// Replace type with NRLib::Distribution if changing to stochastic. 
	// If stochastic, the Get-functions must change return values, 
	// if we want to give access to the whole distribution.
	// The return value should then be a const reference.
	// Example: NRLib::Distribution & const GetBulkModulus() const {return k_;}
	// Sampled values from the distributions would typically also be requested.
	// Example: double GetBulkModulusSample() const;
	double k_;              // Bulk modulus
	double rho_;            // Density
  std::vector<double*> elastics_; // Points to k_ and rho_;

};

#endif
