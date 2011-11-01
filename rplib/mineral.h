#ifndef MINERAL_H
#define MINERAL_H

#include <string>
#include <vector>


class Mineral {
public:

  Mineral(std::string name, double k, double g, double rho);
  Mineral(const Mineral & rhs);

  ~Mineral();

  std::string                GetName()         const {return name_;}
  double                     GetBulkModulus()  const {return k_;}
  double                     GetShearModulus() const {return g_;}
  double                     GetDensity()      const {return rho_;}
  const std::vector<double*> GetElastics()     const {return elastics_;}

private:
	std::string name_;		// Mineral identifier.

	// Elastic properties assumed to be deterministic. 
	// Replace type with NRLib::Distribution if changing to stochastic. 
	// If stochastic, the Get-functions must change return values, 
	// if we want to give access to the whole distribution.
	// The return value should then be a const reference.
	// Example: NRLib::Distribution & const GetBulkModulus() const {return k_;}
	// Sampled values from the distributions would typically also be requested.
	// Example: double GetBulkModulusSample() const;
	double k_;                      // Bulk modulus
	double g_;                      // Shear modulus
	double rho_;                    // Density
  std::vector<double*> elastics_; // Points to k_, g_, and rho_;

};

#endif
