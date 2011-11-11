#ifndef MINERAL_H
#define MINERAL_H

#include <string>
#include <vector>


class Mineral {
public:

  Mineral(const std::string& name, double k, double g, double rho);
  Mineral(const Mineral & rhs);
  ~Mineral();
  Mineral& operator=(const Mineral& rhs);

  const std::string         & GetName()         const {return name_;}
  double                      GetBulkModulus()  const {return elastics_[0];}
  double                      GetShearModulus() const {return elastics_[1];}
  double                      GetDensity()      const {return elastics_[2];}
  const std::vector<double> & GetElastics()     const {return elastics_;}

private:
	std::string name_;		// Mineral identifier.

	// Elastic properties assumed to be deterministic. 
	// Replace vector type double with NRLib::Distribution if changing to stochastic. 
	// If stochastic, the Get-functions must change return values, 
	// if we want to give access to the whole distribution.
	// The return value should then be a const reference.
	// Example: NRLib::Distribution & const GetBulkModulus() const {return elastics_[0];}
	// Sampled values from the distributions would typically also be requested.
	// Example: double GetBulkModulusSample() const;
	
  std::vector<double> elastics_; // elastics_[0] = bulk modulus, elastics_[1] = shear modulus, elastics_[2] = density.

};

#endif
