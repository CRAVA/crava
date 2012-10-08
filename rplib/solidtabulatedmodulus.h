#ifndef RPLIB_SOLID_TABULATED_MODULUS_H
#define RPLIB_SOLID_TABULATED_MODULUS_H

#include "rplib/solid.h"

class SolidTabulatedModulus : public Solid {
public:

  SolidTabulatedModulus(const double                    k,
                        const double                    mu,
                        const double                    rho,
                        std::vector<double>             u);


  // Copy constructor
  SolidTabulatedModulus(const SolidTabulatedModulus & rhs);

  virtual ~SolidTabulatedModulus(){}

  // Assignment operator, not yet implemented.
  /*SolidTabulatedModulus           & operator=(const SolidTabulatedModulus& rhs);*/

  virtual Solid                     * Clone() const;

  virtual void                        ComputeElasticParams(double & k, double & mu, double & rho) const;

  virtual Solid                     * Evolve(const std::vector<int>             & /*delta_time*/,
                                             const std::vector< const Solid * > & /*fluid*/) const;

};

#endif
