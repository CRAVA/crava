#ifndef RPLIB_DRYROCK_TABULATED_MODULUS_H
#define RPLIB_DRYROCK_TABULATED_MODULUS_H

#include "rplib/dryrock.h"

class DryRockTabulatedModulus : public DryRock {
public:

  DryRockTabulatedModulus(const double                    k,
                          const double                    mu,
                          const double                    rho,
                          const double                    total_porosity,
                          const double                    mineral_moduli_k,
                          std::vector<double>             u);


  // Copy constructor
  DryRockTabulatedModulus(const DryRockTabulatedModulus & rhs);

  virtual ~DryRockTabulatedModulus(){}

  // Assignment operator, not yet implemented.
  /*DryRockTabulatedModulus           & operator=(const DryRockTabulatedModulus& rhs);*/

  virtual DryRock                     * Clone() const;

  virtual void                          ComputeElasticParams(double & k, double & mu, double & rho) const;

  virtual DryRock                     * Evolve(const std::vector<int>               & /*delta_time*/,
                                               const std::vector< const DryRock * > & /*fluid*/) const;

};

#endif
