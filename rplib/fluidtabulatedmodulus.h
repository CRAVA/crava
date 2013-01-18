#ifndef RPLIB_FLUID_TABULATED_MODULUS_H
#define RPLIB_FLUID_TABULATED_MODULUS_H

#include "rplib/fluid.h"

class FluidTabulatedModulus : public Fluid {
public:

  FluidTabulatedModulus(const double                    k,
                        const double                    rho,
                        std::vector<double>             u);

  // Copy constructor
  FluidTabulatedModulus(const FluidTabulatedModulus & rhs);

  virtual ~FluidTabulatedModulus();

  virtual Fluid * Clone() const;


};

#endif
