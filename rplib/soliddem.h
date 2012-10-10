#ifndef RPLIB_SOLIDDEM_H
#define RPLIB_SOLIDDEM_H


#include "rplib/solid.h"

#include "rplib/demmodelling.h"

#include <vector>

class SolidDEM : public Solid {
public:

  SolidDEM(const Solid                         * solid,
           const Solid                         * solid_inc,
           const std::vector<double>           & inclusion_spectrum,
           const std::vector<double>           & inclusion_concentration,
           double                                porosity,
           const std::vector<double>           & u);

  SolidDEM();

  virtual         ~SolidDEM();

  // Assignment operator.
  SolidDEM      & operator=(const SolidDEM& rhs);

  virtual Solid * Clone() const;

  const   Solid * GetSolidHost()      const {return solid_;}
  const   Solid * GetSolidInclusion() const {return solid_inc_;}

  virtual Solid * Evolve(const std::vector<int>               & delta_time,
                         const std::vector< const Solid * >   & solids) const;

  virtual void    SetPorosity(double porosity);

private:
  //Copy constructor for getting base class variables , used by Clone:
  SolidDEM(const SolidDEM & rhs) : Solid(rhs) {}

  // Calculate elastic and seismic parameters, to be
  // used whenever new information is sent to class.
  void ComputeElasticParams();

  Solid                               * solid_; // Owned and deleted by this class.
  Solid                               * solid_inc_; // Owned and deleted by this class.
  std::vector<double>                   inclusion_spectrum_;
  std::vector<double>                   inclusion_concentration_;
  double                                porosity_;
};

#endif
