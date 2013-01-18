#ifndef RPLIB_SOLIDDEM_H
#define RPLIB_SOLIDDEM_H


#include "rplib/solid.h"

#include "rplib/demmodelling.h"

#include <vector>

class SolidDEM : public Solid {
public:

  SolidDEM(const Solid                         * solid,
           const std::vector< Solid* >         & solid_inc,
           const std::vector<double>           & inclusion_spectrum,
           const std::vector<double>           & inclusion_concentration, // the first element is the concentration of the host
           const std::vector<double>           & u);

  SolidDEM();

  virtual         ~SolidDEM();

  // Assignment operator.
  SolidDEM                            & operator=(const SolidDEM& rhs);

  virtual Solid                       * Clone() const;

  const Solid                         * GetSolidHost()      const {return solid_;}

  const std::vector<Solid*>           & GetSolidInclusion() const {return solid_inc_;}
  const Solid                         * GetSolidInclusion(size_t i) const { return solid_inc_[i]; } // no error checking on valid index range

private:
  //Copy constructor for getting base class variables , used by Clone:
  SolidDEM(const SolidDEM & rhs) : Solid(rhs) {}

  // Calculate elastic and seismic parameters, to be
  // used whenever new information is sent to class.
  void                                  ComputeElasticParams();
  void                                  Clone(const std::vector< Solid* > & solid_in);
  void                                  DeleteInclusion();


  Solid                               * solid_;                 // Owned and deleted by this class.
  std::vector< Solid* >                 solid_inc_;             // Owned and deleted by this class.
  std::vector<double>                   inclusion_spectrum_;
  std::vector<double>                   inclusion_concentration_; // the first element is the concentration of the host
};

#endif
