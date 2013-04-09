#ifndef RPLIB_ROCK_H
#define RPLIB_ROCK_H

#include <vector>


// Abstract rock class.
class Rock {
public:
                                    Rock();
  virtual                           ~Rock();

  virtual Rock                    * Clone()                                                      const = 0;

  void                              GetSeismicParams(double & vp, double & vs, double & rho)     const {
                                      vp = vp_; vs = vs_; rho = rho_;
                                    }

  // Important:SetPorosity is not a "simple" set function. It triggers a recalculation of a lot of the member variables.
  // See children classes for examples.
  virtual void                      SetPorosity(double porosity)                                        = 0;

  const std::vector<double>       & GetU()                                                        const { return u_; }


protected:
  double                            vp_;
  double                            vs_;
  double                            rho_;
  std::vector<double>               u_;
};

#endif
