#ifndef RPLIB_ROCK_H
#define RPLIB_ROCK_H

#include <assert.h>
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

                                    // Input parameters:
                                    //      delta_time : the set of previous and present incremental time steps
                                    //      rock : the set of previous rock samples
                                    // Recommended in implementation: assert(delta_time.size() == rock.size() + 1);
  virtual                           Rock * Evolve(const std::vector<int>         & delta_time,
                                                  const std::vector< Rock * >    & rock)          const = 0;

  virtual void                      SetPorosity(double porosity)                                        = 0;

  const std::vector<double>       & GetU()                                                        const { return u_; }

protected:
  double                            vp_;
  double                            vs_;
  double                            rho_;
  std::vector<double>               u_;

};

#endif
