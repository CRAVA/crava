#ifndef RPLIB_SOLID_H
#define RPLIB_SOLID_H

#include <vector>

// Abstract solid class.
class Solid {
public:

                             Solid(){}
  virtual                    ~Solid() {}

  virtual Solid *            Clone()                                                  const = 0;

  void                       GetElasticParams(double & k, double & mu, double & rho)  const
                             {
                                   k = k_; mu  = mu_; rho = rho_;
                             }

  const std::vector<double>& GetU()                                                   const { return u_; }

protected:
  double                      k_;
  double                      mu_;
  double                      rho_;
  std::vector<double>         u_;

};

#endif
