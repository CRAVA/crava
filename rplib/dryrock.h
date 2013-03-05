#ifndef RPLIB_DRY_ROCK_H
#define RPLIB_DRY_ROCK_H

#include <vector>

// Abstract DryRock class.
class DryRock {
public:

  DryRock(){}
  virtual ~DryRock() {}

  virtual DryRock                   * Clone()                                                                             const = 0;

  void                                GetElasticParams(double & k, double & mu, double & rho)  const {
                                        k = k_; mu  = mu_; rho = rho_;
                                      }

  double                              GetMineralModuliK() const   { return mineral_moduli_k_; }
  double                              GetTotalPorosity()  const   { return total_porosity_; }
  virtual void                        SetTotalPorosity(double porosity)  { total_porosity_ = porosity; }

  const std::vector<double>&          GetU()                                                                              const { return u_; }

protected:
  double                      k_;
  double                      mu_;
  double                      rho_;

  double                      total_porosity_;
  double                      mineral_moduli_k_;

  std::vector<double>         u_;

};

#endif
