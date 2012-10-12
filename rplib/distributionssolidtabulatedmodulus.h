#ifndef DISTRIBUTIONS_SOLID_TABULATED_MODULUS_H
#define DISTRIBUTIONS_SOLID_TABULATED_MODULUS_H

#include "rplib/distributionssolid.h"
#include "rplib/solidtabulatedmodulus.h"
#include "rplib/tabulated.h"

class DistributionWithTrend;

class DistributionsSolidTabulatedModulus : public DistributionsSolid {
public:

  DistributionsSolidTabulatedModulus(const DistributionWithTrend         * distr_k,
                                     const DistributionWithTrend         * distr_mu,
                                     const DistributionWithTrend         * distr_rho,
                                     const double                          corr_k_mu,
                                     const double                          corr_k_rho,
                                     const double                          corr_mu_rho);

  virtual                               ~DistributionsSolidTabulatedModulus();

  virtual Solid                       * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                          HasDistribution() const;

  virtual std::vector<bool>             HasTrend() const;

protected:
  virtual Solid *                       UpdateSample(const std::vector< double > & /*corr*/,
                                                     const Solid                 & /*solid*/) const;

private:

  Solid                               * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend         * distr_k_;         // Pointer to external object.
  const DistributionWithTrend         * distr_mu_;        // Pointer to external object.
  const DistributionWithTrend         * distr_rho_;       // Pointer to external object.
  const double                          corr_k_mu_;
  const double                          corr_k_rho_;
  const double                          corr_mu_rho_;
  Tabulated                           * tabulated_;
};

#endif
