#ifndef RPLIB_DISTRIBUTIONS_ROCK_TABULATED_MODULUS_H
#define RPLIB_DISTRIBUTIONS_ROCK_TABULATED_MODULUS_H

#include "rplib/rock.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"

class DistributionsRockTabulatedModulus : public DistributionsRock {
public:

  DistributionsRockTabulatedModulus(const DistributionWithTrend * bulk_modulus,
                                    const DistributionWithTrend * shear_modulus,
                                    const DistributionWithTrend * density,
                                    double                        corr_bulk_shear,
                                    double                        corr_bulk_density,
                                    double                        corr_shear_density);

  virtual ~DistributionsRockTabulatedModulus();

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                     * UpdateSample(const std::vector<double> & /*corr*/,
                                                  const Rock                & /*rock*/)       const { return NULL; }

  virtual std::vector<double>        GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>      GetCovariance(const std::vector<double> & trend_params)  const;

  virtual Pdf3D                    * GeneratePdf() const;

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const;

  virtual bool                       GetIsOkForBounding() const;

private:

  Rock                             * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend * bulk_modulus_;
  const DistributionWithTrend * shear_modulus_;
  const DistributionWithTrend * density_;
  double                        corr_bulk_shear_;
  double                        corr_bulk_density_;
  double                        corr_shear_density_;
  Tabulated                   * tabulated_;
  bool                          has_distribution_;
  std::vector<bool>             has_trend_;
};

#endif
