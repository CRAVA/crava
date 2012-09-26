#ifndef RPLIB_DISTRIBUTIONS_ROCK_TABULATED_VELOCITY_H
#define RPLIB_DISTRIBUTIONS_ROCK_TABULATED_VELOCITY_H

#include "rplib/rock.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"

// Abstract class for holding all t = 0 distribution functions for rock parameters.
// One derived class for each rock model, the latter specified in a parallel, derived Rock class.
// The class must be able to produce an object of the specific Rock class.
class DistributionsRockTabulatedVelocity : public DistributionsRock {
public:

  //NB: Class is not completed
  DistributionsRockTabulatedVelocity(const DistributionWithTrend * vp,
                                     const DistributionWithTrend * vs,
                                     const DistributionWithTrend * density,
                                     double                        corr_vp_vs,
                                     double                        corr_vp_density,
                                     double                        corr_vs_density);

  virtual ~DistributionsRockTabulatedVelocity();

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual std::vector<double>        GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>      GetCovariance(const std::vector<double> & trend_params) const;

  virtual Pdf3D                    * GeneratePdf() const;

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const; //dummy function that needs to be implemented

private:
  const DistributionWithTrend * vp_;
  const DistributionWithTrend * vs_;
  const DistributionWithTrend * density_;
  double                        corr_vp_vs_;
  double                        corr_vp_density_;
  double                        corr_vs_density_;
  Tabulated                     tabulated_;
};

#endif
