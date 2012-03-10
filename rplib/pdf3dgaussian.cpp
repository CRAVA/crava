#include "rplib/pdf3dgaussian.h"
#include "rplib/trinormalwith2dtrend.h"

Pdf3DGaussian::Pdf3DGaussian(TriNormalWith2DTrend * tri_normal_distr)
{
  tri_normal_distr_ = tri_normal_distr;
}

Pdf3DGaussian::~Pdf3DGaussian()
{
}

double
Pdf3DGaussian::density(const double & vp,
                       const double & vs,
                       const double & rho,
                       const double & s1,
                       const double & s2) const
{
  double prob;

  tri_normal_distr_->CalculatePDF(s1, s2, vp, vs, rho, prob);

  return(prob);
}
