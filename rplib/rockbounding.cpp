#include "rplib/rockbounding.h"
#include "rplib/demmodelling.h"

RockBounding::RockBounding(const Rock          * upper_rock,
                           const Rock          * lower_rock,
                           double                porosity,
                           double                K_weight,
                           double                M_weight,
                           std::vector<double>   u)
: Rock()
{
  upper_rock_ = upper_rock->Clone();
  lower_rock_ = lower_rock->Clone();
  porosity_   = porosity;
  K_weight_   = K_weight;
  M_weight_   = M_weight;
  u_          = u;               // u contains correlated samples used in quantiles of (porosity, K_weight, M_weight)

  ComputeSeismicVariables();
}

RockBounding::RockBounding()
: Rock()
{
}

RockBounding::RockBounding(const RockBounding & old_rock)
: Rock(old_rock)
{
  porosity_ = old_rock.porosity_;
  K_weight_ = old_rock.K_weight_;
  M_weight_ = old_rock.M_weight_;
  vp_       = old_rock.vp_;
  vs_       = old_rock.vs_;
  rho_      = old_rock.rho_;
  u_        = old_rock.u_;

  delete upper_rock_;
  delete lower_rock_;

  upper_rock_ = old_rock.upper_rock_->Clone();
  lower_rock_ = old_rock.lower_rock_->Clone();

}


RockBounding::~RockBounding()
{
  delete upper_rock_;
  delete lower_rock_;
}


Rock *
RockBounding::Clone() const
{
  return new RockBounding(*this);
}


Rock *
RockBounding::Evolve(const std::vector<int>         & /*delta_time*/,
                     const std::vector< Rock * >    & /*rock*/) const
{
  return new RockBounding(*this);
}

void
RockBounding::SetPorosity(double porosity)
{
  porosity_ = porosity;

  ComputeSeismicVariables();
}

void
RockBounding::ComputeSeismicVariables()
{

  upper_rock_->SetPorosity(porosity_);
  lower_rock_->SetPorosity(porosity_);

  double vp_upper;
  double vs_upper;
  double density_upper;
  upper_rock_->GetSeismicParams(vp_upper, vs_upper, density_upper);

  double K_upper;
  double G_upper;
  DEMTools::CalcElasticParamsFromSeismicParams(vp_upper, vs_upper, density_upper, K_upper, G_upper);

  double M_upper = K_upper + 4.0/3.0 * G_upper;

  double vp_lower;
  double vs_lower;
  double density_lower;
  lower_rock_->GetSeismicParams(vp_lower, vs_lower, density_lower);

  double K_lower;
  double G_lower;
  DEMTools::CalcElasticParamsFromSeismicParams(vp_lower, vs_lower, density_lower, K_lower, G_lower);

  double M_lower = K_lower + 4.0/3.0 * G_lower;

  double K = K_weight_ * K_upper       + (1-K_weight_) * K_lower;
  double M = M_weight_ * M_upper       + (1-M_weight_) * M_lower;
  rho_     = 0.5       * density_upper +  0.5          * density_lower;

  double G  = (M - K) * 3.0/4.0;

  DEMTools::CalcSeismicParamsFromElasticParams(K, G, rho_, vp_, vs_);
}
