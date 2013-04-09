#include "rplib/dem.h"

#include "nrlib/exception/exception.hpp"

#include "rplib/orddiffeqsolver.h"

#include <numeric>
#include <cmath>

static DEM* global_dem;

static std::vector<double> WrapperGEQDEMYPrime(std::vector<double>&       y,
                                               double                     t) {
  return global_dem->GEQDEMYPrime(y, t);
}


DEM::DEM(const std::vector<double>&       bulk_modulus,
         const std::vector<double>&       shear_modulus,
         const std::vector<double>&       aspect_ratio,
         std::vector<double>&             concentration,
         double                           bulk_modulus_bg,
         double                           shear_modulus_bg) :
  bulk_modulus_bg_(bulk_modulus_bg),
  shear_modulus_bg_(shear_modulus_bg),
  bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  aspect_ratio_(aspect_ratio),
  concentration_(concentration) {

  global_dem = this;

}

DEM::~DEM() {

}

void
DEM::CalcEffectiveModulus(double&                    effective_bulk_modulus,
                          double&                    effective_shear_modulus) {

  effective_bulk_modulus = effective_shear_modulus = 0;

  /*
  In case of multiple inclusions replacing all of the host, the DEM
  algorithm gives wrong result. This is avoided by leaving a very small
  part of the host left, and replacing it with 0.999999 of the
  inclusion material.
  */
  // calculate sum of concentration_
  std::vector<double>::iterator it = concentration_.begin();
  double sum_conc = 0;
  for (; it < concentration_.end(); ++it)
    sum_conc += *it;

  if (concentration_.size() > 1 && sum_conc == 1.0) {
    double repair_factor = 0.999999;
    it = concentration_.begin();
    for (; it < concentration_.end(); ++it)
      *it *= repair_factor;

    //update sum after rescaling
    sum_conc *= repair_factor;
  }

  //double phic_ = 1.0; //Not used

  if (sum_conc == 1) {
    //warning("Assumes all inclusions have the same elastic moduli");
    effective_bulk_modulus  = bulk_modulus_.back();
    effective_shear_modulus = shear_modulus_.back();

  }
  else {
    std::vector<double>           y0;
    std::vector<double>                 tout;
    std::vector< std::vector<double> >  yout;

    double tfinal = sum_conc;
    y0.resize(2);
    y0[0] = bulk_modulus_bg_;
    y0[1] = shear_modulus_bg_;

    OrdDiffEqSolver::
    Ode45(&WrapperGEQDEMYPrime,
          0.0,
          tfinal,
          y0,
          tout,
          yout,
          1e-5);

    if (!yout.empty())
      effective_bulk_modulus  = (yout.back())[0];
    if (!yout.empty())
      effective_shear_modulus = (yout.back())[1];

  }

}

std::vector<double>
DEM::GEQDEMYPrime(std::vector<double>&       y,
                  double                     t) {

  size_t ninclusions = aspect_ratio_.size();

  double krhs = 0;
  double murhs = 0;
  double sum_conc = std::accumulate(concentration_.begin(), concentration_.end(), 0.0);

  std::vector<double> yprime(2, 0.0);

  for (size_t index = 0; index < ninclusions; index++) {
    double k2 = bulk_modulus_[index];
    double mu2 = shear_modulus_[index];
    double asp = aspect_ratio_[index];

    double conc = concentration_[index]/sum_conc;

    //double krc = bulk_modulus_bg_*k2/((1 - phic_)*k2 + phic_*bulk_modulus_bg_); //Not used
    //double murc = shear_modulus_bg_*mu2/((1 - phic_)*mu2 + phic_*shear_modulus_bg_); //Not used

    double ka = k2;
    double mua = mu2;


    double k = y[0];
    double mu = y[1];

    // truncation
    if (asp == 1.0)
      asp = 0.99;

    //******* P and Q *****************
    double theta = 0.0;
    double fn = 0.0;

    if (asp < 1.0) {
      theta = (asp/(pow((1 - asp*asp), 3.0/2.0)))*(acos(asp) - asp*sqrt(1 - asp*asp));
      fn = ((asp*asp)/(1 - asp*asp))*(3*theta -2);
    }
    else {
      //theta = (asp/(pow((asp*asp - 1), 3.0/2.0)))*(asp*sqrt(asp*asp-1)-acosh(asp)); //bug in original code???
      fn=((asp*asp)/(asp*asp - 1))*(2-3*theta);
      throw NRLib::Exception("DEM: asp > 1 not supported.");
    }

    double nu = (3*k - 2*mu)/(2*(3*k + mu));
    double r = (1 - 2*nu)/(2*(1 - nu));
    double a = mua/mu - 1;
    double b = (1.0/3.0)*(ka/k - mua/mu);

    double f1a = 1 + a*((3.0/2.0)*(fn + theta)-r*((3.0/2.0)*fn + (5.0/2.0)*theta - (4.0/3.0)));

    double f2a = 1 + a*(1+(3.0/2.0)*(fn + theta)-(r/2)*(3*fn + 5*theta)) + b*(3-4*r);
    f2a += (a/2)*(a + 3*b)*(3-4*r)*(fn + theta-r*(fn-theta + 2*theta*theta));

    double f3a = 1 + a*(1-(fn+(3.0/2.0)*theta) + r*(fn+theta));

    double f4a = 1 + (a/4)*(fn+3*theta - r*(fn - theta));

    double f5a = a*(-fn + r*(fn+theta-(4.0/3.0))) + b*theta*(3 - 4*r);

    double f6a = 1 + a*(1 + fn - r*(fn + theta)) + b*(1-theta)*(3 - 4*r);

    double f7a = 2 + (a/4)*(3*fn + 9*theta - r*(3*fn + 5*theta)) + b*theta*(3 - 4*r);

    double f8a = a*(1 - 2*r + (fn/2)*(r - 1)+(theta/2)*(5*r - 3)) + b*(1 - theta)*(3 - 4*r);

    double f9a = a*((r - 1)*fn - r*theta) + b*theta*(3 - 4*r);

    double pa = 3*f1a/f2a;
    double qa = (2/f3a) + (1/f4a) +((f4a*f5a + f6a*f7a - f8a*f9a)/(f2a*f4a));

    pa = pa/3.0;
    qa = qa/5.0;

    krhs += conc*(ka - k)*pa;
    murhs += conc*(mua - mu)*qa;
  }

  yprime[0] = krhs/(1 - t);
  yprime[1] = murhs/(1 - t);

  return yprime;

}




