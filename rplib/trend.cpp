#include "rplib/trend.h"

#include <cassert>
#include <cmath>

#include "nrlib/exception/exception.hpp"

Trend::Trend() {
}

Trend::~Trend() {

}


ConstantTrend::ConstantTrend(double trend) : trend_(trend) {

}

ConstantTrend::~ConstantTrend() {

}


Trend1D::Trend1D(const std::vector<double>& trend, double s1_min, double s1_max) :
  trend_(trend), s1_min_(s1_min), s1_max_(s1_max) {

  if(!(s1_max_ > s1_min_ && !trend_.empty()))
    throw NRLib::Exception("Trend1D: Input parameters are not valid.");

  inv_s1_inc_ = (trend_.size() - 1)/(s1_max_ - s1_min_);
}

Trend1D::~Trend1D() {

}

double
Trend1D::GetValue(double s1, double /*s2*/) const {
  //Linear interpolation method
  double i = (s1 - s1_min_)*inv_s1_inc_;
  int i1   = static_cast<int>(floor(i));

  double val = 0;
  if (i1 < 0)
    val = trend_.front();
  else if (i1+1 > trend_.size() - 1)
    val = trend_.back();
  else {
    double t = i - i1;
    val = t*trend_[i1+1] + (1.0 - t)*trend_[i1];
  }

  return val;

}

Trend2D::Trend2D(const std::vector<double>& trend,
                 double s1_min, double s1_max, int ns1,
                 double s2_min, double s2_max, int ns2) :
 trend_(trend), s1_min_(s1_min), s1_max_(s1_max), ns1_(ns1),
 s2_min_(s2_min), s2_max_(s2_max), ns2_(ns2)
{
  if(!(s1_max_ > s1_min_ && !trend_.empty() &&
      s2_max_ > s2_min_ && trend_.size() == ns1*ns2))
    throw NRLib::Exception("Trend2D: Input parameters are not valid.");

  inv_s1_inc_ = (ns1_ - 1)/(s1_max_ - s1_min_);

  inv_s2_inc_ = (ns2_ - 1)/(s2_max_ - s2_min_);


}

double
Trend2D::GetValue(double s1, double s2) const {

  //Bilinear interpolation method
  double i = (s1 - s1_min_)*inv_s1_inc_;
  int i1   = static_cast<int>(floor(i));

  double j = (s2 - s2_min_)*inv_s2_inc_;
  int j1   = static_cast<int>(floor(j));

  double t1 = i - i1;
  double t2 = j - j1;

  double val = 0;
  if (i1 < 0) {
    if (j1 < 0)
      val = trend_[0];
    else if (j1 + 1 > ns2_ - 1)
      val = trend_[0 + (ns2_ - 1)*ns1_];
    else
      val = t2*trend_[0 + (j1+1)*ns1_] + (1.0 - t2)*trend_[0 + j1*ns1_];
  }
  else if (i1 + 1 > ns1_ - 1) {
    if (j1 < 0)
      val = trend_[ns1_ - 1 + 0];
    else if (j1 + 1 > ns2_ - 1)
      val = trend_[ns1_ - 1 + (ns2_ - 1)*ns1_];
    else
      val = t2*trend_[ns1_ - 1 + (j1+1)*ns1_] + (1.0 - t2)*trend_[ns1_ - 1 + j1*ns1_];
  }
  else if (j1 < 0)// i1 is inside
    val = t1*trend_[i1+1 + 0] + (1.0 - t1)*trend_[i1 + 0];
  else if (j1 + 1 > ns2_ - 1)
    val = t1*trend_[i1+1 + (ns2_-1)*ns1_] + (1.0 - t1)*trend_[i1 + (ns2_-1)*ns1_];
  else { //i1 AND j1 are inside
    double val_s1_1 = t1*trend_[i1+1 + j1*ns1_] + (1.0 - t1)*trend_[i1 + j1*ns1_];
    double val_s1_2 = t1*trend_[i1+1 + (j1+1)*ns1_] + (1.0 - t1)*trend_[i1 + (j1+1)*ns1_];

    val = t2*val_s1_2 + (1-t2)*val_s1_1;
  }

  return val;

}

Trend2D::~Trend2D() {

}

