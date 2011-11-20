#ifndef RPLIB_TREND_H
#define RPLIB_TREND_H

#include <vector>

class Trend {
public:
  Trend();
  virtual ~Trend();
  virtual double GetValue(double /*s1*/, double /*s2*/) const = 0;

};

class ConstantTrend : public Trend {
public:
  ConstantTrend(double trend);
  virtual ~ConstantTrend();
  virtual double GetValue(double /*s1*/, double /*s2*/) const {return trend_;}

private:
  double trend_;
  

};

class Trend1D : public Trend {
public:
  Trend1D(const std::vector<double>& trend, double s1_min, double s1_max); //If error, throws NRLib::Exception
  virtual ~Trend1D();
  virtual double GetValue(double s1, double /*s2*/) const;

private:
  const std::vector<double>&  trend_;
  double                      s1_min_, s1_max_;
  double                      inv_s1_inc_;

};

class Trend2D : public Trend {
public:
  Trend2D(const std::vector<double>& trend, 
          double s1_min, double s1_max, int ns1,
          double s2_min, double s2_max, int ns2); //If error, throws NRLib::Exception
  virtual ~Trend2D();
  virtual double GetValue(double s1, double s2) const;

private:
  const std::vector<double>&    trend_;
  int                           ns1_, ns2_;

  double                        s1_min_, s1_max_;
  double                        inv_s1_inc_;

  double                        s2_min_, s2_max_;
  double                        inv_s2_inc_;

};
#endif
