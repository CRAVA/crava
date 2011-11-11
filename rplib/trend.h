#ifndef RPLIB_TREND_H
#define RPLIB_TREND_H

#include <vector>

class Trend {
public:
  Trend();
  virtual ~Trend();
  virtual double GetValue(int /*s1*/, int /*s2*/) const = 0;

};

class ConstantTrend : public Trend {
public:
  ConstantTrend(double trend);
  virtual ~ConstantTrend();
  virtual double GetValue(int /*s1*/, int /*s2*/) const {return trend_;}

private:
  double trend_;
  

};

class Trend1D : public Trend {
public:
  Trend1D(const std::vector<double>& trend);
  virtual ~Trend1D();
  virtual double GetValue(int s1, int /*s2*/) const {return trend_[s1];}

private:
  std::vector<double> trend_;
};

class Trend2D : public Trend {
public:
  Trend2D(const std::vector<double>& trend, int ns1);
  virtual ~Trend2D();
  virtual double GetValue(int s1, int s2) const {return trend_[s1 + s2*ns1_];}

private:
  std::vector<double> trend_;
  int                 ns1_;

};
#endif
