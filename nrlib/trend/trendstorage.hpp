#ifndef NRLIB_TREND_STORAGE_HPP
#define NRLIB_TREND_STORAGE_HPP

#include "trend.hpp"
#include <stdio.h>
#include <string>

class TrendStorage {
public:
  TrendStorage();

  virtual ~TrendStorage();

  virtual Trend * GenerateTrend(const std::string & /*path*/,
                                std::string       & /*errTxt*/) const = 0;
};

class TrendConstantStorage : public TrendStorage {
public:
  TrendConstantStorage(const double value);

  virtual ~TrendConstantStorage();

  virtual Trend * GenerateTrend(const std::string & /*path*/,
                                std::string       & /*errTxt*/) const;

private:
  double mean_value_;
};

class Trend1DStorage : public TrendStorage {
public:
  Trend1DStorage(const std::string file_name,
                 const std::string reference_parameter);

  virtual ~Trend1DStorage();

  virtual Trend * GenerateTrend(const std::string & path,
                                std::string       & errTxt) const;

private:
  void             readTrend1D(const std::string   & file_name,
                               std::string         & errText,
                               std::vector<double> & trend1d,
                               double              & s_min,
                               double              & s_max)     const;

  int              getTrend1DFileFormat(const std::string & file_name,
                                        std::string       & errText)   const;

  std::string file_name_;
  std::string reference_parameter_;
};

class Trend2DStorage : public TrendStorage {
public:
  Trend2DStorage(const std::string file_name,
                 const std::string reference_parameter1,
                 const std::string reference_parameter2);

  virtual ~Trend2DStorage();

  virtual Trend * GenerateTrend(const std::string & /*path*/,
                                std::string       & errTxt)   const;

private:
  std::string file_name_;
  std::string reference_parameter_one_;
  std::string reference_parameter_two_;
};

#endif
