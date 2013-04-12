// $Id: trendstorage.hpp 1118 2012-11-13 12:58:17Z pdahle $
#ifndef NRLIB_TRENDSTORAGE_HPP
#define NRLIB_TRENDSTORAGE_HPP

#include "trend.hpp"
#include <stdio.h>
#include <string>

namespace NRLib {
class TrendStorage {
public:
  TrendStorage();

  virtual ~TrendStorage();

  virtual TrendStorage * Clone() const = 0;

  virtual Trend * GenerateTrend(const std::string                       & /*path*/,
                                const std::vector<std::string>          & /*trend_cube_parameters*/,
                                const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                std::string                             & /*errTxt*/) const = 0;
};

//-------------------------------------------------------------------//

class TrendConstantStorage : public TrendStorage {
public:
  TrendConstantStorage(const double value);
  TrendConstantStorage();
  TrendConstantStorage(const TrendConstantStorage & trend_storage);

  virtual ~TrendConstantStorage();

  virtual TrendStorage * Clone()   const { return new TrendConstantStorage(*this) ;}
  virtual double         GetMean() const { return mean_value_                     ;}

  virtual Trend * GenerateTrend(const std::string                       & /*path*/,
                                const std::vector<std::string>          & /*trend_cube_parameters*/,
                                const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                std::string                             & /*errTxt*/) const;

private:
  double mean_value_;
};

//-------------------------------------------------------------------//

class Trend1DStorage : public TrendStorage {
public:
  Trend1DStorage(const std::string file_name,
                 const std::string reference_parameter);
  Trend1DStorage(const Trend1DStorage & trend_storage);

  Trend1DStorage();

  virtual ~Trend1DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend1DStorage(*this) ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                std::string                             & errTxt) const;

private:

  std::string file_name_;
  std::string reference_parameter_;
};

//-------------------------------------------------------------------//

class Trend2DStorage : public TrendStorage {
public:
  Trend2DStorage(const std::string file_name,
                 const std::string reference_parameter1,
                 const std::string reference_parameter2);

  Trend2DStorage(const Trend2DStorage & trend_storage);

  virtual ~Trend2DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend2DStorage(*this) ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                std::string                             & errTxt)   const;

private:
  std::string file_name_;
  std::string reference_parameter_one_;
  std::string reference_parameter_two_;
};
}
#endif
