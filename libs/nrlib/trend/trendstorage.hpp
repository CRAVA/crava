// $Id: trendstorage.hpp 1166 2013-05-03 11:34:58Z ulvmoen $
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
                                const std::vector<std::vector<float> >  & /*blocked_logs*/,
                                std::string                             & /*errTxt*/) const = 0;

  #define RMISSING -99999.000
  #define IMISSING -99999
};

//-------------------------------------------------------------------//

class TrendConstantStorage : public TrendStorage {
public:
  TrendConstantStorage(const double & value,
                       const bool   & estimate);

  TrendConstantStorage(const TrendConstantStorage & trend_storage);

  TrendConstantStorage();

  virtual ~TrendConstantStorage();

  virtual TrendStorage * Clone()   const { return new TrendConstantStorage(*this) ;}
  virtual double         GetMean() const { return mean_value_                     ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<float> >  & blocked_logs,
                                std::string                             & errTxt) const;

private:
  double mean_value_;
  bool   estimate_;
};

//-------------------------------------------------------------------//

class Trend1DStorage : public TrendStorage {
public:
  Trend1DStorage(const std::string & file_name,
                 const std::string & reference_parameter,
                 const bool        & estimate);

  Trend1DStorage(const Trend1DStorage & trend_storage);

  Trend1DStorage();

  virtual ~Trend1DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend1DStorage(*this) ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<float> >  & blocked_logs,
                                std::string                             & errTxt) const;

private:

  void            Estimate1DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                                  std::vector<double>                     & trend) const;

  void            SmoothTrendWithLocalLinearRegression(std::vector<double>    & trend,
                                                       const std::vector<int> & count,
                                                       const int              & iWells) const;
  std::string file_name_;
  std::string reference_parameter_;
  bool        estimate_;
};

//-------------------------------------------------------------------//

class Trend2DStorage : public TrendStorage {
public:
  Trend2DStorage(const std::string & file_name,
                 const std::string & reference_parameter1,
                 const std::string & reference_parameter2,
                 const bool        & estimate);

  Trend2DStorage(const Trend2DStorage & trend_storage);

  virtual ~Trend2DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend2DStorage(*this) ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<float> >  & blocked_logs,
                                std::string                             & errTxt)   const;

private:
  std::string file_name_;
  std::string reference_parameter_one_;
  std::string reference_parameter_two_;
  bool        estimate_;
};

}
#endif
