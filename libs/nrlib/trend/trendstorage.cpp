// $Id: trendstorage.cpp 1035 2012-06-25 11:02:08Z ulvmoen $
#include "trendstorage.hpp"
#include "trendkit.hpp"
#include "trend.hpp"

#include "../surface/regularsurface.hpp"
#include "../iotools/fileio.hpp"

#include <fstream>

using namespace NRLib;

TrendStorage::TrendStorage()
{
}

TrendStorage::~TrendStorage()
{
}

//-------------------------------------------------------------------//

TrendConstantStorage::TrendConstantStorage(const double value)
: mean_value_(value)
{
}

TrendConstantStorage::TrendConstantStorage()
{
}

TrendConstantStorage::TrendConstantStorage(const TrendConstantStorage & trend_storage)
: mean_value_(trend_storage.mean_value_)
{
}

TrendConstantStorage::~TrendConstantStorage()
{
}

Trend *
TrendConstantStorage::GenerateTrend(const std::string                       & /*path*/,
                                    const std::vector<std::string>          & /*trend_cube_parameters*/,
                                    const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                    std::string                             & /*errTxt*/) const
{
  Trend * trend = new TrendConstant(mean_value_);
  return trend;
}

//-------------------------------------------------------------------//

Trend1DStorage::Trend1DStorage(const std::string file_name,
                               const std::string reference_parameter)
: file_name_(file_name),
  reference_parameter_(reference_parameter)
{
}

Trend1DStorage::Trend1DStorage()
{
}

Trend1DStorage::Trend1DStorage(const Trend1DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_(trend_storage.reference_parameter_)
{
}

Trend1DStorage::~Trend1DStorage()
{
}

Trend *
Trend1DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              std::string                             & errTxt) const
{
  int reference = 0;
  for(int i=0; i<static_cast<int>(trend_cube_parameters.size()); i++) {
    if(reference_parameter_ == trend_cube_parameters[i])
      reference = i+1;
  }

  if(reference == 0) {
    errTxt += "The reference parameter of the 1D trend in "+file_name_+" \n"
              "is not the same as the parameter names of the trend cubes\n";
    return(0);
  }

  std::string file_name = path + file_name_;

  int file_format = GetTrend1DFileFormat(file_name, errTxt);

  std::vector<double> trend_values;
  double              s_min;
  double              dz;

  if(file_format < 0) {
    errTxt += "Invalid 1D trend file\n";
    return(0);
  }
  else {
    ReadTrend1DJason(file_name,errTxt,trend_values,s_min,dz);

    double s_max          = s_min + dz*static_cast<int>(trend_values.size());
    int    n_cube_samples = static_cast<int>(trend_cube_sampling[reference-1].size());

    std::vector<double> resampled_trend(static_cast<int>(n_cube_samples));;

    if(trend_cube_sampling[reference-1][0] < s_min) {
      errTxt += "The mimimum value of the 1D trend in "+file_name_+" \n"
                "is lower than the minimum value of "+trend_cube_parameters[reference-1]+"\n";
    }
    else if(trend_cube_sampling[reference-1][n_cube_samples-1] > s_max) {
      errTxt += "The maximum value of the 1D trend in "+file_name_+" \n"
                "is higher than the maximum value of "+trend_cube_parameters[reference-1]+"\n";
    }
    else {

      std::vector<double> trend_sampling(static_cast<int>(trend_values.size()));

      for(int i=0; i<static_cast<int>(trend_values.size()); i++)
        trend_sampling[i] = s_min + i*dz;

      ResampleTrend1D(trend_sampling,
                      trend_values,
                      trend_cube_sampling[reference-1],
                      resampled_trend);

    }

    double increment = trend_cube_sampling[reference-1][1]-trend_cube_sampling[reference-1][0];

    Trend * trend = new Trend1D(resampled_trend, reference, increment);

    return trend;
  }
}


//-------------------------------------------------------------------//

Trend2DStorage::Trend2DStorage(const std::string file_name,
                               const std::string reference_parameter1,
                               const std::string reference_parameter2)
: file_name_(file_name),
  reference_parameter_one_(reference_parameter1),
  reference_parameter_two_(reference_parameter2)
{
}

Trend2DStorage::Trend2DStorage(const Trend2DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_one_(trend_storage.reference_parameter_one_),
  reference_parameter_two_(trend_storage.reference_parameter_two_)
{
}

Trend2DStorage::~Trend2DStorage()
{
}

Trend *
Trend2DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              std::string                             & errTxt) const
{
  int reference1 = 0;
  int reference2 = 0;
  for(int i=0; i<static_cast<int>(trend_cube_parameters.size()); i++) {
    if(reference_parameter_one_ == trend_cube_parameters[i])
      reference1 = i+1;
    else if(reference_parameter_two_ == trend_cube_parameters[i])
      reference2 = i+1;
  }

  if(reference1 == 0 || reference2 == 0)
    errTxt += "The reference parameters of the 2D trend in "+file_name_+" \n"
              "is not the same as the parameter names of the trend cubes\n";

  std::string file_name = path + file_name_;
  RegularSurface<double> surface(file_name);

  RegularSurface<double> resampled_surface;
  bool transpose;
  if(reference1 == 1)
    transpose = false;
  else {
    transpose = true;
    reference1 = 1;
    reference2 = 2;
  }
  resampled_surface = ResampleTrend2D(surface, trend_cube_sampling[0], trend_cube_sampling[1], transpose);

  double increment1 = trend_cube_sampling[0][1] - trend_cube_sampling[0][0];
  double increment2 = trend_cube_sampling[1][1] - trend_cube_sampling[1][0];

  Trend * trend = new Trend2D(resampled_surface, reference1, reference2, increment1, increment2);

  return trend;

}

