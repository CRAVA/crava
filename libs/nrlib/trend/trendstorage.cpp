// $Id: trendstorage.cpp 1166 2013-05-03 11:34:58Z ulvmoen $
#include "trendstorage.hpp"
#include "trendkit.hpp"
#include "trend.hpp"

#include "../surface/regularsurface.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/logkit.hpp"

#include <fstream>

using namespace NRLib;

TrendStorage::TrendStorage()
{
}

TrendStorage::~TrendStorage()
{
}

//-------------------------------------------------------------------//

TrendConstantStorage::TrendConstantStorage(const double & value,
                                           const bool   & estimate)
: mean_value_(value),
  estimate_(estimate)
{
}

TrendConstantStorage::TrendConstantStorage()
{
}

TrendConstantStorage::TrendConstantStorage(const TrendConstantStorage & trend_storage)
: mean_value_(trend_storage.mean_value_),
  estimate_(trend_storage.estimate_)
{
}

TrendConstantStorage::~TrendConstantStorage()
{
}

Trend *
TrendConstantStorage::GenerateTrend(const std::string                       & /*path*/,
                                    const std::vector<std::string>          & /*trend_cube_parameters*/,
                                    const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                    const std::vector<std::vector<float> >  & /*blocked_logs*/,
                                    std::string                             & errTxt) const
{
  if(estimate_ == true)
    errTxt += "Estimation of value in rock physics models has not been implemented yet\n";

  Trend * trend = new TrendConstant(mean_value_);
  return trend;
}

//-------------------------------------------------------------------//

Trend1DStorage::Trend1DStorage(const std::string & file_name,
                               const std::string & reference_parameter,
                               const bool        & estimate)
: file_name_(file_name),
  reference_parameter_(reference_parameter),
  estimate_(estimate)
{
}

Trend1DStorage::Trend1DStorage()
{
}

Trend1DStorage::Trend1DStorage(const Trend1DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_(trend_storage.reference_parameter_),
  estimate_(trend_storage.estimate_)
{
}

Trend1DStorage::~Trend1DStorage()
{
}

Trend *
Trend1DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              const std::vector<std::vector<float> >  & blocked_logs,
                              std::string                             & errTxt) const
{
  Trend * trend = NULL;

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

  int n_cube_samples = static_cast<int>(trend_cube_sampling[reference-1].size());
  double increment   = trend_cube_sampling[reference-1][1]-trend_cube_sampling[reference-1][0];

  std::vector<double> resampled_trend(static_cast<int>(n_cube_samples));;
  std::vector<double> trend_values;
  std::vector<double> trend_sampling;
  double              s_min;
  double              dz;

  if(estimate_ == true) {

    Estimate1DTrend(blocked_logs, trend_values);

    s_min = trend_cube_sampling[reference-1][0];
    double s_max = trend_cube_sampling[reference-1][n_cube_samples-1];
    dz = (s_max - s_min + 1)/(static_cast<int>(trend_values.size()));

    trend_sampling.resize(trend_values.size());

    for(int i=0; i<static_cast<int>(trend_values.size()); i++)
      trend_sampling[i] = s_min + i*dz;

  }
  else {

    std::string file_name = path + file_name_;

    int file_format = GetTrend1DFileFormat(file_name, errTxt);

    if(file_format < 0) {
      errTxt += "Invalid 1D trend file\n";
      return(0);
    }
    else {
      ReadTrend1DJason(file_name,errTxt,trend_values,s_min,dz);

      double s_max          = s_min + dz*static_cast<int>(trend_values.size());

      if(trend_cube_sampling[reference-1][0] < s_min) {
        errTxt += "The mimimum value of the 1D trend in "+file_name_+" \n"
                  " is lower than the minimum value of "+trend_cube_parameters[reference-1]+"\n";
      }
      else if(trend_cube_sampling[reference-1][n_cube_samples-1] > s_max) {
        errTxt += "The maximum value of the 1D trend in "+file_name_+" \n"
                  " is higher than the maximum value of "+trend_cube_parameters[reference-1]+"\n";
      }
      else {

        trend_sampling.resize(trend_values.size());

        for(int i=0; i<static_cast<int>(trend_values.size()); i++)
          trend_sampling[i] = s_min + i*dz;

      }
    }
  }

  if(errTxt == "") {
    ResampleTrend1D(trend_sampling,
                    trend_values,
                    trend_cube_sampling[reference-1],
                    resampled_trend);

    trend = new Trend1D(resampled_trend, reference, increment);

  }

  return trend;
}

//-------------------------------------------------------------------------------
void
Trend1DStorage::Estimate1DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                                std::vector<double>                     & trend) const
{
  int nWells  = static_cast<int>(blocked_logs.size());

  if(nWells > 0) {
    int nBlocks = static_cast<int>(blocked_logs[0].size());

    trend.resize(nBlocks, 0);

    std::vector<int> count(nBlocks, 0);

    int iWells = 0;
    for (int w = 0 ; w < nWells ; w++) {
      if(blocked_logs[w].size() > 0) {
        for (int k = 0 ; k < nBlocks ; k++) {
          if(blocked_logs[w][k] != RMISSING) {
            trend[k] += exp(blocked_logs[w][k]);
            count[k]++;
          }
        }
        iWells++;
      }
    }

    for (int k = 0 ; k < nBlocks ; k++) {
      if (count[k] > 0) {
        trend[k] = trend[k]/count[k];
      }
    }

    SmoothTrendWithLocalLinearRegression(trend, count, iWells);
  }
}
//-------------------------------------------------------------------------------
void
Trend1DStorage::SmoothTrendWithLocalLinearRegression(std::vector<double>    & trend,
                                                     const std::vector<int> & count,
                                                     const int              & iWells) const
{
  // Copied from Background

  //
  // 1. Center-parts of scatter plots
  //
  // In the center parts of the scatter plots the average value should be
  // accepted as a trend value if the number of data points behind each
  // average is fraction * iWells, where fraction is the acceptance
  // fraction, typically larger than one.
  //
  // Sometimes we have only one, two, or three wells available. In the
  // case of two and three wells, the logs for these may differ considerably,
  // in which case the trend need to be stabilised by requiring a minimum
  // number of data points behind each trend value. The denser the log is
  // sampled, the more prone it is to numerical instabilities. This
  // minimum value is therefore linked to the sampling density 'dz'.
  //
  // Finally, we must require a definite minimum number of data points
  // that should be behind each trend value.
  //
  // nDataMin = max(min_req_points, max(min_time_sample, fraction * iWells))
  //
  // 2. End-parts of scatter plot
  //
  // If the blocked wells do not contain any values for the lower layers of
  // the simbox, the end part of the scatter plot needs special attention.
  //
  // We must possibly require a larger min_req_points to avoid an
  // "arbitrary" end behaviour.
  //

  float fraction   = 5.0f;                      // Require minimum 5*iWells
  int   nTimeLimit = static_cast<int>(50.0/4); // The smaller sampling density, the more values are needed.
  int   nLowLimit  = 10;                        // Require minimum 10
  int   nDataMin   = std::max(nLowLimit, std::max(nTimeLimit, int(fraction * iWells)));
  int   nz         = static_cast<int>(trend.size());
  float min_value  = 1300.0f;
  float max_value  = 7000.0f;

  bool  use_weights = true;
  bool  errorMid    = false;
  bool  errorHead   = false;
  bool  errorTrail  = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  std::vector<double> mean = trend;

  //
  // Find first non-missing value
  //
  int firstNonmissing = 0;
  for (int k = 0 ; k < nz ; k++) {
    if (trend[k] > 0.0f) {
      firstNonmissing = k;
      break;
    }
  }

  //
  // Find last non-missing value
  //
  int lastNonmissing = nz - 1;
  for (int k = nz - 1 ; k > 0 ; k--) {
    if (trend[k] > 0.0f) {
      lastNonmissing = k;
      break;
    }
  }

  std::vector<double> x(nz);  // Time indices
  std::vector<double> y(nz);  // Log values
  std::vector<double> w(nz);  // Weights (number of data behind each avg.)

  for (int k = 0 ; k < nz ; k++) {
    int nCurDataMin = nDataMin;
    if (k < firstNonmissing || k > lastNonmissing) {
      nCurDataMin *= 2;
    }

    int n = 0;
    int nData = 0;
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<double>(count[k]);
      x[0]   = static_cast<double>(k);
      y[0]   = trend[k];
      nData += count[k];
      n++;
    }

    //
    // 2. Add local data points to get 'nCurDataMin' points behind each trend
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (nData < nCurDataMin) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]   = static_cast<double>(count[k - i]);
        x[n]   = static_cast<double>(k - i);
        y[n]   = mean [k - i];
        nData += count[k - i];
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]   = static_cast<double>(count[k + i]);
        x[n]   = static_cast<double>(k + i);
        y[n]   = mean [k + i];
        nData += count[k + i];
        n++;
      }
      if (k-i < 0 && k+i >= nz) { // We will never find enough data
        break;
      }
    }

    //
    // Calculate normalised weights
    //
    if (use_weights)
      for (i = 0 ; i < n ; i++)
        w[i] /= nData;
    else
      for (i = 0 ; i < n ; i++)
        w[i] = 1.0f/n;

    //
    // We need at least two points to make a line.
    //
    if (n > 1) {
      //
      // Estimate local regression line: y = bx + a
      //
      double Sx  = x[0]*w[0];
      double Sy  = y[0]*w[0];
      double Sxx = x[0]*w[0]*x[0];
      double Sxy = x[0]*w[0]*y[0];
      for (i = 1 ; i < n ; i++) {
        Sx  += x[i]*w[i];
        Sy  += y[i]*w[i];
        Sxx += x[i]*w[i]*x[i];
        Sxy += x[i]*w[i]*y[i];
      }
      double b = (Sxy - Sx*Sy)/(Sxx - Sx*Sx);
      double a = (Sy - b*Sx);
      //
      // Estimate value of regression line at requested point.
      //
      double value = a + b*static_cast<double>(k);
      if (value < min_value || value > max_value) {
        if (k < firstNonmissing)
          errorHead = true;
        else if (k > lastNonmissing)
          errorTrail = true;
        else {
          errorMid   = true;
          break;
        }
      }
      trend[k] = value;
    }
    else {
      trend[k] = y[0];
    }
  }

  if (errorMid) {
    // Big problem ...
    LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend using local linear\n");
    LogKit::LogFormatted(LogKit::Low,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::Low,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::Low,"          2) There are too many layers in grid compared to well logs available.\n");
    double sum = 0.0f;
    int nData = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum   += mean[k]*count[k];
          nData += count[k];
        }
        else {
          sum += mean[k];
          nData += 1;
        }
      }
    }
    double global_mean = sum/nData;
    for (int k = 0 ; k < nz ; k++) {
      trend[k] = global_mean;
    }
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean\n");
  }
  else {
    if (errorHead) {
      // Fix first part of trend containing missing-values.
      double firstValue = trend[firstNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [0-%d] where the log is undefined. The first\n",firstNonmissing-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(firstValue),firstNonmissing);
      for (int k = 0 ; k < firstNonmissing ; k++) {
        trend[k] = firstValue;
      }
    }
    if (errorTrail) {
      // Fix last part of trend containing missing-values.
      double lastValue = trend[lastNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",lastNonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(lastValue),lastNonmissing);
      for (int k = lastNonmissing + 1 ; k < nz ; k++) {
        trend[k] = lastValue;
      }
    }
  }
}

//-------------------------------------------------------------------//

Trend2DStorage::Trend2DStorage(const std::string & file_name,
                               const std::string & reference_parameter1,
                               const std::string & reference_parameter2,
                               const bool        & estimate)
: file_name_(file_name),
  reference_parameter_one_(reference_parameter1),
  reference_parameter_two_(reference_parameter2),
  estimate_(estimate)
{
}

Trend2DStorage::Trend2DStorage(const Trend2DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_one_(trend_storage.reference_parameter_one_),
  reference_parameter_two_(trend_storage.reference_parameter_two_),
  estimate_(trend_storage.estimate_)
{
}

Trend2DStorage::~Trend2DStorage()
{
}

Trend *
Trend2DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              const std::vector<std::vector<float> >  & /*blocked_logs*/,
                              std::string                             & errTxt) const
{
  Trend * trend = NULL;

  if(estimate_ == true)
    errTxt += "Estimation of 1D trend in rock physics models has not been implemented yet\n";

  else {

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
                " is not the same as the parameter names of the trend cubes\n";

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

    trend = new Trend2D(resampled_surface, reference1, reference2, increment1, increment2);
  }

  return trend;

}

