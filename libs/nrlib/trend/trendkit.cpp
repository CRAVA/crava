// $Id: trendkit.cpp 1075 2012-09-19 13:42:16Z georgsen $
#include "trendkit.hpp"
#include "../iotools/logkit.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/regularsurface.hpp"
#include <fstream>

namespace NRLib {

void EstimateConstantTrend(const std::vector<std::vector<float> >  & blocked_logs,
                           double                                  & trend)
{
  int nWells  = static_cast<int>(blocked_logs.size());

  if(nWells > 0) {
    int nBlocks = static_cast<int>(blocked_logs[0].size());

    trend = 0;

    int count = 0;

    for (int w = 0 ; w < nWells ; w++) {
      if(blocked_logs[w].size() > 0) {
        for (int k = 0 ; k < nBlocks ; k++) {
          if(blocked_logs[w][k] != RMISSING) {
            trend += exp(blocked_logs[w][k]);
            count++;
          }
        }
      }
    }


    if (count > 0)
      trend = trend/count;
  }
}

//-------------------------------------------------------------------------------

void ReadTrend1DJason(const std::string   & file_name,
                      std::string         & errText,
                      std::vector<double> & trend1d,
                      double              & s_min,
                      double              & dz)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  bool lineIsComment = true;
  int  line          = 0;
  int  thisLine      = 0;

  while( lineIsComment == true) {

    if(CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }

    ReadNextToken(file,dummyStr,line);
    if (line == thisLine)
      DiscardRestOfLine(file,line,false);
    thisLine = line;
    if((dummyStr[0]!='*') &  (dummyStr[0]!='"')) {
      lineIsComment = false;
    }
  }

  s_min = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file))  {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);
  thisLine = line;

  dz = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file)) {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);
  thisLine = line;

  int nz = ParseType<int>(dummyStr);

  trend1d.resize(nz);

  for(int i=0; i<nz; i++) {
    if (CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }
    ReadNextToken(file,dummyStr,line);

    trend1d[i] = ParseType<double>(dummyStr);
  }
  file.close();
}

//----------------------------------------------------//

void ReadTrend1DPlainAscii(const std::string   & file_name,
                           std::string         & /*errText*/,
                           std::vector<double> & trend1d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  int line = 0;

  while(CheckEndOfFile(file)==false){
    ReadNextToken(file,dummyStr,line);
    trend1d.push_back(ParseType<double>(dummyStr));
  }
  file.close();

}

//----------------------------------------------------//

int
GetTrend1DFileFormat(const std::string & file_name,
                     std::string       & errText)
{
  std::string   dummyStr;
  std::string   targetString;
  std::ifstream file;

  // test for jason file format
  OpenRead(file,file_name);

  int  fileformat    = -1;
  int  line          = 0;
  int  thisLine      = 0;
  bool lineIsComment = true;
  bool commentFound = false;

  while (lineIsComment == true) {
    ReadNextToken(file,dummyStr,line);
    if (CheckEndOfFile(file)) {
      errText += "End of trend file "+file_name+" is premature\n";
      return 0;
    }
    else {
      if (thisLine == line) {
        DiscardRestOfLine(file,line,false);
        thisLine = line;
      }
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
        lineIsComment = false;
      else
        commentFound = true;
    }
  }
  file.close();
  if (IsNumber(dummyStr) && commentFound == true) // not convertable number
    fileformat = 0; //Same file format as WAVELET::JASON
  else
    fileformat = 1; // plain ascii file
  return fileformat;
}

//----------------------------------------------------//

void InterpolateTrend1DValues(const double & xi,
                              const double & xi1,
                              const double & fxi,
                              const double & fxi1,
                              const double & yj,
                              double       & fyj)
{
  double t;

  t   = (yj-xi1)/(xi-xi1);

  fyj = fxi*t + fxi1*(1-t);

}

//----------------------------------------------------//

void ResampleTrend1D(const std::vector<double> & x,
                     const std::vector<double> & fx,
                     const std::vector<double> & y,
                     std::vector<double>       & fy)
{
  // Resample fx with sampling x into fy with sampling y

  int nx = static_cast<int>(x.size());
  int ny = static_cast<int>(y.size());

  int i=0;
  int j=0;

  while(i<nx-1) {

    while(j<ny && y[j]>=x[i] && y[j]<x[i+1]) {

      InterpolateTrend1DValues(x[i],
                               x[i+1],
                               fx[i],
                               fx[i+1],
                               y[j],
                               fy[j]);
      j++;
    }

    i++;
  }
}
//-------------------------------------------------------------------------------
void
Estimate1DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                std::vector<double>                     & trend)
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
SmoothTrendWithLocalLinearRegression(std::vector<double>    & trend,
                                     const std::vector<int> & count,
                                     const int              & iWells)
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

//----------------------------------------------------//

RegularSurface<double>
ResampleTrend2D(const RegularSurface<double> & surface,
                const std::vector<double>    & x,
                const std::vector<double>    & y,
                const bool                   & transpose)
{
  int length_x = static_cast<int>(x.size());
  int length_y = static_cast<int>(y.size());

  Grid2D<double> resampled_grid(length_x, length_y);

  for(int i=0; i<length_x; i++) {
    for(int j=0; j<length_y; j++) {
      if(transpose)
        resampled_grid(i,j) = surface.GetZ(y[j],x[i]);
      else
        resampled_grid(i,j) = surface.GetZ(x[i],y[j]);
    }
  }

  double x0 = x[0];
  double y0 = y[0];

  RegularSurface<double> resampled_surface(x0, y0, length_x, length_y, resampled_grid);

  return(resampled_surface);
}

//----------------------------------------------------//

void ReadTrend2DPlainAscii(const std::string     & file_name,
                           std::string           & err_txt,
                           NRLib::Grid2D<double> & trend2d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  int line = 0;

  int ni = ReadNext<int>(file, line);
  int nj = ReadNext<int>(file, line);
  trend2d.Resize(ni, nj);

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      if (!CheckEndOfFile(file))
        trend2d(i,j) = ReadNext<double>(file, line);
      else
        err_txt += "Premature end of file for 2D trend: " + file_name + ".\n";
    }
  }

  file.close();
}


}

