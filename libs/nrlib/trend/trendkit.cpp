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

void ReadTrend1D(const std::string        & file_name,
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

void WriteTrend1D(const std::string         & filename,
                  const std::vector<double> & s,
                  const std::vector<double> & trend)
{
  std::ofstream file;
  OpenWrite(file, filename);

  file << std::fixed
       << std::setprecision(6)
       << s.front()   << "\n"
       << s[1] - s[0] << "\n"
       << s.back()    << "\n";

  for (size_t i = 0; i < trend.size(); i++) {
    file << trend[i] << "\n";
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
void Estimate1DTrend(const std::vector<std::vector<float> >  & blocked_logs,
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
void Estimate2DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                     const std::vector<std::vector<double> > & trend_cube_sampling,
                     const std::vector<std::vector<double> > & s1,
                     const std::vector<std::vector<double> > & s2,
                     std::vector<std::vector<double> >       & trend,
                     std::string                             & errTxt)
{
  double              scale = 1.0;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> w;
  std::vector<double> x0 = trend_cube_sampling[0];
  std::vector<double> y0 = trend_cube_sampling[1];
  std::vector<double> x0_regridded;
  std::vector<double> y0_regridded;

  double              bandwidth_x;
  double              bandwidth_y;

  bool valid_dataset = PreprocessData2D(blocked_logs,
                                        trend_cube_sampling,
                                        s1,
                                        s2,
                                        scale,
                                        x,
                                        y,
                                        z,
                                        w,
                                        x0_regridded,
                                        y0_regridded,
                                        bandwidth_x,
                                        bandwidth_y,
                                        errTxt);

  if (valid_dataset) {

    size_t                            x0_regridded_n = x0_regridded.size();
    size_t                            y0_regridded_n = y0_regridded.size();

    std::vector<std::vector<double> > z0_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));
    std::vector<std::vector<double> > w0_regridded_dummy(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));

    bool complete_surface = false;
    bool stop             = false;

    size_t l = 0;
    while (!stop) {
      LocalLinearRegression2DSurface(x,
                                     y,
                                     z,
                                     w,
                                     bandwidth_x,
                                     bandwidth_y,
                                     x0_regridded,
                                     y0_regridded,
                                     z0_regridded,
                                     w0_regridded_dummy,
                                     complete_surface);

      bandwidth_x = (1 + 0.02)*bandwidth_x;
      bandwidth_y = (1 + 0.02)*bandwidth_y;

      if (complete_surface) {
        stop = true;
      } else if (l > 1e4) {
        errTxt +=  "Error: Unable to compute complete trend surface (type 3). \n";
        errTxt +=  "       The interpolated region is too large compared to the support of the data. \n";
        errTxt +=  "       Reduce the size of the estiamted region or use a low dimensional method. \n";
      }
      l = l + 1;
    }
    BilinearInterpolation(x0_regridded, y0_regridded, z0_regridded, x0, y0, trend);
  }
}
//-------------------------------------------------------------------------------
bool PreprocessData2D(const std::vector<std::vector<float> >  & blocked_logs,
                      const std::vector<std::vector<double> > & trend_cube_sampling,
                      const std::vector<std::vector<double> > & s1,
                      const std::vector<std::vector<double> > & s2,
                      const double                            & scale,
                      std::vector<double>                     & x,
                      std::vector<double>                     & y,
                      std::vector<double>                     & z,
                      std::vector<double>                     & w,
                      std::vector<double>                     & x0_regridded,
                      std::vector<double>                     & y0_regridded,
                      double                                  & bandwidth_x,
                      double                                  & bandwidth_y,
                      std::string                             & errTxt)
{
  size_t n_wells   = blocked_logs.size();
  size_t n_samples = 0;

  if (n_wells > 0) {
    for (size_t i = 0; i < n_wells; i++) {
      size_t n_samples_i = blocked_logs[i].size();

      if (n_samples_i > 0) {
        n_samples = n_samples + n_samples_i;

        for (size_t j = 0; j < n_samples_i; j++) {
          double x_j     = s1[i][j];
          double y_j     = s2[i][j];
          double z_j_log = blocked_logs[i][j];

          if (x_j != RMISSING && x_j != RMISSING && z_j_log != RMISSING) {
            x.push_back(x_j);
            y.push_back(y_j);
            z.push_back(std::exp(z_j_log));
            w.push_back(1.0);
          }
        }
      }
    }
  }

  double              delta     = 1e-5;
  double              upp       = 3.0;
  double              large     = 1.0;

  size_t              n_low     = 10;      // see 'weight_total_limit' in KernelSmoother2DSurface() and LocalLinearRegression2DSurface()
  size_t              n_small   = 48;      //const std::string         & filename,
  size_t              n_upp     = 128*128; // see 'n_max' in MakeNewGridResolution()

  double              x_min     = *std::min_element(x.begin(), x.end());
  double              x_max     = *std::max_element(x.begin(), x.end());
  double              x_delta   = std::abs(x_max - x_min);
  double              x_upp     = upp*x_delta;
  double              x_large   = large*x_delta;

  double              y_min     = *std::min_element(y.begin(), y.end());
  double              y_max     = *std::max_element(y.begin(), y.end());
  double              y_delta   = std::abs(y_max - y_min);
  double              y_upp     = upp*y_delta;
  double              y_large   = large*y_delta;

  std::vector<double> x0        = trend_cube_sampling[0];
  std::vector<double> y0        = trend_cube_sampling[1];

  if (!CheckIfVectorIsSorted(x0) || !CheckIfVectorIsSorted(y0)) {
    errTxt += "Error: At least one of the inputs are not sorted. \n";
    return(false);
  }

  if (x0.front() > x_min - x_upp && x0.back() < x_max + x_upp && y0.front() > y_min - y_upp && y0.back() < y_max + y_upp) {
    if (x0.front() < x_min - x_large || x0.back() > x_max + x_large || y0.front() < y_min - y_large || y0.back() > y_max + y_large) {
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The defined region is large compared to the support of the data. This can result in unstable estimates. \n");
      LogKit::LogFormatted(LogKit::Low,"            Consider using a low-dimensional method. \n");
    }

    bandwidth_x = CalculateBandwidth(x, scale*std::pow(2, -0.5), 0.2);
    bandwidth_y = CalculateBandwidth(y, scale*std::pow(2, -0.5), 0.2);

    if (bandwidth_x > delta && bandwidth_y > delta) {
      if (n_samples > n_low) {
        if (n_samples < n_small) {
          LogKit::LogFormatted(LogKit::Low,"\nWARNING : The sample size is relatively small. This can result in unstable estimates. \n");
          LogKit::LogFormatted(LogKit::Low,"            Consider using a low-dimensional method. \n");
        }

        MakeNewGridResolution(bandwidth_x, x0, x0_regridded);
        MakeNewGridResolution(bandwidth_y, y0, y0_regridded);

        if (n_samples > n_upp) {

          MakeBinnedDataset(x0_regridded, y0_regridded, x, y, z, w);

          if (z.size() < n_low) {
            errTxt += "Error: Unable to compute complete trend surface (type 2). \n";
            errTxt += "       The interpolated region is too large compared to the support of the data. \n";
            errTxt += "       Reduce the size of the estiamted region or use a low dimensional method. \n";
            return(false);
          }
        }
      } else {
        errTxt += "Error: Too few observations. \n";
        errTxt += "       The algorithm requires at least 11 observations. \n";
        return(false);
      }
    } else {
      errTxt += "Error: The spread in the data is too low to provide stable estimates. \n ";
      errTxt += "       Use a low dimensional method. \n";
      return(false);
    }
  } else {
     errTxt += "Error: Unable to compute complete trend surface (type 1). \n";
     errTxt += "       The interpolated region is too large compared to the support of the data. \n";
     errTxt += "       Reduce the size of the estiamted region or use a low dimensional method. \n";
    return(false);
  }
  return(true);
}
//-------------------------------------------------------------------------------
bool CheckIfVectorIsSorted(const std::vector<double> & x) {

  bool is_sorted = true;

  for (size_t i = 0; i < x.size() - 1; i++) {
    if (x[i] > x[i + 1]) {
      is_sorted = false;
    }
  }
  return(is_sorted);
}
//-------------------------------------------------------------------------------
double CalculateBandwidth(const std::vector<double> & x,
                          const double              & scale,
                          const double              & power)
{
  size_t n               = x.size();

  double one_over_root_2 = 0.7071068;
  double sd_x            = std::pow(CalculateVariance(x), 0.5);

  return(scale*one_over_root_2*sd_x*std::pow(n, -power));
}
//-------------------------------------------------------------------------------
double CalculateEffectiveSampleSize(const std::vector<double> & x,
                                    const double              & bandwidth)
{
  size_t n     = x.size();
  double x_min = *std::min_element(x.begin(), x.end());
  double x_max = *std::max_element(x.begin(), x.end());

  return(2.506628*bandwidth*n/(x_max - x_min));
}
//-------------------------------------------------------------------------------
double CalculateVariance(const std::vector<double> & x)
{
  int n         = 0;

  double sum_x  = 0.0;
  double sum_x2 = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    if (x[i] != RMISSING) {
      sum_x  = sum_x  + x[i];
      sum_x2 = sum_x2 + x[i]*x[i];
      n      = n + 1;
    } else {
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : TRemoved missing values in variance estimation.\n");
    }
  }
  double var_x  = sum_x2/(n - 1) - sum_x*sum_x/(n*(n - 1));

  return(var_x);
}
//-------------------------------------------------------------------------------
void MakeNewGridResolution(const double              & bandwidth_x,
                           const std::vector<double> & x0,
                           std::vector<double>       & x0_regridded)
{
  double              x0_min           = *std::min_element(x0.begin(), x0.end());
  double              x0_max           = *std::max_element(x0.begin(), x0.end());

  size_t              n_max            = 128;
  size_t              x0_n_regridded   = std::min(static_cast<size_t>(std::ceil(7*(x0_max - x0_min)/(2*bandwidth_x))), n_max);

  double              x0_regridded_inc = (x0_max - x0_min)/(x0_n_regridded - 1);
  x0_regridded.resize(x0_n_regridded);

  for (size_t i = 0; i < x0_n_regridded; i++) {
    x0_regridded[i] = x0_min + x0_regridded_inc*i;
  }
}
//-------------------------------------------------------------------------------
void MakeBinnedDataset(const std::vector<double>         & x0,
                       const std::vector<double>         & y0,
                       std::vector<double>               & x,
                       std::vector<double>               & y,
                       std::vector<double>               & z,
                       std::vector<double>               & w)
{
  double delta = 1e-5;

  size_t n     = x0.size();
  size_t m     = y0.size();

  std::vector<std::vector<double> > z0(n, std::vector<double>(m, 0.0));
  std::vector<std::vector<double> > w0(n, std::vector<double>(m, 0.0));

  BilinearBinning(x, y, z, x0, y0, z0, w0);

  std::vector<double> x_new(n*m, RMISSING);
  std::vector<double> y_new(n*m, RMISSING);
  std::vector<double> z_new(n*m, RMISSING);
  std::vector<double> w_new(n*m, RMISSING);

  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      if (w0[i][j] > 0.0 + delta) {
        x_new[k] = x0[i];
        y_new[k] = y0[j];
        z_new[k] = z0[i][j]/w0[i][j];
        w_new[k] = w0[i][j];

        k = k + 1;
      }
    }
  }
  if (k > 0) {
    x_new.resize(k);
    y_new.resize(k);
    z_new.resize(k);
    w_new.resize(k);
  } else {
    x_new.resize(1);
    y_new.resize(1);
    z_new.resize(1);
    w_new.resize(1);
  }
  x = x_new;
  y = y_new;
  z = z_new;
  w = w_new;
}
//-------------------------------------------------------------------------------
void BilinearBinning(const std::vector<double>               & x,
                     const std::vector<double>               & y,
                     const std::vector<double>               & z,
                     const std::vector<double>               & x0,
                     const std::vector<double>               & y0,
                     std::vector<std::vector<double> >       & z0,
                     std::vector<std::vector<double> >       & w0)
{
  double delta    = 1e-5;

  double x0_first = x0.front();
  double x0_last  = x0.back();
  double x0_inc   = x0[1] - x0[0];
  double x0_scale = 1/x0_inc;

  double y0_first = y0.front();
  double y0_last  = y0.back();
  double y0_inc   = y0[1] - y0[0];
  double y0_scale = 1/y0_inc;

  for (size_t i = 0; i < z.size(); i++) {

    if ((x0_first - delta < x[i]) && (x[i] < x0_last + delta) && (y0_first - delta < y[i]) && (y[i] < y0_last + delta)) {

      size_t k;
      size_t l;
      double w_y;
      double w_x;

      if (x[i] < x0_first + delta) {
        k   = 0;
        w_x = 0.0;
      } else if (x[i] > x0_last - delta) {
        k   = x0.size() - 2;
        w_x = 1.0;
      } else {
        k   = static_cast<size_t>(std::floor((x[i] - x0_first)*x0_scale));
        w_x = (x[i] - x0_first - k*x0_inc)*x0_scale;
      }

      if (y[i] < y0_first + delta) {
        l   = 0;
        w_y = 0.0;
      } else if (y[i] > y0_last - delta) {
        l   = y0.size() - 2;
        w_y = 1.0;
      } else {
        l   = static_cast<size_t>(std::floor((y[i] - y0_first)*y0_scale));
        w_y = (y[i] - y0_first - l*y0_inc)*y0_scale;
      }

      double w00 = (1 - w_x)*(1 - w_y);
      double w10 = w_x*(1 - w_y);
      double w01 = (1 - w_x)*w_y;
      double w11 = w_x*w_y;

      z0[k    ][l    ] = z0[k    ][l    ] + z[i]*w00;
      z0[k + 1][l    ] = z0[k + 1][l    ] + z[i]*w10;
      z0[k    ][l + 1] = z0[k    ][l + 1] + z[i]*w01;
      z0[k + 1][l + 1] = z0[k + 1][l + 1] + z[i]*w11;

      w0[k    ][l    ] = w0[k    ][l    ] + w00;
      w0[k + 1][l    ] = w0[k + 1][l    ] + w10;
      w0[k    ][l + 1] = w0[k    ][l + 1] + w01;
      w0[k + 1][l + 1] = w0[k + 1][l + 1] + w11;
    }
  }
}
//-------------------------------------------------------------------------------
void LocalLinearRegression2DSurface(const std::vector<double>         & x,
                                    const std::vector<double>         & y,
                                    const std::vector<double>         & z,
                                    const std::vector<double>         & w,
                                    const double                      & bandwidth_x,
                                    const double                      & bandwidth_y,
                                    const std::vector<double>         & x0,
                                    const std::vector<double>         & y0,
                                    std::vector<std::vector<double> > & z0,
                                    std::vector<std::vector<double> > & w0,
                                    bool                              & complete_surface)
{
  double               epsilon               =  1e-4;
  double               threshold             =  std::pow(std::log(1e5), 0.5);
  double               one_over_bandwidth_x  =  1.0/bandwidth_x;
  double               one_over_bandwidth_y  =  1.0/bandwidth_y;
  double               weight_total_limit    = 10.0;
  double               delta_limit           =  2*std::log(1e5);

  size_t               nSamples              = z.size();

  std::vector<double> weights(nSamples, RMISSING);
  std::vector<size_t> index(nSamples, nSamples + 1);

  std::vector<double> x_order(nSamples, 0);
  std::vector<size_t> x_backward(nSamples, -1);
  std::vector<size_t> x_forward(nSamples, -1);
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  std::vector<double> y_order(nSamples, 0);
  std::vector<size_t> y_backward(nSamples, -1);
  std::vector<size_t> y_forward(nSamples, -1);
  SortOrderAndRank(y, y_order, y_backward, y_forward);


  complete_surface = true;

  for (size_t i = 0; i < x0.size(); i++) {

    double x_low        = x0[i] - threshold*bandwidth_x;
    double x_upp        = x0[i] + threshold*bandwidth_x;

    size_t x_low_k      = FindLowerBoundInSortedVector(x_order, x_low);
    size_t x_upp_k      = FindUpperBoundInSortedVector(x_order, x_upp);

    for (size_t j = 0; j < y0.size(); j++) {

      if (z0[i][j] == RMISSING) {

        double weight_total = 0.0;

        double y_low        = y0[j] - threshold*bandwidth_y;
        double y_upp        = y0[j] + threshold*bandwidth_y;

        size_t y_low_l      = FindLowerBoundInSortedVector(y_order, y_low);
        size_t y_upp_l      = FindUpperBoundInSortedVector(y_order, y_upp);

        size_t K = 0;
        for (size_t k = x_low_k; k < x_upp_k + 1; k++) {

          size_t l = x_backward[k];

          if (y_low_l <= y_forward[l] && y_forward[l] <= y_upp_l) {
            double x_delta = (x0[i] - x[l])*one_over_bandwidth_x;
            double y_delta = (y0[j] - y[l])*one_over_bandwidth_y;
            double delta   = std::pow(x_delta, 2) + std::pow(y_delta, 2);

            if (delta < delta_limit && z[l] != RMISSING) {
              double weight = w[l]*std::exp(-0.5*delta);

              weight_total  = weight_total + weight;
              weights[K]    = weight;
              index[K]      = l;
              K             = K + 1;
            }
          }
        }
        w0[i][j] = weight_total;

        if (K > 0 && weight_total > weight_total_limit) {
          // exact calculation of the inverse.
          // A_inv = B/det_A, where A = X^{t}WX and A_inv = (X^{t}WX)^{-1}
          double a11 = 0.0;
          double a12 = 0.0;
          double a13 = 0.0;
          double a22 = 0.0;
          double a23 = 0.0;
          double a33 = 0.0;

          // C = X^{t}Wz
          double c1  = 0.0;
          double c2  = 0.0;
          double c3  = 0.0;

          double weight;

          size_t m;

          for (size_t k = 0; k < K; k++) {
            m      = index[k];
            weight = weights[k];

            a11    = a11 +                          weight;
            a12    = a12 + x[m]                    *weight;
            a13    = a13 +           y[m]          *weight;
            a22    = a22 + x[m]*x[m]               *weight;
            a23    = a23 + x[m]     *y[m]          *weight;
            a33    = a33 +           y[m]*y[m]     *weight;

            c1     = c1  +                     z[m]*weight;
            c2     = c2  + x[m]               *z[m]*weight;
            c3     = c3  +           y[m]     *z[m]*weight;
          }

          // add a smll number to the diagonal to make the matrix nonsingular.
          //a11        = a11*(1 + epsilon);
          a22          = a22*(1 + epsilon);
          a33          = a33*(1 + epsilon);

          double b11   =  (a33*a22 - a23*a23);
          double b12   = -(a33*a12 - a23*a13);
          double b13   =  (a23*a12 - a22*a13);
          double b22   =  (a33*a11 - a13*a13);
          double b23   = -(a23*a11 - a12*a13);
          double b33   =  (a22*a11 - a12*a12);

          double det_A = a11*b11 + a12*b12 + a13*b13;

          if (det_A > epsilon) {
            double tmp = 0.0;

            tmp        = tmp +       (b11*c1 + b12*c2 + b13*c3)/det_A;
            tmp        = tmp + x0[i]*(b12*c1 + b22*c2 + b23*c3)/det_A;
            tmp        = tmp + y0[j]*(b13*c1 + b23*c2 + b33*c3)/det_A;

            z0[i][j]   = tmp;

          } else {
            z0[i][j]         = RMISSING;
            complete_surface = false;
          }
        } else {
          z0[i][j]         = RMISSING;
          complete_surface = false;
        }
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------------------
size_t FindLowerBoundInSortedVector(const std::vector<double> & x,
                                    const double              & x_0)
{
  size_t n   = x.size();
  size_t low = 0;
  size_t upp = n - 1;

  if (x_0 <= x[0]) {
    return(0);
  } else if (x_0 >= x[n - 1]) {
    return(n - 1);
  } else {
    size_t j = std::floor((upp - low)*0.5);
    while (upp - low > 1) {
      while (x[j] <= x_0 && upp - low > 1) {
        low = j;
        j   = low + std::floor((upp - low)*0.5);
      }
      while (x_0 < x[j] && upp - low > 1) {
        upp = j;
        j   = low + std::floor((upp - low)*0.5);
      }
    }
    return(j);
  }
}
//--------------------------------------------------------------------------------------------------------------
size_t FindUpperBoundInSortedVector(const std::vector<double> & x,
                                    const double              & x_0)
{
  size_t n   = x.size();
  size_t low = 0;
  size_t upp = n - 1;

  if (x_0 <= x[0]) {
    return(0);
  } else if (x_0 >= x[n - 1]) {
    return(n - 1);
  } else {
    size_t j = std::floor((upp - low)*0.5);
    while (upp - low > 1) {
      while (x[j] < x_0 && upp - low > 1) {
        low = j;
        j   = low + std::floor((upp - low)*0.5);
      }
      while (x_0 <= x[j] && upp - low > 1) {
        upp = j;
        j   = low + std::floor((upp - low)*0.5);
      }
    }
    return(std::min(j + 1, n - 1));
  }
}
//-------------------------------------------------------------------------------
void SortOrderAndRank(const std::vector<double> & x,
                      std::vector<double>       & x_sort,
                      std::vector<size_t>       & x_order,
                      std::vector<size_t>       & x_rank)
{
  size_t n = x.size();

  std::vector<std::pair<double, size_t> > x_order_tmp(n);
  std::vector<std::pair<size_t, size_t> > x_rank_tmp(n);

  for (size_t i = 0; i < n; i++) {
    x_order_tmp[i].first  = x[i];
    x_order_tmp[i].second = i;
  }

  std::sort(x_order_tmp.begin(), x_order_tmp.end());

  for (size_t i = 0; i < n; i++) {
    x_sort[i]            = x_order_tmp[i].first;
    x_order[i]           = x_order_tmp[i].second;

    x_rank_tmp[i].first  = x_order_tmp[i].second;
    x_rank_tmp[i].second = i;
  }

  std::sort(x_rank_tmp.begin(), x_rank_tmp.end());

  for (size_t i = 0; i < n; i++) {
    x_rank[i] = x_rank_tmp[i].second;
  }
}
//-------------------------------------------------------------------------------
void BilinearInterpolation(const std::vector<double>               & x,
                           const std::vector<double>               & y,
                           const std::vector<std::vector<double> > & z,
                           const std::vector<double>               & x0,
                           const std::vector<double>               & y0,
                           std::vector<std::vector<double> >       & z0)
{
  for (size_t i = 0; i < x0.size(); i++) {
    for (size_t j = 0; j < y0.size(); j++) {
      z0[i][j] = Interpolate(x, y, z, x0[i], y0[j]);
    }
  }
}
//-------------------------------------------------------------------------------
double Interpolate(const std::vector<double>               & x,
                   const std::vector<double>               & y,
                   const std::vector<std::vector<double> > & z,
                   const double                            & x0,
                   const double                            & y0)
{
  double delta   = 1e-5;

  double x_first = x.front();
  double x_last  = x.back();
  double x_inc   = x[1] - x[0];
  double x_scale = 1/x_inc;

  double y_first = y.front();
  double y_last  = y.back();
  double y_inc   = y[1] - y[0];
  double y_scale = 1/y_inc;

  if ((x_first - delta < x0) && (x0 < x_last + delta) && (y_first - delta < y0) && (y0 < y_last + delta)) {

    size_t k;
    size_t l;
    double w_y;
    double w_x;

    if (x0 < x_first + delta) {
      k   = 0;
      w_x = 0.0;
    } else if (x0 > x_last - delta) {
      k   = x.size() - 2;
      w_x = 1.0;
    } else {
      k   = static_cast<size_t>(std::floor((x0 - x_first)*x_scale));
      w_x = (x0 - x_first - k*x_inc)*x_scale;
    }

    if (y0 < y_first + delta) {
      l   = 0;
      w_y = 0.0;
    } else if (y0 > y_last - delta) {
      l   = y.size() - 2;
      w_y = 1.0;
    } else {
      l   = static_cast<size_t>(std::floor((y0 - y_first)*y_scale));
      w_y = (y0 - y_first - l*y_inc)*y_scale;
    }

    double v00 = z[k    ][l    ];
    double v10 = z[k + 1][l    ];
    double v01 = z[k    ][l + 1];
    double v11 = z[k + 1][l + 1];

    if (v00 != RMISSING && v10 != RMISSING && v01 != RMISSING && v11 != RMISSING) {
      return(v00*(1 - w_x)*(1 - w_y) + v10*w_x*(1 - w_y) + v01*(1 - w_x)*w_y + v11*w_x*w_y);
    } else {
      return(RMISSING);
    }
  } else {
    return(RMISSING);
  }
}
//-------------------------------------------------------------------------------
void Estimate2DVariance(const std::vector<std::vector<float> >  & blocked_logs,
                        const std::vector<std::vector<double> > & trend_cube_sampling,
                        const std::vector<std::vector<double> > & s1,
                        const std::vector<std::vector<double> > & s2,
                        const std::vector<std::vector<double> > & trend,
                        std::vector<std::vector<double> >       & var,
                        std::string                             & errTxt)
{
  double              scale = 1.0;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> w;
  std::vector<double> x0 = trend_cube_sampling[0];
  std::vector<double> y0 = trend_cube_sampling[1];
  std::vector<double> x0_regridded;
  std::vector<double> y0_regridded;

  double              bandwidth_x;
  double              bandwidth_y;

  bool valid_dataset = PreprocessData2D(blocked_logs,
                                        trend_cube_sampling,
                                        s1,
                                        s2,
                                        scale,
                                        x,
                                        y,
                                        z,
                                        w,
                                        x0_regridded,
                                        y0_regridded,
                                        bandwidth_x,
                                        bandwidth_y,
                                        errTxt);

  if (valid_dataset) {

    double              z_var_weight  = 10.0;
    size_t              nSamples      = z.size();

    /* -- estiamte global variance --*/
    std::vector<double> z_z_mean_squared(nSamples, RMISSING);
    double              z_var = 0.0;
    int                 z_n   = 0;
    for (size_t i = 0; i < nSamples; i++) {
      double z_mean_i = Interpolate(x0, y0, trend, x[i], y[i]);
      if (z_mean_i != RMISSING && z[i] != RMISSING) {
        z_z_mean_squared[i]  = (z[i] - z_mean_i)*(z[i] - z_mean_i);
        z_var                = z_var + z_z_mean_squared[i];
        z_n                  = z_n + 1;
      }
    }
    z_var = z_var/(z_n - 1);
    /* ------------------------------*/

    size_t                            x0_regridded_n = x0_regridded.size();
    size_t                            y0_regridded_n = y0_regridded.size();

    std::vector<std::vector<double> > z0_weighted_var_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));
    std::vector<std::vector<double> > z0_var_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));
    std::vector<std::vector<double> > w0_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));

    bool complete_surface_dummy = false;
    KernelSmoother2DSurface(x,
                            y,
                            z_z_mean_squared,
                            w,
                            bandwidth_x,
                            bandwidth_y,
                            x0_regridded,
                            y0_regridded,
                            z0_var_regridded,
                            w0_regridded,
                            complete_surface_dummy);

    for (size_t i = 0; i < x0_regridded_n; i++) {
      for (size_t j = 0; j < y0_regridded_n; j++) {
        double numerator   = z_var_weight*z_var;
        double denominator = z_var_weight;

        if (z0_var_regridded[i][j] != RMISSING) {
          numerator   = numerator   + w0_regridded[i][j]*z0_var_regridded[i][j];
          denominator = denominator + w0_regridded[i][j];
        }
        z0_weighted_var_regridded[i][j] = numerator/denominator;
      }
    }
    BilinearInterpolation(x0_regridded, y0_regridded, z0_weighted_var_regridded, x0, y0, var);
  }
}
//-------------------------------------------------------------------------------
void KernelSmoother2DSurface(const std::vector<double>         & x,
                             const std::vector<double>         & y,
                             const std::vector<double>         & z,
                             const std::vector<double>         & w,
                             const double                      & bandwidth_x,
                             const double                      & bandwidth_y,
                             const std::vector<double>         & x0,
                             const std::vector<double>         & y0,
                             std::vector<std::vector<double> > & z0,
                             std::vector<std::vector<double> > & w0,
                             bool                              & complete_surface)
{
  double               threshold             =  std::pow(std::log(1e5), 0.5);
  double               one_over_bandwidth_x  =  1.0/bandwidth_x;
  double               one_over_bandwidth_y  =  1.0/bandwidth_y;
  double               weight_total_limit    = 10.0;
  double               delta_limit           =  2*std::log(1e5);

  size_t               nSamples              = z.size();

  std::vector<double> weights(nSamples, RMISSING);
  std::vector<size_t> index(nSamples, nSamples + 1);

  std::vector<double> x_order(nSamples, 0);
  std::vector<size_t> x_backward(nSamples, -1);
  std::vector<size_t> x_forward(nSamples, -1);
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  std::vector<double> y_order(nSamples, 0);
  std::vector<size_t> y_backward(nSamples, -1);
  std::vector<size_t> y_forward(nSamples, -1);
  SortOrderAndRank(y, y_order, y_backward, y_forward);


  complete_surface = true;

  for (size_t i = 0; i < x0.size(); i++) {

    double x_low        = x0[i] - threshold*bandwidth_x;
    double x_upp        = x0[i] + threshold*bandwidth_x;

    size_t x_low_k      = FindLowerBoundInSortedVector(x_order, x_low);
    size_t x_upp_k      = FindUpperBoundInSortedVector(x_order, x_upp);

    for (size_t j = 0; j < y0.size(); j++) {

      if (z0[i][j] == RMISSING) {

        double weight_total = 0.0;
        double z0_hat       = 0.0;

        double y_low        = y0[j] - threshold*bandwidth_y;
        double y_upp        = y0[j] + threshold*bandwidth_y;

        size_t y_low_l      = FindLowerBoundInSortedVector(y_order, y_low);
        size_t y_upp_l      = FindUpperBoundInSortedVector(y_order, y_upp);

        for (size_t k = x_low_k; k < x_upp_k + 1; k++) {

          size_t l = x_backward[k];

          if (y_low_l <= y_forward[l] && y_forward[l] <= y_upp_l) {
            double x_delta = (x0[i] - x[l])*one_over_bandwidth_x;
            double y_delta = (y0[j] - y[l])*one_over_bandwidth_y;
            double delta   = std::pow(x_delta, 2) + std::pow(y_delta, 2);

            if (delta < delta_limit && z[l] != RMISSING) {
              double weight = w[l]*std::exp(-0.5*delta);

              weight_total  = weight_total + weight;
              z0_hat        = z0_hat       + weight*z[l];
            }
          }
        }
        w0[i][j] = weight_total;

        if (weight_total > weight_total_limit) {
          z0[i][j] = z0_hat/weight_total;
        } else {
          z0[i][j]         = RMISSING;
          complete_surface = false;
        }
      }
    }
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
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean: %.2f\n\n",global_mean);
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

