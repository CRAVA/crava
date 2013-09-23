// $Id: trendkit.hpp 1075 2012-09-19 13:42:16Z georgsen $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NRLIB_TRENDKIT_HPP
#define NRLIB_TRENDKIT_HPP
#include <string>
#include <vector>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../surface/regularsurface.hpp"
#include "trend.hpp"

namespace NRLib {

#define RMISSING -99999.000
#define IMISSING -99999

  void                    EstimateConstantTrend(const std::vector<std::vector<float> >  & blocked_logs,
                                                double                                  & trend);

  int                     GetTrend1DFileFormat(const std::string & file_name,
                                               std::string       & errText);

  void                    ReadTrend1D(const std::string   & file_name,
                                      std::string         & errText,
                                      std::vector<double> & trend1d,
                                      double              & s_min,
                                      double              & dz);

  void                    WriteTrend1D(const std::string         & file_name,
                                       const std::vector<double> & s1,
                                       const std::vector<double> & trend);

  void                    ReadTrend1DPlainAscii(const std::string   & file_name,
                                                std::string         & /*errText*/,
                                                std::vector<double> & trend1d);

  void                    InterpolateTrend1DValues(const double & xi,
                                                   const double & xi1,
                                                   const double & fxi,
                                                   const double & fxi1,
                                                   const double & yj,
                                                   double       & fyj);

  void                   ResampleTrend1D(const std::vector<double> & x,
                                         const std::vector<double> & fx,
                                         const std::vector<double> & y,
                                         std::vector<double>       & fy);

  void                   Estimate1DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                                         std::vector<double>                     & trend);

  void                   Estimate2DTrend(const std::vector<std::vector<float> >  & blocked_logs,
                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                         const std::vector<std::vector<double> > & s1,
                                         const std::vector<std::vector<double> > & s2,
                                         std::vector<std::vector<double> >       & trend,
                                         std::string                             & errTxt);

  bool                   PreprocessData2D(const std::vector<std::vector<float> >  & blocked_logs,
                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                          const std::vector<std::vector<double> > & s1,
                                          const std::vector<std::vector<double> > & s2,
                                          const double                              scale,
                                          std::vector<double>                     & x,
                                          std::vector<double>                     & y,
                                          std::vector<double>                     & z,
                                          std::vector<double>                     & w,
                                          std::vector<double>                     & x0,
                                          std::vector<double>                     & y0,
                                          double                                  & bandwidth_x,
                                          double                                  & bandwidth_y,
                                          std::string                             & errTxt);

  bool                   CheckIfVectorIsSorted(const std::vector<double> & x);

  double                 CalculateBandwidth(const std::vector<double> & x,
                                            const double                scale,
                                            const double                power);

  double                 CalculateVariance(const std::vector<double> & x);

  void                   MakeNewGridResolution(const double                bandwidth_x,
                                               const std::vector<double> & x0,
                                               std::vector<double>       & x0_regridded);


  void                   MakeBinnedDataset(const std::vector<double>         & x0,
                                           const std::vector<double>         & y0,
                                           std::vector<double>               & x,
                                           std::vector<double>               & y,
                                           std::vector<double>               & z,
                                           std::vector<double>               & w);

  void                   BilinearBinning(const std::vector<double>               & x,
                                         const std::vector<double>               & y,
                                         const std::vector<double>               & z,
                                         const std::vector<double>               & x0,
                                         const std::vector<double>               & y0,
                                         std::vector<std::vector<double> >       & z0,
                                         std::vector<std::vector<double> >       & w0);

  void                   LocalLinearRegression2DSurface(const std::vector<double>         & x,
                                                        const std::vector<double>         & y,
                                                        const std::vector<double>         & z,
                                                        const std::vector<double>         & w,
                                                        const double                        bandwidth_x,
                                                        const double                        bandwidth_y,
                                                        const std::vector<double>         & x0,
                                                        const std::vector<double>         & y0,
                                                        std::vector<std::vector<double> > & z0,
                                                        std::vector<std::vector<double> > & w0,
                                                        bool                              & complete_surface);

  size_t                 FindLowerBoundInSortedVector(const std::vector<double> & x,
                                                      const double                x_0);

  size_t                 FindUpperBoundInSortedVector(const std::vector<double> & x,
                                                      const double                x_0);

  void                   SortOrderAndRank(const std::vector<double> & x,
                                          std::vector<double>       & x_sort,
                                          std::vector<size_t>       & x_order,
                                          std::vector<size_t>       & x_rank);

  void                   BilinearInterpolation(const std::vector<double>               & x,
                                               const std::vector<double>               & y,
                                               const std::vector<std::vector<double> > & z,
                                               const std::vector<double>               & x0,
                                               const std::vector<double>               & y0,
                                               std::vector<std::vector<double> >       & z0);

  double                 Interpolate(const std::vector<double>               & x,
                                     const std::vector<double>               & y,
                                     const std::vector<std::vector<double> > & z,
                                     const double                              x0,
                                     const double                              y0);

  void                   Estimate2DVariance(const std::vector<std::vector<float> >  & blocked_logs,
                                            const std::vector<std::vector<double> > & trend_cube_sampling,
                                            const std::vector<std::vector<double> > & s1,
                                            const std::vector<std::vector<double> > & s2,
                                            const std::vector<std::vector<double> > & trend,
                                            std::vector<std::vector<double> >       & var,
                                            std::string                             & errTxt);

  void                   KernelSmoother2DSurface(const std::vector<double>         & x,
                                                 const std::vector<double>         & y,
                                                 const std::vector<double>         & z,
                                                 const std::vector<double>         & w,
                                                 const double                        bandwidth_x,
                                                 const double                        bandwidth_y,
                                                 const std::vector<double>         & x0,
                                                 const std::vector<double>         & y0,
                                                 std::vector<std::vector<double> > & z0,
                                                 std::vector<std::vector<double> > & w0,
                                                 bool                              & complete_surface);

  void                   SmoothTrendWithLocalLinearRegression(std::vector<double>    & trend,
                                                              const std::vector<int> & count,
                                                              const int              & iWells);

  void                   ReadTrend2DPlainAscii(const std::string     & file_name,
                                               std::string           & /*errText*/,
                                               NRLib::Grid2D<double> & trend2d);

  RegularSurface<double> ResampleTrend2D(const RegularSurface<double> & surface,
                                         const std::vector<double>    & x,
                                         const std::vector<double>    & y,
                                         const bool                   & transpose);


}
#endif
