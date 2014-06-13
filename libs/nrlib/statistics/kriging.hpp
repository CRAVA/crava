// $Id: kriging.hpp 1249 2014-02-26 10:52:16Z vigsnes $

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

#ifndef NRLIB_STATISTICS_KRIGING_HPP
#define NRLIB_STATISTICS_KRIGING_HPP

#include <cstdlib>
#include <vector>

#include "../flens/nrlib_flens.hpp"
#include "../variogram/variogram.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/variogram/covgrid2d.hpp"
#include "nrlib/statistics/krigingdata2d.hpp"


namespace NRLib {
class Kriging
{
public:
  static void Krig1D(std::vector<double> &field,
              const std::vector<bool>    &is_known,
              const std::vector<double>  &obs,
              double                      dx,
              const Variogram            &vario);

  static void Krig2D(Grid2D<double>      & trend,
                     const KrigingData2D & krigingData,
                     const CovGrid2D     & cov);
private:
  static void SetBlock(const Grid2D<double>       & trend,
                       const KrigingData2D        & kriging_data,
                       const double               & Rx,
                       const double               & Ry,
                       const double               & dx,
                       const double               & dy,
                       const size_t               & n_threads,
                       const size_t               & max_data_in_range,
                       const size_t               & min_data_in_block,
                       const size_t               & min_data_in_range,
                       std::vector<KrigingData2D> & kriging_data_blocks,
                       std::vector<size_t>        & block_x,
                       std::vector<size_t>        & block_y,
                       size_t                     & nxb,
                       size_t                     & nyb,
                       size_t                     & n_blocks_x,
                       size_t                     & n_blocks_y);

  static void FindOptimalDataInBlock(const size_t               & max_data_in_range,
                                     const size_t               & min_data_in_block,
                                     const size_t               & min_data_in_range,
                                     const Grid2D<double>       & trend,
                                     const KrigingData2D        & kriging_data,
                                     const double               & Rx,
                                     const double               & Ry,
                                     const double               & dx,
                                     const double               & dy,
                                     const size_t               & n_blocks_x,
                                     const size_t               & n_blocks_y,
                                     std::vector<KrigingData2D> & kriging_data_blocks,
                                     std::vector<size_t>        & block_x,
                                     std::vector<size_t>        & block_y,
                                     const bool                 & all_data_within_one_range = false);

  static void FindOptimalBlock(const Grid2D<double> & trend,
                               const double         & Rx,
                               const double         & Ry,
                               const double         & dx,
                               const double         & dy,
                               const size_t         & n_data,
                               const size_t         & n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y);

  static size_t CountDataInBlock(const KrigingData2D & kriging_data,
                                 size_t xmin, size_t xmax, size_t ymin, size_t ymax);

  static void AddDataToBlock(const KrigingData2D & kriging_data,
                            KrigingData2D       & kriging_data_block,
                            size_t xmin, size_t xmax, size_t ymin, size_t ymax);

  static double OverallTime(size_t n,
                            size_t Ns,
                            size_t nx,
                            size_t ny);

  static void Krig2DBlock(const Grid2D<double> & trend_orig,
                          const KrigingData2D  & kriging_data,
                          const CovGrid2D      & cov,
                          const size_t         & block_x,
                          const size_t         & block_y,
                          const size_t         & nxb,
                          const size_t         & nyb,
                          const size_t         & n_blocks_x,
                          const size_t         & n_blocks_y,
                          Grid2D<double>       & trend,
                          Grid2D<double>       & filled);

  static void FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj);

  static void FillKrigingVector(NRLib::Vector          & k,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj,
                                int i,
                                int j);
};
}
#endif

