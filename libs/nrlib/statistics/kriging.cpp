// $Id: kriging.cpp 1249 2014-02-26 10:52:16Z vigsnes $

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

#include <vector>


#include "kriging.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/variogram/covgrid2d.hpp"
#include "nrlib/statistics/krigingdata2d.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"

////for debug print
//#include <nrlib/surface/regularsurface.hpp>
//#include <nrlib/surface/surfaceio.hpp>

#ifdef PARALLEL
#include <omp.h>
#endif

using namespace NRLib;

void Kriging::Krig1D(std::vector<double>       &field,
                     const std::vector<bool>   &is_known,
                     const std::vector<double> &obs,
                     double                     dx,
                     const Variogram           &vario)
{
  assert(field.size() == is_known.size());
  assert(obs.size() <= field.size());
  int n_obs = static_cast<int>(obs.size());

  NRLib::Vector residual(static_cast<int> (n_obs));
  int j = 0;
  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true){
      residual(j) = (obs[j]-field[i]);
      j++;
    }

  }

  NRLib::SymmetricMatrix K(n_obs);
  NRLib::Vector k_vec(n_obs);
  std::vector<double> x_known;
  std::vector<double> x_unknown;
  NRLib::Vector K_res(n_obs);

  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true)
      x_known.push_back(i*dx);
    else
      x_unknown.push_back(i*dx);
  }


  for(int i = 0; i < n_obs; i++)
    for(int j = 0; j <= i; j++)
      K(j,i) = vario.GetCorr(x_known[i]-x_known[j]);

  CholeskyInvert(K);
  K_res = K * residual;
  j = 0;
  size_t kk = 0;
  for(size_t k = 0; k < field.size(); k++){
    if(is_known[k] == false){
      for(int i = 0; i < n_obs; i++)
        k_vec(i) = vario.GetCorr(x_known[i]-x_unknown[j]);

      double prod = k_vec * K_res;
      field[k] += prod;
      j++;
    }
    else{
      field[k] = obs[kk];
      kk++;
    }
  }

}


//----------------------------------------------------------
void Kriging::Krig2D(Grid2D<double>      & trend,
                     const KrigingData2D & kriging_data,
                     const CovGrid2D     & cov)
//----------------------------------------------------------
{
  double dx = cov.GetDX();
  double dy = cov.GetDY();
  double Rx = cov.GetRangeX();
  double Ry = cov.GetRangeY();

  // blocking
  size_t nxb, nyb, n_blocks_x, n_blocks_y;
  std::vector<size_t> block_x, block_y;
  std::vector<KrigingData2D> kriging_data_blocks;

  // select constants
  size_t max_data_in_range = 200;
  size_t min_data_in_block = 15;
  size_t min_data_in_range = 25;
  int n_threads            = 1;

#ifdef PARALLEL
  n_threads = omp_get_num_procs();
#endif

  SetBlock(trend,
           kriging_data,
           Rx, Ry, dx, dy,
           n_threads,
           max_data_in_range,
           min_data_in_block,
           min_data_in_range,
           kriging_data_blocks,
           block_x, block_y,
           nxb,
           nyb,
           n_blocks_x,
           n_blocks_y);

  // kriging in parallel
  int nx = static_cast<int>(trend.GetNI());
  int ny = static_cast<int>(trend.GetNJ());
  Grid2D<double> filled(nx,ny,0);
  Grid2D<double> trend_orig(trend);

#ifdef PARALLEL
  int  chunk_size = 1;
#pragma omp parallel
#pragma omp master
  //int n_processors = omp_get_num_procs();
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
  for (size_t i = 0; i < kriging_data_blocks.size(); ++i){
    Krig2DBlock(trend_orig,
                kriging_data_blocks[i],
                cov,
                block_x[i],
                block_y[i],
                nxb,
                nyb,
                n_blocks_x,
                n_blocks_y,
                trend,
                filled);
  }
}

//------------------------------------------
bool compare(const std::pair<int, int> & i1,
             const std::pair<int, int> & i2)
//------------------------------------------
{
  return (i1.second > i2.second);
}

//-----------------------------------------------------------------------
void Kriging::SetBlock(const Grid2D<double>       & trend,
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
                       size_t                     & n_blocks_y)
//-----------------------------------------------------------------------
{

  FindOptimalBlock(trend,
                   Rx,
                   Ry,
                   dx,
                   dy,
                   kriging_data.GetNumberOfData(),
                   n_threads,
                   n_blocks_x,
                   n_blocks_y); // returns n_block_x and n_block_y

  size_t nx = trend.GetNI();
  size_t ny = trend.GetNJ();

  // number of grid cells in block in x and y direction
  nxb = static_cast<size_t>(floor(nx/static_cast<double>(n_blocks_x)));
  nyb = static_cast<size_t>(floor(ny/static_cast<double>(n_blocks_y)));

  size_t n_blocks = n_blocks_x * n_blocks_y;
  kriging_data_blocks.resize(n_blocks);
  block_x.resize(n_blocks);
  block_y.resize(n_blocks);

  FindOptimalDataInBlock(max_data_in_range,
                         min_data_in_block,
                         min_data_in_range,
                         trend,
                         kriging_data,
                         Rx,
                         Ry,
                         dx,
                         dy,
                         n_blocks_x,
                         n_blocks_y,
                         kriging_data_blocks,
                         block_x,
                         block_y,             // return kriging_data_blocks, block_x and block_y
                         false);              // if true, keep data within one range - not try to optimize.
}


//----------------------------------------------------------------------
void Kriging::FindOptimalBlock(const Grid2D<double> & trend,
                               const double         & Rx,
                               const double         & Ry,
                               const double         & dx,
                               const double         & dy,
                               const size_t         & n_data,
                               const size_t         & n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y)
//----------------------------------------------------------------------
{
  // optimize time with respect to PETROSIM document

  const size_t nx      = trend.GetNI();
  const size_t ny      = trend.GetNJ();
  const double xlength = nx*dx;
  const double ylength = ny*dy;

  const double Vcell   = dx*dy;
  const double V       = xlength * ylength;
  const double R       = Rx * Ry;
  double       mintime = std::numeric_limits<double>::infinity();

  size_t max_points = 200;
  n_blocks_x = 10;
  n_blocks_y = 10;
  double factor_x = 1.0;
  double factor_y = 1.0;

  if (n_data > max_points && n_data > 0){

    size_t min_grid_cells = 5;                                 // minimum size of blocks
    size_t min_blocks     = 5;                                 // minimum number of blocks in each direction
    double Smin           = std::max((min_grid_cells * dx / Rx),(min_grid_cells * dy / Ry));
    double Smax           = std::min((xlength / (min_blocks * Rx)),(ylength / (min_blocks * Ry)));
    double minatS         = Smax;
    size_t NS             = 200;

    for (size_t k = 0; k <= NS; ++k){                          // interval for S = [Smin, Smax] is divided into NS equal parts
      double       S       = Smin + k * (Smax - Smin) / NS;
      double       v       = std::max(R*S*S, Vcell);           // volume of segment
      double       vd      = R*(2 + S)*(2 + S);                // volume of neighbourhood
      size_t       Ns      = static_cast<size_t>(V/v);
      size_t       nd      = static_cast<size_t>(n_data * vd / V);
      double       totTime = OverallTime(nd, Ns, nx, ny);
      if (totTime <= mintime){
        n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(S*Rx)));
        n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(S*Ry)));
        if ((n_blocks_x * n_blocks_y) >= n_threads){
          minatS  = S;
          mintime = totTime;
        }
      }
    }
    n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(minatS*Rx)));
    n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(minatS*Ry)));
  }
}


//------------------------------------
double Kriging::OverallTime(size_t n,
                            size_t Ns,
                            size_t nx,
                            size_t ny)
//------------------------------------
{
  // Relative times (i.e. divided by the time it takes to assemble one k vector)
  // Relative times esimtated based on 6 tests
  double       V     = static_cast<double>(nx * ny);
  double       tChol = 0.00032;
  double       tK    = 0.75088;
  double       tSol  = 0.12702;

  double       Tchol = tChol * Ns * n*n*n;
  double       TK    = tK    * Ns * n*n;
  double       Tsol  = tSol  * Ns * n*n;
  double       Tk    = 1.0   * V  * n;

  double       T     = Tchol + TK + Tsol + Tk;

  return T;
}


//-------------------------------------------------------------------------------------------
void Kriging::FindOptimalDataInBlock(const size_t               & max_data_in_range,
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
                                     std::vector<size_t>        & block_index_x,
                                     std::vector<size_t>        & block_index_y,
                                     const bool                 & all_data_within_one_range)
//-------------------------------------------------------------------------------------------
{
  // find number of data in each block and block+range
  size_t nx  = trend.GetNI();
  size_t ny  = trend.GetNJ();
  size_t nxb = static_cast<size_t>(floor(nx/static_cast<double>(n_blocks_x)));
  size_t nyb = static_cast<size_t>(floor(ny/static_cast<double>(n_blocks_y)));
  int    rx  = static_cast<int>(Rx/dx);
  int    ry  = static_cast<int>(Ry/dy);

  size_t n_blocks = n_blocks_x * n_blocks_y;
  std::vector<size_t> n_data_blocks(n_blocks);
  std::vector<size_t> n_data_blocks_range(n_blocks);

  size_t xmin,   xmax,   ymin,   ymax;
  size_t xmin_r, xmax_r, ymin_r, ymax_r;

  size_t index;
  for (size_t i = 0; i < n_blocks_x; ++i){
    xmin   = i*nxb;
    xmax   = std::min(static_cast<int>((i+1)*nxb), static_cast<int>(nx));
    xmin_r = std::max((static_cast<int>(i*nxb)-rx), 0);
    xmax_r = std::min((static_cast<int>((i+1)*nxb)+rx), static_cast<int>(nx));
    for (size_t j = 0; j < n_blocks_y; ++j){
      index  = i*n_blocks_y + j;
      ymin   = j*nyb;
      ymax   = std::min(static_cast<int>((j+1)*nyb), static_cast<int>(ny));
      ymin_r = std::max((static_cast<int>(j*nyb)-ry), 0);
      ymax_r = std::min((static_cast<int>((j+1)*nyb)+ry), static_cast<int>(ny));
      block_index_x[index] = i;
      block_index_y[index] = j;
      n_data_blocks[index] = CountDataInBlock(kriging_data,
                                              xmin,
                                              xmax,
                                              ymin,
                                              ymax);
      n_data_blocks_range[index] = CountDataInBlock(kriging_data,
                                                    xmin_r,
                                                    xmax_r,
                                                    ymin_r,
                                                    ymax_r);
    }
  }

  //loop over blocks:
  //  if # data in block+range < max_data_in_range,
  //    block ok, include one range
  //  else
  //    loop over each neighbour block (8):
  //      if # data < min_data_in_block
  //      include out to one range outside this neighbour block

  size_t index_block;
  if (n_blocks_x > 1 || n_blocks_y > 1) {
    for (size_t i = 0; i < n_blocks_x; ++i) {
      for (size_t j = 0; j < n_blocks_y; ++j) {
        index = i*n_blocks_y + j;
        //---------------------------------------------------------------------
        //          IF TOO MANY DATA IN BLOCK+RANGE
        //---------------------------------------------------------------------
        if (n_data_blocks_range[index] > max_data_in_range // too many data?
            && rx > nxb && ry > nyb                        // range too small compared to blocks?
            && all_data_within_one_range == false){        // option to ALWAYS use one range overall
          //-----------------------------------------
          // y = j
          //-----------------------------------------
          ymin = j*nyb;
          if (j == n_blocks_y - 1)
            ymax = ny;
          else
            ymax = (j+1)*nyb;
          // x = i-1
          if (i > 0){
            index_block = (i-1)*n_blocks_y + j;
            xmax = i*nxb;
            if (n_data_blocks[index_block] < min_data_in_block)
              xmin = std::max((static_cast<int>(i*nxb)-rx), 0);
            else
              xmin = (i-1)*nxb;
            AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
          }
          // x = i
          index_block = i*n_blocks_y + j;
          xmin = i*nxb;
          if (i == n_blocks_x - 1)
            xmax = nx;
          else
            xmax = (i+1)*nxb;
          AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
          // x = i+1
          if ((i+1) < n_blocks_x){
            index_block = (i+1)*n_blocks_y + j;
            xmin = (i+1)*nxb;
            if (n_data_blocks[index_block] < min_data_in_block){
              if ((i+1) == (n_blocks_x-1))
                xmax = nx;
              else
                xmax = std::min((static_cast<int>((i+1)*nxb)+rx), static_cast<int>(nx));
            }
            else{
              if ((i+1) == (n_blocks_x-1))
                xmax = nx;
              else
                xmax = (i+2)*nxb;
            }
            AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
          }
          //-----------------------------------------
          //y = j+1
          //-----------------------------------------
          if (j < (n_blocks_y-1)){
            ymin = (j+1)*nyb;
            // x = i-1
            if (i > 0){
              index_block = (i-1)*n_blocks_y + j+1;
              xmax = i*nxb;
              if (n_data_blocks[index_block] < min_data_in_block){
                if ((j+1) == (n_blocks_y-1))
                  ymax = ny;
                else
                  ymax = std::min((static_cast<int>((j+1)*nyb)+ry), static_cast<int>(ny));
                xmin = std::max((static_cast<int>(i*nxb)-rx), 0);
              }
              else {
                if ((j+1) == (n_blocks_y-1))
                  ymax = ny;
                else
                  ymax = (j+2)*nyb;
                xmin = (i-1)*nxb;
              }
              AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            }
            // x = i
            index_block = i*n_blocks_y + j+1;
            xmin = i*nxb;
            if (i == n_blocks_x - 1)
              xmax = nx;
            else
              xmax = (i+1)*nxb;
            if (n_data_blocks[index_block] < min_data_in_block){
              if ((j+1) == (n_blocks_y-1))
                ymax = ny;
              else
                ymax = std::min((static_cast<int>((j+1)*nyb)+ry), static_cast<int>(ny));
            }
            else{
              if ((j+1) == (n_blocks_y-1))
                ymax = ny;
              else
                ymax = (j+2)*nyb;
            }
            AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            // x = i+1
            if ((i+1) < n_blocks_x){
              index_block = (i+1)*n_blocks_y + j+1;
              xmin = (i+1)*nxb;
              if (n_data_blocks[index_block] < min_data_in_block){
                if ((j+1) == (n_blocks_y-1))
                  ymax = ny;
                else
                  ymax = std::min((static_cast<int>((j+1)*nyb)+ry), static_cast<int>(ny));
                if ((i+1) == (n_blocks_x-1))
                  xmax = nx;
                else
                  xmax = std::min((static_cast<int>((i+1)*nxb)+rx), static_cast<int>(nx));
              }
              else{
                if ((j+1) == (n_blocks_y-1))
                  ymax = ny;
                else
                  ymax = (j+2)*nyb;
                if ((i+1) == (n_blocks_x-1))
                  xmax = nx;
                else
                  xmax = (i+2)*nxb;
              }
              AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            }
          }
          //-----------------------------------------
          // y = j-1
          //-----------------------------------------
          if (j > 0){
            ymax = j*nyb;
            // x = i-1
            if (i > 0){
              index_block = (i-1)*n_blocks_y + j-1;
              xmax = i*nxb;
              if (n_data_blocks[index_block] < min_data_in_block){
                ymin = std::max((static_cast<int>(j*nyb)-ry), 0);
                xmin = std::max((static_cast<int>(i*nxb)-rx), 0);
              }
              else {
                ymin = (j-1)*nyb;
                xmin = (i-1)*nxb;
              }
              AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            }
            // x = i
            index_block = i*n_blocks_y + j-1;
            xmin = i*nxb;
            if (i == n_blocks_x - 1)
              xmax = nx;
            else
              xmax = (i+1)*nxb;
            if (n_data_blocks[index_block] < min_data_in_block){
              ymin = std::max((static_cast<int>(j*nyb)-ry), 0);
            }
            else
              ymin = (j-1)*nyb;
            AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            // x = i+1
            if ((i+1) < n_blocks_x){
              index_block = (i+1)*n_blocks_y + j-1;
              xmin = (i+1)*nxb;
              if (n_data_blocks[index_block] < min_data_in_block){
                ymin = std::max((static_cast<int>(j*nyb)-ry), 0);
                if ((i+1) == (n_blocks_x-1))
                  xmax = nx;
                else
                  xmax = std::min((static_cast<int>((i+1)*nxb)+rx), static_cast<int>(nx));
              }
              else{
                ymin = (j-1)*nyb;
                if ((i+1) == (n_blocks_x-1))
                  xmax = nx;
                else
                  xmax = (i+2)*nxb;
              }
              AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin, xmax, ymin, ymax);
            }
          }
        }
        //---------------------------------------------------------------------
        // if NOT too many data in block+range OR use all data within one range
        //---------------------------------------------------------------------
        else{
          double rf;
          //if too little data in block+range, increase to 1.5 range
          if (n_data_blocks_range[index] < min_data_in_range)
            rf = 1.5; //corresponding to 0.01
          else
            rf = 1;   //corresponding to 0.05
          xmin_r = std::max(static_cast<int>(i * nxb - (rx * rf)), 0);
          xmax_r = std::min(static_cast<int>((i + 1) * nxb + (rx * rf)), static_cast<int>(nx));
          ymin_r = std::max(static_cast<int>(j * nyb - (ry * rf)), 0);
          ymax_r = std::min(static_cast<int>((j + 1) * nyb + (ry * rf)), static_cast<int>(ny));
          AddDataToBlock(kriging_data, kriging_data_blocks[index], xmin_r, xmax_r, ymin_r, ymax_r);
        }
      }
    }
  }
  else // only one block
    kriging_data_blocks[0] = kriging_data;

  // sort blocks after largest number of data, for parallelizing
  typedef std::pair<int, int> Pair;
  std::vector<Pair> nums;

  for (size_t i = 0 ; i < n_blocks ; ++i) {
    kriging_data_blocks[i].FindMeanValues();
    nums.push_back(Pair(i, kriging_data_blocks[i].GetNumberOfData()));
  }

  std::sort(nums.begin(), nums.end(), compare);

  std::vector<KrigingData2D>  kriging_data_blocks_sort;
  std::vector<size_t>         block_x_sort;
  std::vector<size_t>         block_y_sort;

  for (size_t i = 0 ; i < n_blocks ; ++i) {
    kriging_data_blocks_sort.push_back( kriging_data_blocks[nums[i].first] );
    block_x_sort.push_back( block_index_x[nums[i].first] );
    block_y_sort.push_back( block_index_y[nums[i].first] );
  }
  kriging_data_blocks = kriging_data_blocks_sort;
  block_index_x = block_x_sort;
  block_index_y = block_y_sort;
}

//----------------------------------------------------------------------
size_t Kriging::CountDataInBlock(const KrigingData2D & kriging_data,
                                 size_t xmin, size_t xmax, size_t ymin, size_t ymax)
//----------------------------------------------------------------------
{
  size_t index_i, index_j;
  size_t count = 0;
  for (int k = 0; k < kriging_data.GetNumberOfData(); ++k){
    index_i = kriging_data.GetIndexI(k);
    index_j = kriging_data.GetIndexJ(k);
    if (index_i <= xmax)
      if (index_i >= xmin)
        if (index_j <= ymax)
          if (index_j >= ymin)
            count++;
  }
  return count;
}

//------------------------------------------------------------------------------
void Kriging::AddDataToBlock(const KrigingData2D & kriging_data,
                             KrigingData2D       & kriging_data_block,
                             size_t xmin, size_t xmax, size_t ymin, size_t ymax)
//------------------------------------------------------------------------------
{
  size_t index_i, index_j;
  for (int k = 0; k < kriging_data.GetNumberOfData(); ++k){
    index_i = static_cast<size_t>(kriging_data.GetIndexI(k));
    index_j = static_cast<size_t>(kriging_data.GetIndexJ(k));
    if (index_i <= xmax)
      if (index_i >= xmin)
        if (index_j <= ymax)
          if (index_j >= ymin)
            kriging_data_block.AddData(static_cast<int>(index_i), static_cast<int>(index_j), kriging_data.GetData(k));
  }
}


//----------------------------------------------------------------
void Kriging::Krig2DBlock(const Grid2D<double> & trend_orig,
                          const KrigingData2D  & kriging_data,
                          const CovGrid2D      & cov,
                          const size_t         & block_x,
                          const size_t         & block_y,
                          const size_t         & nxb,
                          const size_t         & nyb,
                          const size_t         & n_blocks_x,
                          const size_t         & n_blocks_y,
                          Grid2D<double>       & trend,
                          Grid2D<double>       & filled)
//----------------------------------------------------------------
{
  //
  // This routine returns z(x) = m(x) + k(x)K^{-1}(d - m).
  //
  int                      md     = kriging_data.GetNumberOfData();
  const std::vector<int> & indexi = kriging_data.GetIndexI();
  const std::vector<int> & indexj = kriging_data.GetIndexJ();
  std::vector<float>       data   = kriging_data.GetData();   // Take an editable copy

  int                      nx     = static_cast<int>(trend.GetNI());
  int                      ny     = static_cast<int>(trend.GetNJ());


  // find min and max of block
  int xmin = std::max(static_cast<int>(block_x*nxb), 0);
  int xmax;
  if (block_x == (n_blocks_x-1))
    xmax = nx;
  else
    xmax = std::min(static_cast<int>((block_x+1)*nxb),nx);

  int ymin = std::max(static_cast<int>(block_y*nyb), 0);
  int ymax;
  if (block_y == (n_blocks_y-1))
    ymax = ny;
  else
    ymax = std::min(static_cast<int>((block_y+1)*nyb),ny);


  if (md > 0 && md <= nx*ny) {

    SymmetricMatrix K(md);
    Vector residual(md);
    Vector k(md);
    Vector x(md);
    Vector Kk(md);

    // subtract trend from original trend grid
    for (int i = 0 ; i < md ; i++)
      residual(i) = data[i] - static_cast<float>(trend_orig(indexi[i],indexj[i]));

    FillKrigingMatrix(K, cov, indexi, indexj);

    CholeskySolve(K, residual, x);

    for(size_t i=0;i<md;i++) {
      if (indexi[i] < xmax)
        if (indexi[i] >= xmin)
          if (indexj[i] < ymax)
            if (indexj[i] >= ymin) {
               trend(indexi[i],indexj[i]) += residual(static_cast<int>(i));
               filled(indexi[i],indexj[i]) = 1.0;
            }
    }

    for (int i = xmin ; i < xmax ; i++) {
      for (int j = ymin ; j < ymax ; j++) {
        if(!(filled(i,j) > 0.0)){ // if this is not a datapoint
          FillKrigingVector(k, cov, indexi, indexj, i, j);
          trend(i,j) += k * x;
        }
      }
    }
  }
}

//---------------------------------------------------------------
void Kriging::FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj)
//---------------------------------------------------------------
{
  int n = K.dim();
  //Grid2D<double> K_test(n, n, 0); //for debug print

  for(int i=0 ; i < n  ; i++) {
    for(int j=0 ; j <= i ; j++) {
      int deltai = indexi[i] - indexi[j];
      int deltaj = indexj[i] - indexj[j];
      K(j,i) = static_cast<double>(cov.GetCov(deltai,deltaj));
      //K_test(i,j) = static_cast<double>(cov.GetCov(deltai,deltaj)); //for debug print
    }
  }
  //RegularSurface<double> K_surface(0,0,n,n,K_test); //for debug print
  //std::string filename = "K_surface.irap";
  //WriteIrapClassicAsciiSurf(K_surface, 0.0, filename);
}

//--------------------------------------------------------------
void Kriging::FillKrigingVector(NRLib::Vector          & k,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj,
                                int i, int j)
//--------------------------------------------------------------
{
  for(int ii=0 ; ii < k.length() ; ii++) {
    int deltai = indexi[ii] - i;
    int deltaj = indexj[ii] - j;
    k(ii) = static_cast<double>(cov.GetCov(deltai,deltaj));
  }
}
