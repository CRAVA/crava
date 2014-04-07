/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/seismicstorage.h"
#include "src/definitions.h"
//#include "src/commondata.h"
#include "src/fftgrid.h"

SeismicStorage::SeismicStorage()
{
}

SeismicStorage::SeismicStorage(std::string   file_name,
                               int           seismic_type,
                               float         angle,
                               NRLib::SegY * segy)
 : file_name_(file_name),
   angle_(angle),
   seismic_type_(seismic_type),
   segy_(segy)
   //storm_grid_(NULL),
   //fft_grid_(NULL)
{
}

SeismicStorage::SeismicStorage(std::string     file_name,
                               int             seismic_type,
                               float           angle,
                               StormContGrid * storm_grid)
 : file_name_(file_name),
   angle_(angle),
   seismic_type_(seismic_type),
   storm_grid_(storm_grid)
   //segy_(NULL),
   //fft_grid_(NULL)
{
}

SeismicStorage::SeismicStorage(std::string     file_name,
                               int             seismic_type,
                               float           angle,
                               FFTGrid       * fft_grid)
 : file_name_(file_name),
   angle_(angle),
   seismic_type_(seismic_type),
   fft_grid_(fft_grid)
   //segy_(NULL),
   //storm_grid_(NULL)
{
}

SeismicStorage::~SeismicStorage()
{
  //if (segy_ != NULL)
  //  delete segy_;
  //if (storm_grid_ != NULL)
  //  delete storm_grid_;
  //if (fft_grid_ != NULL)
  //  delete fft_grid_;
}

//std::vector<float>
//SeismicStorage::GetTraceData(int index) const
//{
//  //Return a vector(float) of trace data for trace index.
//
//  std::vector<float> trace_data;
//
//  if (seismic_type_ == SEGY) {
//    NRLib::SegYTrace * segy_tmp = segy_->getTrace(index);
//
//    if (segy_tmp != NULL) {
//      size_t start = segy_tmp->GetStart();
//      size_t end = segy_tmp->GetEnd();
//      for (size_t i = start; i < end; i++)
//        trace_data.push_back(segy_tmp->GetValue(i));
//    }
//  }
//  else if (seismic_type_ == STORM || seismic_type_ == SGRI) {
//
//    double x_tmp = 0.0;
//    double y_tmp = 0.0;
//    double z_tmp = 0.0;
//
//    size_t i;
//    size_t j;
//    size_t k;
//
//    storm_grid_->GetIJK(index, i, j, k);
//
//    for (size_t kk = 0; kk < storm_grid_->GetNK(); kk++) {
//      storm_grid_->FindCenterOfCell(i, j, kk, x_tmp, y_tmp, z_tmp);
//      trace_data.push_back(storm_grid_->GetValueClosestInZ(x_tmp, y_tmp, z_tmp));
//    }
//  }
//  else { //FFTGrid
//    double x_tmp = 0.0;
//    double y_tmp = 0.0;
//    double z_tmp = 0.0;
//
//    size_t i;
//    size_t j;
//    size_t k;
//
//    fft_grid_->getRealTrace(i,j);
//
//  }
//
//  return trace_data;
//}

int SeismicStorage::GetNx() const {

  int nx  = 0;

  if (seismic_type_ == SEGY)
    nx = segy_->GetGeometry()->GetNx();
  else if (seismic_type_ == STORM || seismic_type_ == SGRI)
    nx = storm_grid_->GetNI();
  else if (seismic_type_ == FFTGRID)
    nx = fft_grid_->getNxp();

  return nx;
}

int SeismicStorage::GetNy() const {

  int ny  = 0;

  if (seismic_type_ == SEGY)
    ny = segy_->GetGeometry()->GetNy();
  else if (seismic_type_ == STORM || seismic_type_ == SGRI)
    ny = storm_grid_->GetNJ();
  else if (seismic_type_ == FFTGRID)
    ny = fft_grid_->getNyp();

  return ny;
}

int SeismicStorage::GetNz() const {

  int nz  = 0;

  if (seismic_type_ == SEGY)
    nz = segy_->GetNz();
  else if (seismic_type_ == STORM || seismic_type_ == SGRI)
    nz = storm_grid_->GetNK();
  else if (seismic_type_ == FFTGRID)
    nz = fft_grid_->getNzp();

  return nz;
}

void
SeismicStorage::GetSparseTraceData(std::vector<std::vector<float> > & trace_data,
                                   std::vector<float>               & trace_length,
                                   int                                n)
{
  //Return a vector(float) of trace data with n elements.
  //SEGY: n traces. STORM/SGRI: n grid-elements.

  if (seismic_type_ == SEGY) {
    int n_traces = segy_->GetNTraces();

    if (n > n_traces)
      n = n_traces;

    trace_length.resize(n, 1);
    trace_data.resize(n);

    for (int i = 0; i < n; i++) {

      int trace_index = i*(static_cast<int>(n_traces / n));

      NRLib::SegYTrace * segy_tmp = segy_->getTrace(trace_index);

      if (segy_tmp != NULL) {
        size_t start = segy_tmp->GetStart();
        size_t end = segy_tmp->GetEnd();
        for (size_t j = start; j < end; j++)
          trace_data[i].push_back(segy_tmp->GetValue(j));

      }
    }
  }
  else if (seismic_type_ == STORM || seismic_type_ == SGRI) {

    double x_tmp = 0.0;
    double y_tmp = 0.0;
    double z_tmp = 0.0;

    size_t index_i     = 0;
    size_t index_j     = 0;
    int    trace_index = 0;

    int n_elements = static_cast<int>(std::sqrt(static_cast<double>(n)));

    for (int i = 0; i < n_elements; i++) {

      index_i = static_cast<size_t>(i*storm_grid_->GetNI()/n_elements);
      if (index_i >= storm_grid_->GetNI())
        index_i = storm_grid_->GetNI() -1;

      for (int j = 0; j < n_elements; j++) {

        index_j = static_cast<size_t>(j*storm_grid_->GetNJ()/n_elements);
        if (index_j >= storm_grid_->GetNJ())
          index_j = storm_grid_->GetNJ()-1;

        for (size_t k = 0; k < storm_grid_->GetNK(); k++) {
          storm_grid_->FindCenterOfCell(index_i, index_j, k, x_tmp, y_tmp, z_tmp);
          trace_data[trace_index].push_back(storm_grid_->GetValueClosestInZ(x_tmp, y_tmp, z_tmp));
        }

        //Store length
        double top = storm_grid_->GetTopSurface().GetZ(x_tmp, y_tmp);
        double bot = storm_grid_->GetBotSurface().GetZ(x_tmp, y_tmp);
        trace_length.push_back(std::abs(static_cast<float>(bot-top)));

        trace_index++;
      }
    }
  }
  else { //FFTGrid

    int index_i     = 0;
    int index_j     = 0;
    int trace_index = 0;

    int n_elements = static_cast<int>(std::sqrt(static_cast<double>(n)));

    trace_length.resize(n, 1);

    for (int i = 0; i < n_elements; i++) {

      index_i = i*fft_grid_->getNx()/n_elements;
      if (index_i >= fft_grid_->getNx())
        index_i = fft_grid_->getNx() - 1;

      for (int j = 0; j < n_elements; j++) {

        index_j = j*fft_grid_->getNy()/n_elements;
        if (index_j >= fft_grid_->getNy())
          index_j = fft_grid_->getNy()-1;

        trace_data[trace_index] = fft_grid_->getRealTrace2(index_i, index_j);

        trace_index++;

      }
    }
  }
}

std::vector<float>
SeismicStorage::GetRealTrace(const Simbox * estimation_simbox,
                             int i,
                             int j) const
{
  std::vector<float> value;
  double x_tmp; double y_tmp; double z_tmp;

  if (seismic_type_ == SEGY) {
    for (size_t k = 0; k < segy_->GetNz(); k++) {
      estimation_simbox->getCoord(i, j, k, x_tmp, y_tmp, z_tmp);
      value.push_back(segy_->GetValue(x_tmp, y_tmp, z_tmp)); ///H Is this correct? Want to get det same value as FFTGrid::getRealTrace2(int i, int j)
    }
  }
  else if (seismic_type_ == STORM || seismic_type_ == SGRI) {
    for (size_t k = 0; k < storm_grid_->GetNK(); k++) {
      storm_grid_->FindCenterOfCell(i, j, k, x_tmp, y_tmp, z_tmp);
      value.push_back(storm_grid_->GetValueClosestInZ(x_tmp, y_tmp, z_tmp));
    }
  }
  else { //FFTGrid
    value = fft_grid_->getRealTrace2(i, j);
  }

  return value;
}

float
SeismicStorage::GetRealTraceValue(const Simbox * estimation_simbox,
                                  int i,
                                  int j,
                                  int k) const
{
  float value;
  double x_tmp; double y_tmp; double z_tmp;

  if (seismic_type_ == SEGY) {
    estimation_simbox->getCoord(i, j, k, x_tmp, y_tmp, z_tmp);
    value = segy_->GetValue(x_tmp, y_tmp, z_tmp);
  }
  else if (seismic_type_ == STORM || seismic_type_ == SGRI) {
    storm_grid_->FindCenterOfCell(i, j, k, x_tmp, y_tmp, z_tmp);
    value = storm_grid_->GetValueClosestInZ(x_tmp, y_tmp, z_tmp);
  }
  else {
    value = fft_grid_->getRealValue(i, j, k);
  }

  return value;
}


void
SeismicStorage::FindSimbox(const Simbox & full_inversion_simbox,
                           double         lz_limit,
                           Simbox       & seismic_simbox,
                           std::string  & err_txt) const
{
  switch(seismic_type_) {
    case SEGY: {
      const NRLib::SegyGeometry * geometry = segy_->GetGeometry();
      seismic_simbox.setArea(geometry, err_txt);
      seismic_simbox.setILXL(geometry);
      double z_top  = static_cast<double>(segy_->GetTop());
      double z_base = static_cast<double>(z_top+segy_->GetDz()*segy_->GetNz());
      double x_min, x_max, y_min, y_max;
      seismic_simbox.getMinAndMaxXY(x_min, x_max, y_min, y_max);
      NRLib::RegularSurface<double> top(x_min, y_min, x_max-x_min, y_max-y_min, 2, 2, z_top);
      NRLib::RegularSurface<double> base(x_min, y_min, x_max-x_min, y_max-y_min, 2, 2, z_base);
      seismic_simbox.setDepth(top, base, segy_->GetNz(), true);
      seismic_simbox.calculateDz(lz_limit, err_txt);
      seismic_simbox.SetNoPadding();
      break;
    }
    case STORM:
    case SGRI:
      seismic_simbox.setArea(storm_grid_,
                             static_cast<int>(storm_grid_->GetNI()),
                             static_cast<int>(storm_grid_->GetNJ()),
                             err_txt);
      seismic_simbox.setDepth(storm_grid_->GetTopSurface(), storm_grid_->GetBotSurface(), storm_grid_->GetNK(), true);
      seismic_simbox.calculateDz(lz_limit, err_txt);
      seismic_simbox.SetNoPadding();
      break;
    default:
      seismic_simbox.setArea(&full_inversion_simbox,
                             static_cast<int>(full_inversion_simbox.getnx()),
                             static_cast<int>(full_inversion_simbox.getny()),
                             err_txt);
      seismic_simbox.CopyAllPadding(full_inversion_simbox, lz_limit, err_txt);
  }
}
