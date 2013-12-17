/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SEISMICSTORAGE_H
#define SEISMICSTORAGE_H

#include <math.h>
#include <string>

#include "nrlib/segy/segy.hpp"
#include "src/definitions.h"
#include "src/simbox.h"

class FFTGrid;

class SeismicStorage
{
public:
  SeismicStorage();

  SeismicStorage(std::string   filename,
                 int           seismic_type,
                 float         angle,
                 NRLib::SegY * segy);

  SeismicStorage(std::string     filename,
                 int             seismic_type,
                 float           angle,
                 StormContGrid * storm_grid);

  SeismicStorage(std::string     file_name,
                 int             seismic_type,
                 float           angle,
                 FFTGrid       * fft_grid);

  ~SeismicStorage();

  enum          seismicTypes{SEGY,
                             STORM,
                             SGRI,
                             FFTGRID, //When reading from crava-format
                             UNKNOWN};

  //GET FUNCTIONS
  std::string     GetFileName()    const { return file_name_    ;}
  int             GetSeismicType() const { return seismic_type_ ;}
  float           GetAngle()       const { return angle_        ;}
  NRLib::SegY *   GetSegY()        const { return segy_         ;}
  StormContGrid * GetStorm()       const { return storm_grid_   ;}
  FFTGrid       * GetFFTGrid()     const { return fft_grid_     ;}

  void            GetSparseTraceData(std::vector<std::vector<float> > & trace_data,
                                     std::vector<float>               & trace_length,
                                     int                                n);

  //std::vector<float> GetTraceData(int index) const;

  std::vector<float> GetRealTrace(const Simbox * estimation_simbox,
                                  int i,
                                  int j) const;

  float GetRealTraceValue(const Simbox * estimation_simbox,
                          int i,
                          int j,
                          int k) const;

  int GetNx() const;
  int GetNy() const;
  int GetNz() const;

  //SET FUNCTIONS
  void            SetFileName(std::string file_name)      { file_name_    = file_name ; }
  void            SetSeismicType(int seismic_type)        { seismic_type_ = seismic_type ;}
  void            SetAngle(float angle)                   { angle_        = angle ;}
  void            SetSegY(NRLib::SegY * segy)             { segy_         = segy ;}
  void            SetStorm(StormContGrid * stormgrid)     { storm_grid_   = stormgrid ;}



private:
  std::string   file_name_;
  float         angle_;
  int           seismic_type_; /// Se seismicTypes above

  NRLib::SegY   * segy_;
  StormContGrid * storm_grid_;
  FFTGrid       * fft_grid_;

};

#endif
