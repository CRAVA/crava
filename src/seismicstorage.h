/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SEISMICSTORAGE_H
#define SEISMICSTORAGE_H

#include <math.h>
#include <string>

//#include "src/commondata.h"
#include "nrlib/segy/segy.hpp"
#include "src/definitions.h"

class SeismicStorage
{
public:
  SeismicStorage();

  SeismicStorage(std::string   filename,
                 int           seismic_type,
                 float         angle,
                 NRLib::SegY * segy);

  SeismicStorage(std::string   filename,
                 int           seismic_type,
                 float         angle,
                 StormContGrid * storm_grid);

  ~SeismicStorage();

  enum          seismicTypes{SEGY,
                             STORM,
                             SGRI,
                             UNKNOWN};

  //GET FUNCTIONS
  std::string   getFileName()    const { return file_name_ ; }
  int           getSeismicType() const { return seismic_type_ ; }
  float         getAngle()       const { return angle_; }
  NRLib::SegY * getSegY()        const { return segy_; }
  StormContGrid * getStorm()       const { return storm_grid_; }

  //SET FUNCTIONS
  void          setFileName(std::string file_name)      { file_name_ = file_name ; }
  void          setSeismicType(int seismic_type)        { seismic_type_ = seismic_type ;}
  void          setAngle(float angle)                   { angle_ = angle ;}
  void          setSegY(NRLib::SegY * segy)             { segy_ = segy ;}
  void          setStorm(StormContGrid * stormgrid)     { storm_grid_ = stormgrid ;}

private:
  std::string   file_name_;
  float         angle_;
  int           seismic_type_; /// Se seismicTypes above

  NRLib::SegY * segy_;
  StormContGrid * storm_grid_;

};

#endif
