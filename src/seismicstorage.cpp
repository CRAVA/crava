/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

//#include "src/commondata.h"
#include "src/seismicstorage.h"
#include "src/definitions.h"

SeismicStorage::SeismicStorage()
{
}

SeismicStorage::SeismicStorage(std::string   filename,
                               int           seismic_type,
                               float         angle,
                               NRLib::SegY * segy)
 : file_name_(filename),
   angle_(angle),
   seismic_type_(seismic_type),
   segy_(segy)
{
}

SeismicStorage::SeismicStorage(std::string   filename,
                               int           seismic_type,
                               float         angle,
                               StormContGrid * storm_grid)
 : file_name_(filename),
   angle_(angle),
   seismic_type_(seismic_type),
   storm_grid_(storm_grid)
{
}


SeismicStorage::~SeismicStorage()
{
}

//void
//SeismicStorage::SeismicStorageSegy(std::string   filename,
//                                   int           seismic_type,
//                                   float         angle,
//                                   NRLib::SegY * segy)
//
//{
//  setFileName(filename);
//  setAngle(angle);
//  setSeismicType(seismic_type);
//  setSegY(segy);
//}

//void
//SeismicStorage::SeismicStorageStorm(std::string   filename,
//                                    int           seismic_type,
//                                    float         angle,
//                                    StormContGrid storm_grid)
//
//{
//  setFileName(filename);
//  setAngle(angle);
//  setSeismicType(seismic_type);
//  setStorm(storm_grid);
//}




