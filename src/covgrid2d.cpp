/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdarg.h>
#include <ostream>

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "src/covgrid2d.h"

CovGrid2D::CovGrid2D(Vario * vario,
                     int      nx,
                     int      ny,
                     double   dx,
                     double   dy)
  : cov_(2*nx*ny),
    nx_(nx),
    ny_(ny),
    dx_(dx),
    dy_(dy)
{
  for(int i=0;i<nx;i++)
    for(int j=-ny;j<ny;j++)
    {
      int index = i*2*ny + ny + j;
      cov_[index] = vario->corr(i*float(dx), j*float(dy));
    }
}

float
CovGrid2D::getCov(int deltai, int deltaj) const
{
  if(deltai<0)
  {
    deltai = -deltai;
    deltaj = -deltaj;
  }
  int index = deltai*2*ny_ + deltaj + ny_;
  return cov_[index];
}

void
CovGrid2D::writeToFile(const std::string & name) const
{
  //
  // Write grid using an ASCII Irap Classic surface format
  //
  std::ofstream file;
  NRLib::OpenWrite(file, name);
  file.precision(14);
  file << 2*nx_-1  << " "
       << 2*ny_-1  << " "
       << dx_      << " "
       << dy_      << "\n"
       << -dx_*nx_ << " "
       <<  dx_*nx_ << " "
       << -dy_*ny_ << " "
       <<  dy_*ny_ << "\n";

  for(int j=-ny_+1;j<ny_;j++)
    for(int i=-nx_+1;i<nx_;i++)
      file << getCov(i,j) << " ";
  file << std::endl;
  file.close();
}
