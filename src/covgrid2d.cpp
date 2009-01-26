#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdarg.h>
#include <ostream>
#include "src/covgrid2d.h"

CovGrid2D::CovGrid2D(Vario *vario, int nx, int ny, float dx, float dy)
{
  nx_ = nx;
  ny_ = ny;
  dx_ = dx;
  dy_ = dy;
  cov_.resize(2*nx_*ny_);
  int index;
  int i,j;
  for(i=0;i<nx_;i++)
    for(j=-ny_;j<ny_;j++)
    {
      //index = i+(j+ny_)*nx_;
      index = i*2*ny_+ny_+j;
      cov_[index] = vario->corr(i*dx, j*dy);
    }


}

float
CovGrid2D::getCov(int deltai, int deltaj)
{
  if(deltai<0)
  {
    deltai = -deltai;
    deltaj = -deltaj;
  }
  int index = deltai*2*ny_+deltaj+ny_;
  return cov_[index];

}

void
CovGrid2D::writeToFile()
{
  std::string filename = "correlation.storm";
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
  //if (!file) {
  //  throw new IOError("Error opening " + filename + " for writing.");
 // }

  file.precision(14);

  file << 2*nx_-1 << " " << 2*ny_-1 << " " 
       << dx_ << " " << dy_ << "\n"
       << -dx_*nx_ << " " << dx_*nx_ << " "
       << -dy_*ny_ << " " << dy_*ny_ << "\n";

  int i, j;
  
    for(j=-ny_+1;j<ny_;j++)
for(i=-nx_+1;i<nx_;i++)
      file << getCov(i,j) << " ";
  file.close();


}
