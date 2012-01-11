#include <iostream>
#include <iomanip>

#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/krigingdata2d.h"
#include "src/modelsettings.h"


KrigingData2D::KrigingData2D(int nData)
  : data_(0),   // Do not reserve space here (gives trouble with .push_back())
    indexI_(0),
    indexJ_(0),
    count_(0)
{
  //
  // Using .reserve() we set aside space for vectors, but such that
  // .push_back() used in addData() adds the 0'th element first.
  //
  data_.reserve(nData);
  indexI_.reserve(nData);
  indexJ_.reserve(nData);
  count_.reserve(nData);
}

//---------------------------------------------------------------------
KrigingData2D::~KrigingData2D(void)
{
}

//---------------------------------------------------------------------
void
KrigingData2D::addData(int   i,
                       int   j,
                       float value)
{
  if (value != RMISSING) // Do not add missing values
  {
    int ij = gotBlock(i,j);
    if (ij == -1) // Data in new location
    {
      indexI_.push_back(i);
      indexJ_.push_back(j);
      count_.push_back(1);
      data_.push_back(value);
    }
    else // Add more data to same location
    {
      count_[ij] += 1;
      data_[ij] += value;
    }
  }
}

//---------------------------------------------------------------------
int
KrigingData2D::gotBlock(int i, int j) const
{
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    if (indexI_[k]==i && indexJ_[k]==j) // Do we already have data in block (i,j)?
      return int(k);
  return -1;
}

//---------------------------------------------------------------------
void
KrigingData2D::findMeanValues(void)
{
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    data_[k] /= count_[k];
}

//---------------------------------------------------------------------
void
KrigingData2D::writeToFile(const std::string & fileName)
{
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);

  file << "    i     j     value\n";
  file << "---------------------\n";
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    file  <<std::fixed << std::setprecision(2)
          << std::setw(5)  << indexI_[k] << " "
          << std::setw(5)  << indexJ_[k] << " "
          << std::setw(10) << data_[k]   << std::endl;
  file.close();
}

