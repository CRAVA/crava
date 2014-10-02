/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/fftfilegrid.h"
#include "src/simbox.h"
#include "src/io.h"

FFTFileGrid::FFTFileGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp) :
FFTGrid(nx, ny, nz, nxp, nyp, nzp)
{
  genFileName();
  accMode_=NONE;
}

FFTFileGrid::FFTFileGrid(FFTFileGrid  * fftGrid, bool expTrans) :
FFTGrid()
{
  float value;
  int   i,j,k;
  genFileName();

  cubetype_       = fftGrid->cubetype_;
  theta_          = fftGrid->theta_;
  nx_             = fftGrid->nx_;
  ny_             = fftGrid->ny_;
  nz_             = fftGrid->nz_;
  nxp_            = fftGrid->nxp_;
  nyp_            = fftGrid->nyp_;
  nzp_            = fftGrid->nzp_;

  cnxp_           = nxp_/2+1;
  rnxp_           = 2*(cnxp_);

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = 0;
  counterForSet_  = 0;
  istransformed_  = fftGrid->istransformed_;
  fNameIn_        = "";
  accMode_        = NONE;

  setAccessMode(WRITE);
  fftGrid->setAccessMode(READ);

  if(istransformed_ == false){
    createRealGrid();
    for(k=0;k<nzp_;k++) {
      for(j=0;j<nyp_;j++) {
        for(i=0;i<rnxp_;i++) {
          value=fftGrid->getNextReal();
          if (expTrans)
            setNextReal(exp(value));
          else
            setNextReal(value);
        }
      }
    }
  }
  else{
    createComplexGrid();
    for(int k=0;k<nzp_;k++) {
      for(int j=0;j<nyp_;j++) {
        for(int i=0;i<cnxp_;i++) {
          fftw_complex value = fftGrid->getNextComplex();
          setNextComplex(value);
        }
      }
    }
  }
  endAccess();
  fftGrid->endAccess();
}



FFTFileGrid::~FFTFileGrid()
{
  endAccess();
  if(fNameIn_ != "")
  {
    remove(fNameIn_.c_str());
  }
  remove(fNameOut_.c_str());
}


void
FFTFileGrid::setAccessMode(int mode)
{
  assert(accMode_ == NONE);
  switch(mode)
  {
  case READ:
    NRLib::OpenRead(inFile_,fNameIn_,std::ios::in | std::ios::binary);
    break;
  case WRITE:
    NRLib::OpenWrite(outFile_,fNameOut_,std::ios::out | std::ios::binary);
    break;
  case READANDWRITE:
    NRLib::OpenRead(inFile_,fNameIn_,std::ios::in | std::ios::binary);
    NRLib::OpenWrite(outFile_,fNameOut_,std::ios::out | std::ios::binary);
    break;
  case RANDOMACCESS:
    modified_ = 0;
    load();
    break;
  }
  accMode_ = mode;
}

void
FFTFileGrid::endAccess()
{
  std::string tmp("");
  switch(accMode_)
  {
  case READ:
    inFile_.close();
    break;
  case READANDWRITE:
    inFile_.close(); //Intentional fallthrough to WRITE
  case WRITE:
    outFile_.close();
    tmp = fNameIn_;
    fNameIn_ = fNameOut_;
    if(tmp != "")
      fNameOut_ = tmp;
    else
      fNameOut_ = fNameIn_+"b";
    break;
  case RANDOMACCESS:
    if(modified_ != 0)
      save();
    else
      unload();
    break;
  }
  accMode_ = NONE;
}

void
FFTFileGrid::createRealGrid(bool /*add*/)
{
  istransformed_ = false;
}

void
FFTFileGrid::createComplexGrid()
{
  istransformed_ = true;
}


fftw_complex
FFTFileGrid::getNextComplex()
{
  assert(istransformed_==true);
  assert(accMode_ == READ || accMode_ == READANDWRITE);
  fftw_complex cVal;
  char * buffer = reinterpret_cast<char *>(&cVal);
  inFile_.read(buffer,sizeof(fftw_complex));
  return(cVal);
}



float
FFTFileGrid::getNextReal()
{
  assert(istransformed_ == false);
  assert(accMode_ == READ || accMode_ == READANDWRITE);
  float rVal;
  char * buffer = reinterpret_cast<char *>(&rVal);
  inFile_.read(buffer,sizeof(float));
  return float(rVal);
}


float
FFTFileGrid::getRealValue(int i, int j, int k, bool extSimbox)
{
  // i index in x direction
  // j index in y direction
  // k index in z direction
  assert(istransformed_==false);
  assert(accMode_ == RANDOMACCESS);

 bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
    ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));

  float value;
 if( inSimbox && notMissing )
  { // if index in simbox
  int index=i+rnxp_*j+k*rnxp_*nyp_;
  value = static_cast<float>(rvalue_[index]);
  }
  else
  {
    value = RMISSING;
  }

  return( value );
 // assert(index<rsize_);
  //return(static_cast<float>(rvalue_[index]));
}

float
FFTFileGrid::getRealValueInterpolated(int i, int j, float kindex)
{

  // i index in x direction
  // j index in y direction
  // k index in z direction, float, should interpolate
  assert(istransformed_==false);
  assert(accMode_ == RANDOMACCESS);

  return(FFTGrid::getRealValueInterpolated(i, j, kindex));
}

int
FFTFileGrid::setRealValue(int i, int j, int k, float value, bool extSimbox)
{
  // i index in x direction
  // j index in y direction
  // k index in z direction
  assert(istransformed_== false);
  assert(accMode_ == RANDOMACCESS);
  modified_ = 1;
  return(FFTGrid::setRealValue(i, j, k, value, extSimbox));
}


int
FFTFileGrid::SetNextComplex(std::complex<double> & value)
{
  assert(istransformed_==true);
  assert(accMode_ == READANDWRITE || accMode_ == WRITE);
  fftw_complex tmp;
  tmp.re = static_cast<fftw_real>(value.real());
  tmp.im = static_cast<fftw_real>(value.imag());
  char * buffer = reinterpret_cast<char *>(&tmp);
  outFile_.write(buffer,sizeof(fftw_complex));
  return(0);
}


int
FFTFileGrid::setNextComplex(fftw_complex value)
{
  assert(istransformed_==true);
  assert(accMode_ == READANDWRITE || accMode_ == WRITE);
  char * buffer = reinterpret_cast<char *>(&value);
  outFile_.write(buffer,sizeof(fftw_complex));
  return(0);
}


int
FFTFileGrid::setNextReal(float  value)
{
  assert(istransformed_== false);
  assert(accMode_ == READANDWRITE || accMode_ == WRITE);
  char * buffer = reinterpret_cast<char *>(&value);
  outFile_.write(buffer,sizeof(float));
  return(0);
}

float
FFTFileGrid::getFirstRealValue(){
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);

  float value;

  if(accMode_ != RANDOMACCESS) {
    setAccessMode(READ);
    value =  getNextReal();
    endAccess();
  }
  else
    value = static_cast<float>(rvalue_[0]);

  return value;
}

int
FFTFileGrid::square()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::square();
  if(accMode_ != RANDOMACCESS)
    save();
  return(0);
}

int
FFTFileGrid::expTransf()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::expTransf();
  if(accMode_ != RANDOMACCESS)
    save();
  return(0);
}

int
FFTFileGrid::logTransf()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::logTransf();
  if(accMode_ != RANDOMACCESS)
    save();
  return(0);
}

int
FFTFileGrid::collapseAndAdd(float * grid)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::collapseAndAdd(grid);
  if(accMode_ != RANDOMACCESS)
    save();
  return(0);
}

void
FFTFileGrid::fftInPlace()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::fftInPlace();
  if(accMode_ != RANDOMACCESS)
    save();
}


void
FFTFileGrid::invFFTInPlace()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::invFFTInPlace();
  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::multiplyByScalar(float scalar)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::multiplyByScalar(scalar);
  if(accMode_ != RANDOMACCESS)
    save();
}


void
FFTFileGrid::add(FFTGrid * fftGrid)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  assert(nxp_==fftGrid->getNxp());
  fftGrid->setAccessMode(READ);

  if(istransformed_==true)
  {
    int i;
    fftw_complex value;
    for(i=0;i<csize_;i++)
    {
      value = fftGrid->getNextComplex();
      cvalue_[i].re += value.re;
      cvalue_[i].im += value.im;
    }
  }
  else
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] += fftGrid->getNextReal();
    }
  }
  fftGrid->endAccess();

  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::addScalar(float scalar)
{
 assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::addScalar(scalar);
  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::subtract(FFTGrid * fftGrid)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  assert(nxp_==fftGrid->getNxp());
  fftGrid->setAccessMode(READ);

  if(istransformed_==true)
  {
    int i;
    fftw_complex value;
    for(i=0;i<csize_;i++)
    {
      value = fftGrid->getNextComplex();
      cvalue_[i].re -= value.re;
      cvalue_[i].im -= value.im;
    }
  }
  else
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] -= fftGrid->getNextReal();
    }
  }
  fftGrid->endAccess();

  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::changeSign()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;

  if(istransformed_==true)
  {
    int i;
    for(i=0;i<csize_;i++)
    {

      cvalue_[i].re = -cvalue_[i].re;
      cvalue_[i].im = -cvalue_[i].im;
    }
  }
  else
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] = -rvalue_[i];
    }
  }


  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::multiply(FFTGrid * fftGrid)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  assert(nxp_==fftGrid->getNxp());
  fftGrid->setAccessMode(READ);

  if(istransformed_==true)
  {
    int i;
    fftw_complex value;
    for(i=0;i<csize_;i++)
    {
      value = fftGrid->getNextComplex();
      cvalue_[i].re *= value.re;
      cvalue_[i].im *= value.im;
    }
  }
  else
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] *= fftGrid->getNextReal();
    }
  }
  fftGrid->endAccess();

  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::conjugate()
{
  assert(istransformed_);
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;

  int i;
  for(i=0;i<csize_;i++)
  {
    cvalue_[i].im = -cvalue_[i].im;
  }

  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::fillInComplexNoise(RandomGen * ranGen)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::fillInComplexNoise(ranGen);
  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::writeFile(const std::string & fileName,
                       const std::string & subDir,
                       const Simbox      * simbox,
                       const std::string   sgriLabel,
                       const float         z0,
                       const GridMapping * depthMap,
                       //const GridMapping * timeMap,
                       const TraceHeaderFormat & thf,
                       bool                             padding,
                       bool                             scientific_format,
                       const std::vector<std::string> & headerText)

{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeFile(fileName, subDir, simbox, sgriLabel, z0, depthMap, thf,padding, scientific_format, headerText);
  if(accMode_ != RANDOMACCESS)
    unload();
}

void
FFTFileGrid::writeStormFile(const std::string & fileName, const Simbox * simbox, bool ascii,
                            bool padding, bool flat, bool scientific_format)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeStormFile(fileName, simbox, ascii, padding, flat, scientific_format);
  if(accMode_ != RANDOMACCESS)
    unload();
}


int
FFTFileGrid::writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0, const TraceHeaderFormat &thf,
                           const std::vector<std::string> & headerText)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  int ok = FFTGrid::writeSegyFile(fileName, simbox, z0, thf, headerText);
  if(accMode_ != RANDOMACCESS)
    unload();
  return(ok);
}

int
FFTFileGrid::writeSgriFile(const std::string & fileName, const Simbox * simbox, const std::string label)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  int ok = FFTGrid::writeSgriFile(fileName, simbox, label);
  if(accMode_ != RANDOMACCESS)
    unload();
  return(ok);
}

void
FFTFileGrid::writeCravaFile(const std::string & fileName, const Simbox * simbox)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeCravaFile(fileName,simbox);
  if(accMode_ != RANDOMACCESS)
    unload();
}


void
FFTFileGrid::readCravaFile(const std::string & fileName, std::string & error, bool nopadding)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::readCravaFile(fileName, error, nopadding);
  modified_ = 1;
  if(accMode_ != RANDOMACCESS)
    save();
}

void
FFTFileGrid::load()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(!istransformed_)
    FFTGrid::createRealGrid();
  else
    FFTGrid::createComplexGrid();
  if(fNameIn_ != "") //Something has been saved.
  {
    int i, nRead = 1;
    NRLib::OpenRead(inFile_,fNameIn_,std::ios::in | std::ios::binary);
    //Real/complex does not matter in next line, since same meory is used.
    for(i=0;((i<rsize_) && (nRead == 1));i++)
    {
      char * buffer = reinterpret_cast<char *>(&(rvalue_[i]));
      inFile_.read(buffer,sizeof(float));
    }
    inFile_.close();
  }
}

void
FFTFileGrid::save()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  NRLib::OpenWrite(outFile_,fNameOut_,std::ios::out | std::ios::binary);
  //Real/complex does not matter in next line, since same meory is used.
  int i, nWritten = 1;
  for(i=0;((i<rsize_) && (nWritten == 1));i++)
  {
    char * buffer = reinterpret_cast<char *>(&(rvalue_[i]));
    outFile_.write(buffer,sizeof(float));
  }
  outFile_.close();
  unload();
  std::string tmp = fNameIn_;
  fNameIn_ = fNameOut_;
  if(tmp != "")
    fNameOut_ = tmp;
  else
    fNameOut_ = fNameIn_+"b";
}

void
FFTFileGrid::unload()
{
  fftw_free(rvalue_); // changed
  nGrids_ = nGrids_ - 1;
// LogKit::LogFormatted(LogKit::Error,"\nFFTFileGrid unload: nGrids_ = %d\n",nGrids_);
  rvalue_ = NULL;
  cvalue_ = NULL;
}

void
FFTFileGrid::genFileName()
{
  fNameIn_  = "";
  // Store tmp file in top directory.
  std::string baseName = IO::PrefixTmpGrids() + NRLib::ToString(gNum);
  std::string fileName = IO::makeFullFileName(IO::PathToTmpFiles(), baseName);
  fNameOut_ = fileName;
  gNum++;
}

void
FFTFileGrid::writeResampledStormCube(const GridMapping * gridmapping,
                                     const std::string & fileName,
                                     const Simbox      * simbox,
                                     const int           format)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeResampledStormCube(gridmapping, fileName, simbox, format);
  if(accMode_ != RANDOMACCESS)
    unload();
}

void
FFTFileGrid::getRealTrace(float * value, int i, int j)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  FFTGrid::getRealTrace(value, i,j);
  if(accMode_ != RANDOMACCESS)
    save();
}

int
FFTFileGrid::setRealTrace(int i, int j, float *value)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  else
    modified_ = 1;
  int notok = FFTGrid::setRealTrace(i,j, value);
  if(accMode_ != RANDOMACCESS)
    save();
  return(notok);
}


int FFTFileGrid::gNum = 0; //Starting value
