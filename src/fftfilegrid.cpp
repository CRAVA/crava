#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "nrlib/iotools/logkit.hpp"
#include "lib/lib_misc.h"

#include "src/fftfilegrid.h"
#include "src/simbox.h"
#include "src/model.h"

FFTFileGrid::FFTFileGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp) :
FFTGrid(nx, ny, nz, nxp, nyp, nzp)
{
  genFileName();
  accMode_=NONE;
}

FFTFileGrid::FFTFileGrid(FFTFileGrid  * fftGrid) :
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
  rnxp_	          = 2*(cnxp_);	   

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = 0; 
  counterForSet_  = 0;
  istransformed_  = false;		
  createRealGrid();
  accMode_=NONE;

  setAccessMode(WRITE);
  fftGrid->setAccessMode(READ);
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)   
      {
        value=fftGrid->getNextReal();
        setNextReal(value);
      }// k,j,i
      endAccess();
      fftGrid->endAccess();
}



FFTFileGrid::~FFTFileGrid()
{
  endAccess();
  if(fNameIn_ != NULL)
  {
    remove(fNameIn_);
    delete [] fNameIn_;
  }
  remove(fNameOut_);
  delete [] fNameOut_;
}


void
FFTFileGrid::setAccessMode(int mode)
{
  assert(accMode_ == NONE);
  switch(mode)
  {
  case READ:
    inFile_ = fopen(fNameIn_, "rb");
    break;
  case WRITE:
    outFile_ = fopen(fNameOut_,"wb");
    break;
  case READANDWRITE:
    inFile_ = fopen(fNameIn_, "rb");
    outFile_ = fopen(fNameOut_,"wb");
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
  char * tmp;
  switch(accMode_)
  {
  case READ:
    fclose(inFile_);
    break;
  case READANDWRITE:
    fclose(inFile_); //Intentional fallthrough to WRITE
  case WRITE:
    fclose(outFile_);
    tmp = fNameIn_;
    fNameIn_ = fNameOut_;
    if(tmp != NULL)
      fNameOut_ = tmp;
    else
    {
      fNameOut_ = new char[50];
      sprintf(fNameOut_, "%sb",fNameIn_);
    }
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
FFTFileGrid::createRealGrid()
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
  fread(&cVal, sizeof(fftw_complex), 1, inFile_);
  return(cVal);
}



float 
FFTFileGrid::getNextReal()
{
  assert(istransformed_ == false);
  assert(accMode_ == READ || accMode_ == READANDWRITE);
  float rVal;
  fread(&rVal, sizeof(float), 1, inFile_);
  return float(rVal);
} 


float        
FFTFileGrid::getRealValue(int i, int j, int k)
{ 
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 
  assert(istransformed_==false);
  assert(accMode_ == RANDOMACCESS);

  int index=i+rnxp_*j+k*rnxp_*nyp_;

  assert(index<rsize_); 
  return(static_cast<float>(rvalue_[index]));
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
FFTFileGrid::setRealValue(int i, int j, int k, float value)
{ 
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 
  assert(istransformed_== false);
  assert(accMode_ == RANDOMACCESS);

  return(FFTGrid::setRealValue(i, j, k, value));
}


int 
FFTFileGrid::setNextComplex(fftw_complex value)
{
  assert(istransformed_==true);
  assert(accMode_ == READANDWRITE || accMode_ == WRITE);
  fwrite(&value, sizeof(fftw_complex), 1, outFile_);
  return(0);  
}


int   
FFTFileGrid::setNextReal(float  value)
{    
  assert(istransformed_== false);
  assert(accMode_ == READANDWRITE || accMode_ == WRITE);
  fwrite(&value, sizeof(float), 1, outFile_);
  return(0);  
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
FFTFileGrid::writeFile(const std::string & fileName, const Simbox * simbox, 
                       const std::string sgriLabel, float z0, 
                       GridMapping * depthMap, GridMapping * timeMap)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeFile(fileName, simbox, sgriLabel, z0, depthMap, timeMap);
  if(accMode_ != RANDOMACCESS)
    unload();
}

void
FFTFileGrid::writeStormFile(const std::string & fileName, const Simbox * simbox, bool expTrans, bool ascii, bool padding, bool flat)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeStormFile(fileName, simbox, expTrans, ascii, padding, flat);
  if(accMode_ != RANDOMACCESS)
    unload();
}


int
FFTFileGrid::writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  int ok = FFTGrid::writeSegyFile(fileName, simbox, z0);
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
FFTFileGrid::load()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(!istransformed_)
    FFTGrid::createRealGrid();
  else
    FFTGrid::createComplexGrid();
  if(fNameIn_ != NULL) //Something has been saved.
  {
    int i, nRead = 1;
    inFile_ = fopen(fNameIn_,"rb");
    //Real/complex does not matter in next line, since same meory is used.
    for(i=0;((i<rsize_) && (nRead == 1));i++)
      nRead = fread(&(rvalue_[i]), sizeof(float), 1, inFile_);
    fclose(inFile_);
  }
}

void 
FFTFileGrid::save()
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  outFile_ = fopen(fNameOut_,"wb");
  //Real/complex does not matter in next line, since same meory is used.
  int i, nWritten = 1;
  for(i=0;((i<rsize_) && (nWritten == 1));i++)
    nWritten = fwrite(&(rvalue_[i]), sizeof(float), 1, outFile_);
  fclose(outFile_);
  unload();
  char * tmp = fNameIn_;
  fNameIn_ = fNameOut_;
  if(tmp != NULL)
    fNameOut_ = tmp;
  else
  {
    fNameOut_ = new char[strlen(fNameIn_)+2];
    sprintf(fNameOut_, "%sb",fNameIn_);
  }
}

void
FFTFileGrid::unload()
{
  fftw_free(rvalue_); // changed
  rvalue_ = NULL;
  cvalue_ = NULL;
}

void 
FFTFileGrid::genFileName()
{
  fNameIn_  = NULL;
  std::string fileName = ModelSettings::makeFullFileName("tmpgrid" + NRLib2::ToString(gNum));
  fNameOut_ = new char[MAX_STRING];
  strcpy(fNameOut_,fileName.c_str());
  gNum++;
}

void
FFTFileGrid::writeResampledStormCube(GridMapping       * gridmapping, 
                                     const std::string & fileName, 
                                     const Simbox      * simbox,
                                     const int           format,
                                     bool                expTrans)
{
  assert(accMode_ == NONE || accMode_ == RANDOMACCESS);
  if(accMode_ != RANDOMACCESS)
    load();
  FFTGrid::writeResampledStormCube(gridmapping, fileName, simbox, format, expTrans);
  if(accMode_ != RANDOMACCESS)
    unload();
}
  

int FFTFileGrid::gNum = 0; //Starting value
