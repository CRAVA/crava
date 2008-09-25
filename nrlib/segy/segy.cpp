#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "segy.hpp"
#include "commonheaders.hpp"
#include "traceheader.hpp"

#include "../iotools/logkit.hpp"
#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/surface.hpp"
#include "../iotools/stringtools.hpp"

const float segyRMISSING = -99999.000;

using namespace NRLib2;

SegY::SegY(const std::string       & fileName,
           float                     z0, 
           const TraceHeaderFormat & traceHeaderFormat)
  : traceHeaderFormat_(traceHeaderFormat)
{
  rmissing_    = segyRMISSING;
  fileName_    = fileName;
  singleTrace_ = true;

  /// \todo Replace with safe open function.
  file_.open(fileName.c_str(), std::ios::in | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }

  // EBCDIC header
  char* junk=new char[3200];
  file_.read(junk,3200);
  delete [] junk;
  binaryHeader_ = new BinaryHeader(file_);
  nz_ = binaryHeader_->getHns();
  dz_ = static_cast<float>(binaryHeader_->getHdt()/1000);
  z0_ = z0;
  if(binaryHeader_->getFormat() == 3)
    datasize_ = 2;
  else
    datasize_ = 4;
  if(binaryHeader_->getFormat() == 4)
  {
    throw FileFormatError("Can not read SegY-file \"" + fileName + "\" that use floating point with gain representation.");
    return;
  }
  geometry_ = NULL;
}

SegY::SegY(const std::string       & fileName, 
           float                     z0, 
           int                       nz, 
           float                     dz,
           const TextualHeader     & ebcdicHeader, 
           const TraceHeaderFormat & traceHeaderFormat)
  : traceHeaderFormat_(traceHeaderFormat)
{
  rmissing_ = segyRMISSING;
  geometry_ = NULL;

  /// \todo Replace with safe open function.
  file_.open(fileName.c_str(), std::ios::out | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }
  else
  {
    z0_ = z0;
    dz_ = dz;
    nz_ = nz;
    writeMainHeader(ebcdicHeader);
  }
}

SegY::~SegY()
{
  if(geometry_!=NULL)
    delete geometry_;
  geometry_ = NULL;
  if(binaryHeader_!=NULL)
    delete binaryHeader_;
  binaryHeader_ = NULL;
  int i;
  for(i=0;i<nTraces_;i++)
    delete traces_[i];
  file_.close();
}

const SegYTrace *
SegY::getNextTrace(double zPad, Volume * volume, bool onlyVolume )
{
  bool outsideSurface = false;
  bool duplicateHeader; // Only needed for memory allocations in readAllTraces()
  if(file_.eof()==true || singleTrace_ == false)
    return NULL;
  else 
    return readTrace(volume, 
                     zPad,
                     duplicateHeader,
                     onlyVolume,
                     outsideSurface);  
}

void 
SegY::readAllTraces(Volume *volume, double zPad, bool onlyVolume)
{
  std::ios::pos_type fSize = findFileSize(fileName_);

  singleTrace_ = false;
  nTraces_     = static_cast<int>(ceil( (static_cast<double>(fSize)-3600.0)/
                                         static_cast<double>(datasize_*nz_+240.0)) );
  traces_.resize(nTraces_);

  LogKit::LogMessage(LogKit::LOW,"\nReading SEGY file " );
  LogKit::LogMessage(LogKit::LOW, fileName_);

  bool outsideSurface = false;
  bool duplicateHeader; // Needed for memory allocations.
  traces_[0] = readTrace(volume, 
                         zPad, 
                         duplicateHeader,
                         onlyVolume,
                         outsideSurface);

  double writeInterval = 0.02;
  double nextWrite = writeInterval;
  LogKit::LogMessage(LogKit::LOW,"\n  0%        20%      40%       60%       80%       100%");
  LogKit::LogMessage(LogKit::LOW,"\n  |    |    |    |    |    |    |    |    |    |    |  ");
  LogKit::LogMessage(LogKit::LOW,"\n  ^");
  int traceSize = datasize_*nz_+240;
  long long bytesRead = 3600+traceSize;
  for(unsigned int i=1 ; i < static_cast<unsigned int>(nTraces_) ; i++)
  {
    double percentDone = bytesRead/static_cast<double>(fSize);
    if(percentDone > nextWrite)
    {
      LogKit::LogMessage(LogKit::LOW,"^");
      nextWrite+=writeInterval;
    }

    traces_[i] = readTrace(volume, 
                           zPad, 
                           duplicateHeader,
                           onlyVolume,
                           outsideSurface);
    bytesRead += traceSize;
    if (duplicateHeader)
      bytesRead += 3600;

    if(file_.eof()==true)
      break;
  }
  LogKit::LogMessage(LogKit::LOW,"^\n");
 
  int count = 0;
  for(unsigned int i=1 ; i<traces_.size() ; i++)
    if (traces_[i] != NULL)
      count++;
  if (count == 0) 
  {
    std::string text;
    text += "No valid traces found. The specified time surfaces do not cover any part of the\n";
    text += " seismic data. Do you need to bypass the coordinate scaling in the SegY header?";
    throw Exception(text);
  }
}

SegYTrace *
SegY::readTrace(Volume * volume, 
                double   zPad, 
                bool   & duplicateHeader, 
                bool     onlyVolume,
                bool   & outsideSurface)
{
  TraceHeader * traceHeader = new TraceHeader(traceHeaderFormat_);

  duplicateHeader = readHeader(traceHeader);

   float x, y;
   if (traceHeaderFormat_.getCoordSys() == TraceHeaderFormat::UTM) {
     x = traceHeader->getUtmx();
     y = traceHeader->getUtmy();
   }
   else if (traceHeaderFormat_.getCoordSys() == TraceHeaderFormat::ILXL) {
     x = static_cast<float>(traceHeader->getInline());
     y = static_cast<float>(traceHeader->getCrossline());
   }
   else {
    throw Exception("Invalid coordinate system number ("
                    +ToString(traceHeaderFormat_.getCoordSys())+")");
   }
  
  int j0 = 0;
  int j1 = nz_-1;
  float zTop, zBot;
  if(volume != NULL)
  {
    if(volume->isInside(x,y) == 0 && onlyVolume == true)
    {
      readDummyTrace(file_,binaryHeader_->getFormat(),nz_);
      delete traceHeader;
      return(NULL);
    }
    try {
      zTop = static_cast<float>(volume->GetTopSurface().GetZ(x,y));
    }
    catch (NRLib2::Exception & e) {
      outsideSurface = true;
      readDummyTrace(file_,binaryHeader_->getFormat(),nz_);
      return(NULL);
    }

    if(volume->GetTopSurface().IsMissing(zTop))
    {
      readDummyTrace(file_,binaryHeader_->getFormat(),nz_);
      delete traceHeader;
      return(NULL);
    }

    try {
      zBot = static_cast<float>(volume->GetBotSurface().GetZ(x,y));
    }
    catch (NRLib2::Exception & e) {
      outsideSurface = true;
      readDummyTrace(file_,binaryHeader_->getFormat(),nz_);
      return(NULL);
    }

    if(volume->GetBotSurface().IsMissing(zBot))
    {
      readDummyTrace(file_,binaryHeader_->getFormat(),nz_);
      delete traceHeader;
      return(NULL);
    }
  }
  else
  {
    zTop = z0_;
    zBot = nz_*dz_+z0_;
  }

  if(static_cast<int>((zTop - z0_)/dz_) < 0) {
    std::string text;
    text+= "A part of the top time surface reaches above the seismic gather. (seismic\n";
    text+= " start time = "+ToString(z0_)+ " while surface z-value = "+ToString(zTop);
    text+= "). Include more seismic\n data or lower the top surface.";
    throw Exception(text); 
  } 
  //NBNB-PAL: Tried nz_ - 1 below, but that failed  (incorrectly
  //NBNB-PAL: I think) in the facies case of test suite
  if(static_cast<int>((zBot - z0_)/dz_) + 1 > nz_) {
    std::string text;
    text+= "A part of the base time surface reaches below the seismic gather. Include\n";
    text+=" more seismic data or highten the base surface.";
    throw Exception(text);
  } 

  float pad = static_cast<float>(0.5*zPad*(zBot - zTop));
  j0 = static_cast<int>((zTop - pad - z0_)/dz_);
  j1 = static_cast<int>((zBot + pad - z0_)/dz_); //Add 1 since int truncates.

  if(j0 < 0)
    j0 = 0;
  if(j1 > nz_-1)
    j1 = nz_-1;

  if(j0>j1)
    throw Exception(" Lower horizon above SegY region or upper horizon below SegY region");

  SegYTrace * trace = NULL;
  if(file_.eof()==false)
  {
    // Copy elements from j0 til j1.
    trace = new SegYTrace(file_, j0, j1, 
                          binaryHeader_->getFormat(), 
                          traceHeader->getUtmx(), 
                          traceHeader->getUtmy(), 
                          traceHeader->getInline(), 
                          traceHeader->getCrossline(), 
                          nz_ );
  }
  delete traceHeader;
  return trace;
}

bool
SegY::readHeader(TraceHeader * header)
{
  bool duplicateHeader;
  header->read(file_,binaryHeader_->getLino());
  if(header->getDt()/1000 != dz_)
    dz_ = static_cast<float>(header->getDt()/1000.0);

  switch(header->getStatus()) {
  case -2:
    throw Exception("Incorrect header encountered in SegY file"+fileName_);
    break;
  case -1:
    binaryHeader_->Update(file_);
    header->read(file_, binaryHeader_->getLino());
    duplicateHeader = true;  // Duplicate header found
    break;
  default:
    duplicateHeader = false;
    break;
  }
  return duplicateHeader;
}

void 
SegY::setGeometry(const SegyGeometry * geometry)
{
  geometry_ = new SegyGeometry(geometry);
  nTraces_  = geometry_->getNx()*geometry_->getNy();
  traces_.resize(nTraces_);
  for(int i=0 ; i < nTraces_ ; i++)
  {
    traces_[i] = NULL;
  }
}

void SegY::createRegularGrid()
{
  geometry_ = new SegyGeometry(traces_);
  nTraces_  = static_cast<int>(traces_.size());

}
std::vector<float> 
SegY::getAllValues(void)
{
  int i,nTot = 0;
  //int nTraces = nx_*ny_;
  for(i=0;i<nTraces_;i++)
    if(traces_[i] != NULL)
      nTot += traces_[i]->getEnd() - traces_[i]->getStart() + 1;

  std::vector<float> result(nTot);
  int k, kS, kE, oInd = 0;
  for(i=0;i<nTraces_;i++)
    if(traces_[i] != NULL)
    {
      kS = traces_[i]->getStart();
      kE = traces_[i]->getEnd();
      for(k = kS;k<=kE;k++)
      {
        result[oInd] = traces_[i]->getValue(k);
        oInd++;
      }
    }
    return(result);
}

float 
SegY::getValue(double x, double y, double z, int outsideMode)
{
  int i, j;
  float xind,yind;
  float value;
  double x0 = geometry_->getX0()+0.5*geometry_->getDx()*geometry_->getCosRot()-0.5*geometry_->getDy()*geometry_->getSinRot();
  double y0 = geometry_->getY0()+0.5*geometry_->getDy()*geometry_->getCosRot()+0.5*geometry_->getDx()*geometry_->getSinRot();
  double sx =  (x-x0)*geometry_->getCosRot() + (y-y0)*geometry_->getSinRot() + 0.5*geometry_->getDx();
  double sy = -(x-x0)*geometry_->getSinRot() + (y-y0)*geometry_->getCosRot() + 0.5*geometry_->getDy();
  if(geometry_!=NULL)
  {
    int ok = geometry_->returnIndex(static_cast<float>(x),static_cast<float>(y),xind,yind);

    i = static_cast<int>(xind);
    j = static_cast<int>(yind);
    int nx = geometry_->getNx();
    int ny = geometry_->getNy();

    int index;

    index = j*nx+i; //NBNB er dette rett??

    if(traces_[index]!=0 && ok==1 && z>=z0_ && z<=z0_+nz_*dz_)
    {
      int zind = static_cast<int>(floor((z-z0_)/dz_));  //NBNB   irap grid rounding different

      float v1 = traces_[index]->getValue(zind);
      if(v1 == rmissing_ && outsideMode == CLOSEST)
      {
        zind = traces_[index]->getLegalIndex(zind);
        v1 = traces_[index]->getValue(zind);
        if(traces_[index]->getValue(zind-1) == rmissing_)
          z = z0_+zind*dz_;          // Want edge value, hence 0/1 dz_ added
        else                         // (0.5 would give center of cell).
          z = z0_+(zind+0.99f)*dz_;   
      }
      if(v1 != rmissing_)
      {
        // Computes interpolated value ax^2+by^2+cz^2+dx+ey+fz+g.
        // abcdefg estimated from closest point and its closest neighbours.
        int maxInd = nx*ny-1;
        float v0, v2, a, b, c, d, e, f, g;

        // Along x:
        v0 = rmissing_;
        v2 = rmissing_;
        if(index-1 >= 0 && traces_[index-1] != NULL)
          v0 = traces_[index-1]->getValue(zind);
        if(index+1 <= maxInd && traces_[index+1] != NULL)
          v2 = traces_[index+1]->getValue(zind);
        if(v0 == rmissing_)
        {
          a = 0;
          if(v2 == rmissing_)
            d = 0;
          else
            d = v2 - v1; // Using unit coordinates in each direction.
        } 
        else if(v2 == rmissing_)
        {
          a = 0;
          d = v1 - v0;
        }
        else
        {
          a = (v2+v0-2*v1)/2.0f;
          d = (v2-v0)/2.0f;
        }
        // Along y:
        v0 = rmissing_;
        v2 = rmissing_;
        int tmpInd = index-nx;
        if(tmpInd <= maxInd && tmpInd >= 0 && traces_[tmpInd] != NULL)
          v0 = traces_[tmpInd]->getValue(zind);
        tmpInd = index +nx;
        if(tmpInd <= maxInd && tmpInd >= 0 && traces_[tmpInd] != NULL)
          v2 = traces_[tmpInd]->getValue(zind);
        if(v0 == rmissing_)
        {
          b = 0;
          if(v2 == rmissing_)
            e = 0;
          else
            e = v2 - v1; // Using unit coordinates in each direction.
        }
        else if(v2 == rmissing_)
        {
          b = 0;
          e = v1 - v0;
        }
        else
        {
          b = (v2+v0-2*v1)/2.0f;
          e = (v2-v0)/2.0f;
        }
        //Along z:
        v0 = traces_[index]->getValue(zind-1);
        v2 = traces_[index]->getValue(zind+1);
        if(v0 == rmissing_)
        {
          c = 0;
          if(v2 == rmissing_)
            f = 0;
          else
            f = v2 - v1; //Using unit coordinates in each direction.
        }
        else if(v2 == rmissing_)
        {
          c = 0;
          f = v1 - v0;
        }
        else
        {
          c = (v2+v0-2*v1)/2.0f;
          f = (v2-v0)/2.0f;
        }
        g = v1;
        double dx = geometry_->getDx();
        double dy = geometry_->getDy();
        float ux = static_cast<float>(sx/dx) - static_cast<float>(floor(sx/dx) + 0.5);
        float uy = static_cast<float>(sy/dy) - static_cast<float>(floor(sy/dy) + 0.5);
        float uz = static_cast<float>((z-z0_)/dz_) - static_cast<float>(floor((z-z0_)/dz_) + 0.5);
        value = a*ux*ux+b*uy*uy+c*uz*uz+d*ux+e*uy+f*uz+g;
      }
      else
      {
        if(outsideMode == ZERO)
          value = 0;
        else
          value = rmissing_;
      }
    }
    else
    {
      if(outsideMode == ZERO)
        value = 0;
      else
        value = rmissing_;
    }
  }
  else
    value = rmissing_;

  return(value);

}
void
SegY::writeMainHeader(const TextualHeader& ebcdicHeader)
{
  assert(file_);
  ebcdicHeader.write(file_);
  binaryHeader_ = new BinaryHeader();
  binaryHeader_->write(file_, geometry_, dz_, nz_);
}

void
SegY::storeTrace(float x, float y, std::vector<float> data, Volume *volume, float topVal,float baseVal)
{ 
  assert(file_);
  assert(geometry_>0);
  TraceHeader header(traceHeaderFormat_);
  header.setNSamples(nz_);
  header.setDt(static_cast<unsigned short>(dz_*1000));

  header.setUtmx(x);
  header.setUtmy(y);
  int IL,XL;
  int ok = geometry_->returnILXL(IL,XL,x,y);
  if(ok==1)
  {
    header.setInline(IL);
    header.setCrossline(XL);

    // header.write(file_);
    double ztop = volume->GetTopSurface().GetZ(x,y);
    std::vector<float>  trace(nz_);
    int k;
    if(volume->GetTopSurface().IsMissing(ztop))
    {
      for(k=0;k<nz_;k++)
        trace[k] = 0;
    }
    else
    {
      int firstData = static_cast<int>((ztop-z0_)/dz_);
      for(k=0;k<firstData;k++)
        trace[k] = topVal; //data[0];
      for(k=firstData;k<(firstData+static_cast<int>(data.size()));k++)
        trace[k] = data[k-firstData];
      for(k=(firstData+static_cast<int>(data.size()));k<nz_;k++)
        trace[k] = baseVal; //data[simbox_->getnz()-1];
    }

    int index;
    int j0 = 0;
    int j1 = nz_-1;
    int i,j;
    geometry_->findIJfromXY(x,y,i,j);
    index = i+geometry_->getNx()*j;
    traces_[index] = new SegYTrace(trace, j0, j1, x, y, IL, XL);
  }
  else
    throw Exception(" Coordinates are outside grid.");
}

void
SegY::writeTrace(float x, float y, std::vector<float> data, const Volume *volume, float topVal,float baseVal)
{ 
  assert(file_);
  assert(geometry_>0);
  TraceHeader header(traceHeaderFormat_);
  header.setNSamples(nz_);
  header.setDt(static_cast<unsigned short>(dz_*1000));

  header.setUtmx(x);
  header.setUtmy(y);
  int IL,XL;
  int ok = geometry_->returnILXL(IL,XL,x,y);
  if(ok==1)
  {
    header.setInline(IL);
    header.setCrossline(XL);
    header.write(file_);
    double ztop = volume->GetTopSurface().GetZ(x,y);

    std::vector<float> trace;
    trace.resize(nz_);
    int k;
    if(volume->GetTopSurface().IsMissing(ztop))
    {
      for(k=0;k<nz_;k++)
        trace[k] = 0;
    }
    else
    {
      int firstData = static_cast<int>((ztop-z0_)/dz_);
      for(k=0;k<firstData;k++)
        trace[k] = topVal; //data[0];
      for(k=firstData;k<(firstData+static_cast<int>(data.size()));k++)
        trace[k] = data[k-firstData];
      for(k=(firstData+static_cast<int>(data.size()));k<nz_;k++)
        trace[k] = baseVal; //data[simbox_->getnz()-1];
    }
    WriteBinaryIbmFloatArray(file_,trace.begin(),trace.end());
  }
  else
    throw Exception("Coordinates are outside grid.");
}

void 
SegY::writeTrace(TraceHeader * traceHeader, std::vector<float> data,Volume *volume, float topVal, float baseVal)
{
  traceHeader->setNSamples(nz_);
  traceHeader->write(file_);

  double z = z0_;
  int   nz = nz_;
  if(volume != NULL)
  {
    z = volume->GetTopSurface().GetZ(traceHeader->getUtmx(),traceHeader->getUtmy());
    nz = static_cast<int>(data.size());
  }
  std::vector<float> trace(nz_);
  int k;
  if(volume->GetTopSurface().IsMissing(z))
  {
    for(k=0;k<nz_;k++)
      trace[k] = 0;
  }
  else
  {
    int firstData = static_cast<int>((z-z0_)/dz_);
    for(k=0;k<firstData;k++)
      trace[k] = topVal; //data[0];
    for(k=firstData;k<firstData+nz;k++)
      trace[k] = data[k-firstData];
    for(k=firstData+nz;k<nz_;k++)
      trace[k] = baseVal; //data[simbox_->getnz()-1];
  }
  WriteBinaryIbmFloatArray(file_,trace.begin(),trace.end());
}

void 
SegY::WriteAllTracesToFile()
{
  int i,k;
  std::vector<float>  trace(nz_);

  float x, y;
  int IL, XL;
  int j, index;
  int il0,xl0;
  int nxl,nil;
  if(geometry_->getILxflag()==true)
  {
    nxl = geometry_->getNy();
    nil = geometry_->getNx();
  }
  else
  {
    nxl = geometry_->getNx();
    nil = geometry_->getNy();
  }
  if(geometry_->getXLstep()>0)
    xl0 = geometry_->getCrossLine0();
  else
    xl0 = geometry_->getCrossLine0()+(nxl-1)*geometry_->getXLstep();
  if(geometry_->getILstep()>0)
    il0 = geometry_->getInLine0();
  else
    il0 = geometry_->getInLine0()+(nil-1)*geometry_->getILstep();
  int ilend, xlend;
  if(geometry_->getILxflag()==true)
  {
    ilend = il0+geometry_->getILstep()*(geometry_->getNx()-1);
    xlend = xl0 + geometry_->getXLstep()*(geometry_->getNy()-1);
  }
  else
  {
    ilend = il0+geometry_->getILstep()*(geometry_->getNy()-1);
    xlend = xl0 + geometry_->getXLstep()*(geometry_->getNx()-1);
  }
  for(IL=il0;IL<ilend;IL++)
  {
    for(XL=xl0;XL<xlend;XL++)
    {
      geometry_->findIJFromILXL(IL,XL,i,j);
      index = i+j*geometry_->getNx();
      if(traces_[index]!=NULL)
      {
        for(k=0;k<nz_;k++)
          trace[k] = traces_[index]->getValue(k);

        // write traceheader
        x = traces_[index]->getX();
        y = traces_[index]->getY();
        TraceHeader header(traceHeaderFormat_);
        header.setNSamples(nz_);
        header.setDt(static_cast<unsigned short>(dz_*1000));
        header.setUtmx(x);
        header.setUtmy(y);
        header.setInline(IL);
        header.setCrossline(XL);
        header.write(file_);
        WriteBinaryIbmFloatArray(file_,trace.begin(),trace.end());
      }
    }
  }
}

std::ios::pos_type
SegY::findFileSize(const std::string& fileName)
{
  std::ios::pos_type length;
  std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }

  file.seekg(0, std::ios_base::end);
  length = file.tellg(); // ftell64bitWrapper(file);

  return length;
}

int 
SegY::findNumberOfTraces(const std::string       & fileName, 
                         const TraceHeaderFormat & traceHeaderFormat)
{
  float dummy_z0 = 0.0f;
  SegY segy(fileName, dummy_z0, traceHeaderFormat);
  int size = segy.findNumberOfTraces();
  return size;
}

int 
SegY::findNumberOfTraces(void)
{
  std::ios::pos_type fSize = findFileSize(fileName_);
  nTraces_ = static_cast<int>(ceil( (static_cast<double>(fSize)-3600.0)/
                                     static_cast<double>(datasize_*nz_+240.0)));

  TraceHeader * traceHeader = new TraceHeader(traceHeaderFormat_);
  readHeader(traceHeader);

  char * buffer = new char[nz_*datasize_];
  file_.read(buffer,(nz_*datasize_));

  int ntr = 1;
  for(int i=1 ; i<nTraces_ ; i++)
  {   
    if(file_.eof()==false)
    {
      ntr++;
      readHeader(traceHeader);
      file_.read(buffer,(nz_*datasize_));
    }
  }
  delete buffer;
  delete traceHeader;
  return(ntr);
}

SegyGeometry *
SegY::findGridGeometry(const std::string       & fileName, 
                       const TraceHeaderFormat & traceHeaderFormat)
{
  float dummy_z0 = 0.0f;
  SegY segy(fileName, dummy_z0, traceHeaderFormat);
  SegyGeometry * geometry = segy.findGridGeometry();
  return geometry;
}

SegyGeometry * 
SegY::findGridGeometry(void)
{
  TraceHeader * traceHeader = new TraceHeader(traceHeaderFormat_);

  std::ios::pos_type fSize = findFileSize(fileName_);
  nTraces_ = static_cast<int>(ceil((static_cast<double>(fSize)-3600.0)/
                                    static_cast<double>(datasize_*nz_+240.0)));
  int i;
  int tenPercent;
  if(nTraces_/10>1)
    tenPercent =nTraces_/10;
  else
    tenPercent = 1;
  char * buffer = new char[nz_*datasize_];
  int minil, maxil, minxl, maxxl;
  int IL,XL;
  int deltaIL, deltaXL;
  int dIL,dXL, maxdIL,maxdXL;
  float dxIL, dyIL, dxXL, dyXL;
  dxIL = rmissing_;
  dyIL = rmissing_;
  dxXL = rmissing_;
  dyXL = rmissing_;
  int IL0,XL0;
  float x, y, x0,y0;
  maxdIL = 0;
  maxdXL = 0;
  readHeader(traceHeader);
  IL0 = traceHeader->getInline();
  XL0 = traceHeader->getCrossline();
  x0  = traceHeader->getUtmx();
  y0  = traceHeader->getUtmy();
  minil = IL0;
  maxil = IL0;
  minxl = XL0;
  maxxl = XL0;
  file_.read(buffer,(nz_*datasize_));
  deltaIL = 1000000000;
  deltaXL = 1000000000;

  for(i=1;i<nTraces_;i++)
  {   
    if(file_.eof()==false)
    {
      readHeader(traceHeader);
      IL = traceHeader->getInline();
      XL = traceHeader->getCrossline();
      if(IL<minil )
        minil = IL;
      if(IL>maxil)
        maxil = IL;
      if(XL<minxl)
        minxl = XL;
      if(XL>maxxl)
        maxxl = XL;
      dIL = IL-IL0;
      dXL = XL-XL0;
      if(XL==XL0 && dIL>maxdIL)
      {
        x = traceHeader->getUtmx();
        y = traceHeader->getUtmy();
        dxIL = (x-x0)/dIL;
        dyIL = (y-y0)/dIL;
        maxdIL = dIL;
      }
      if(IL==IL0 && dXL>maxdXL)
      {
        x = traceHeader->getUtmx();
        y = traceHeader->getUtmy();
        dxXL = (x-x0)/dXL;
        dyXL = (y-y0)/dXL;
        maxdXL = dXL;
      }
      if(dIL>0 && dIL<deltaIL)
        deltaIL = dIL;
      if(dXL>0 && dXL<deltaXL)
        deltaXL = dXL;
      file_.read(buffer,(nz_*datasize_));
    }
  }
  float lx0 = x0-XL0*dxXL-IL0*dxIL;
  float ly0 = y0-XL0*dyXL-IL0*dyIL;

  delete traceHeader;

  // Find the four corners

  std::vector<double> cornerx(4);
  std::vector<double> cornery(4);
  std::vector<double> cornerxl(4);
  std::vector<double> corneril(4);
  cornerx[0] = lx0+(minxl-0.5*deltaXL)*dxXL+(minil-0.5*deltaIL)*dxIL;
  cornery[0] = ly0+(minxl-0.5*deltaXL)*dyXL+(minil-0.5*deltaIL)*dyIL;
  cornerx[1] = lx0+(minxl-0.5*deltaXL)*dxXL+(maxil+0.5*deltaIL)*dxIL;
  cornery[1] = ly0+(minxl-0.5*deltaXL)*dyXL+(maxil+0.5*deltaIL)*dyIL;
  cornerx[2] = lx0+(maxxl+0.5*deltaXL)*dxXL+(minil-0.5*deltaIL)*dxIL;
  cornery[2] = ly0+(maxxl+0.5*deltaXL)*dyXL+(minil-0.5*deltaIL)*dyIL;
  cornerx[3] = lx0+(maxxl+0.5*deltaXL)*dxXL+(maxil+0.5*deltaIL)*dxIL;
  cornery[3] = ly0+(maxxl+0.5*deltaXL)*dyXL+(maxil+0.5*deltaIL)*dyIL;

  cornerxl[0] = minxl-0.5*deltaXL;
  corneril[0] =minil-0.5*deltaIL;
  cornerxl[1] = minxl-0.5*deltaXL;
  corneril[1] =  maxil+0.5*deltaIL;
  cornerxl[2] = maxxl+0.5*deltaXL;
  corneril[2] = minil-0.5*deltaIL;
  cornerxl[3] = maxxl+0.5*deltaXL;
  corneril[3] = maxil+0.5*deltaIL;

  double miny, maxx, minx, maxy;
  miny = cornery[0];
  maxx = cornerx[0];
  minx = cornerx[0];
  maxy = cornery[0];
  int index = 0;
  int index2 = 0;
  int index3 = 0;
  int index4 = 0;
  for(i=1;i<4;i++)
  {
    if(cornery[i]<miny)
    {
      miny = cornery[i];
      index = i;
    }
    if(cornerx[i]>maxx)
    {
      maxx = cornerx[i];
      index2 = i;
    }
    if(cornerx[i]<minx)
    {
      minx = cornerx[i];
      index3 = i;
    }
    if(cornery[i]>maxy)
    {
      maxy = cornery[i];
      index4 = i;
    }
  }
  double rot = atan((cornery[index]-cornery[index2])/(cornerx[index]-cornerx[index2]));

  double xl,yl, dx, dy;
  double x00,y00;
  long nx, ny;
  if(fabs(rot)<=0.25*M_PI)
  {

    x00 = cornerx[index];
    y00 = cornery[index];
    xl = sqrt((cornery[index]-cornery[index2])*(cornery[index]-cornery[index2])+(cornerx[index]-cornerx[index2])*(cornerx[index]-cornerx[index2]));

    nx = long(abs(cornerxl[index]-cornerxl[index2])/deltaXL);
    if(nx==0)
      nx = long(abs(corneril[index]-corneril[index2])/deltaIL);
    dx = xl/nx;
    if(index3==index)
      index3 = index4;
    yl = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));

    ny = long(abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if(ny==0)
      ny = long(abs(corneril[index]-corneril[index3])/deltaIL);
    dy = yl/ny;
    if(index==2 || index==3)
      deltaXL = -deltaXL;
    if(index==1 || index==3)
      deltaIL = -deltaIL;
  }
  else
  {
    x00 = cornerx[index3];
    y00 = cornery[index3];
    if(index3==index)
      index3 = index2;
    rot = atan((cornery[index]-cornery[index3])/(cornerx[index]-cornerx[index3]));
    xl = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));

    nx = long(abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if(nx==0)
      nx = long(abs(corneril[index]-corneril[index3])/deltaIL);
    dx = xl/nx;
    yl = sqrt((cornery[index4]-cornery[index3])*(cornery[index4]-cornery[index3])+(cornerx[index4]-cornerx[index3])*(cornerx[index4]-cornerx[index3]));

    ny = long(abs(cornerxl[index4]-cornerxl[index3])/deltaXL);
    if(ny==0)
      ny = long(abs(corneril[index4]-corneril[index3])/deltaIL);
    dy = yl/ny;
    if(index3==2 || index3==3)
      deltaXL = -deltaXL;
    if(index3==1 || index3==3)
      deltaIL = -deltaIL;
  }
  delete buffer;

  IL0 = static_cast<int>((dyXL*(x00-lx0)-dxXL*(y00-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL+0.5);
  XL0 = static_cast<int>((dyIL*(x00-lx0)-dxIL*(y00-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL+0.5);
  x   = static_cast<float>(x00+dx*cos(rot));
  y   = static_cast<float>(y00+dx*sin(rot));

  int IL10 = static_cast<int>((dyXL*(x-lx0)-dxXL*(y-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL+0.5);
  int XL10 = static_cast<int>((dyIL*(x-lx0)-dxIL*(y-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL+0.5);
  int ilStep, xlStep;
  ilStep = IL10-IL0;
  xlStep = XL10-XL0;
  bool ILxflag = true;
  if(XL0 != XL10)
  {
    xlStep = XL10-XL0;
    ILxflag = false;
  }
  if(IL0 != IL10)
  {
    ilStep = IL10-IL0;
    ILxflag = true;
  }
  x = static_cast<float>(x00-dy*sin(rot));
  y = static_cast<float>(y00+dy*cos(rot));

  int IL01 = static_cast<int>((dyXL*(x-lx0)-dxXL*(y-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL+0.5);
  int XL01 = static_cast<int>((dyIL*(x-lx0)-dxIL*(y-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL+0.5);
  if(ILxflag == false)
    ilStep = IL01-IL0;
  else
    xlStep = XL01-XL0;

  SegyGeometry *geometry = new SegyGeometry(x00,y00,dx,dy,nx,ny,IL0,XL0,ilStep,xlStep,ILxflag,rot);

  return(geometry);
}

void SegY::readDummyTrace(std::fstream & file, int format, int nz)
{
std::vector<float> predata;
  predata.resize(nz);
  
  if(format==1)
  {
    //IBM
    ReadBinaryIbmFloatArray(file, predata.begin(), nz);
   
  }
  else if(format==2)
  {
    std::vector<int> b;
    ReadBinaryIntArray(file, b.begin(), nz);

  }
  else if(format==3)
  {
    std::vector<short> b;
    ReadBinaryShortArray(file, b.begin(), nz);
 
  }
  else if(format==5)
  {
    ReadBinaryFloatArray(file, predata.begin(), nz);

  }
  else
    throw FileFormatError("Bad format");

}

//===============================================================================
//  Class SegYTrace
//===============================================================================

SegYTrace::SegYTrace(std::fstream & file, int jStart, int jEnd, int format, float x, float y, int inLine, int crossLine, int nz)
{
  rmissing_ = segyRMISSING;
  jStart_ = jStart;
  jEnd_ = jEnd;
  x_ = x;
  y_ = y;
  inLine_ = inLine;
  crossLine_ = crossLine;
  int nData = jEnd - jStart + 1;
  int i;
  std::vector<float> predata;
  predata.resize(nz);
  data_.resize(nData);
  if(format==1)
  {
    //IBM
    ReadBinaryIbmFloatArray(file, predata.begin(), nz);
    for(i=0;i<nData;i++)
      data_[i] = predata[jStart+i];
  }
  else if(format==2)
  {
    std::vector<int> b;
    ReadBinaryIntArray(file, b.begin(), nz);
    for(i=0;i<nData;i++)
      data_[i] =  static_cast<float> (b[jStart+i]);
  }
  else if(format==3)
  {
    std::vector<short> b;
    ReadBinaryShortArray(file, b.begin(), nz);
    for(i=0;i<nData;i++)
      data_[i] =  static_cast<float> (b[jStart+i]);
  }
  else if(format==5)
  {
    ReadBinaryFloatArray(file, predata.begin(), nz);
    for(i=0;i<nData;i++)
      data_[i] = predata[jStart+i];
  }
  else
    throw FileFormatError("Bad format");
}

SegYTrace::SegYTrace(std::vector<float> indata, int jStart, int jEnd, float x, float y, int inLine, int crossLine)
{
  jStart_ = jStart;
  jEnd_ = jEnd;
  x_ = x;
  y_ = y;
  inLine_ = inLine;
  crossLine_ = crossLine;
  int nData = jEnd - jStart + 1;
  int i;
  data_.resize(nData);
  for(i=0;i<nData;i++)
    data_[i] = indata[jStart+i];
}

SegYTrace::~SegYTrace()
{
}

float
SegYTrace::getValue(int j)const
{
  float value;
  if(j < jStart_ || j > jEnd_)
    value = rmissing_;
  else
    value = data_[j-jStart_];
  return(value);
}

int
SegYTrace::getLegalIndex(int index)
{
  if(index < jStart_)
    return(jStart_);
  else if(index > jEnd_)
    return(jEnd_);
  else
    return(index);
}

SegyGeometry::SegyGeometry(std::vector<SegYTrace *> &traces)
{
  int ntraces = int(traces.size());
  int minil, maxil;
  int minxl, maxxl;
  int il,xl;
  int i;
  int ii = 0;
  while(traces[ii]==0)
    ii++;
  minil = traces[ii]->getInline();
  maxil = minil;
  minxl = traces[ii]->getCrossline();
  maxxl = minxl;

  for(i=ii+1;i<ntraces;i++)
  {
    if(traces[i]!=0)
    {
      il = traces[i]->getInline();
      xl = traces[i]->getCrossline();
      if(il<minil)
        minil = il;
      if(il>maxil)
        maxil = il;
      if(xl<minxl)
        minxl = xl;
      if(xl>maxxl)
        maxxl = xl;
    }
  }

  xl = traces[ii]->getCrossline();
  il = traces[ii]->getInline();
  int index = 0;
  int diff, maxdiff, mindiff;
  maxdiff = 0;
  mindiff = 0;
  int first = 1;
  for(i=ii+1;i<ntraces;i++)
  {
    if(traces[i]!=0)
    {
      if(traces[i]->getCrossline()==xl)
      {
        diff = abs(traces[i]->getInline()-il);
        if(diff>maxdiff)
        {
          maxdiff = diff;
          index = i;
        }
        if(first==1 || (diff<mindiff && diff>0))
        {
          mindiff = diff;
          first = 0;
        }
      }
    }
  }
  int deltaIL = mindiff;
  float dxIL = (traces[ii]->getX()-traces[index]->getX())/(il-traces[index]->getInline());
  float dyIL = (traces[ii]->getY()-traces[index]->getY())/(il-traces[index]->getInline());

  maxdiff = 0;
  first = 1;
  for(i=ii+1;i<ntraces;i++)
  {
    if(traces[i]!=0)
    {
      if(traces[i]->getInline()==il)
      {
        diff = abs(traces[i]->getCrossline()-xl);
        if(diff>maxdiff)
        {
          maxdiff = diff;
          index = i;
        }
        if(first==1 || (diff<mindiff && diff>0))
        {
          mindiff = diff;
          first = 0;
        }
      }
    }
  }

  int deltaXL = mindiff;
  float dxXL = (traces[ii]->getX()-traces[index]->getX())/(xl-traces[index]->getCrossline());
  float dyXL = (traces[ii]->getY()-traces[index]->getY())/(xl-traces[index]->getCrossline());

  float lx0 = traces[ii]->getX()-traces[ii]->getCrossline()*dxXL-traces[ii]->getInline()*dxIL;
  float ly0 = traces[ii]->getY()-traces[ii]->getCrossline()*dyXL-traces[ii]->getInline()*dyIL;

  std::vector<double> cornerx(4);
  std::vector<double> cornery(4);
  std::vector<double> cornerxl(4);
  std::vector<double> corneril(4);

  cornerx[0] = lx0+(minxl-0.5*deltaXL)*dxXL+(minil-0.5*deltaIL)*dxIL;
  cornery[0] = ly0+(minxl-0.5*deltaXL)*dyXL+(minil-0.5*deltaIL)*dyIL;
  cornerx[1] = lx0+(minxl-0.5*deltaXL)*dxXL+(maxil+0.5*deltaIL)*dxIL;
  cornery[1] = ly0+(minxl-0.5*deltaXL)*dyXL+(maxil+0.5*deltaIL)*dyIL;
  cornerx[2] = lx0+(maxxl+0.5*deltaXL)*dxXL+(minil-0.5*deltaIL)*dxIL;
  cornery[2] = ly0+(maxxl+0.5*deltaXL)*dyXL+(minil-0.5*deltaIL)*dyIL;
  cornerx[3] = lx0+(maxxl+0.5*deltaXL)*dxXL+(maxil+0.5*deltaIL)*dxIL;
  cornery[3] = ly0+(maxxl+0.5*deltaXL)*dyXL+(maxil+0.5*deltaIL)*dyIL;

  cornerxl[0] = minxl-0.5*deltaXL;
  corneril[0] =minil-0.5*deltaIL;
  cornerxl[1] = minxl-0.5*deltaXL;
  corneril[1] =  maxil+0.5*deltaIL;
  cornerxl[2] = maxxl+0.5*deltaXL;
  corneril[2] = minil-0.5*deltaIL;
  cornerxl[3] = maxxl+0.5*deltaXL;
  corneril[3] = maxil+0.5*deltaIL;
  double miny, maxx, minx, maxy;
  miny = cornery[0];
  maxx = cornerx[0];
  minx = cornerx[0];
  maxy = cornery[0];
  index = 0;
  int index2 = 0;
  int index3 = 0;
  int index4 = 0;
  for(i=1;i<4;i++)
  {
    if(cornery[i]<miny)
    {
      miny = cornery[i];
      index = i;
    }
    if(cornerx[i]>maxx)
    {
      maxx = cornerx[i];
      index2 = i;
    }
    if(cornerx[i]<minx)
    {
      minx = cornerx[i];
      index3 = i;
    }
    if(cornery[i]>maxy)
    {
      maxy = cornery[i];
      index4 = i;
    }
  }
  rot_ = atan((cornery[index]-cornery[index2])/(cornerx[index]-cornerx[index2]));

  double lx,ly;
  if(fabs(rot_)<=0.25*M_PI)
  {
    x0_ = cornerx[index];
    y0_= cornery[index];
    lx = sqrt((cornery[index]-cornery[index2])*(cornery[index]-cornery[index2])+(cornerx[index]-cornerx[index2])*(cornerx[index]-cornerx[index2]));
    nx_ = long(abs(cornerxl[index]-cornerxl[index2])/deltaXL);
    if(nx_==0)
      nx_ = long(abs(corneril[index]-corneril[index2])/deltaIL);
    dx_ = double(lx/nx_);

    if(index3==index)
      index3 = index4;
    ly = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));
    ny_ = long(abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if(ny_==0)
      ny_ = long(abs(corneril[index]-corneril[index3])/deltaIL);
    dy_ = double(ly/ny_);
    if(index==2 || index==3)
      deltaXL = -deltaXL;
    if(index==1 || index==3)
      deltaIL = -deltaIL;
  }
  else
  {

    x0_ = cornerx[index3];
    y0_= cornery[index3];
    if(index3==index)
      index = index2;
    rot_ = atan((cornery[index]-cornery[index3])/(cornerx[index]-cornerx[index3]));
    lx = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));
    nx_ = long(abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if(nx_==0)
      nx_ = long(abs(corneril[index]-corneril[index3])/deltaIL);
    dx_ = double(lx/nx_);
    if(index3==index4)
      index3 = index;
    ly = sqrt((cornery[index4]-cornery[index3])*(cornery[index4]-cornery[index3])+(cornerx[index4]-cornerx[index3])*(cornerx[index4]-cornerx[index3]));
    ny_ = long(abs(cornerxl[index4]-cornerxl[index3])/deltaXL);
    if(ny_==0)
      ny_ = long(abs(corneril[index4]-corneril[index3])/deltaIL);
    dy_ = double(ly/ny_);
    if(index3==2 || index3==3)
      deltaXL = -deltaXL;
    if(index3==1 || index3==3)
      deltaIL = -deltaIL;
  }
  sinRot_ = sin(rot_);
  cosRot_ = cos(rot_);

  int XL0 = static_cast<int>((dyIL*(x0_-lx0)-dxIL*(y0_-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL + 0.5);
  int IL0 = static_cast<int>((dyXL*(x0_-lx0)-dxXL*(y0_-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL + 0.5);
  inLine0_ = IL0;
  crossLine0_ = XL0;
  float x,y;
  findXYfromIJ(1,0,x,y);
  int IL10 = static_cast<int>((dyXL*(x-lx0)-dxXL*(y-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL+0.5);
  int XL10 = static_cast<int>((dyIL*(x-lx0)-dxIL*(y-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL+0.5);
  if(XL0 != XL10)
  {
    xlStep_ = XL10-XL0;
    ILxflag_ = false;
  }
  if(IL0 != IL10)
  {
    ilStep_ = IL10-IL0;
    ILxflag_ = true;
  }
  findXYfromIJ(0,1,x,y);
  int IL01 = static_cast<int>((dyXL*(x-lx0)-dxXL*(y-ly0))/(dxIL*dyXL-dyIL*dxXL)+0.5*deltaIL+0.5);
  int XL01 = static_cast<int>((dyIL*(x-lx0)-dxIL*(y-ly0))/(dxXL*dyIL-dyXL*dxIL)+0.5*deltaXL+0.5);
  if(ILxflag_ == false)
    ilStep_ = IL01-IL0;
  else
    xlStep_ = XL01-XL0;

  int nTraces = nx_*ny_;
  std::vector<SegYTrace*> tracestmp;
  tracestmp.resize(nTraces);

  int j, k;
  for(k=0;k<nTraces;k++)
    tracestmp[k] = 0;

  for(k=0;k<ntraces;k++)
  {
    if(traces[k]!=0)
    {
      x = traces[k]->getX();
      y = traces[k]->getY();
      findIJfromXY(x,y,i,j);
      tracestmp[i+nx_*j] = traces[k];
    }
  }

  traces.resize(nTraces);
  for(k=0;k<nTraces;k++)
    traces[k] = tracestmp[k];
}

SegyGeometry::SegyGeometry(double x0, double y0, double dx, double dy, int nx, int ny, 
                           int IL0, int XL0, int ilStep, int xlStep, bool ILxflag, double rot)
{
  x0_         = x0;
  y0_         = y0;
  dx_         = dx;
  dy_         = dy;
  nx_         = nx;
  ny_         = ny;
  inLine0_    = IL0;
  crossLine0_ = XL0;
  ilStep_     = ilStep;
  xlStep_     = xlStep;
  ILxflag_    = ILxflag;
  cosRot_     = cos(rot);
  sinRot_     = sin(rot);
  rot_        = rot;
}

SegyGeometry::SegyGeometry(const SegyGeometry *geometry)
{
  x0_         = geometry->x0_; 
  y0_         = geometry->y0_;
  dx_         = geometry->dx_;
  dy_         = geometry->dy_;
  nx_         = geometry->nx_;
  ny_         = geometry->ny_;
  inLine0_    = geometry->inLine0_;
  crossLine0_ = geometry->crossLine0_;
  ilStep_     = geometry->ilStep_; 
  xlStep_     = geometry->xlStep_;
  ILxflag_    = geometry->ILxflag_;
  sinRot_     = geometry->sinRot_;
  cosRot_     = geometry->cosRot_;
  rot_        = geometry->rot_;
}

SegyGeometry::~SegyGeometry()
{
}

void 
SegyGeometry::findXYfromIJ(int i, int j, float &x, float &y)const
{
  x = static_cast<float>(x0_+i*dx_*cosRot_-j*dy_*sinRot_);
  y = static_cast<float>(y0_+i*dx_*sinRot_+j*dy_*cosRot_);
}

void 
SegyGeometry::findIJfromXY(float x, float y, int &i, int &j)
{
  double teller1 = (x-x0_)*cosRot_+(y-y0_)*sinRot_;
  i = static_cast<int>(teller1/dx_);
  double teller2 = -(x-x0_)*sinRot_+(y-y0_)*cosRot_;
  j = static_cast<int>(teller2/dy_);
}

int 
SegyGeometry::returnIndex(float x, float y, float &xind, float  &yind)
{
  xind = static_cast<float>(( (x-x0_)*cosRot_+(y-y0_)*sinRot_)/dx_);
  yind = static_cast<float>((-(x-x0_)*sinRot_+(y-y0_)*cosRot_)/dy_);
  if(xind>=0 && xind<nx_ && yind>=0 && yind<ny_)
    return 1;
  else
    return 0;
}

int 
SegyGeometry::returnILXL(int &IL, int &XL, float x, float y)
{
  int i, j;
  findIJfromXY(x,y,i,j);
  if(ILxflag_==false)
  {
    XL  = crossLine0_+i*xlStep_;
    IL = inLine0_+j*ilStep_;
  }
  else
  {
    XL = crossLine0_+j*xlStep_;
    IL = inLine0_+i*ilStep_;
  }
  if(i>=0 && i<nx_ && j>=0 && j<ny_)
  {
    return 1;
  }
  else
    return 0;
}

void 
SegyGeometry::findIJFromILXL(int IL, int XL, int &i, int &j)
{
  if(ILxflag_==false)
  {
    i = (XL-crossLine0_)/xlStep_;
    j = (IL-inLine0_)/ilStep_;
  }
  else
  {
    i = (IL-inLine0_)/ilStep_;
    j = (XL-crossLine0_)/xlStep_;
  }
}

void SegyGeometry::writeGeometry() const
{
  double geoangle = -rot_*180/(M_PI);
  if (geoangle < 0)
    geoangle += 360.0;

  LogKit::LogFormatted(LogKit::HIGH,"\n                        x0           y0            lx        ly         dx      dy     azimuth\n");
  LogKit::LogFormatted(LogKit::HIGH,"----------------------------------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::HIGH,"Seismic area:  %11.2f  %11.2f    %7.2f %7.2f    %8.3f\n", x0_, y0_, nx_*dx_, ny_*dy_, dx_, dy_, geoangle);
}
