// $Id: segy.cpp 1199 2013-10-02 08:24:02Z anner $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>

#include "segy.hpp"
#include "commonheaders.hpp"
#include "traceheader.hpp"

#include "../iotools/logkit.hpp"
#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/surface.hpp"
#include "../iotools/stringtools.hpp"


const float segyRMISSING = -99999.000;

using namespace NRLib;


SegY::SegY(const std::string       & fileName,
           float                     z0,
           const TraceHeaderFormat & traceHeaderFormat)
  : trace_header_format_(traceHeaderFormat)
{
  rmissing_    = segyRMISSING;
  file_name_    = fileName;
  single_trace_ = true;

  /// \todo Replace with safe open function.
 // file_.open(fileName.c_str(), std::ios::in | std::ios::binary);

  OpenRead(file_, fileName, std::ios::in | std::ios::binary);

  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }

  // EBCDIC header
  char* junk=new char[3200];
  file_.read(junk,3200);
  delete [] junk;
  binary_header_ = new BinaryHeader(file_);
  nz_ = binary_header_->GetHns();
  dz_ = static_cast<float>(binary_header_->GetHdt()/1000);
  z0_ = z0;
  if (binary_header_->GetFormat() == 3)
    datasize_ = 2;
  else
    datasize_ = 4;
  if (binary_header_->GetFormat() == 4)
  {
    delete binary_header_;
    throw FileFormatError("Can not read SegY-file \"" + fileName + "\" that use floating point with gain representation.");
    return;
  }
  geometry_ = NULL;

  unsigned long long fSize = FindFileSize(file_name_);
  n_traces_     = static_cast<int>(ceil( (static_cast<double>(fSize)-3600.0)/
                                         static_cast<double>(datasize_*nz_+240.0)) );
}

//-------------------------------------------------------------------------
SegY::SegY(const std::string               & fileName,
           float                             z0,
           std::vector<TraceHeaderFormat*>   thf,
           bool                              searchStandardFormats)
{
  rmissing_     = segyRMISSING;
  file_name_    = fileName;
  single_trace_ = true;

  /// \todo Replace with safe open function.
 // file_.open(fileName.c_str(), std::ios::in | std::ios::binary);
  OpenRead(file_, fileName, std::ios::in | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }

  // EBCDIC header
  char* junk=new char[3200];
  file_.read(junk,3200);
  delete [] junk;
  binary_header_ = new BinaryHeader(file_);
  nz_ = binary_header_->GetHns();
  dz_ = static_cast<float>(binary_header_->GetHdt()/1000);
  z0_ = z0;
  if (binary_header_->GetFormat() == 3)
    datasize_ = 2;
  else
    datasize_ = 4;
  if (binary_header_->GetFormat() == 4)
  {
    delete binary_header_;
    throw FileFormatError("Can not read SegY-file \"" + fileName + "\" that use floating point with gain representation.");
    return;
  }
  geometry_ = NULL;

  //Find which trace header to use

  std::vector<TraceHeaderFormat*> stdList(TraceHeaderFormat::GetListOfStandardHeaders());

  if (searchStandardFormats)
  {
    for (size_t i = 0; i < stdList.size(); i++)
      thf.push_back(stdList[i]); // thf.insert() gives an UMR in Purify
  }

  // Check that all formats are only defined once.
  size_t nthf = thf.size();
  std::vector<bool> remove(nthf);
  size_t i;
  for (i = 0; i < nthf; i++)
    remove[i] = 0;
  for (i = 0; i < nthf; i++)
    for (size_t j = i + 1; j < nthf ; j++)
      if (thf[i]->IsDifferent(*thf[j])==0)
        remove[j] = 1;
  for (size_t j = nthf - 1; j >= 1; j--)
    if (remove[j] == 1)
      thf.erase(thf.begin()+j);
  if(remove[0] == 1)
    thf.erase(thf.begin());

  size_t ntraceheaders = thf.size();
  bool ok = 0;
  i = 0;

  while(ok==false && i<ntraceheaders)
  {
    ok = TraceHeaderOK(file_, thf[i]);
    i++;
    if(!ok)
    {
      file_.close();
      file_.clear();
      OpenRead(file_, fileName, std::ios::in | std::ios::binary);
      if (!file_) {
        throw new IOError("Error opening " + fileName);
      }
    junk = new char[3600];
    file_.read(junk,3600);
    delete [] junk;
    }


  }
  if(ok==false)
  {
  throw FileFormatError("Can not read SegY-file \"" + fileName + "\". Unknown traceheader format.");
    return;
  }
  else
  {
    TraceHeaderFormat traceHeaderFormat = *thf[i-1];
    trace_header_format_ = traceHeaderFormat;
    LogKit::LogMessage(LogKit::Medium,"\nSeismic data of type \'" + thf[i-1]->GetFormatName() + "\' detected\n");
    thf[i-1]->WriteValues();
  }
  //
  // Delete the trace header formats given in the standard list. The other
  // THFs must be deleted in the program calling this constructor.
  //
  for (unsigned int i = 0; i < stdList.size(); i++)
    delete stdList[i];

 // for (i = 0; i < ntraceheaders; i++)
 // {
//    delete traceHeader[i];
//    delete traceHeader2[i];
//  }

  // file_.seekg(0, std::ios::beg);
  file_.close();
  file_.clear();

 // file_.open(fileName.c_str(), std::ios::in | std::ios::binary);
  OpenRead(file_, fileName, std::ios::in | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }

  junk = new char[3600];
  file_.read(junk,3600);
  delete [] junk;

  unsigned long long fSize = FindFileSize(file_name_);
  n_traces_     = static_cast<int>(ceil( (static_cast<double>(fSize)-3600.0)/
                                          static_cast<double>(datasize_*nz_+240.0)) );
}


bool SegY::CompareTraces(TraceHeader *header1, TraceHeader *header2, int &delta, int &deltail, int &deltaxl)
{
  double deltax  = std::abs(header1->GetUtmx() - header2->GetUtmx());
  double deltay  = std::abs(header1->GetUtmy() - header2->GetUtmy());
  deltail = abs(header1->GetInline() - header2->GetInline());
  deltaxl = abs(header1->GetCrossline() - header2->GetCrossline());
  double deltaxSq = deltax*deltax;
  double deltaySq = deltay*deltay;
  double utmDist = sqrt(deltaxSq+deltaySq);
  int   deltaxlSq = deltaxl*deltaxl;
  int   deltailSq = deltail*deltail;

  double k;
  if (deltail == 0 && deltaxl == 0)
    return false;
  if (deltax == 0.0 && deltay == 0)
    return false;
  bool ok = false;
  delta = 0;
  if (deltail==0)
  {
    k = utmDist/deltaxl;
    double help = fmod(k,double(6.25));
      if (k/6.25<8.5 && k/6.25>0.5 && (help<=1.0 || (6.25-help)<=1.0))
      {
        ok = true;
        delta = 1;
      }

  }
  else if (deltaxl==0)
  {
    k = utmDist/deltail;
    double help = fmod(k,double(6.25));
      if (k/6.25<8.5 && k/6.25>0.5 && (help<=1.0 || (6.25-help)<=1.0))
      {
        ok = true;
        delta = 2;
      }
  }
  else
  {
    // solution of a^2dXL^2+b^2dIL^2=dx^2+dy^2, a and b multiples of 6.25.
    int i,j;
    double prodSum = utmDist/6.25;
    double prodSumSqRelXL = prodSum*prodSum/deltaxlSq;
    double deltailxlRelSq = deltailSq/deltaxlSq;
    i = static_cast<int>(floor(0.5+prodSum/static_cast<float>(deltail)));
    while(i > 0 && !ok) {
      int j0 = static_cast<int>(floor(sqrt(prodSumSqRelXL-i*i*deltailxlRelSq)));
      for(j=j0;j<=j0+1;j++) {
        k = utmDist/sqrt(static_cast<float>(i*i*deltailSq+j*j*deltaxlSq));
        double help = fmod(k,6.25);
        if (k/6.25<8.5 && k/6.25>0.5 && (help<=1.0 || (6.25-help)<=1.0))
          ok = true;
      }
      i--;
    }
  }

  return ok;

}

SegY::SegY(const std::string       & fileName,
           float                     z0,
           size_t                    nz,
           float                     dz,
           const TextualHeader     & ebcdicHeader,
           const TraceHeaderFormat & traceHeaderFormat)
  : trace_header_format_(traceHeaderFormat)
{
  rmissing_ = segyRMISSING;
  geometry_ = NULL;
  binary_header_ = NULL;

  /// \todo Replace with safe open function.
  //file_.open(fileName.c_str(), std::ios::out | std::ios::binary);
  OpenWrite(file_,fileName, std::ios::out | std::ios::binary);
  if (!file_) {
    throw new IOError("Error opening " + fileName);
  }
  else
  {
    z0_ = z0;
    dz_ = dz;
    nz_ = nz;
    WriteMainHeader(ebcdicHeader);
  }
}

SegY::SegY(const StormContGrid     * storm_grid,
           const SegyGeometry      * geometry,
           float                     z0,
           int                       nz,
           const std::string       & file_name,
           bool                      write_to_file,
           const TraceHeaderFormat & trace_header_format,
           bool                      is_seismic)
{
  rmissing_      = segyRMISSING;
  geometry_      = NULL;
  binary_header_ = NULL;

  int i,k,j;
  TextualHeader header = TextualHeader::standardHeader();
  int nx = static_cast<int>(storm_grid->GetNI());
  int ny = static_cast<int>(storm_grid->GetNJ());
  float dz = float(floor((storm_grid->GetLZ()/storm_grid->GetNK())));
  //float z0 = 0.0;
  //int nz   = int(ceil((storm_grid->GetZMax())/dz));

  //SegY segyout(file_name,0,nz,dz,header);
  trace_header_format_ = trace_header_format;
  z0_ = z0;
  nz_ = nz;
  dz_ = dz;

  if (write_to_file) {
    OpenWrite(file_, file_name, std::ios::out | std::ios::binary);
    WriteMainHeader(header);
  }

  if(geometry == NULL) {
    SegyGeometry dummy_geometry(storm_grid->GetXMin(), storm_grid->GetYMin(), storm_grid->GetDX(), storm_grid->GetDY(),
                                nx, ny, storm_grid->GetAngle());
    SetGeometry(&dummy_geometry);
  }
  else
    SetGeometry(geometry);


  std::vector<float> data_vec;
  data_vec.resize(nz);
  double z_shift = 0;
  if(is_seismic == false)
    z_shift = 0.5*dz_;
  float x, y, xt, yt, z;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      xt = float((i+0.5)*geometry_->GetDx());
      yt = float((j+0.5)*geometry_->GetDy());
      x  = float(geometry_->GetX0()+xt*geometry_->GetCosRot()-yt*geometry_->GetSinRot());
      y  = float(geometry_->GetY0()+yt*geometry_->GetCosRot()+xt*geometry_->GetSinRot());

      double z_bot = storm_grid->GetBotSurface().GetZ(x,y);
      double z_top = storm_grid->GetTopSurface().GetZ(x,y);

      if (!storm_grid->GetTopSurface().IsMissing(z_top) && !storm_grid->GetBotSurface().IsMissing(z_bot)) {
        z_bot         -= z0;
        z_top         -= z0;
        int first_data = static_cast<int>(floor((z_top+z_shift)/dz));
        int end_data   = static_cast<int>(floor((z_bot-z_shift)/dz));

        if (end_data > nz) {
          printf("Internal warning: SEGY-grid too small (%d, %d needed). Truncating data.\n", nz, end_data);
          end_data = nz;
        }
        for (k = 0; k < first_data; k++) {
          data_vec[k] = 0.0;
        }

        for (k = first_data; k < end_data; k++) {
          z           = z0 + k*dz + static_cast<float>(z_shift);
          data_vec[k] = float(storm_grid->GetValueZInterpolated(x,y,z));
        }
        for (k = end_data; k < nz; k++) {
          data_vec[k] = 0.0;
        }
      }
      else {
        for (k = 0; k < nz; k++) {
          data_vec[k] = 0.0;
        }
      }

      StoreTrace(x, y, data_vec, NULL);
    }

  }

  if (write_to_file)
    WriteAllTracesToFile();

}

SegY::~SegY()
{
  if (geometry_!=NULL)
    delete geometry_;
  geometry_ = NULL;
  if (binary_header_!=NULL)
    delete binary_header_;
  binary_header_ = NULL;
  if (traces_.size() > 0) {
    for (size_t i = 0; i < n_traces_; i++)
      delete traces_[i];
  }
  file_.close();
}

SegYTrace *
SegY::GetNextTrace(double zPad, const Volume * volume, bool onlyVolume )
{
  bool outsideSurface = false;
  bool duplicateHeader; // Only needed for memory allocations in readAllTraces()
  if (single_trace_ == false) {
    return NULL;
  }
  try {
    return ReadTrace(volume,
                     zPad,
                     duplicateHeader,
                     onlyVolume,
                     outsideSurface);
  }
  catch (EndOfFile& ) {
    return NULL;
  }
}

void
SegY::GetTraceData(int IL, int XL, std::vector<float> & result, const Volume * volume)
{
  if(geometry_ == NULL)
    throw Exception("Can not find data without geometry set.\n");

  float x,y;
  geometry_->FindXYFromILXL(IL, XL, x, y);
  GetTraceData(x, y, result, volume);
}

void
SegY::GetTraceData(float x, float y, std::vector<float> & result, const Volume * volume)
{
  if(geometry_ == NULL)
    throw Exception("Can not find data without geometry set.\n");

  std::streampos pos = GetFilePos(x, y);

  bool outside = false;
  float z_top = -1; //Indicates inactive.
  float z_bot = -1; //Indicates inactive.
  if(volume != NULL) {
    try {
      z_top = static_cast<float>(volume->GetTopSurface().GetZ(x,y));
      z_bot = static_cast<float>(volume->GetBotSurface().GetZ(x,y));
    }
    catch (NRLib::Exception & ) {
      outside = true;
      result.resize(0);
    }
  }

  if(outside == false)
    GetTraceData(pos, result, z_top, z_bot);
}

const TraceHeader &
SegY::GetTraceHeader(int IL, int XL) const
{
  if(geometry_ == NULL)
    throw Exception("Can not find file position without geometry set.\n");
  size_t i, j;
  geometry_->FindIndex(IL, XL, i, j);
  size_t index = j*geometry_->GetNx()+i;
  if(traces_[index] == NULL)
    throw Exception("Trace is not defined.\n");

  return(traces_[index]->GetTraceHeader());
}

const TraceHeader &
SegY::GetTraceHeader(float x, float y) const
{
  if(geometry_ == NULL)
    throw Exception("Can not find file position without geometry set.\n");
  size_t i, j;
  geometry_->FindIndex(x, y, i, j);
  size_t index = j*geometry_->GetNx()+i;
  if(traces_[index] == NULL)
    throw Exception("Trace is not defined.\n");

  return(traces_[index]->GetTraceHeader());
}

std::streampos
SegY::GetFilePos(int IL, int XL) const
{
  if(geometry_ == NULL)
    throw Exception("Can not find file position without geometry set.\n");
  else if(single_trace_ == false)
    throw Exception("Can not find file position for memory-stored grid.\n");

  size_t i, j;
  geometry_->FindIndex(IL, XL, i, j);
  size_t index = j*geometry_->GetNx()+i;
  if(traces_[index] != NULL)
    return(traces_[index]->GetFilePos());
  else
    return(0);
}

std::streampos
SegY::GetFilePos(float x, float y) const
{
  if(geometry_ == NULL)
    throw Exception("Can not find file position without geometry set.\n");
  else if(single_trace_ == false)
    throw Exception("Can not find file position for memory-stored grid.\n");

  size_t i, j;
  geometry_->FindIndex(x, y, i, j);
  size_t index = j*geometry_->GetNx()+i;
  if(traces_[index] != NULL)
    return(traces_[index]->GetFilePos());
  else
    return(0);
}

void
SegY::GetTraceData(std::streampos pos, std::vector<float> & result, float z_top, float z_bot)
{
  assert(file_);

  size_t j0, j1;

  if(z_top >= 0)
    j0 = static_cast<size_t>((z_top - z0_)/dz_);
  else
    j0 = 0;

  if(z_bot > 0)
    j1 = static_cast<size_t>((z_bot - z0_)/dz_);
  else
    j1 = nz_-1;

  if (j0 < 0)
    j0 = 0;
  if (j1 > nz_ - 1)
    j1 = nz_ - 1;

  if (j0 > j1)
    throw Exception(" Lower horizon above SegY region or upper horizon below SegY region");

  FILE * seek_file = fopen(file_name_.c_str(),"rb");
  Seek(seek_file, pos, SEEK_SET);
  char * buffer = new char[nz_*datasize_];
  size_t n_read = fread(buffer, datasize_, nz_, seek_file);
  if (n_read < nz_)
    throw Exception("Failed to read from SEGY-file or unexpected end of file.");

  result.resize(nz_);
  switch(binary_header_->GetFormat()) {
    case 1:
      for (size_t i = j0; i <= j1; ++i)
        ParseIBMFloatBE(&buffer[4*i], result[i]);
      break;
    case 2: {
        int tmp;
        for (size_t i = j0; i <= j1; ++i) {
          ParseInt32BE(&buffer[4*i], tmp);
          result[i] = static_cast<float>(tmp);
        }
      }
      break;
    case 3: {
        short tmp;
        for (size_t i = j0; i <= j1; ++i) {
          ParseInt16BE(&buffer[2*i], tmp);
          result[i] = static_cast<float>(tmp);
        }
      }
      break;
    case 5:
      for (size_t i = j0; i <= j1; ++i)
        ParseIEEEFloatBE(&buffer[4*i], result[i]);
      break;
    default:
      assert(0); //We should catch this much earlier.
  }
  delete [] buffer;
  fclose(seek_file);
}

void
SegY::ReadAllTraces(const Volume * volume,
                    double         zPad,
                    bool           onlyVolume,
                    bool           relative_padding)
{
  single_trace_ = false;
  traces_.resize(n_traces_);

  LogKit::LogMessage(LogKit::Low,"\nReading SEGY file " );
  LogKit::LogMessage(LogKit::Low, file_name_);

  bool outsideSurface = false;
  bool duplicateHeader; // Needed for memory allocations.
  double outsideTopBot[6];
  double outsideTopMax[6]; //Largest lack of data top
  double outsideBotMax[6]; //Largest lack of data bot

  traces_[0] = ReadTrace(volume,
                         zPad,
                         duplicateHeader,
                         onlyVolume,
                         outsideSurface,
                         true,
                         outsideTopBot,
                         relative_padding);
  int k;
  for (k=0;k<6;k++) {
    outsideTopMax[k] = outsideTopBot[k];
    outsideBotMax[k] = outsideTopBot[k];
  }
  double writeInterval = 0.02;
  double nextWrite = writeInterval;
  LogKit::LogMessage(LogKit::Low,"\n  0%        20%      40%       60%       80%       100%");
  LogKit::LogMessage(LogKit::Low,"\n  |    |    |    |    |    |    |    |    |    |    |  ");
  LogKit::LogMessage(LogKit::Low,"\n  ^");
  size_t traceSize = datasize_ * nz_ + 240;
  size_t fSize = 3600 + n_traces_ * traceSize;
  long long bytesRead = 3600+traceSize;
  for (unsigned int i=1 ; i < static_cast<unsigned int>(n_traces_) ; i++)
  {
    double percentDone = bytesRead/static_cast<double>(fSize);
    if (percentDone > nextWrite)
    {
      LogKit::LogMessage(LogKit::Low,"^");
      nextWrite+=writeInterval;
    }

    try {
      traces_[i] = ReadTrace(volume,
                             zPad,
                             duplicateHeader,
                             onlyVolume,
                             outsideSurface,
                             false,
                             outsideTopBot,
                             relative_padding);
    }
    catch (EndOfFile& ) {
      break;
    }

    if (outsideTopBot[0] > outsideTopMax[0])
      for (k=0;k<6;k++)
        outsideTopMax[k] = outsideTopBot[k];
    if (outsideTopBot[1] > outsideBotMax[1])
      for (k=0;k<6;k++)
        outsideBotMax[k] = outsideTopBot[k];

    bytesRead += traceSize;
    if (duplicateHeader)
      bytesRead += 3600;
  }
  LogKit::LogMessage(LogKit::Low,"^\n");
  n_traces_ = traces_.size();

  if (outsideTopBot[0] > outsideTopMax[0])
    for (k=0;k<6;k++)
      outsideTopMax[k] = outsideTopBot[k];
  if (outsideTopBot[1] > outsideBotMax[1])
    for (k=0;k<6;k++)
      outsideBotMax[k] = outsideTopBot[k];

  try {
    CheckTopBotError(outsideTopMax, outsideBotMax); //Throws exception if > 0.
  }
  catch (NRLib::Exception & e) {
    throw Exception(e.what());
  }

  int count = 0;
  for (unsigned int i=1 ; i<traces_.size() ; i++)
    if (traces_[i] != NULL)
      count++;
  if (count == 0)
  {
    std::string text;
    text += " No valid traces found. The specified time surfaces do not cover any part of the\n";
    text += " seismic data. The reason can be that you have given an incorrect SEGY trace header\n";
    text += " format, or that you need to bypass the coordinate scaling.";
    throw Exception(text);
  }
}

void
SegY::CheckTopBotError(const double * tE, const double * bE)
{
  std::string text = "";
  if (tE[0] > 0.0) {
    text+= "There is a region between the top surface of the inversion interval and the seismic data volume \n";
    text+= "with no seismic data. The largest gap is for the seismic trace at position ("+ToString(tE[2],0)+","+ToString(tE[3],0)+") \n";
    text+= "where the surface z-value = "+ToString(tE[4],2)+" while the seismic start time = "+ToString(z0_,2)+". Please reduce \n";
    text+= "the start time of your seismic data or lower the top surface "+ToString(tE[0],2)+"ms.\n";
  }
  if (bE[1] > 0.0) {
    double sMax = bE[5] - bE[1];
    text+= "There is a region between the base surface of the inversion interval and the seismic data volume \n";
    text+= "with no seismic data. The largest gap is for the seismic trace at position ("+ToString(bE[2], 0)+","+ToString(bE[3], 0)+") \n";
    text+= "where the surface z-value = "+ToString(bE[5], 2)+" while the seismic end time = "+ToString(sMax, 2)+". Please include \n";
    text+= "more seismic data or heighten the base surface "+ToString(bE[1], 2)+"ms.\n";
  }
  if (text != "")
    throw Exception(text);
}


SegYTrace *
SegY::ReadTrace(const Volume * volume,
                double         zPad,
                bool         & duplicateHeader,
                bool           onlyVolume,
                bool         & outsideSurface,
                bool           writevalues,
                double       * outsideTopBot,
                bool           relative_padding)
{
  TraceHeader traceHeader(trace_header_format_);

  duplicateHeader = ReadHeader(traceHeader);
  if (writevalues == 1)
    traceHeader.WriteValues();

  if (outsideTopBot != NULL) {
    outsideTopBot[0] = 0; // > 0 indicates top error
    outsideTopBot[1] = 0; // > 0 indicates bot error
  }
  double x, y;
  if (trace_header_format_.GetCoordSys() == TraceHeaderFormat::UTM) {
    x = traceHeader.GetUtmx();
    y = traceHeader.GetUtmy();
  }
  else if (trace_header_format_.GetCoordSys() == TraceHeaderFormat::ILXL) {
    x = static_cast<double>(traceHeader.GetInline());
    y = static_cast<double>(traceHeader.GetCrossline());
  }
  else {
   throw Exception("Invalid coordinate system number ("
                   +ToString(trace_header_format_.GetCoordSys())+")");
  }

  size_t j0 = 0;
  size_t j1 = nz_-1;
  float zTop, zBot;
  if (volume != NULL)
  {
    if (onlyVolume && !volume->IsInside(x,y))
    {
      ReadDummyTrace(file_,binary_header_->GetFormat(),nz_);
      return(NULL);
    }

    try {
      zTop = static_cast<float>(volume->GetTopSurface().GetZ(x,y));
    }
    catch (NRLib::Exception & ) {
      outsideSurface = true;
      ReadDummyTrace(file_,binary_header_->GetFormat(),nz_);
      return(NULL);
    }

    try {
      zBot = static_cast<float>(volume->GetBotSurface().GetZ(x,y));
    }
    catch (NRLib::Exception & ) {
      outsideSurface = true;
      ReadDummyTrace(file_,binary_header_->GetFormat(),nz_);
      return(NULL);
    }

    if (volume->GetTopSurface().IsMissing(zTop) || volume->GetBotSurface().IsMissing(zBot))
    {
      ReadDummyTrace(file_,binary_header_->GetFormat(),nz_);
      return(NULL);
    }
  }
  else {
    zTop = z0_;
    zBot = z0_ + nz_*dz_;
  }

  if (zTop < z0_) {
    if (outsideTopBot == NULL) {
      std::string text;
      text+= "For the trace at position ("+ToString(x,0)+","+ToString(y,0)+") the top time surface";
      text+= " reaches below\n the seismic gather. Seismic start time = "+ToString(z0_,2);
      text+= " while surface z-value = "+ToString(zTop,2)+".\n Include more seismic data or lower";
      text+= " the top surface "+ToString(z0_-zTop,2)+"ms";
      throw Exception(text);
    }
    else {
      outsideTopBot[0] = z0_-zTop;
      outsideTopBot[2] = x;
      outsideTopBot[3] = y;
      outsideTopBot[4] = zTop;
    }
  }
  double sBot = z0_ + nz_*dz_;
  if (zBot > sBot) {
    if (outsideTopBot == NULL) {
      std::string text;
      text+= "For the trace at position ("+ToString(x,0)+","+ToString(y,0)+") the base time surface";
      text+= " reaches below\n the seismic gather. Seismic end time = "+ToString(sBot,2);
      text+= " while surface z-value = "+ToString(zBot,2)+".\n Include more seismic data or heighten";
      text+= " the base surface "+ToString(zBot-sBot,2)+"ms";
      throw Exception(text);
    }
    else {
      outsideTopBot[1] = zBot-sBot;
      outsideTopBot[2] = x;
      outsideTopBot[3] = y;
      outsideTopBot[5] = zBot;
    }
  }
  if (outsideTopBot != NULL && (outsideTopBot[0] > 0.0 || outsideTopBot[1] > 0.0)) {
    ReadDummyTrace(file_,binary_header_->GetFormat(),nz_);
    return(NULL);
  }

  float pad;
  if (relative_padding)
    pad = static_cast<float>(0.5*zPad*(zBot - zTop));
  else
    pad = static_cast<float>(0.5 * zPad);

  int j0_temp = static_cast<int>((zTop - pad - z0_)/dz_); // Use <int> to handle negative numbers
  j1 = static_cast<size_t>((zBot + pad - z0_)/dz_);

  if (j0_temp < 0)
    j0 = 0;
  else
    j0 = static_cast<size_t>(j0_temp);
  if (j1 > nz_ - 1)
    j1 = nz_ - 1;

  if (j0 > j1)
    throw Exception(" Lower horizon above SegY region or upper horizon below SegY region");

  SegYTrace * trace = NULL;
  if (file_.eof() == false)
  {
    // Copy elements from j0 til j1.
    trace = new SegYTrace(file_, j0, j1,
                          binary_header_->GetFormat(), nz_,
                          &traceHeader);
  }
  return trace;
}

bool
SegY::ReadHeader(TraceHeader & header)
{
  bool duplicateHeader;
  header.Read(file_,binary_header_->GetLino());
  switch(header.GetStatus()) {
  case -1:
    binary_header_->Update(file_);
    header.Read(file_, binary_header_->GetLino());
    duplicateHeader = true;  // Duplicate header found
    break;
  default:
    duplicateHeader = false;
    break;
  }
  if (header.GetDt()/1000 != dz_) {
    if(dz_ == 0)
      dz_ = static_cast<float>(header.GetDt()/1000.0);
    else if(header.GetDt() > 0) {
      std::string error = "Different sampling densities given.";
      error += "Initial sampling density of "+ToString(dz_)+" ms changed to " +
        ToString(header.GetDt()/1000.0) + " in trace with XL " +
        ToString(header.GetCrossline()) + " and inline " +
        ToString(header.GetInline()) + ".\n";
      throw(Exception(error));
    }
  }
  return duplicateHeader;
}

void
SegY::SetGeometry(const SegyGeometry * geometry)
{
  geometry_ = new SegyGeometry(geometry);
  n_traces_  = geometry_->GetNx()*geometry_->GetNy();
  traces_.resize(n_traces_);
  for (size_t i = 0 ; i < n_traces_ ; i++)
  {
    traces_[i] = NULL;
  }
}

void SegY::CreateRegularGrid()
{
  geometry_  = new SegyGeometry(traces_);
  n_traces_  = static_cast<int>(traces_.size());
}

std::vector<float>
SegY::GetAllValues(void)
{
  size_t i,nTot = 0;
  //int nTraces = nx_*ny_;
  for (i = 0; i < n_traces_; i++)
    if (traces_[i] != NULL)
      nTot += traces_[i]->GetEnd() - traces_[i]->GetStart() + 1;

  std::vector<float> result(nTot);
  size_t k, kS, kE, oInd = 0;
  for (i = 0; i < n_traces_; i++)
    if (traces_[i] != NULL)
    {
      kS = traces_[i]->GetStart();
      kE = traces_[i]->GetEnd();
      for (k = kS; k <= kE; k++)
      {
        result[oInd] = traces_[i]->GetValue(k);
        oInd++;
      }
    }
    return(result);
}

size_t
SegY::FindNumberOfSamplesInLongestTrace(void) const
{
  size_t max_length = 0;
  for (size_t i=0 ; i < traces_.size() ; i++) {
    if (traces_[i] != NULL) {
      size_t istart = traces_[i]->GetStart();
      size_t iend   = traces_[i]->GetEnd();
      max_length = std::max(max_length, iend - istart + 1);
    }
  }
  return max_length;
}


void
SegY::GetNearestTrace(std::vector<float> & trace_data,
                      bool               & missing,
                      float              & z0_data,
                      float                x,
                      float                y) const
{
  size_t i = geometry_->FindIndex(x, y);

  if (traces_[i] != NULL) {
    trace_data = traces_[i]->GetTrace();
    // NBNB: The 0.5f below is a shift we have introduced when reading
    // in seismic data to get data values in centre of grid cells rather
    // than on their borders. This choice and its implications need to
    // be looked at more carefully.
    z0_data    = z0_ + (traces_[i]->GetStart() + 0.5f)*dz_;
    missing    = false;
  }
  else {
    missing = true;
  }
}


float
SegY::GetValue(double x, double y, double z, int outsideMode) const
{
  if(geometry_ == NULL)
    throw Exception("Geometry is not defined.\n");

  int i, j;
  float xind,yind;
  float value;
  double x0 = geometry_->GetX0()+0.5*geometry_->GetDx()*geometry_->GetCosRot()-0.5*geometry_->GetDy()*geometry_->GetSinRot();
  double y0 = geometry_->GetY0()+0.5*geometry_->GetDy()*geometry_->GetCosRot()+0.5*geometry_->GetDx()*geometry_->GetSinRot();
  double sx =  (x-x0)*geometry_->GetCosRot() + (y-y0)*geometry_->GetSinRot() + 0.5*geometry_->GetDx();
  double sy = -(x-x0)*geometry_->GetSinRot() + (y-y0)*geometry_->GetCosRot() + 0.5*geometry_->GetDy();
  if (geometry_!=NULL)
  {
    int ok = geometry_->FindContIndex(static_cast<float>(x),static_cast<float>(y),xind,yind);

    i = static_cast<int>(xind);
    j = static_cast<int>(yind);
    size_t nx = geometry_->GetNx();
    size_t ny = geometry_->GetNy();

    size_t index;

    index = j*nx+i; //NBNB er dette rett??

    if (ok==1 && traces_[index]!=0 && z>=z0_ && z<=z0_+nz_*dz_)
    {
      size_t zind = static_cast<size_t>(floor((z-z0_)/dz_));  //NBNB   irap grid rounding different

      float v1 = traces_[index]->GetValue(zind);
      if (v1 == rmissing_ && outsideMode == CLOSEST)
      {
        zind = traces_[index]->GetLegalIndex(zind);
        v1 = traces_[index]->GetValue(zind);
        if (traces_[index]->GetValue(zind-1) == rmissing_)
          z = z0_+zind*dz_;          // Want edge value, hence 0/1 dz_ added
        else                         // (0.5 would give center of cell).
          z = z0_+(zind+0.99f)*dz_;
      }
      if (v1 != rmissing_)
      {
        // Computes interpolated value ax^2+by^2+cz^2+dx+ey+fz+g.
        // abcdefg estimated from closest point and its closest neighbours.
        size_t maxInd = nx*ny - 1;
        float v0, v2, a, b, c, d, e, f, g;

        // Along x:
        v0 = rmissing_;
        v2 = rmissing_;
        if (index >= 1 && traces_[index-1] != NULL)
          v0 = traces_[index-1]->GetValue(zind);
        if (index+1 <= maxInd && traces_[index+1] != NULL)
          v2 = traces_[index+1]->GetValue(zind);
        if (v0 == rmissing_)
        {
          a = 0;
          if (v2 == rmissing_)
            d = 0;
          else
            d = v2 - v1; // Using unit coordinates in each direction.
        }
        else if (v2 == rmissing_)
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
        size_t tmpInd;
        if(index >= nx) {
          tmpInd = index - nx;
          if (tmpInd <= maxInd && traces_[tmpInd] != NULL)
            v0 = traces_[tmpInd]->GetValue(zind);
        }
        tmpInd = index + nx;
        if (tmpInd <= maxInd && traces_[tmpInd] != NULL)
          v2 = traces_[tmpInd]->GetValue(zind);
        if (v0 == rmissing_)
        {
          b = 0;
          if (v2 == rmissing_)
            e = 0;
          else
            e = v2 - v1; // Using unit coordinates in each direction.
        }
        else if (v2 == rmissing_)
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
        v0 = traces_[index]->GetValue(zind-1);
        v2 = traces_[index]->GetValue(zind+1);
        if (v0 == rmissing_)
        {
          c = 0;
          if (v2 == rmissing_)
            f = 0;
          else
            f = v2 - v1; //Using unit coordinates in each direction.
        }
        else if (v2 == rmissing_)
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
        double dx = geometry_->GetDx();
        double dy = geometry_->GetDy();
        float ux = static_cast<float>(sx/dx) - static_cast<float>(floor(sx/dx) + 0.5);
        float uy = static_cast<float>(sy/dy) - static_cast<float>(floor(sy/dy) + 0.5);
        float uz = static_cast<float>((z-z0_)/dz_) - static_cast<float>(floor((z-z0_)/dz_) + 0.5);
        value = a*ux*ux+b*uy*uy+c*uz*uz+d*ux+e*uy+f*uz+g;
      }
      else
      {
        if (outsideMode == ZERO)
          value = 0;
        else
          value = rmissing_;
      }
    }
    else
    {
      if (outsideMode == ZERO)
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
SegY::WriteMainHeader(const TextualHeader& ebcdicHeader)
{
  assert(file_);
  ebcdicHeader.Write(file_);
  if (binary_header_ != NULL)
    delete binary_header_;
  binary_header_ = new BinaryHeader();
  binary_header_->Write(file_, geometry_, dz_, nz_);
}

void
SegY::StoreTrace(double x, double y, const std::vector<float> data, const Volume *volume, float topVal,float baseVal)
{
  assert(file_);
  assert(geometry_ != 0);
 // TraceHeader header(trace_header_format_);
 // header.SetNSamples(nz_);
 // header.SetDt(static_cast<unsigned short>(dz_*1000));
 // header.SetScalCo(scalcoinitial);
 // header.SetUtmx(x);
 // header.SetUtmy(y);
  int IL,XL;
  bool ok = geometry_->IsInside(x,y);
  if (ok==true)
  {
    geometry_->FindILXL(x,y,IL,XL);
  //  header.SetInline(IL);
  //  header.SetCrossline(XL);

    // header.write(file_);
    double ztop;
    if (volume != NULL)
      ztop = volume->GetTopSurface().GetZ(x,y);
    else
      ztop = z0_;
    std::vector<float>  trace(nz_);
    size_t k;
    if (volume != NULL && volume->GetTopSurface().IsMissing(ztop))
    {
      for (k = 0; k < nz_; k++)
        trace[k] = 0;
    }
    else
    {
      size_t firstData = static_cast<size_t>((ztop-z0_) / dz_);
      for (k = 0; k < firstData; k++)
        trace[k] = topVal; //data[0];
      for (k = firstData; k < firstData + data.size(); k++)
        trace[k] = data[k - firstData];
      for (k = firstData + data.size(); k < nz_; k++)
        trace[k] = baseVal; //data[simbox_->getnz()-1];
    }

    size_t index;
    size_t j0 = 0;
    size_t j1 = nz_-1;
    size_t i,j;
    geometry_->FindIndex(x, y, i, j);
    index = i + geometry_->GetNx() * j;
    traces_[index] = new SegYTrace(trace, j0, j1, x, y, IL, XL);
    traces_[index]->SetTableIndex(index);
  }
  else
    throw Exception(" Coordinates are outside grid.");
}


void
SegY::WriteTrace(double x, double y, const std::vector<float> data, const Volume *volume, float topVal,float baseVal, short scalcoinitial)
{
  assert(file_);
  assert(geometry_ != 0);
  TraceHeader header(trace_header_format_);
  header.SetNSamples(nz_);
  header.SetDt(static_cast<unsigned short>(dz_*1000));
  header.SetScalCo(scalcoinitial);

  header.SetUtmx(x);
  header.SetUtmy(y);
  int IL,XL;
  bool ok = geometry_->IsInside(x,y);
  if (ok==true)
  {
    geometry_->FindILXL(x,y,IL,XL);
    header.SetInline(IL);
    header.SetCrossline(XL);
    header.Write(file_);
    double ztop = z0_;
    int nData = static_cast<int>(data.size());
    if (volume != NULL)
      ztop = volume->GetTopSurface().GetZ(x,y);

    std::vector<float> trace;
    trace.resize(nz_);
    size_t k;
    if (volume != NULL && volume->GetTopSurface().IsMissing(ztop))
    {
      for (k = 0; k < nz_; k++)
        trace[k] = 0;
    }
    else
    {
      size_t firstData = static_cast<size_t>((ztop-z0_)/dz_);
      for (k = 0 ; k < firstData; k++)
        trace[k] = topVal; //data[0];
      for ( ; k < firstData + nData; k++)
        trace[k] = data[k-firstData];
      for ( ; k < nz_; k++)
        trace[k] = baseVal; //data[simbox_->getnz()-1];
    }
    WriteBinaryIbmFloatArray(file_,trace.begin(),trace.end());
  }
  else
    throw Exception("Coordinates are outside grid.");
}

void
SegY::WriteTrace(const TraceHeader & origTraceHeader,
                 const std::vector<float> & data,
                 const Volume * volume, float topVal, float baseVal)
{
  TraceHeader traceHeader(origTraceHeader);

  traceHeader.SetNSamples(nz_);
  traceHeader.Dump(file_, true);

  double z = z0_;
  size_t nz = nz_;
  if (volume != NULL)
  {
    z = volume->GetTopSurface().GetZ(traceHeader.GetUtmx(),traceHeader.GetUtmy());
    nz = static_cast<int>(data.size());
  }

  std::vector<float> trace(nz_);
  size_t k;

  if (volume != NULL && volume->GetTopSurface().IsMissing(z))
  {
    for (k = 0; k < nz_; k++)
      trace[k] = 0;
  }
  else
  {
    size_t firstData = static_cast<size_t>((z-z0_)/dz_);
    for (k = 0; k < firstData; k++)
      trace[k] = topVal; //data[0];
    for (k = firstData; k < firstData + nz; k++)
      trace[k] = data[k-firstData];
    for (k = firstData + nz; k < nz_; k++)
      trace[k] = baseVal; //data[simbox_->getnz()-1];
  }
  WriteBinaryIbmFloatArray(file_, trace.begin(), trace.end());
}


//Sorting functions for WriteAllTracesToFile
bool SortILXL(const SegYTrace * t1, const SegYTrace * t2)
{
  if (t1 == NULL)
    return(false);
  if (t2 == NULL)
    return(true);
  if (t1->GetInline() < t2->GetInline())
    return(true);
  if (t1->GetInline() == t2->GetInline()) {
    if (t1->GetCrossline() < t2->GetCrossline())
      return(true);
    else
      if (t1->GetCrossline() == t2->GetCrossline() &&
            t1->GetTableIndex() < t2->GetTableIndex())
            return(true);
  }
  return(false);
}

bool SortIndex(const SegYTrace * t1, const SegYTrace * t2)
{
  if (t2 == NULL) //Note: This msut come first, since NULL,NULL must give false.
    return(false);
  if (t1 == NULL)
    return(true); //Move NULL to start.

  if (t1->GetTableIndex() < t2->GetTableIndex())
    return(true);
  else
    return(false);
}

void
SegY::WriteAllTracesToFile(short scalcoinitial)
{
  size_t i, k;
  std::vector<float>  trace(nz_);
  double x,y;

  std::sort(traces_.begin(), traces_.end(), SortILXL);

  for (i = 0; i < traces_.size(); i++)
  {
    if (traces_[i] != NULL)
    {
      for (k = 0; k < nz_; k++)
        trace[k] = traces_[i]->GetValue(k);

      // write traceheader
      x = traces_[i]->GetX();
      y = traces_[i]->GetY();
      TraceHeader header(trace_header_format_);
      header.SetNSamples(nz_);
      header.SetDt(static_cast<unsigned short>(dz_*1000));
      header.SetScalCo(scalcoinitial);
      header.SetUtmx(static_cast<double>(x));
      header.SetUtmy(static_cast<double>(y));
      header.SetInline(traces_[i]->GetInline());
      header.SetCrossline(traces_[i]->GetCrossline());
      header.Write(file_);
      WriteBinaryIbmFloatArray(file_,trace.begin(),trace.end());
    }
  }
  sort(traces_.begin(), traces_.end(), SortIndex);
  //Traces are sorted, but NULL traces are all at beginning, instead of at correct location.
  for (i = 0; i < traces_.size(); i++) {
    if (traces_[i] != NULL) {
      size_t index = traces_[i]->GetTableIndex();
      if(i > index) {
        traces_[index] = traces_[i];
        traces_[i] = NULL;
      }
      else
        i = traces_.size(); //No more NULL traces.
    }
  }
}

size_t
SegY::FindNumberOfTraces(const std::string       & fileName,
                         const TraceHeaderFormat * traceHeaderFormat)
{
  float dummy_z0 = 0.0f;
  if (traceHeaderFormat!=NULL)
  {
    SegY segy(fileName, dummy_z0, (*traceHeaderFormat));
    size_t size = segy.FindNumberOfTraces();
    return size;
  }
  else
  {
    SegY segy(fileName, dummy_z0);
    size_t size = segy.FindNumberOfTraces();
    return size;
  }
}

size_t
SegY::FindNumberOfTraces(void)
{
  unsigned long long fSize = FindFileSize(file_name_);
  n_traces_ = static_cast<size_t>(ceil( (static_cast<double>(fSize)-3600.0)/
                                        static_cast<double>(datasize_*nz_+240.0)));

  TraceHeader traceHeader(trace_header_format_);
  ReadHeader(traceHeader);

  char * buffer = new char[nz_ * datasize_];
  file_.read(buffer, static_cast<std::streamsize>(nz_ * datasize_));

  size_t ntr = 1;
  for (size_t i = 1; i < n_traces_; i++)
  {
    if (file_.eof() == false)
    {
      ntr++;
      ReadHeader(traceHeader);
      file_.read(buffer, static_cast<std::streamsize>(nz_ * datasize_));
    }
  }
  delete [] buffer;
  return(ntr);
}

TraceHeaderFormat
SegY::FindTraceHeaderFormat(const std::string & fileName)
{
  float z0 = 0.0f;
  SegY segy(fileName,z0);
  TraceHeaderFormat thf = segy.GetTraceHeaderFormat();
  return thf;
}

SegyGeometry *
SegY::FindGridGeometry(const std::string       & fileName,
                       const TraceHeaderFormat *traceHeaderFormat)
{
  float dummy_z0 = 0.0f;
  if (traceHeaderFormat!=NULL)
  {
    SegY segy(fileName, dummy_z0, (*traceHeaderFormat));
    SegyGeometry * geometry = segy.FindGridGeometry(); // returns a new SegyGeometry object
    return geometry;
  }
  else
  {
    SegY segy(fileName, dummy_z0);
    SegyGeometry * geometry = segy.FindGridGeometry(); // returns a new SegyGeometry object
    return geometry;
  }
}

void
SegY::FindAndSetGridGeometry(bool only_ilxl, bool keep_header)
{
  geometry_ = FindGridGeometry(only_ilxl, keep_header);
}

SegyGeometry *
SegY::FindGridGeometry(bool only_ilxl, bool keep_header)
{
  if(geometry_ != NULL)
    return(new SegyGeometry(geometry_));

  if(file_.tellg() != static_cast<std::streampos>(3600))
    throw(Exception("Can not find SegY geometry for a file where traces have already been read.\n"));

  TraceHeader traceHeader(trace_header_format_);

  std::streampos pos  = 3840;
  std::streampos step = static_cast<std::streampos>(nz_*datasize_+240);

  size_t i;
  traces_.resize(n_traces_);
  char * buffer = new char[nz_*datasize_];
  for (i = 0; i < n_traces_; i++)
  {
    try {
      if (file_.eof()==false)
      {
        bool extra_header = ReadHeader(traceHeader);
        traces_[i] = new SegYTrace(traceHeader,keep_header);
        file_.read(buffer, static_cast<std::streamsize>(nz_*datasize_));
        if(only_ilxl == true)
          traces_[i]->RemoveXY();
        traces_[i]->SetFilePos(pos);
        pos += step;
        if(extra_header == true)
          pos += 3600;
      }
    }
    catch(Exception & e) {
      throw(Exception("In trace number " + ToString(i) + ":\n" + e.what()));
    }
  }
  delete [] buffer;


  SetBogusILXLUndefined(traces_);

  SegyGeometry * geometry = new SegyGeometry(traces_);
  n_traces_  = static_cast<int>(traces_.size());
  return(geometry);
}

void
SegY::SetBogusILXLUndefined(std::vector<NRLib::SegYTrace*> & traces)
{
  //--------------------------------------------------------------------
  // Try to identify incorrect IL-XL numbers
  //
  // 1. Tabulate all IL and XL numbers
  // 2. Sort in increasing order
  // 3. List the IL and XL steps and count how often they occur
  // 4. Set the most frequent steps as correct
  // 5. Identify and traces with XL IL numbers not matching step
  // 6. Remove bogus traces
  //--------------------------------------------------------------------

  std::list<int> ilList;
  std::list<int> xlList;
  //
  // 1. Tabulate all IL and XL numbers
  //
  for (size_t i = 0; i < traces.size() ; i++) {
    int IL = traces[i]->GetInline();
    int XL = traces[i]->GetCrossline();
    if (find(ilList.begin(), ilList.end(), IL) == ilList.end())
      ilList.push_back(IL);
    if (find(xlList.begin(), xlList.end(), XL) == xlList.end())
      xlList.push_back(XL);
  }

  //
  // 2. Sort in increasing order
  //
  ilList.sort();
  xlList.sort();

  // Move to vector to simplify step handling
  std::vector<int> ils;
  std::vector<int> xls;
  std::copy(ilList.begin (), ilList.end (), std::back_inserter(ils));
  std::copy(xlList.begin (), xlList.end (), std::back_inserter(xls));

  //
  // 3. List the IL and XL steps and count how often they occur
  //
  std::map<int,int> ilstep;
  std::map<int,int> xlstep;

  size_t i = 1;
  while(i < ils.size()) {
    int step = ils[i] - ils[i-1];
    if (ilstep.find(step) == ilstep.end())
      ilstep[step] = 1;
    else
      ilstep[step]++;
    i++;
  }
  i = 1;
  while(i < xls.size()) {
    int step = xls[i] - xls[i-1];
    if (xlstep.find(step) == xlstep.end())
      xlstep[step] = 1;
    else
      xlstep[step]++;
    i++;
  }

  if (ilstep.size() > 1 || xlstep.size() > 1) { // Problem identified
    //
    // 4. Set the most frequent steps as correct
    //
    int step  = 0;
    int value = 0;
    std::map<int,int>::iterator m;
    for (m = ilstep.begin() ; m != ilstep.end() ; m++ ) {
      if (m->second > value) {
        step  = m->first;
        value = m->second;
      }
    }
    int ILstep = step;

    step  = 0;
    value = 0;
    for (m = xlstep.begin() ; m != xlstep.end() ; m++ ) {
      if (m->second > value) {
        step  = m->first;
        value = m->second;
      }
    }
    int XLstep = step;

    //
    // 5. Identify and traces with IL and XL numbers not matching step
    //
    std::list<int> ilWrong;
    std::list<int> xlWrong;

    i = 1;
    while (i + 1 < ils.size()) {
      int step0 = ils[i    ] - ils[i - 1];
      int step1 = ils[i + 1] - ils[i    ];
      if (step0 != ILstep)
        ilWrong.push_back(ils[i - 1]);
      else if (step1 != ILstep)
        ilWrong.push_back(ils[i + 1]);
      i++;
    }
    i = 1;
    while (i + 1 < xls.size()) {
      int step0 = xls[i    ] - xls[i - 1];
      int step1 = xls[i + 1] - xls[i    ];
      if (step0 != XLstep)
        xlWrong.push_back(xls[i - 1]);
      else if (step1 != XLstep)
        xlWrong.push_back(xls[i + 1]);
      i++;
    }

    //
    // 6. Remove bogus traces
    //
    std::vector<int> removedID;
    std::vector<int> removedIL;
    std::vector<int> removedXL;
    for (i = 0; i < traces.size() ; i++) {
      int IL = traces[i]->GetInline();
      int XL = traces[i]->GetCrossline();
      if (find(ilWrong.begin(), ilWrong.end(), IL) != ilWrong.end() ||
          find(xlWrong.begin(), xlWrong.end(), XL) != xlWrong.end()) {
           traces[i]->SetILXLUndefined();
           removedID.push_back(static_cast<int>(i));
           removedIL.push_back(IL);
           removedXL.push_back(XL);
        }
    }

    LogKit::LogFormatted(LogKit::Warning,"\nWarning: There are more than one step size for IL and/or XL in the data\n\n");
    LogKit::LogFormatted(LogKit::Warning,"          Step     Count\n");
    LogKit::LogFormatted(LogKit::Warning,"------------------------\n");
    for (m = ilstep.begin() ; m != ilstep.end() ; m++ )
      LogKit::LogFormatted(LogKit::Warning,"IL    %8d  %8d\n", m->first, m->second);
    for (m = xlstep.begin() ; m != xlstep.end() ; m++ )
      LogKit::LogFormatted(LogKit::Warning,"XL    %8d  %8d\n", m->first, m->second);
    LogKit::LogFormatted(LogKit::Warning,"\nThe corrects steps are assumed to be IL: %d and XL: %d\n",ILstep, XLstep);
    LogKit::LogFormatted(LogKit::High,"\nThe following IL and XL numbers are treated as incorrect\n");
    std::list<int>::iterator it;
    for (it = ilWrong.begin() ; it != ilWrong.end() ; it++ )
      LogKit::LogFormatted(LogKit::High,"  IL  %8d\n", (*it));
    for (it = xlWrong.begin() ; it != xlWrong.end() ; it++ )
      LogKit::LogFormatted(LogKit::High,"  XL  %8d\n", (*it));
    LogKit::LogFormatted(LogKit::Warning,"\nWe have found %d traces with bogus IL and XL numbers.\n",removedID.size());
    LogKit::LogFormatted(LogKit::High,"\nThese are (counting from zero):\n");
    for (i = 0 ; i < removedID.size() ; i++ )
      LogKit::LogFormatted(LogKit::High,"%10d :  (IL, XL) = (%5d,%5d)\n", removedID[i], removedIL[i], removedXL[i]);
    LogKit::LogFormatted(LogKit::High,"These IL and XL values will be set undefined.\n");
  }
}


void SegY::ReadDummyTrace(std::fstream & file, int format, size_t nz)
{
std::vector<float> predata;
  predata.resize(nz);

  if (format==1)
  {
    //IBM
    ReadBinaryIbmFloatArray(file, predata.begin(), nz);

  }
  else if (format==2)
  {
    std::vector<int> b(nz);
    ReadBinaryIntArray(file, b.begin(), nz);

  }
  else if (format==3)
  {
    std::vector<short> b(nz);
    ReadBinaryShortArray(file, b.begin(), nz);

  }
  else if (format==5)
  {
    ReadBinaryFloatArray(file, predata.begin(), nz);

  }
  else
    throw FileFormatError("Bad format");

}

bool
SegY::TraceHeaderOK(std::fstream &file, const TraceHeaderFormat *headerFormat)
{
  char buffer[240];
  memset(buffer, 0, 240);
  file.read(buffer, 240);

  std::string dummy;
  dummy.assign(buffer, 240);

  TraceHeader *t1 = new TraceHeader(*headerFormat);
  std::istringstream header1(dummy);
  t1->Read(header1, binary_header_->GetLino());
  if (t1->GetCrossline() < 0 || t1->GetInline() < 0
       || t1->GetUtmx() < 0.0 || t1->GetUtmy() < 0.0)
  {
    delete t1;
    return false;
  }
  ReadDummyTrace(file, binary_header_->GetFormat(), nz_);
  file.read(buffer, 240);
  std::string dummy2;
  dummy2.assign(buffer, 240);
  TraceHeader *t2 = new TraceHeader(*headerFormat);
  std::istringstream header2(dummy2);
  t2->Read(header2, binary_header_->GetLino());
  if (t2->GetCrossline() < 0 || t2->GetInline() < 0
       || t2->GetUtmx() < 0.0 || t2->GetUtmy() < 0.0)
  {
    delete t1;
    delete t2;
    return false;
  }
  int delta = -1;
  int deltail1, deltaxl1;
  bool ok = CompareTraces(t1,t2,delta, deltail1, deltaxl1); // Check that distance between t1 and t2 is multiplum of 6.25
  if(ok==false)
  {
    delete t1;
    delete t2;
    return false;
  }

  bool contin = true;
  std::string dummy3;
  TraceHeader *t3 = NULL;
  int delta2 = -1;
  int deltail2, deltaxl2;
  // Find point number 3, not on line with the other
  while(contin== true && !file.eof())
  {
    ReadDummyTrace(file, binary_header_->GetFormat(), nz_);
    file.read(buffer, 240);
    dummy3.assign(buffer, 240);
    std::istringstream header3(dummy3);
    if(t3!=NULL)
      delete t3;
    t3 = new TraceHeader(*headerFormat);
    t3->Read(header3, binary_header_->GetLino());
    if (t3->GetCrossline() < 0 || t3->GetInline() < 0
      || t3->GetUtmx() < 0.0 || t3->GetUtmy() < 0.0)
    {
      delete t1;
      delete t2;
      delete t3;
      return false;
    }
    ok = CompareTraces(t2,t3,delta2, deltail2, deltaxl2);
    if(ok==false)
    {
      delete t1;
      delete t2;
      delete t3;
      return false;
    }
    if(delta==1 && delta2!=1) // 1 og 2 p[ samme inline, 3 p[ annen
    {
      contin = false;
    }
    else if(delta==2 && delta2!=2) // 1 og 2 p[ samme crossline, 3 p[ annen
    {
      contin = false;
    }
    else if(delta==0) // 1 og 2 har ulik inline og crossline. Finn t3 som ikke ligger p[ linje
    {
      if(delta2==0) // check if straight line
      {
        double s1 = (t1->GetUtmy()-t2->GetUtmy())/(t1->GetUtmx()-t2->GetUtmx());
        double s2 = (t2->GetUtmy()-t3->GetUtmy())/(t2->GetUtmx()-t3->GetUtmx());
        if(fabs(s1-s2)>0.01)
          contin = false;
      }
      else
        contin = false;
    }
    else if(deltaxl1 !=0 && deltaxl1*deltaxl2<=0)
    {
      if(deltail2==0)
      {
        delete t1;
        delete t2;
        delete t3;
        return false;
      }
    }
    else if(deltail1 !=0 && deltail1*deltail2<=0)
    {
      if(deltaxl2==0)
      {
        delete t1;
        delete t2;
        delete t3;
        return false;
      }

    }
  }
  if(file.eof()==true && contin==true)// a line
  {
    delete t1;
    delete t2;
    delete t3;
    return true;
  }

//Given t1, t2 and t3, calculate dxXL, ...
  double dxXL, dyXL, dxIL, dyIL;

  FindDeltaILXL(t1,t2,t3,dxIL,dxXL,true);
  FindDeltaILXL(t1,t2,t3,dyIL,dyXL,false);

  std::string dummy4;
  TraceHeader *t4 = NULL;
  contin = true;

  // Find a fourth point for control. Three points decides all values.
  // A fourth can confirm that everything is ok.
  while(contin== true && !file.eof())
  {
    ReadDummyTrace(file, binary_header_->GetFormat(), nz_);
    file.read(buffer, 240);
    dummy4.assign(buffer, 240);
    std::istringstream header4(dummy4);
    if(t4!=NULL)
      delete t4;
    t4 = new TraceHeader(*headerFormat);
    t4->Read(header4, binary_header_->GetLino());
    if (t4->GetCrossline() < 0 || t4->GetInline() < 0
      || t4->GetUtmx() < 0.0 || t4->GetUtmy() < 0.0)
    {
      delete t1;
      delete t2;
      delete t3;
      delete t4;
      return false;
    }
   // Check that t4 is not on line with any other line
    double s1 = (t1->GetUtmy()-t2->GetUtmy())/(t1->GetUtmx()-t2->GetUtmx());
    double s2 = (t2->GetUtmy()-t3->GetUtmy())/(t2->GetUtmx()-t3->GetUtmx());
    double s3 = (t1->GetUtmy()-t3->GetUtmy())/(t1->GetUtmx()-t3->GetUtmx());
    double s4 = (t2->GetUtmy()-t4->GetUtmy())/(t2->GetUtmx()-t4->GetUtmx());
    double s5 = (t3->GetUtmy()-t4->GetUtmy())/(t3->GetUtmx()-t4->GetUtmx());
    if(fabs(s1-s4)>0.01 && fabs(s2-s5)>0.01 && fabs(s3-s5)>0.01)
      contin = false;
  }
  double utmx4 = t1->GetUtmx() + (t4->GetCrossline()-t1->GetCrossline())*dxXL + (t4->GetInline()-t1->GetInline())*dxIL;
  if(fabs(utmx4-t4->GetUtmx())> 5.0)
  {
    delete t1;
    delete t2;
    delete t3;
    delete t4;
    return false;
  }
  double utmy4 = t1->GetUtmy() + (t4->GetCrossline()-t1->GetCrossline())*dyXL + (t4->GetInline()-t1->GetInline())*dyIL;
  if(fabs(utmy4-t4->GetUtmy())> 5.0)
  {
    delete t1;
    delete t2;
    delete t3;
    delete t4;
    return false;
  }
  delete t1;
  delete t2;
  delete t3;
  delete t4;

  return true;
}


//solves the equations x1 +(xl2-xl1)dxl + (il2-il1)dil = x2
//                     x1 + (xl3-xl1)dxl + (il3-il1)dil = x3
void SegY::FindDeltaILXL(TraceHeader *t1, TraceHeader *t2, TraceHeader *t3, double &dil, double &dxl, bool x)
{
  double x1, x2, x3;
  int il1, il2, il3, xl1, xl2, xl3;
  if(x==true)
  {
    x1 = t1->GetUtmx();
    x2 = t2->GetUtmx();
    x3 = t3->GetUtmx();
  }
  else
  {
    x1 = t1->GetUtmy();
    x2 = t2->GetUtmy();
    x3 = t3->GetUtmy();
  }
  il1 = t1->GetInline();
  il2 = t2->GetInline();
  il3 = t3->GetInline();
  xl1 = t1->GetCrossline();
  xl2 = t2->GetCrossline();
  xl3 = t3->GetCrossline();

  double teller = x3-x1;
  double nevner = ((xl3-xl1)*(x2-x1)-(il2-il1))/(xl2-xl1)+(il3-il1);
  dil = teller/nevner;
  dxl = ((x2-x1)-(il2-il1)*dil)/(xl2-xl1);
}
