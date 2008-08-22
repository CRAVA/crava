// $Id$ 

#include "traceheader.hpp"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>

#include "../iotools/fileio.hpp"
const float RMISSING = -99999.000;
const int IMISSING = -99999;

using namespace std;
using namespace NRLib2;

enum standardLoc {
  SCALCO_LOC    = 71,
  SX_LOC        = 73,
  SY_LOC        = 77,
  NS_LOC        = 115,
  DT_LOC        = 117,
  INLINE_LOC    = 9, // 189,
  CROSSLINE_LOC = 21 // 193
};

TraceHeaderFormat::TraceHeaderFormat(int headerformat)
{
  if(headerformat==SEISWORKS)
  {
   scalCoLoc_ = SCALCO_LOC;
    utmxLoc_ = SX_LOC;
    utmyLoc_ = SY_LOC;
    inlineLoc_ = INLINE_LOC;
    crosslineLoc_ = CROSSLINE_LOC;
    coordSys_ = UTM;
  }
  else if(headerformat==IESX)
  {
    scalCoLoc_ = SCALCO_LOC;
    utmxLoc_ = SX_LOC;
    utmyLoc_ = SY_LOC;
    inlineLoc_ = 221;
    crosslineLoc_ = CROSSLINE_LOC;
    coordSys_ = UTM;


  }
}

TraceHeaderFormat::TraceHeaderFormat(int scalCoLoc,
                                     int utmxLoc,
                                     int utmyLoc,
                                     int inlineLoc,
                                     int crosslineLoc,
                                     coordSys_t coordSys)
  : scalCoLoc_(scalCoLoc),
    utmxLoc_(utmxLoc),
    utmyLoc_(utmyLoc),
    inlineLoc_(inlineLoc),
    crosslineLoc_(crosslineLoc),
    coordSys_(coordSys)
{}



const char* TraceHeaderFormat::toString() const
{
  static char output[80];
  sprintf(output, "Coord used: %s Location in trace header: "
    "UTM-X: %d UTM-Y: %d IL: %d, XL: %d",
    (coordSys_ == UTM ? "UTM" : "IL/XL"), 
    utmxLoc_, utmyLoc_, inlineLoc_, crosslineLoc_);
  return output;
}

void TraceHeaderFormat::bypassCoordinateScaling()
{
  scalCoLoc_ = -1;
}
TraceHeader::TraceHeader(const TraceHeaderFormat& format)
  : format_(format),
    status_(0),
    scalCo_(1)
{
  memset(buffer_, 0, 240);
  rmissing_ = RMISSING;
  imissing_ = IMISSING;
}


void TraceHeader::read(std::istream& inFile, int lineNo)
{
 // if (std::fread(buffer_, 1, 240, inFile) != 240) {
  //  status_ = -2;
  //  return;
 // }
  if(!(inFile.read(buffer_,240))) {
    status_ = -2;
    return;
  }
  
  if(buffer_[0] == 'Ã' && buffer_[1] == '@' && buffer_[2] == 'ñ' 
     && buffer_[80] == 'Ã' && buffer_[160] == 'Ã')
  {
    // This is not a trace header, but the start of an EDBDIC-header.
    // Set file pointer at end of EDBDIC header.
   // fseek(inFile, 2960, SEEK_CUR);
    inFile.seekg(2960);
    status_ = -1;
    return;
  }

  inFile.seekg(-240,ios_base::cur);

  short dummy;
  int i = 0;
  while(i<240) 
  {   
    if(i==(format_.getScalCoLoc()-1))
    {
      scalcoinitial_ = ReadBinaryShort(inFile);
      i=i+2;
    }
    else if(i==(format_.getUtmxLoc()-1))
    {
      utmx_ = float(ReadBinaryInt(inFile));
      i=i+4;
    }
    else if(i==(format_.getUtmyLoc()-1))
    {
      utmy_ = float(ReadBinaryInt(inFile));
      i=i+4;
    }
    else if(i==(NS_LOC-1))
    {
      ns_ = ReadBinaryShort(inFile);
      i=i+2;
    }
    else if(i==(DT_LOC-1))
    {
      dt_ = ReadBinaryShort(inFile);
      i=i+2;
    }
    else if(i==(format_.getInlineLoc()-1))
    {
      inline_ = ReadBinaryInt(inFile);
      i=i+4;
    }
    else if(i==(format_.getCrosslineLoc()-1))
    {
      crossline_ = ReadBinaryInt(inFile);
      i=i+4;
    }
    else
    {
      dummy = ReadBinaryShort(inFile);
      i=i+2;
    }

  }


 // swapBuffer();

  if (format_.getScalCoLoc() > 0) {
    switch (scalcoinitial_) {
      case -10000: 
        scalCo_ = 0.0001f;
        break;
      case -1000:
        scalCo_ = 0.001f;
        break;
      case -100:
        scalCo_ = 0.01f;
        break;
      case -10:
        scalCo_ = 0.1f;
        break;
      case 10:
      case 100:
      case 1000:
      case 10000:
        scalCo_ = static_cast<float>(scalcoinitial_);
    }
  }
useBinaryInline = false;
  // if (format.useInlineFromBinHeader()) {
  if (lineNo > 0 && getInline() == 0 && useBinaryInline==true) { 
    setInline(lineNo);
  }
}

int TraceHeader::write(std::ostream& outFile) 
{
  int errCode = 0;

 // write on correct locations. What to write between?
  int i = 0;

  while(i<240) 
  {   
    if(i==(format_.getScalCoLoc()-1))
    {
      WriteBinaryShort(outFile, scalcoinitial_);
      i = i+2;
    }
    else if(i==(format_.getUtmxLoc()-1))
    {
      WriteBinaryInt(outFile, int(utmx_));
      i=i+4;
    }
    else if(i==(format_.getUtmyLoc()-1))
    {
      WriteBinaryInt(outFile, int(utmy_));
      i=i+4;
    }
    else if(i==(NS_LOC-1))
    {
      WriteBinaryShort(outFile, ns_);
      i=i+2;
    }
    else if(i==(DT_LOC-1))
    {
      WriteBinaryShort(outFile, dt_);
      i=i+2;
    }
    else if(i==(format_.getInlineLoc()-1))
    {
      WriteBinaryInt(outFile, inline_);
      i=i+4;
    }
    else if(i==(format_.getCrosslineLoc()-1))
    {
      WriteBinaryInt(outFile, crossline_);
      i=i+4;
    }
    else
    {
      outFile.write(&(buffer_[i]), 2);
      i=i+2;
    }
  }
  

  //swapBuffer();
  //int nWritten = fwrite(buffer_, 1, 240, outFile);
  //if (nWritten != 240) {
  //  errCode = -1;
  //}
  // Swap back to host byte order.
 // swapBuffer();

  return errCode;
}


float TraceHeader::getUtmx() const 
{
  int loc = format_.getUtmxLoc();
  if (loc < 0) {
    return rmissing_;
  }
 // return scalCo_ * float(getInt32(loc));
  return(scalCo_*utmx_);
}

void TraceHeader::setUtmx(float utmx) 
{
  int loc = format_.getUtmxLoc();
  if (loc > 0) {
  //  setInt32(loc, int(utmx/scalCo_));
    utmx_ = utmx;
  }
}

float TraceHeader::getUtmy() const 
{
  int loc = format_.getUtmyLoc();
  if (loc < 0) {
    return RMISSING;
  }
 // return scalCo_ * float(getInt32(loc));
  return(scalCo_*utmy_);
}

void TraceHeader::setUtmy(float utmy) 
{
  int loc = format_.getUtmyLoc();
  if (loc > 0) {
 //   setInt32(loc, int(utmy/scalCo_));
    utmy_ = utmy;
  }
}


int TraceHeader::getInline() const 
{
  int loc = format_.getInlineLoc();
  if (loc < 0) {
    return imissing_;
  }
//  return getInt32(loc);
  return inline_;
}


void TraceHeader::setInline(int inLine) 
{
  int loc = format_.getInlineLoc();  
  if (loc > 0) {
  //  setInt32(loc, inLine);
    inline_ = inLine;
  }
}


int TraceHeader::getCrossline() const 
{
  int loc = format_.getCrosslineLoc();
  if (loc < 0) {
    return IMISSING;
  }
//  return getInt32(loc);
  return crossline_;
}


void TraceHeader::setCrossline(int crossLine) 
{
  int loc = format_.getCrosslineLoc();  
  if (loc > 0) {
  //  setInt32(loc, crossLine);
    crossline_ = crossLine;
  }
}


void TraceHeader::setNSamples(int ns) 
{
 // setInt16(NS_LOC, short(ns));
  ns_ = short(ns);
}


void TraceHeader::setDt(int dt) 
{
//  setInt16(DT_LOC, short(dt));
  dt_ = short(dt);
}


short TraceHeader::getScalCo() const 
{
  int loc = format_.getScalCoLoc();
  if (loc < 0) {
    return 0;
  }
 // return getInt16(loc);
  return scalcoinitial_;
}


int TraceHeader::getInt32(int loc) const
{
  return *reinterpret_cast<const int*>(&buffer_[loc - 1]);  
}


void TraceHeader::setInt32(int loc, int val) 
{
  int* buff = reinterpret_cast<int*>(&buffer_[loc - 1]);
  *buff = val;
}


short TraceHeader::getInt16(int loc) const
{
  return *reinterpret_cast<const short*>(&buffer_[loc - 1]);
}


void TraceHeader::setInt16(int loc, short val) 
{
  short* buff = reinterpret_cast<short*>(&buffer_[loc - 1]);
  *buff = val;
}


