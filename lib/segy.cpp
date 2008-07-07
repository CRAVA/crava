#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

#include "lib/global_def.h"
#include "lib/segy.h"
#include "lib/segylib.h"
#include "lib/lib_misc.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/simbox.h"
#include "src/fftgrid.h"

using namespace NRLib2;

SegY::SegY(char * fileName, const Simbox * simbox, float zPad, float z0)
{
  //  long int timestart, timeend;
  //  time(&timestart);
  simbox_  = simbox;
  zPad_    = static_cast<double>(zPad);
  error_   = 0;
  errMsg_  = NULL;
  outFile_ = NULL;
  unsigned long fSize = fileSize(fileName);
  FILE * file = fopen(fileName,"rb");
  char* junk=new char[3200];
  fread(junk, 1,3200, file); 
  char* binaryHeader= new char[400];
  BinaryHeader * sHead  = (BinaryHeader *)  binaryHeader;
  readBinaryHeader(file, binaryHeader);
  format_ = sHead->format;
  if(format_ == 3)
    datasize_ = 2;
  else
    datasize_ = 4;
  if(format_ == 4)
  {
    error_ = 1;
    errMsg_ = new char[120];
    LogKit::LogFormatted(LogKit::LOW,"failed\n");
    sprintf(errMsg_,
      "Can not read SegY-file \'%s\' that use floating point with gain representation.\n",
      fileName);
    return;
  }

  char* header= new char[240];
  SeisworksHeader * sTrace = (SeisworksHeader *) header;
  readHeader(file, header, binaryHeader);
  nz_ = sTrace->ns;
  dz_ = static_cast<double>(sTrace->dt/1000.0);
  x0_ = static_cast<double>(sTrace->sx);
  y0_ = static_cast<double>(sTrace->sy);
  if(z0_ == RMISSING) //Note: z0_ is common, and hence static.
    z0_ = static_cast<double>(z0);
  else
    assert(fabs(z0_-z0) < 0.01);
  inLine0_ = sTrace->fldr;
  crossLine0_ = sTrace->cdp;

  nx_ = ny_ = 0;
  int actL = 0, curL = 0;
  double x1 = 0.0 , y1 = 0.0 , x2 = 0.0 , y2 = 0.0 , prevx = 0.0 , prevy = 0.0 ;
  int nTraces = static_cast<int>(ceil((static_cast<double>(fSize)-3600.0)/static_cast<double>(datasize_*nz_+240.0)));
  //nTraces = 1161347;
  LogKit::LogFormatted(LogKit::LOW,"\nReading SEGY file %s.\n", fileName); 
  char*  buffer = new char[nz_*4];
  traces_ = new SegYTrace *[nTraces];
  traces_[0] = readTrace(file, buffer, x0_, y0_);
  yDir_ = 1;
  int i, flag = -1;
  int count = 0;
  int activeTraces = 0;
  int nHeaders = 0;
  int monitorSize = MAXIM(1,static_cast<int>(nTraces*0.02f));
  for(i=1;i<nTraces;i++)
  {
    nHeaders += readHeader(file, header, binaryHeader);
    traces_[i] = readTrace(file, buffer, static_cast<double>(sTrace->sx), static_cast<double>(sTrace->sy));
    if(traces_[i] != NULL)
      activeTraces++;
    if(nz_ != sTrace->ns || dz_ != sTrace->dt/1000.0)
    {
      error_ = 1;
      errMsg_ = new char[120];
      LogKit::LogFormatted(LogKit::LOW,"failed\n");
      sprintf(errMsg_,"Heading information about ns or dt changed in file %s\n",
        fileName);
      break;
    }
    if(flag == -1)
    {
      if(sTrace->cdp - crossLine0_ > 0)
      {
        flag = 0;
        curL = crossLine0_;
      }
      else 
      {
        flag = 1;
        curL = inLine0_;
      }
    }
    if(flag < 2)
    {
      count++;
      if(flag == 0)
        actL = sTrace->cdp;
      else
        actL = sTrace->fldr;
      if(curL < actL)
      {
        curL = actL;
        prevx = static_cast<double>(sTrace->sx);
        prevy = static_cast<double>(sTrace->sy);
      }
      else
      {
        if(nHeaders > 0)
          nTraces = count*(fSize/(3600+count*(datasize_*nz_+240)));
        LogKit::LogFormatted(LogKit::LOW,"%d traces (%d x %d), %d samples per trace:\n",
                         nTraces, count, nTraces/count, nz_);
        if(flag == 0)
          nx_ = count; //Associate x with crossline
        else
          nx_ = nTraces/count;
        ny_ = nTraces/nx_;
        x1 = prevx;
        y1 = prevy;
        flag += 2;
        printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
        printf("\n  |    |    |    |    |    |    |    |    |    |    |  ");
        printf("\n  ^");
      }
    }
    else if((i % count) == 0)
    {
      if(i+count == nTraces)
      {
        x2 = static_cast<double>(sTrace->sx);
        y2 = static_cast<double>(sTrace->sy);
      }
    }
    if ((i + 1) % monitorSize == 0) 
    { 
      printf("^");
      fflush(stdout);
    }
  }
  double deltaX = x1 - x0_;
  double deltaY = y1 - y0_;
  double xl = sqrt(deltaX*deltaX+deltaY*deltaY);
  dx_ = xl/static_cast<double>(nx_-1);
  cosRot_ = deltaX/xl;
  sinRot_ = deltaY/xl;
  deltaX = x2 - x0_;
  deltaY = y2 - y0_;
  double yl = sqrt(deltaX*deltaX+deltaY*deltaY);
  dy_ = yl/static_cast<double>(ny_-1);
  fclose(file);

  xlStep_ = (sTrace->cdp - crossLine0_)/(nx_-1);
  ilStep_ = (sTrace->fldr - inLine0_)/(ny_-1);
  if((sinRot_ > 0 && x2>x0_) || (sinRot_ < 0 && x2<x0_)) 
  {
    ilStep_ *= -1;
    yDir_ = -1;
    x0_ = x2;
    y0_ = y2;
  }

  if(simbox_->status() == Simbox::BOXOK)
  {
    double rot = acos(cosRot_);
    if(sinRot_< 0 )
      rot = -rot;
    if(simbox_->insideRectangle(x0_-0.5*(dx_*cosRot_-dy_*sinRot_),  
                                y0_-0.5*(dy_*cosRot_+dx_*sinRot_), 
                                rot,
                                xl+dx_, 
                                yl+dy_) == 0)
    {
      error_ = 1;
      errMsg_ = new char[MAX_STRING];
      sprintf(errMsg_,"Data in file %s does not cover simulation area.\n",
        fileName);
    }
  }
  delete [] buffer;
  delete [] junk;
  delete [] binaryHeader;
  delete [] header;
  //  time(&timeend);
  // LogKit::LogFormatted(LogKit::LOW,"SEGY read  in %ld seconds \n",timeend-timestart);
}

SegY::SegY(char * fileName, const Simbox * simbox)
{
  error_ = 0;
  errMsg_ = 0;
  traces_ = NULL;
  outFile_ = fopen(fileName,"wb");
  if(outFile_ == NULL)
  {
    errMsg_ = new char[MAX_STRING];
    error_ = 1;
    sprintf(errMsg_, "Could not open file %s for writing segy-file.\n", fileName);
  }
  else
  {
    if(z0_ == RMISSING)
      z0_ = 0;
    simbox_ = simbox;
    double min;
    double max;
    simbox_->getMinMaxZ(min, max);
    dz_ = floor(simbox_->getdz()+0.5);
    if(dz_ < 1.0)
      dz_ = 1.0;
    nz_ = static_cast<int>((max - z0_)/dz_)+1;
    writeMainHeader();
  }
}

SegY::~SegY()
{
  if(errMsg_ != NULL)
    delete [] errMsg_;
  if(outFile_ != 0)
    fclose(outFile_);
  int i;
  if(traces_ != NULL)
  {
    for(i=0;i<nx_*ny_;i++)
      if(traces_[i] != 0)
        delete traces_[i];
    delete [] traces_;
  }
}

float 
SegY::getValue(int i, int j, int k, int outsideMode)
{
  double x,y,z;
  simbox_->getCoord(i, j, k, x, y, z);
  float value;
  double sx =  (x-x0_)*cosRot_ + (y-y0_)*sinRot_ + 0.5*dx_;
  double sy = -(x-x0_)*sinRot_ + (y-y0_)*cosRot_ + 0.5*dy_;

  //printf("\n");
  //printf("x,x0_  y,y0_ = %.5f %.5f %.5f %.5f \n",x,x0_,y,y0_);
  //printf("(x-x0_)  (y-y0_) = %.5f %.5f \n",(x-x0_),(y-y0_));
  //printf("cosRot_ , sinRot_  = %.5f %.5f \n",cosRot_,sinRot_);
  //printf("(x-x0_)*cosRot_, (y-y0_)*sinRot_, 0.5*dx_  = %.5f %.5f %.5f\n",(x-x0_)*cosRot_,(y-y0_)*sinRot_,0.5*dx_);


  int tmpInd, index = static_cast<int>(floor(sx/dx_));
  if(yDir_ == 1)
    index += nx_*static_cast<int>(floor(sy/dy_));
  else
    index += nx_*(ny_-1-static_cast<int>(floor(sy/dy_)));
  if(z == RMISSING || sx < 0 || sx > nx_*dx_ || 
    sy < 0 || sy > ny_*dy_ || traces_[index] == 0)
    value = 0;
  else
  {
    int zind = static_cast<int>(floor((z-z0_)/dz_));
    float v1 = traces_[index]->getValue(zind);
    if(v1 == RMISSING && outsideMode == CLOSEST)
    {
      zind = traces_[index]->getLegalIndex(zind);
      v1 = traces_[index]->getValue(zind);
      if(traces_[index]->getValue(zind-1) == RMISSING)
        z = z0_+zind*dz_;          //Want edge value, hence 0/1 dz_ added
      else                         //(0.5 would give center of cell).
        z = z0_+(zind+0.99)*dz_;   
    }
    if(v1 != RMISSING)
    {
      //Computes interpolated value ax^2+by^2+cz^2+dx+ey+fz+g.
      //abcdefg estimated from closest point and its closest neighbours.
      int maxInd = nx_*ny_-1;
      float v0, v2, a, b, c, d, e, f, g;

      //Along x:
      v0 = RMISSING;
      v2 = RMISSING;
      if(index-1 >= 0 && traces_[index-1] != NULL)
        v0 = traces_[index-1]->getValue(zind);
      if(index+1 <= maxInd && traces_[index+1] != NULL)
        v2 = traces_[index+1]->getValue(zind);
      if(v0 == RMISSING)
      {
        a = 0;
        if(v2 == RMISSING)
          d = 0;
        else
          d = v2 - v1; //Using unit coordinates in each direction.
      }
      else if(v2 == RMISSING)
      {
        a = 0;
        d = v1 - v0;
      }
      else
      {
        a = (v2+v0-2*v1)/2.0f;
        d = (v2-v0)/2.0f;
      }
      //Along y:
      v0 = RMISSING;
      v2 = RMISSING;
      tmpInd = index-yDir_*nx_;
      if(tmpInd <= maxInd && tmpInd >= 0 && traces_[tmpInd] != NULL)
        v0 = traces_[tmpInd]->getValue(zind);
      tmpInd = index+yDir_*nx_;
      if(tmpInd <= maxInd && tmpInd >= 0 && traces_[tmpInd] != NULL)
        v2 = traces_[tmpInd]->getValue(zind);
      if(v0 == RMISSING)
      {
        b = 0;
        if(v2 == RMISSING)
          e = 0;
        else
          e = v2 - v1; //Using unit coordinates in each direction.
      }
      else if(v2 == RMISSING)
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
      if(v0 == RMISSING)
      {
        c = 0;
        if(v2 == RMISSING)
          f = 0;
        else
          f = v2 - v1; //Using unit coordinates in each direction.
      }
      else if(v2 == RMISSING)
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
      float ux = static_cast<float>(sx/dx_ - (floor(sx/dx_) + 0.5));
      float uy = static_cast<float>(sy/dy_ - (floor(sy/dy_) + 0.5));
      float uz = static_cast<float>((z-z0_)/dz_ - (floor((z-z0_)/dz_) + 0.5));
      value = a*ux*ux+b*uy*uy+c*uz*uz+d*ux+e*uy+f*uz+g;
    }
    else
    {
      if(outsideMode == ZERO)
        value = 0;
      else
        value = RMISSING;
    }
  }
  return(value);
}



void 
SegY::readBinaryHeader(FILE * file, char * buffer)
{
  fread(buffer, 1, 400, file);

#ifndef BIGENDIAN
  //The following fixes the byte order for 
  //inline, crossline, x, y, ns, dt on PC
  swap4Bytes(&(buffer[0]),3);
  swap2Bytes(&(buffer[12]),24);
#endif
}


int 
SegY::readHeader(FILE * file, char * buffer, char * b2)
{
  int i, newHead = 0;
  fread(buffer, 1, 240, file);
  if(buffer[0] == '√' && buffer[1] == '@' && buffer[2] == 'Ò' 
    && buffer[80] == '√' && buffer[160] == '√')
  {
    fread(buffer, 1, 160, file);
    for(i=0;i<7;i++)
      fread(b2, 1, 400, file);
    fread(b2, 1, 400, file);
    fread(buffer, 1, 240, file);
#ifndef BIGENDIAN
    //The following fixes the byte order for 
    //inline, crossline, x, y, ns, dt on PC
    swap4Bytes(&(b2[0]),3);
    swap2Bytes(&(b2[12]),24);
#endif
    newHead = 1;
  }

#ifndef BIGENDIAN
  //The following fixes the byte order for 
  //inline, crossline, x, y, ns, dt on PC
  swap4Bytes(&(buffer[8]),1);
  swap4Bytes(&(buffer[20]),1);
  swap2Bytes(&(buffer[70]),1);
  swap4Bytes(&(buffer[72]),2);
  swap2Bytes(&(buffer[114]),2);

  swap4Bytes(&(buffer[220]),1); //NBNB ENI

#endif
  if(((float *)&(buffer[8]))[0] == 0)
    for(i=0;i<4;i++)
      buffer[8+i] = b2[4+i];

#ifndef BYPASS_COORDINATE_SCALING
  double fac = 1;
  SeisworksHeader * sTrace = (SeisworksHeader *) buffer;
  switch(sTrace->scalco) {
    case -10000:
      fac = 0.0001;
      break;
    case -1000:
      fac = 0.001;
      break;
    case -100:
      fac = 0.01;
      break;
    case -10:
      fac = 0.1;
      break;
    case 10:
    case 100:
    case 1000:
    case 10000:
      fac = static_cast<double>(sTrace->scalco);
  }
  sTrace->sx = static_cast<int>(static_cast<double>(sTrace->sx)*fac);
  sTrace->sy = static_cast<int>(static_cast<double>(sTrace->sy)*fac);
#endif

  return(newHead);
}


SegYTrace *
SegY::readTrace(FILE * file, char * buffer, double x, double y)
{
  fread(buffer, datasize_, nz_, file);
  //Something like the following might be included, but be careful with padding.
  //  if(simbox_->isInside(x,y) == 0)
  //    return(NULL);
  double zTop = simbox_->getTop(x, y);
  if(zTop == RMISSING || zTop == WELLMISSING)
    //    zTop = z0_;
    return(NULL);
  double zBot = simbox_->getBot(x, y);
  if(zBot == RMISSING || zBot == WELLMISSING)
    //    zBot = nz_*dz_+z0_;
    return(NULL);

  if(static_cast<int>((zTop - z0_)/dz_) < 0) {
    LogKit::LogFormatted(LogKit::LOW,"\nERROR: A part of the top time surface reaches above the seismic gather. (seismic\n");
    LogKit::LogFormatted(LogKit::LOW,"       start time =%7.1f). Include more seismic data or lower the top surface.\n",z0_);
    LogKit::LogFormatted(LogKit::LOW,"       zTop = %7.2f\n");    
    exit(1);
  } 
  //NBNB-PAL: Tried nz_ - 1 below, but that failed  (incorrectly
  //NBBN-PAL: I think) in the facies case of test suite
  if(static_cast<int>((zBot - z0_)/dz_) + 1 > nz_) { 
    LogKit::LogFormatted(LogKit::LOW,"\nERROR: A part of the base time surface reaches below the seismic gather. Include\n");
    LogKit::LogFormatted(LogKit::LOW,"       more seismic data or highten the base surface.\n");
    exit(1);
  } 

  double pad = 0.5*zPad_*(zBot - zTop);
  int j0 = static_cast<int>((zTop - pad - z0_)/dz_);
  int j1 = static_cast<int>((zBot + pad - z0_)/dz_) + 1; //Add 1 since int truncates.

  if(j0 < 0)
    j0 = 0;
  if(j1 > nz_-1)
    j1 = nz_-1;
  return(new SegYTrace((float *) buffer, j0, j1, format_));
}

void
SegY::writeMainHeader()
{
  assert(outFile_);

  char * ebcdicH = new char[3200];
  ebcdicHeader(ebcdicH);  // writes values to eHead
  ebcdicH[3199]='@';
  fwrite(ebcdicH, 1, 3200, outFile_);
  delete [] ebcdicH;
  char * binaryHeader = new char[400];
  memset(binaryHeader,0,400); // initialize buffer as zero
  BinaryHeader * sHead  = (BinaryHeader *)  binaryHeader;
  sHead->jobid  = 9999; // dummy
  sHead->lino   = 1; // dummy

  sHead->reno   = 1; // dummy
  sHead->ntrpr  = short(simbox_->getnx()); 
  sHead->nart   = 0; // dummy
  sHead->hdt    = short(dz_*1000); 
  sHead->hns    = short(nz_);
  sHead->format = 1; // traces in float format
  sHead->fold   = 1; // dummy
  sHead->tsort  = 4; // dummy
  sHead->mfeet  = 1; // units meter 


#ifndef BIGENDIAN
  swap4Bytes(&(binaryHeader[0]),3);
  swap2Bytes(&(binaryHeader[12]),24);
#endif
  fwrite(binaryHeader, 1, 400, outFile_);
  delete [] binaryHeader;
}

void
SegY::writeTrace(int i, int j, FFTGrid * grid)
{
  char * buffer = new char[240];
  memset(buffer,0,240); // initialize buffer as zero
  assert(outFile_);
  SeisworksHeader * sTrace = (SeisworksHeader *) buffer;
  sTrace->ns = static_cast<unsigned short>(nz_);
  sTrace->dt = static_cast<unsigned short>(dz_*1000);

  double x, y, z;
  simbox_->getCoord(i, j, 0, x, y, z); 
  sTrace->sx = static_cast<int>(x); 
  sTrace->sy = static_cast<int>(y); 

  // begin x-y change (for float definition )
  //float * xypos= new float[2];
  //float z;
  //int* xyposInt= (int*)xypos;
  //simbox_->getCoord(i,j,0,xypos[0],xypos[1],z);
  //pb_Ieee2Ibm(xypos, xypos, 2);
  //sTrace->sx = xyposInt[0];
  //sTrace->sy = xyposInt[1];
  // end x-y change

  sTrace->fldr   = simbox_->getIL0()+j*simbox_->getILStep(); //inline
  sTrace->cdp    = simbox_->getXL0()+i*simbox_->getXLStep(); //crossline
  sTrace->counit = 1;   // units in length (meters or feet)
  sTrace->tatyp  = 3;   // other taper type
  // dummyvalues are assiged to the parameters below
  sTrace->tracl  = 1;
  sTrace->cdpt   = 1;
  sTrace->tracr  = ((i+1)*sTrace->cdpt);
  sTrace->tracf  = j+1;
  sTrace->trid   = 3;
  sTrace->scalel = 1;
  sTrace->scalco = 0;
  //sTrace->laga   = 0;
  //  pb_Ieee2Ibm(&(buffer[72]), (float*) &(buffer[72]), 2);
#ifndef BIGENDIAN
  //The following fixes the byte order for 
  //inline, crossline, x, y, ns, dt on PC
  swap4Bytes(&(buffer[8]),1);
  swap4Bytes(&(buffer[20]),1);
  swap2Bytes(&(buffer[70]),1);
  swap4Bytes(&(buffer[72]),2);
  swap2Bytes(&(buffer[114]),2);
  swap4Bytes(&(buffer[220]),1); //ENI

#endif
  //printf("Writing trace %3d,%3d:\n",i,j);
  int nWritten = fwrite(buffer, 1, 240, outFile_);
  //printf("H: %d\n ",nWritten);
  delete [] buffer;

  float * trace = new float[nz_];
  int k;
  if(z == RMISSING || z == WELLMISSING)
  {
    //printf("Missing trace.\n");
    for(k=0;k<nz_;k++)
      trace[k] = 0;
  }
  else
  {
    int firstData = int((z-z0_)/dz_);
    //printf("Generating trace %d %d %d.\n", firstData, simbox_->getnz(),nz_);
    for(k=0;k<firstData;k++)
      trace[k] = 0; //data[0];
    // NBNB-PAL : change div div to mult
    double fac = dz_/(simbox_->getdz()*simbox_->getRelThick(i,j));
    int endData = firstData + static_cast<int>(simbox_->getnz()/fac);

    if(endData > nz_)
    {

      //printf("\ni,j = %d %d   firstData  %d    simbox_->getnz() fac %d %.5f    endData = %d    nz_ = %d\n",
      //       i,j,firstData,simbox_->getnz(),fac,endData,nz_);
      //printf(" dz_  simbox_->getdz() simbox_->getRelThick(i,j)    %.5f %.5f  %.5f\n",
      //       dz_,simbox_->getdz(),simbox_->getRelThick(i,j)) ;
      //exit(1);
      

      printf("Internal warning: SEGY-grid too small (%d, %d needed). Truncating data.\n", nz_, endData); 
      endData = nz_;
    }
    for(k=firstData;k<endData;k++)
      trace[k] = grid->getRealValue(i,j,static_cast<int>((k-firstData)*fac+0.5));
    for(k=endData;k<nz_;k++)
      trace[k] = 0; //data[simbox_->getnz()-1];
    //printf("Trace generated.\n");
  }

  buffer = (char *) trace;
  pb_Ieee2Ibm(buffer, (float*) buffer, nz_);
#ifndef BIGENDIAN //Rearrange bytes for PC
  swap4Bytes(buffer, nz_);
#endif
  nWritten = fwrite(buffer, 4, nz_, outFile_);
  //printf("T: %d\n",nWritten);  
  delete [] trace;
}


void
SegY::writeTrace(int i, int j, float * data)
{
  char * buffer = new char[240];
  memset(buffer,0,240); // initialize buffer as zero
  assert(outFile_);
  SeisworksHeader * sTrace = (SeisworksHeader *) buffer;
  sTrace->ns = (unsigned short) nz_;
  sTrace->dt = (unsigned short) dz_*1000;

  double x, y, z;
  simbox_->getCoord(i, j, 0, x, y, z); 
  sTrace->sx = static_cast<int>(x); 
  sTrace->sy = static_cast<int>(y); 

  // begin x-y change (for float definition )
  //float * xypos= new float[2];
  //float z;
  //int* xyposInt= (int*)xypos;
  //simbox_->getCoord(i,j,0,xypos[0],xypos[1],z);
  //pb_Ieee2Ibm(xypos, xypos, 2);
  //sTrace->sx = xyposInt[0];
  //sTrace->sy = xyposInt[1];
  // end x-y change

  sTrace->fldr   = simbox_->getIL0()+j*simbox_->getILStep(); //inline
  sTrace->cdp    = simbox_->getXL0()+i*simbox_->getXLStep(); //crossline
  sTrace->counit = 1;   // units in length (meters or feet)
  sTrace->tatyp  = 3;   // other taper type
  // dummyvalues are assiged to the parameters below
  sTrace->tracl  = 1;
  sTrace->cdpt   = 1;
  sTrace->tracr  = ((i+1)*sTrace->cdpt);
  sTrace->tracf  = j+1;
  sTrace->trid   = 3;
  sTrace->scalel = 1;
  sTrace->scalco = 0;
  //sTrace->laga   = 0;
  //  pb_Ieee2Ibm(&(buffer[72]), (float*) &(buffer[72]), 2);
#ifndef BIGENDIAN
  //The following fixes the byte order for 
  //inline, crossline, x, y, ns, dt on PC
  swap4Bytes(&(buffer[8]),1);
  swap4Bytes(&(buffer[20]),1);
  swap2Bytes(&(buffer[70]),1);
  swap4Bytes(&(buffer[72]),2);
  swap2Bytes(&(buffer[114]),2);

  swap4Bytes(&(buffer[220]),1); //NBNB ENI

#endif
  //printf("Writing trace %3d,%3d:\n",i,j);
  int nWritten = fwrite(buffer, 1, 240, outFile_);
  //printf("H: %d\n ",nWritten);
  delete [] buffer;

  float * trace = new float[nz_];
  int k;
  if(z == RMISSING || z == WELLMISSING)
  {
    //printf("Missing trace.\n");
    for(k=0;k<nz_;k++)
      trace[k] = 0;
  }
  else
  {
    int firstData = int((z-z0_)/dz_);
    //printf("Generating trace %d %d %d.\n", firstData, simbox_->getnz(),nz_);
    for(k=0;k<firstData;k++)
      trace[k] = 0; //data[0];
    for(k=firstData;k<firstData+simbox_->getnz();k++)
      trace[k] = data[k-firstData];
    for(k=firstData+simbox_->getnz();k<nz_;k++)
      trace[k] = 0; //data[simbox_->getnz()-1];
    //printf("Trace generated.\n");
  }

  buffer = (char *) trace;
  pb_Ieee2Ibm(buffer, (float*) buffer, nz_);
#ifndef BIGENDIAN //Rearrange bytes for PC
  swap4Bytes(buffer, nz_);
#endif
  nWritten = fwrite(buffer, 4, nz_, outFile_);
  //printf("T: %d\n",nWritten);  
  delete [] trace;
}

long
SegY::fileSize(char * fileName)
{
  long length;
  FILE * file = fopen(fileName,"rb");
  fseek(file, 0, SEEK_END);
  length = ftell(file);
  fclose(file);

  return length;
}


int
SegY::completeTimeSimbox(Simbox * simbox, double lzLimit, char * errText)
{
  double rot = acos(cosRot_);
  if(sinRot_ < 0)
    rot = -rot;

  double ax0 = x0_ - 0.5*(dx_*cosRot_-dy_*sinRot_);
  double ay0 = y0_ - 0.5*(dy_*cosRot_+dx_*sinRot_);
  double lx  = nx_*dx_; 
  double ly  = ny_*dy_; 

  simbox->setArea(ax0, ay0, lx, ly, rot, dx_, dy_);
  simbox->setSeisLines(inLine0_, crossLine0_, ilStep_, xlStep_);
  int error = simbox->checkError(lzLimit, errText);

  LogKit::LogFormatted(LogKit::LOW,"\n\nExtracting inversion grid information from seismic volume:\n");
  LogKit::LogFormatted(LogKit::LOW,"\n                  x0          y0          lx        ly       dx     dy      angle");
  LogKit::LogFormatted(LogKit::LOW,"\n---------------------------------------------------------------------------------");
  LogKit::LogFormatted(LogKit::LOW,"\nArea:    %11.2f %11.2f   %9.2f %9.2f   %6.2f %6.2f   %8.3f  ", ax0, ay0, lx, ly, dx_, dy_, 180.0*rot/PI);
  LogKit::LogFormatted(LogKit::LOW,"\nSeismic: %11.2f %11.2f   %9.2f %9.2f   %6.2f %6.2f   %8.3f\n", x0_, y0_, lx, ly, dx_, dy_, 180.0*rot/PI);
  return(error);
}

void
SegY::completeDepthSimbox(Simbox * depthSimbox)
{
  double rot = acos(cosRot_);
  if(sinRot_ < 0)
    rot = -rot;

  double ax0 = x0_ - 0.5*(dx_*cosRot_-dy_*sinRot_);
  double ay0 = y0_ - 0.5*(dy_*cosRot_+dx_*sinRot_);
  double lx  = nx_*dx_; 
  double ly  = ny_*dy_; 

  depthSimbox->setArea(ax0, ay0, lx, ly, rot, dx_, dy_);
}

void SegY::ebcdicHeader(char * ebcdicH)
{
  sprintf(&(ebcdicH[0]),  "√@Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[80]), "√@Ú@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[160]),"√@Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[240]),"√@Ù@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[320]),"√@ı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[400]),"√@ˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[480]),"√@˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[560]),"√@¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[640]),"√@˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[720]),"√Ò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[800]),"√ÒÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[880]),"√ÒÚ@@@@@@@@@„»…‚@‚≈«Ë@∆…”≈@Ê¡‚@ÊŸ…„„≈’@¬Ë@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[960]),"√ÒÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1040]),"√ÒÙ@@@@@@@@@„»≈@√Ÿ¡Â¡@◊Ÿ÷«Ÿ¡‘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1120]),"√Òı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1200]),"√Òˆ@@@@@@@@@ƒ≈Â≈”÷◊≈ƒ@¬Ë@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1280]),"√Ò˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1360]),"√Ò¯@@@@@@@@@Ÿ¡«’¡Ÿ@»¡‰«≈@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1440]),"√Ò˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1520]),"√Ú@@@@@@@@@¡’’≈@Ÿ¡’ƒ…@‚ËÂ≈Ÿ‚Â≈≈’@@@@@@@@@@@¡’ƒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1600]),"√ÚÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1680]),"√ÚÚ@@@@@@@@@÷ƒƒ@“÷”¬—÷Ÿ’‚≈’@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1760]),"√ÚÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1840]),"√ÚÙ@@@@@@@@@¡„@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[1920]),"√Úı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2000]),"√Úˆ@@@@@@@@@„»≈@’÷ŸÊ≈«…¡’@√÷‘◊‰„…’«@√≈’„≈Ÿ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2080]),"√Ú˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2160]),"√Ú¯@@@@@@@@@ÚÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2240]),"√Ú˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2320]),"√Û@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2400]),"√ÛÒ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2480]),"√ÛÚ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2560]),"√ÛÛ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2640]),"√ÛÙ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2720]),"√Ûı@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2800]),"√Ûˆ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2880]),"√Û˜@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[2960]),"√Û¯@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[3040]),"√Û˘@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  sprintf(&(ebcdicH[3120]),"√Ù@≈’ƒ@≈¬√ƒ…√@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
  ebcdicH[3199]='@';
}


double SegY::z0_ = RMISSING;



SegYTrace::SegYTrace(float * data, int jStart, int jEnd, int format)
{
  int i, nData = jEnd - jStart + 1;  
  if(nData < 0)
    LogKit::LogFormatted(LogKit::LOW,"Oooops...\n");
  data_ = new float[nData];
  if(data_ == NULL)
  {
    LogKit::LogFormatted(LogKit::LOW,"Ran out of memory while reading segy file.\n");
    exit(1);
  }
  if(format < 3)
  {
    for(i=jStart;i<=jEnd;i++)
      data_[i-jStart] = data[i];
    jStart_ = jStart;
    jEnd_ = jEnd;
    char * buffer = (char *) data_;
#ifndef BIGENDIAN //Rearrange bytes for PC
    swap4Bytes(buffer, nData);
#endif
    if(format == 1)
      pb_Ibm2Ieee(buffer, (float*) buffer, nData);
  }
  else
  {
    short int * siData = (short int *) data;
    char * buffer = (char *) &(siData[jStart]);
#ifndef BIGENDIAN //Rearrange bytes for PC
    swap2Bytes(buffer, nData);
#endif
    for(i=jStart;i<=jEnd;i++)
      data_[i-jStart] = (float) siData[i];
    jStart_ = jStart;
    jEnd_ = jEnd;
  }
  //for (int k = 0 ; k < nData ; k++) {
  //  printf("k = %d      data_[k] = %.5f\n",k,data_[k]);
  //} 
}


SegYTrace::~SegYTrace()
{
  delete [] data_;
}

float
SegYTrace::getValue(int j)
{
  float value;
  if(j < jStart_ || j > jEnd_)
    value = RMISSING;
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









