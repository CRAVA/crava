#ifndef SEGY_HPP
#define SEGY_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "traceheader.hpp"
#include "commonheaders.hpp"
#include "../volume/volume.hpp"


const int segyIMISSING = -99999;
class SegYTrace;
class SegyGeometry;
class BinaryHeader;
class TextualHeader;



class SegY{
public:
/**
  Constructor for reading 
  Read only the headers on top of the file
  @param[in] fileName  Name of file to read data from
  @param[in] z0
  */
  SegY(const std::string       & fileName, 
       float                     z0, 
       const TraceHeaderFormat & traceHeaderFormat = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)); 
/**
  Constructor for writing
  @param[in] fileName  Name of file to write data to
  */
 SegY(const std::string       & fileName, 
      float                     z0, 
      int                       nz, 
      float                     dz, 
      const TextualHeader     & ebcdicHeader,
      const TraceHeaderFormat & traceHeaderFormat = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));
  
  ~SegY();

  //>>>Begin read all traces mode
  void                      readAllTraces(NRLib2::Volume * volume, 
                                          double           zPad, 
                                          bool             onlyVolume = false); ///< Read all traces with header
  float                     getValue(double x, 
                                     double y, 
                                     double z, 
                                     int    outsideMode = segyIMISSING);
  
  std::vector<float>        getAllValues(); ///< Return vector with all values. 
  void                      createRegularGrid();
  const SegyGeometry      * getGeometry(void)  { return geometry_ ;} //Only makes sense after command above.
  //<<<End read all trace mode

  //>>>Begin read single trace mode
  const SegYTrace         * getNextTrace(double           zPad = 0, 
                                         NRLib2::Volume * volume = NULL, 
                                         bool             onlyVolume = false);
  //<<<End read single trace mode

  //>>>Begin write mode
  void                      setGeometry(const SegyGeometry * geometry);  
  void                      storeTrace(float               x, 
                                       float               y, 
                                       std::vector<float>  data,
                                       NRLib2::Volume    * volume,
                                       float               topVal=0.0f,
                                       float               baseVal=0.0f);
  void                      writeTrace(TraceHeader       * traceHeader, 
                                       std::vector<float>  data,
                                       NRLib2::Volume    * volume,
                                       float               topVal=0.0f,
                                       float               baseVal=0.0f);    ///< Write single trace to file
  void                      writeTrace(float                  x, 
                                       float                  y, 
                                       std::vector<float>     data,
                                       const NRLib2::Volume * volume,
                                       float                  topVal=0.0f,
                                       float                  baseVal=0.0f); ///< Write single trace to internal memory
  void                      WriteAllTracesToFile(); ///< Use only after writeTrace with x and y as input is used for the whole cube
  //<<<End write mode


  // int checkError(char * errText) 
  //   {if(error_ > 0) strcpy(errText, errMsg_);return(error_);}
  /// Return (possibly upper limit for) number of traces
  
  int                       getNTraces() { return nTraces_ ;} 
  int                       getNz()      { return nz_      ;}
  float                     getDz()      { return dz_      ;}

  enum                      outsideModes{MISSING, ZERO, CLOSEST};

  int                       findNumberOfTraces(void);
  static int                findNumberOfTraces(const std::string       & fileName, 
                                               const TraceHeaderFormat & traceHeaderFormat 
                                               = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));

  SegyGeometry            * findGridGeometry(void);                                 
  static SegyGeometry     * findGridGeometry(const std::string       & fileName, 
                                             const TraceHeaderFormat & traceHeaderFormat 
                                             = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));
private:
  void                      ebcdicHeader(std::string& outstring);               ///<
  bool                      readHeader(TraceHeader * header);                   ///< Trace header
  SegYTrace               * readTrace(NRLib2::Volume * volume, 
                                      double           zPad, 
                                      bool           & duplicateHeader, 
                                      bool             onlyVolume,
                                      bool           & outsideSurface);     ///< Read single trace from file
  void                      writeMainHeader(const TextualHeader& ebcdicHeader); ///< Quasi-dummy at the moment. 
  void                      readDummyTrace(std::fstream & file, int format, int nz);
  double                    interpolate(int xyidx, 
                                        int zidx);
  std::ios::pos_type        findFileSize(const std::string & fileName);

  TraceHeaderFormat         traceHeaderFormat_;

  SegyGeometry            * geometry_;             ///< Parameters to find final index from i and j
  BinaryHeader            * binaryHeader_;         ///<
  
  bool                      singleTrace_;          ///< Read one and one trace
  bool                      simboxOnly_;           ///<
  bool                      checkSimbox_;          ///<

  std::vector<SegYTrace*>   traces_;               ///< All traces
  int                       nTraces_;              ///< WARNING: Counts many different things

  int                       datasize_;             ///< Bytes per datapoint in file.         

  int                       nz_; 

  float                     z0_;                   ///< Top of segy cube
  float                     dz_;                   ///< Sampling density
                                                                            
  std::fstream              file_;                                 
  std::string               fileName_;              
                                                           
  float                     rmissing_;                    
};


class SegYTrace{
public:
  /**
  Constructor
  @param[in] data           Data in trace
  @param[in] jStart         Start index in trace
  @param[in] jEnd           End index in trace
  @param[in] format
  @param[in] x              x coordinate of trace
  @param[in] y              y coordinate of trace
  @param[in] inLine         inline number
  @param[in] crossLine      crossline number
 
  */
  SegYTrace(std::fstream & file, int jStart, int jEnd, int format, float x, float y, int inLine, int crossLine, int nz);
  SegYTrace(std::vector<float> indata, int jStart, int jEnd, float x, float y, int inLine, int crossLine);
  ~SegYTrace();
  /// get trace value at index j
  float getValue(int j) const;
  /// 
  int getLegalIndex(int index);
  /// Get start index
  int getStart() {return(jStart_);}
  /// Get end index
  int getEnd() {return(jEnd_);}
  /// Get x coordinate
  float getX() {return(x_);}
  /// Get y coordinate
  float getY() {return(y_);}
  /// Get inline number
  int getInline() {return(inLine_);}
  /// Get crossline number
  int getCrossline() {return(crossLine_);}

private:
  /// Data in trace
  std::vector<float> data_;
  /// Start index
  int jStart_;
  /// End index
  int jEnd_;
  /// x coord
  float x_;
  /// y coord
  float y_;
  /// inline index
  int inLine_;
  /// crossline index
  int crossLine_;
float rmissing_;
};

class SegyGeometry{
public:
  /**
  Constructor
  @param[in] traces  Traces from segY cube
  */
  SegyGeometry(std::vector<SegYTrace *> &traces);
  SegyGeometry(double x0,double y0,double dx,double dy,
               int nx,int ny,int IL0,int XL0,int ilStep,int xlStep,bool ILxflag,double rot);
  SegyGeometry(const SegyGeometry *geometry);   ///< Copy constructor
  ~SegyGeometry();
  
  int    returnIndex(float x, float y, float &xind, float  &yind); ///< Return grid index for i and j
  int    returnIndex(int IL, int XL, int &i, int &j);
  int    returnILXL(int &IL, int &XL, float x, float y);
  int    getNx()           const { return nx_        ;}            ///< return nx
  int    getNy()           const { return ny_        ;}            ///< return ny
  bool   getILxflag()      const { return ILxflag_   ;}
  // Fyll på med get-funksjonar. etter behov.
  double getDx()           const { return dx_        ;}
  double getDy()           const { return dy_        ;}
  double getX0()           const { return x0_        ;}
  double getY0()           const { return y0_        ;}
  double getlx()           const { return nx_*dx_    ;}
  double getly()           const { return ny_*dy_    ;}
  double getAngle()        const { return rot_       ;}
  double getCosRot()       const { return cosRot_    ;}
  double getSinRot()       const { return sinRot_    ;}
  int    getInLine0()      const { return inLine0_   ;}
  int    getCrossLine0()   const { return crossLine0_;}
  int    getILstep()       const { return ilStep_    ;}
  int    getXLstep()       const { return xlStep_    ;}
  void   findIJFromILXL(int IL, int XL, int &i, int &j);
  void   findIJfromXY(float x, float y, int &i, int &j);
  void   findXYfromIJ(int i, int j, float &x, float &y) const;
  void   writeGeometry() const;

private:
  double x0_, y0_;              ///<
  double dx_, dy_;              ///< Cell increments
  long nx_, ny_;                ///< Grid dimensions  
  int inLine0_, crossLine0_;    ///< Start value for inline and crossline
  int ilStep_, xlStep_;         ///< 
  bool ILxflag_;                ///<
  double sinRot_, cosRot_;      ///<
  double rot_;                  ///<
};

#endif
