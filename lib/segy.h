#ifndef SEGY_H
#define SEGY_H

#include <stdio.h>
#include <string.h>

class Simbox;
class SegYTrace;
class FFTGrid;

class SegY{
public:

  SegY(char * fileName, const Simbox * simbox, float zPad, float z0); //For reading
  SegY(char * fileName, const Simbox * simbox);                       //For writing
  ~SegY();

  float          getValue(int i, int j, int k, int outsideMode = MISSING);
  int            checkError(char * errText) {if(error_ > 0) strcpy(errText, errMsg_);return(error_);}
  void           writeTrace(int i, int j, float * data);
  void           writeTrace(int i, int j, FFTGrid * grid);
  void           getAreaParameters(double * areaParams);
  void           getSeisLines(int * seisLines);
  enum           outsideModes{MISSING, ZERO, CLOSEST};

private:
  int            readHeader(FILE * file, char * buffer, char * b2);
  void           readBinaryHeader( FILE * file, char * buffer);
  SegYTrace    * readTrace(FILE * file, char * buffer, double x, double y);

  void           writeMainHeader(); //Quasi-dummy at the moment. 
  void           ebcdicHeader(char* outstring);
  double         interpolate(int xyidx, int zidx);

  long           fileSize(char * fileName);

  const Simbox * simbox_;
  SegYTrace   ** traces_;
  int            format_;           // Type of input format.
  int            datasize_;         // Bytes per datapoint in file.

  double         x0_;               
  double         y0_;               
  static double  z0_;               // Must be common for all SegY files, static for writing.

  double         zPad_;             // The padding length in each end of a trace.

  double         dx_;               
  double         dy_;               
  double         dz_;

  double         sinRot_;          
  double         cosRot_;          

  int            nx_;
  int            ny_;
  int            nz_;

  int            inLine0_;
  int            crossLine0_;

  int            ilStep_;
  int            xlStep_;

  int            yDir_;             // -1 if inLine is reverse y.

  int            error_;
  char         * errMsg_;
  FILE         * outFile_;
};


class SegYTrace{
public:

  SegYTrace(float * data, int jStart, int jEnd, int format);
  ~SegYTrace();

  float          getValue(int j);
  int            getLegalIndex(int index); //Returns closest legal index.

private:
  float        * data_;
  int            jStart_;
  int            jEnd_;
};

#endif
