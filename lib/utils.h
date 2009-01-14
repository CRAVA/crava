
#ifndef UTILS_H
#define UTILS_H

#include "src/definitions.h"
#include "nrlib/iotools/logkit.hpp"

class Utils
{
public:
  static void    writeHeader(const std::string & text, LogKit::MessageLevels logLevel = LogKit::LOW);
  static void    writeTitler(const char * text);

  static void    copyVector(const float * from,
                            float       * to,
                            int           ndim);
  static void    copyMatrix(const float ** from,
                            float       ** to,
                            int            ndim1,
                            int            ndim2);

  static void    writeVector(float * vector,
                             int     ndim);
  static void    writeVector(double * vector,
                             int      ndim);

  static void    writeMatrix(float ** matrix,
                             int      ndim1,
                             int      ndim2);  
  static void    writeMatrix(double ** matrix,
                             int       ndim1,
                             int       ndim2);  
};

#endif
