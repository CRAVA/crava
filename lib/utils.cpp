#include <iostream>

#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"

using namespace NRLib2;

//------------------------------------------------------------
void    
Utils::writeTitler(const char * text)
{
  const int length = static_cast<int>(strlen(text));
  LogKit::LogFormatted(LogKit::LOW,"\n%s\n",text);
  for (int i = 0 ; i < length ; i++)
    LogKit::LogFormatted(LogKit::LOW,"-");
  LogKit::LogFormatted(LogKit::LOW,"\n");
}  

//------------------------------------------------------------
void
Utils::copyVector(const float * from,
                  float       * to,
                  int           ndim)
{
  for (int i = 0 ; i < ndim ; i++)
    to[i] = from[i];
}

//------------------------------------------------------------
void
Utils::copyMatrix(const float ** from,
                  float       ** to,
                  int            ndim1,
                  int            ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++)
    for (int j = 0 ; j < ndim2 ; j++)
      to[i][j] = from[i][j];
}

//------------------------------------------------------------
void    
Utils::writeVector(float * vector,
                   int     ndim)
{
  for (int i = 0 ; i < ndim ; i++) {
    LogKit::LogFormatted(LogKit::LOW,"%10.6f ",vector[i]);
  }
  LogKit::LogFormatted(LogKit::LOW,"\n");
}

//------------------------------------------------------------
void    
Utils::writeVector(double * vector,
                   int      ndim)
{
  for (int i = 0 ; i < ndim ; i++) {
    LogKit::LogFormatted(LogKit::LOW,"%10.6f ",vector[i]);
  }
  LogKit::LogFormatted(LogKit::LOW,"\n");
}

//------------------------------------------------------------
void    
Utils::writeMatrix(float ** matrix,
                   int      ndim1,
                   int      ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++) {
    for (int j = 0 ; j < ndim2 ; j++) {
      LogKit::LogFormatted(LogKit::LOW,"%10.6f ",matrix[i][j]);
    }
    LogKit::LogFormatted(LogKit::LOW,"\n");
  }
}

//------------------------------------------------------------
void    
Utils::writeMatrix(double ** matrix,
                   int       ndim1,
                   int       ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++) {
    for (int j = 0 ; j < ndim2 ; j++) {
      LogKit::LogFormatted(LogKit::LOW,"%10.6f ",matrix[i][j]);
    }
    LogKit::LogFormatted(LogKit::LOW,"\n");
  }
}





