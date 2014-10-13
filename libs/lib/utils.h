/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef UTILS_H
#define UTILS_H

#include "src/definitions.h"
#include "nrlib/iotools/logkit.hpp"
#include "fftw.h"


class Utils
{
public:
  static void    writeTitler(const std::string & text);

  static void    copyVector(const int * from,
                            int       * to,
                            int         ndim);
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

  static void    writeVectorToFile(const std::string & filename,
                                   float             * vector,
                                   int                 ndim);
  static void    writeVectorToFile(const std::string        & filename,
                                   const std::vector<float> & vector);

  static void    writeMatrix(float ** matrix,
                             int      ndim1,
                             int      ndim2);
  static void    writeMatrix(double ** matrix,
                             int       ndim1,
                             int       ndim2);

  static void    fft(fftw_real    * rAmp,
                     fftw_complex * cAmp,
                     int            nt);

  static void    fftInv(fftw_complex * cAmp,
                        fftw_real    * rAmp,
                        int            nt);

  static  void   readUntilStop(int           pos,
                               std::string & in,
                               std::string & out,
                               std::string   read);

  static int     findEnd(std::string & seek,
                         int           start,
                         std::string & find);

  static void     ShiftTrace(std::vector<fftw_real> & trace, bool shift_up = true);

  static void     ShiftTrace(fftw_real * trace,
                             size_t      n_data,
                             bool        shift_up = true);

};

#endif
