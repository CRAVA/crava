/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef WAVELETFILTER_H
#define WAVELETFILTER_H

#include "src/definitions.h"

class WaveletFilter {
public:
  //Constructor
  WaveletFilter(const std::string & filename,
                int               & errCode,
                std::string       & errText);
  WaveletFilter() {};
  ~WaveletFilter();

  void                    readFile(const std::string  & filename,
                                   int                & errCode,
                                   std::string        & errText);

  const float            getAlpha1(double phi, double psi) const;
  const float            getHalpha(double phi, double psi) const;
  const float            getBeta1(double phi, double psi)  const ;
  const float            getHbeta(double phi, double psi)  const ;
  const bool             hasHalpha() const {return hasHalpha_;}

private:
//  void                    createGrid();

  NRLib::RegularSurface<float> alpha1_;
  NRLib::RegularSurface<float> Halpha_;
  NRLib::RegularSurface<float> beta1_;
  NRLib::RegularSurface<float> Hbeta_;

  bool                    hasHalpha_;

};

#endif

