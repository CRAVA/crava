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

  bool                    readFile(const std::string  & filename,
                                   int                & errCode,
                                   std::string        & errText);

  const float            getAlpha1(double phi, double psi) const {return alpha1_.GetZ(phi,psi)  ;}
  const float            getHalpha(double phi, double psi) const {return Halpha_.GetZ(phi,psi)  ;}
  const float            getBeta1(double phi, double psi)  const {return beta1_.GetZ(phi,psi)   ;}
  const float            getHbeta(double phi, double psi)  const {return Hbeta_.GetZ(phi,psi)   ;}
  const bool              hasHalpha() const {return hasHalpha_;}  

private:
//  void                    createGrid();

  NRLib::RegularSurface<float> alpha1_;
  NRLib::RegularSurface<float> Halpha_;
  NRLib::RegularSurface<float> beta1_;
  NRLib::RegularSurface<float> Hbeta_;

  bool                    hasHalpha_;

};

#endif

