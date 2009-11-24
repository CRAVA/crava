#ifndef WAVELETFILTER_H
#define WAVELETFILTER_H

#include "src/definitions.h"

class WaveletFilter {
public: 
  //Constructor
  WaveletFilter(const std::string & filename,
                int               & errCode,
                char              * errText);
  ~WaveletFilter();

  bool                    readFile(const std::string  & filename,
                                   int                & errCode,
                                   char               * errText);

  const double            getAlpha1(double phi, double psi) const {return alpha1_.GetZ(phi,psi)  ;}
  const double            getHalpha(double phi, double psi) const {return Halpha_.GetZ(phi,psi)  ;}
  const double            getBeta1(double phi, double psi)  const {return beta1_.GetZ(phi,psi)   ;}
  const double            getHbeta(double phi, double psi)  const {return Hbeta_.GetZ(phi,psi)   ;}

private:
  
  Surface                 alpha1_;
  Surface                 Halpha_;
  Surface                 beta1_;
  Surface                 Hbeta_;

};

#endif

