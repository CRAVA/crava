#ifndef CORR_H
#define CORR_H

#include "nrlib/surface/regularsurface.hpp"

class Corr{
public:
  Corr(float   ** pointVar0,
       float   ** Var0, 
       float    * CorrT, 
       int        n, 
       float      dt, 
       NRLib2::RegularSurface<double> * CorrXY);

  ~Corr(void); 
  const float   ** getVar0() const;
  const float    * getCorrT(int &n, float &dt) const;
  
  const NRLib2::RegularSurface<double>*  getCorrXY() const;
  
  const int        getn()  const { return n_;} ;
  const int        getnx() const { return CorrXY_->GetNI();} 
  const int        getny() const { return CorrXY_->GetNJ();} 
  float            getdt() const { return(dt_); }
  void             setVar0(float ** var0);

  void             dumpResult() const;
  void             printVariancesToScreen();

private:
  float         ** Var0_;      // blocked variance this is used in the computations
  float         ** pointVar0_;
  float          * CorrT_;
  int              n_;         // length of CorT
  float            dt_;        // time step for CorT

  NRLib2::RegularSurface<double> * CorrXY_;
};
#endif
