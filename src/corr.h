#ifndef CORR_H
#define CORR_H

#include "lib/irapgrid.h"

class Corr{
public:
  Corr(float   ** pointVar0,
       float   ** Var0, 
       float    * CorrT, 
       int        n, 
       float      dt, 
       irapgrid * CorrXY);
  ~Corr(void); 
  const float   ** getVar0() const;
  const float    * getCorrT(int &n, float &dt) const;
  const irapgrid * getCorrXY() const;
  const int        getn()  const { return n_;} ;
  const int        getnx() const { return CorrXY_->nx ;} 
  const int        getny() const { return CorrXY_->ny ;} 
  float            getdt() const { return(dt_); }
  void             setVar0(float ** var0);

  void             dumpResult() const;
  void             printVariancesToScreen();

private:
  float         ** Var0_;      // blocked variance this is used in the computations
  float         ** pointVar0_;
  float          * CorrT_;
  irapgrid       * CorrXY_;
  int              n_;         // length of CorT
  float            dt_;        // time step for CorT
};
#endif
