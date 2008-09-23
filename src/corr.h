#ifndef CORR_H
#define CORR_H

#include "src/definitions.h"
#include "nrlib/surface/regularsurface.hpp"

class Corr{
public:
  Corr(float   ** pointVar0,
       float   ** Var0, 
       float    * CorrT, 
       int        n, 
       float      dt, 
       Surface  * CorrXY);
  ~Corr(void); 
  float    * getCorrT(int &n, float &dt) const;
  float   ** getVar0()                   const { return Var0_            ;}
  Surface  * getCorrXY()                 const { return CorrXY_          ;}
  int        getn()                      const { return n_               ;}
  int        getnx()                     const { return CorrXY_->GetNI() ;} 
  int        getny()                     const { return CorrXY_->GetNJ() ;} 
  float      getdt()                     const { return dt_              ;}
  void       dumpResult()                const;
  void       setVar0(float ** var0);
  void       printVariancesToScreen();

private:
  float   ** Var0_;      // blocked variance this is used in the computations
  float   ** pointVar0_;
  float    * CorrT_;
  int        n_;         // length of CorT
  float      dt_;        // time step for CorT
  Surface  * CorrXY_;
};
#endif
