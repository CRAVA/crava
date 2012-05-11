#ifndef FILTERWELLLOGS_H
#define FILTERWELLLOGS_H
#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

class Simbox;
class CKrigingAdmin;
class KrigingData;
class Corr;

class FilterWellLogs
{
public:
  FilterWellLogs(const Simbox          * timeSimboxConstThick,
                 const Simbox          * timeSimboxOrigThick,
                 const Corr            * correlations,
                 int                     nzp,
                 int                     nz,
                 std::vector<WellData *> wells,
                 int                     nWells,
                 float                   lowCut,
                 float                   highCut,
                 int                     relative);
  ~FilterWellLogs();

  float       ** getVtAlphaFiltered() const { return vtAlphaFiltered_ ;}
  float       ** getVtBetaFiltered()  const { return vtBetaFiltered_  ;}
  float       ** getVtRhoFiltered()   const { return vtRhoFiltered_   ;}

  float       ** getVtAlpha()         const { return vtAlpha_       ;}
  float       ** getVtBeta()          const { return vtBeta_        ;}
  float       ** getVtRho()           const { return vtRho_         ;}

private:
  void           doFiltering(const Simbox          * timeSimboxConstThick,
                             const Simbox          * timeSimboxOrigThick,
                             std::vector<WellData *> wells,
                             int                     nWells,
                             float                ** sigma0,
                             fftw_real             * postcova,
                             fftw_real             * postcovb,
                             fftw_real             * postcovr,
                             fftw_real             * postcrab,
                             fftw_real             * postcrar,
                             fftw_real             * postcrbr,
                             fftw_real             * corrprior,
                             float                   lowCut,
                             float                   highCut,
                             int                     relative,
                             int                     nz,
                             int                     nzp);
  void           extrapolate(float * log,
                             int     nz) ;
  void           calcFilter(fftw_complex ** sigmaK,
                            fftw_complex ** sigmaE,
                            double       ** F);

  float       ** vtAlphaFiltered_;
  float       ** vtBetaFiltered_;
  float       ** vtRhoFiltered_;

  float       ** vtAlpha_;
  float       ** vtBeta_;
  float       ** vtRho_;

  const int      nWells_;
};
#endif

