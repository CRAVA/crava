#ifndef SEISMICPARAMETERSHOLDER
#define SEISMICPARAMETERSHOLDER

#include <src/fftgrid.h>
#include <src/corr.h>

// A class holding the pointers for the seismic parameters
// for easy parameter transmission of the pointers to the class TimeEvolution.

class SeismicParametersHolder
{
public:
  SeismicParametersHolder(void);

  ~SeismicParametersHolder(void);
  //void getMeanReferenceVector(std::vector<FFTGrid *> mu);
  //void getCovReferenceVector(std::vector<FFTGrid *> sigma);

  void setSeismicParameters(FFTGrid  * muAlpha,
                            FFTGrid  * muBeta,
                            FFTGrid  * muRho,
                            Corr     * correlations);

  void setBackgroundParameters(FFTGrid  * muAlpha,
                               FFTGrid  * muBeta,
                               FFTGrid  * muRho);

  void allocateGrids(const int nx,const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad);

  FFTGrid * GetMuAlpha()        { return muAlpha_        ;}
  FFTGrid * GetMuBeta()         { return muBeta_         ;}
  FFTGrid * GetMuRho()          { return muRho_          ;}
  FFTGrid * GetCovAlpha()       { return covAlpha_       ;}
  FFTGrid * GetCovBeta()        { return covBeta_        ;}
  FFTGrid * GetCovRho()         { return covRho_         ;}
  FFTGrid * GetCrCovAlphaBeta() { return crCovAlphaBeta_ ;}
  FFTGrid * GetCrCovAlphaRho()  { return crCovAlphaRho_  ;}
  FFTGrid * GetCrCovBetaRho()   { return crCovBetaRho_   ;}

  void      invFFT();

private:
    FFTGrid * muAlpha_;
    FFTGrid * muBeta_ ;
    FFTGrid * muRho_  ;
    FFTGrid * covAlpha_;
    FFTGrid * covBeta_ ;
    FFTGrid * covRho_  ;
    FFTGrid * crCovAlphaBeta_;
    FFTGrid * crCovAlphaRho_ ;
    FFTGrid * crCovBetaRho_  ;
};
#endif
