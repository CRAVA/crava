#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

class FFTGrid;
class Crava;
class Simbox;
class RandomGen;
class CKrigingAdmin;
class KrigingData;

class FaciesProb{
public:
	FaciesProb(int filegrid, const float ** sigma0, float *corrprior, Simbox *simbox, const Simbox& osimbox, int nzp, int nz, FFTGrid* bgAlpha, 
             FFTGrid* bgBeta, FFTGrid* bgRho, RandomGen *random, float p_undef);
	~FaciesProb();
  void       makeFaciesProb(int nfac, FFTGrid *alphagrid, FFTGrid *betagrid, FFTGrid *rhogrid);
	FFTGrid*   getFaciesProb(int i){return faciesProb_[i];};
  void       filterWellLogs(WellData **wells, int nwells,
                    fftw_real *postcova,fftw_real *postcovb,fftw_real *postcovr,fftw_real *postcrab,fftw_real *postcrar,
                    fftw_real *postcrbr, float lowCut, float highCut);
private:

  void       calcFilter(fftw_complex **sigmaK, fftw_complex **sigmaE, double **F);
  void       getMinMax();
  void       calculateVariances();
  void       extrapolate(float * log, int nz); 
  Simbox*    simbox_;
  const Simbox   & origSimbox_;
  int        nWells_;
  fftw_real* corrprior_;
  float**    sigma0_;
	float      findDensity(float alpha, float beta, float rho, int facies);
	void       calculateFaciesProb( FFTGrid *alphagrid, FFTGrid *betagrid, FFTGrid *rhogrid);
	int        nFacies_, nbins_, nbinsr_, nobs_;
	int *      nData_;
	float      dalpha_, dbeta_, drho_;
	FFTGrid**  density_;  // Note: normalizes to the number of observations of facies type. 
	float      alphamin_,betamin_,rhomin_, alphamax_,betamax_,rhomax_;
	FFTGrid**   faciesProb_;
  int        fileGrid_;
  int        nzp_, nz_, ndata_;
  FFTGrid*   bgAlpha_;
  FFTGrid*   bgBeta_;
  FFTGrid*   bgRho_;
  float*     alphafiltered_, *betafiltered_, *rhofiltered_;
  int*       facieslog_;
  float      varAlpha_,varBeta_,varRho_;  // variances
  RandomGen* random_;
  float      p_undefined_;
};
#endif


