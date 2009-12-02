#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

class Corr;
class FFTGrid;
class Crava;
class Simbox;
class RandomGen;
class CKrigingAdmin;
class KrigingData;
class ModelSettings;
class FilterWellLogs;
class WellData;
class SpatialWellFilter;

class FaciesProb
{
public:
  FaciesProb(FFTGrid                      * alpha,
             FFTGrid                      * beta,
             FFTGrid                      * rho,
             int                            nFac,
             float                          p_undef, 
             const float                  * priorFacies,
             FFTGrid                     ** priorFaciesCubes,
             const double                ** sigmaEOrig,
             const WellData              ** wells,
             int                            nWells,
             const std::vector<Surface *> & faciesEstimInterval,
             const double                   dz,
             bool                           relative,
             bool                           noVs);
  ~FaciesProb();

  FFTGrid              * getFaciesProb(int i){return faciesProb_[i];};

  void                   calculateConditionalFaciesProb(WellData                    ** wells, 
                                                        int                            nwells, 
                                                        const std::vector<Surface *> & faciesEstimInterval,
                                                        const ModelSettings          * modelSettings,
                                                        const double                   dz);
  void calculateFaciesProbGeomodel(const float *                  priorFacies,
                                             FFTGrid                    ** priorFaciesCubes);

private:
  int             nFacies_;
  FFTGrid      ** faciesProb_;
  FFTGrid      *  faciesProbUndef_;

  void                   makeFaciesProb(int                            nfac, 
                                        FFTGrid                      * postAlpha, 
                                        FFTGrid                      * postBeta, 
                                        FFTGrid                      * postRho,
                                        const double                ** sigmaEOrig, 
                                        const WellData              ** wells, 
                                        int                            nWells,
                                        const std::vector<Surface *> & faciesEstimInterval,
                                        const double                   dz,
                                        bool                           relative,
                                        bool                           noVs,
                                        float                          p_undef,
                                        const float                  * priorFacies,
                                        FFTGrid                     ** priorFaciesCubes);

  std::vector<FFTGrid *> makeFaciesHistAndSetPriorProb(const std::vector<float> & alpha,
                                                       const std::vector<float> & beta,
                                                       const std::vector<float> & rho,
                                                       const std::vector<int>   & facies,
                                                       const Simbox             * volume);

  void                   makeFaciesDens(int nfac, const double  ** sigmaEOrig,
                                        bool                       noVs,  
                                        const std::vector<float> & alphaFiltered,
                                        const std::vector<float> & betaFiltered,
                                        const std::vector<float> & rhoFiltered,
                                        const std::vector<int>   & faciesLog,
                                        std::vector<FFTGrid *>   & density,
                                        Simbox                  ** volume);

  void                   setNeededLogsSpatial(const WellData              ** wells,
                                              int                            nWells,
                                              const std::vector<Surface *> & faciesEstimInterval,
                                              const double                   dz,
                                              bool                           relative,
                                              bool                           noVs,
                                              std::vector<float>           & alphaFiltered,
                                              std::vector<float>           & betaFiltered,
                                              std::vector<float>           & rhoFiltered,
                                              std::vector<int>             & faciesLog);
  
  void                   calculateVariances(const std::vector<float> & alpha,
                                            const std::vector<float> & beta,
                                            const std::vector<float> & rho,
                                            const std::vector<int>   & facies,
                                            float                    & varAlpha,
                                            float                    & varBeta,
                                            float                    & varRho);

  float                  findDensity(float          alpha, 
                                     float          beta, 
                                     float          rho, 
                                     FFTGrid      * density, 
                                     const Simbox * volume);

  void                   calculateFaciesProb(FFTGrid                      * alphagrid, 
                                             FFTGrid                      * betagrid, 
                                             FFTGrid                      * rhogrid,
                                             const std::vector<FFTGrid *> & density, 
                                             const Simbox                 * volume,
                                             float                          p_undef,
                                             const float                  * priorFacies,
                                             FFTGrid                     ** priorFaciesCubes);

  void                   normalizeCubes(FFTGrid **priorFaciesCubes);

};
#endif
