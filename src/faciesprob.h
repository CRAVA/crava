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
             const std::vector<double **> & sigmaEOrig,
             const WellData              ** wells,
             int                            nWells,
             const std::vector<Surface *> & faciesEstimInterval,
             const double                   dz,
             bool                           relative,
             bool                           noVs,
             Crava                         *cravaResult,
             const std::vector<Grid2D *>   & noiseScale);
  ~FaciesProb();

  FFTGrid              * getFaciesProb(int i){return faciesProb_[i];};

  void                   calculateConditionalFaciesProb(WellData                    ** wells, 
                                                        int                            nwells, 
                                                        const std::vector<Surface *> & faciesEstimInterval,
                                                        const ModelSettings          * modelSettings,
                                                        const double                   dz);
  void calculateFaciesProbGeomodel(const float *                  priorFacies,
                                    FFTGrid                    ** priorFaciesCubes);
  
  void                   calculateChiSquareTest(WellData                    ** wells, 
                                                int                            nWells, 
                                                const std::vector<Surface *> & faciesEstimInterval);

  void                   calculateChiSquareAlternativeTest(WellData                    ** wells, 
                                                           int                            nWells, 
                                                           const std::vector<Surface *> & faciesEstimInterval,
                                                           const ModelSettings          * modelSettings);

  void                   writeBWFaciesProb(WellData ** wells, 
                                           int         nWells);

private:
  int             nFacies_;
  FFTGrid      ** faciesProb_;
  FFTGrid      *  faciesProbUndef_;

  void                   makeFaciesProb(int                            nfac, 
                                        FFTGrid                      * postAlpha, 
                                        FFTGrid                      * postBeta, 
                                        FFTGrid                      * postRho,
                                        const std::vector<double **> & sigmaEOrig, 
                                        const WellData              ** wells, 
                                        int                            nWells,
                                        const std::vector<Surface *> & faciesEstimInterval,
                                        const double                   dz,
                                        bool                           relative,
                                        bool                           noVs,
                                        float                          p_undef,
                                        const float                  * priorFacies,
                                        FFTGrid                     ** priorFaciesCubes,
                                        Crava                         * cravaResult,
                                        const std::vector<Grid2D *>   & noiseScale);

  std::vector<FFTGrid *> makeFaciesHistAndSetPriorProb(const std::vector<float> & alpha,
                                                       const std::vector<float> & beta,
                                                       const std::vector<float> & rho,
                                                       const std::vector<int>   & facies,
                                                       const Simbox             * volume);

  void                   makeFaciesDens(int nfac, 
                                        const std::vector<double **> & sigmaEOrig,
                                        bool                           noVs,  
                                        const std::vector<float>     & alphaFiltered,
                                        const std::vector<float>     & betaFiltered,
                                        const std::vector<float>     & rhoFiltered,
                                        const std::vector<int>       & faciesLog,
                                        std::vector<FFTGrid *>       & density,
                                        Simbox                      ** volume,
                                        int                            index, 
                                        double                       **G,
                                        Crava                         *cravaResult,
                                        const std::vector<Grid2D *>   & noiseScale);

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
                                     std::vector<std::vector<FFTGrid*> >       density, 
                                     int facies,
                                     const std::vector<Simbox *> volume,
                                      std::vector<float> t, 
                                      int nAng);

  void                   calculateFaciesProb(FFTGrid                      * alphagrid, 
                                             FFTGrid                      * betagrid, 
                                             FFTGrid                      * rhogrid,
                                             const std::vector<std::vector<FFTGrid *> > & density, 
                                             const std::vector<Simbox *>    volume,
                                             float                          p_undef,
                                             const float                  * priorFacies,
                                             FFTGrid                     ** priorFaciesCubes,
                                             const std::vector<Grid2D *>   & noiseScale);

  void                   normalizeCubes(FFTGrid **priorFaciesCubes);

};
#endif
