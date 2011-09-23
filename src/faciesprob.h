#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

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
             bool                           useFilter,
             const WellData              ** wells,
             int                            nWells,
             const std::vector<Surface *> & faciesEstimInterval,
             const double                   dz,
             bool                           relative,
             bool                           noVs,
             Crava                         *cravaResult,
             const std::vector<Grid2D *>  & noiseScale,
             const ModelSettings          * modelSettings,
             FFTGrid                      * seismicLH);

  ~FaciesProb();

  FFTGrid              * getFaciesProb(int i){return faciesProb_[i];};

  FFTGrid              * getFaciesProbUndef(){return faciesProbUndef_;};

  void                   calculateConditionalFaciesProb(WellData                      ** wells,
                                                        int                              nwells,
                                                        const std::vector<Surface *>   & faciesEstimInterval,
                                                        const std::vector<std::string> & faciesNames,
                                                        const double                     dz);

  void                   calculateFaciesProbGeomodel(const float *                  priorFacies,
                                                     FFTGrid                    ** priorFaciesCubes);

  std::vector<double>    calculateChiSquareTest(WellData                    ** wells,
                                                int                            nWells,
                                                const std::vector<Surface *> & faciesEstimInterval);

  void                   calculateChiSquareAlternativeTest(WellData                    ** wells,
                                                           int                            nWells,
                                                           const std::vector<Surface *> & faciesEstimInterval,
                                                           const ModelSettings          * modelSettings);

  FFTGrid *              createLHCube(FFTGrid     * likelihood,
                                      int           fac,
                                      const float * priorFacies,
                                      FFTGrid    ** priorFaciesCubes);

  void writeBWFaciesProb(WellData ** wells,
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
                                        bool                           useFilter,
                                        const WellData              ** wells,
                                        int                            nWells,
                                        const std::vector<Surface *> & faciesEstimInterval,
                                        const double                   dz,
                                        bool                           relative,
                                        bool                           noVs,
                                        float                          p_undef,
                                        const float                  * priorFacies,
                                        FFTGrid                     ** priorFaciesCubes,
                                        Crava                        * cravaResult,
                                        const std::vector<Grid2D *>  & noiseScale,
                                        const ModelSettings          * modelSettings,
                                        FFTGrid                      * seismicLH);

  std::vector<FFTGrid *> makeFaciesHistAndSetPriorProb(const std::vector<float> & alpha,
                                                       const std::vector<float> & beta,
                                                       const std::vector<float> & rho,
                                                       const std::vector<int>   & facies,
                                                       const Simbox             * volume);

  void                   makeFaciesDens(int nfac,
                                        const std::vector<double **> & sigmaEOrig,
                                        bool                           useFilter,
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
                                              bool                           useFilter,
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

  float                  findDensity(float                                       alpha,
                                     float                                       beta,
                                     float                                       rho,
                                     const std::vector<std::vector<FFTGrid*> > & density,
                                     int                                         facies,
                                     const std::vector<Simbox *>               & volume,
                                     const std::vector<float>                  & t,
                                     int                                         nAng);

  void                   resampleAndWriteDensity(const FFTGrid     * density,
                                                 const std::string & fileName,
                                                 const Simbox      * origVol,
                                                 Simbox            * volume,
                                                 int                 gridNo,
                                                 bool                writeSurface);
  Simbox *               createExpVol(const Simbox * volume);


  void                   calculateFaciesProb(FFTGrid                      * alphagrid,
                                             FFTGrid                      * betagrid,
                                             FFTGrid                      * rhogrid,
                                             const std::vector<std::vector<FFTGrid *> > & density,
                                             const std::vector<Simbox *>    volume,
                                             float                          p_undef,
                                             const float                  * priorFacies,
                                             FFTGrid                     ** priorFaciesCubes,
                                             const std::vector<Grid2D *>   & noiseScale,
                                             FFTGrid                       * seismicLH);

  void                   normalizeCubes(FFTGrid **priorFaciesCubes);

  void                   checkConditionalProbabilities(float                         ** condFaciesProb,
                                                       const std::vector<std::string> & faciesNames,
                                                       const int                        nFacies,
                                                       const std::string              & identifier,
                                                       const bool                       accumulative,
                                                       bool                           & lowProbs);
};
#endif
