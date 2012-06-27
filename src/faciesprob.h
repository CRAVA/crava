#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include <rplib/distributionsrock.h>

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
class ModelGeneral;
class CravaTrend;

class FaciesProb
{
public:
  FaciesProb(FFTGrid                      * alpha,
             FFTGrid                      * beta,
             FFTGrid                      * rho,
             int                            nFac,
             float                          p_undef,
             const std::vector<float>     & priorFacies,
             std::vector<FFTGrid *>         priorFaciesCubes,
             const std::vector<double **> & sigmaEOrig,
             bool                           useFilter,
             std::vector<WellData *>        wells,
             int                            nWells,
             const std::vector<Surface *> & faciesEstimInterval,
             const double                   dz,
             bool                           relative,
             bool                           noVs,
             Crava                         *cravaResult,
             const std::vector<Grid2D *>  & noiseScale,
             const ModelSettings          * modelSettings,
             FFTGrid                      * seismicLH);

  FaciesProb(FFTGrid                                           * alpha,
             FFTGrid                                           * beta,
             FFTGrid                                           * rho,
             int                                                 nFac,
             float                                               p_undef,
             FFTGrid                                           * seismicLH,
             ModelGeneral                                      * modelGeneral,
             const std::map<std::string, DistributionsRock *>    rock_distributions,
             const CravaTrend                                  & trend_cubes);

  ~FaciesProb();

  FFTGrid              * getFaciesProb(int i){return faciesProb_[i];};

  FFTGrid              * getFaciesProbUndef(){return faciesProbUndef_;};

  void                   calculateConditionalFaciesProb(std::vector<WellData *>          wells,
                                                        int                              nwells,
                                                        const std::vector<Surface *>   & faciesEstimInterval,
                                                        const std::vector<std::string> & faciesNames,
                                                        const double                     dz);

  void                   calculateFaciesProbGeomodel(const std::vector<float>           & priorFacies,
                                                     std::vector<FFTGrid *>               priorFaciesCubes);

  std::vector<double>    calculateChiSquareTest(std::vector<WellData *>        wells,
                                                int                            nWells,
                                                const std::vector<Surface *> & faciesEstimInterval);

  void                   calculateChiSquareAlternativeTest(std::vector<WellData *>        wells,
                                                           int                            nWells,
                                                           const std::vector<Surface *> & faciesEstimInterval,
                                                           const ModelSettings          * modelSettings);

  FFTGrid *              createLHCube(FFTGrid                  * likelihood,
                                      int                        fac,
                                      const std::vector<float> & priorFacies,
                                      std::vector<FFTGrid *>     priorFaciesCubes);

  void writeBWFaciesProb(std::vector<WellData *> wells,
                         int                     nWells);

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
                                        std::vector<WellData *>        wells,
                                        int                            nWells,
                                        const std::vector<Surface *> & faciesEstimInterval,
                                        const double                   dz,
                                        bool                           relative,
                                        bool                           noVs,
                                        float                          p_undef,
                                        const std::vector<float>     & priorFacies,
                                        std::vector<FFTGrid *>       & priorFaciesCubes,
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

  void                   setNeededLogsSpatial(std::vector<WellData  *>       wells,
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
                                             const std::vector<float>     & priorFacies,
                                             std::vector<FFTGrid *>       & priorFaciesCubes,
                                             const std::vector<Grid2D *>   & noiseScale,
                                             FFTGrid                       * seismicLH);

  void                   calculateFaciesProbFromRockPhysicsModel(FFTGrid                                          * alphagrid,
                                                                 FFTGrid                                          * betagrid,
                                                                 FFTGrid                                          * rhogrid,
                                                                 float                                              p_undef,
                                                                 FFTGrid                                          * seismicLH,
                                                                 ModelGeneral                                     * modelGeneral,
                                                                 const std::map<std::string, DistributionsRock *> & rock_distributions,
                                                                 const CravaTrend                                 & trend_cubes);

  void                   normalizeCubes(std::vector<FFTGrid *> & priorFaciesCubes);

  void                   checkConditionalProbabilities(float                         ** condFaciesProb,
                                                       const std::vector<std::string> & faciesNames,
                                                       const int                        nFacies,
                                                       const std::string              & identifier,
                                                       const bool                       accumulative,
                                                       bool                           & lowProbs,
                                                       int                            * faciesCount = NULL);
};
#endif
