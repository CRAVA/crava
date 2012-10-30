#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include <rplib/distributionsrock.h>
#include <src/posteriorelasticpdf.h>
#include <src/posteriorelasticpdf3d.h>

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
             FFTGrid                      * seismicLH,
             const std::vector<std::string> facies_names);

  /// Constructor. Creates a PosteriorElasticPDF object and calculates the probabilities.
  /// @param[in] alpha: an FFTGrid of Vp.
  /// @param[in] beta: an FFTGrid of Vs.
  /// @param[in] rho: an FFTGrid of rho.
  /// @param[in] nFac: number of Facies.

  FaciesProb(FFTGrid                                            * alpha,
             FFTGrid                                            * beta,
             FFTGrid                                            * rho,
             int                                                  nFac,
             float                                                p_undef,
             const std::vector<float>                           & priorFacies,
             std::vector<FFTGrid *>                               priorFaciesCubes,
             FFTGrid                                            * seismicLH,
             const std::map<std::string, DistributionsRock *>   & rock_distributions,
             const double                                         trend1_min,
             const double                                         trend1_max,
             const double                                         trend2_min,
             const double                                         trend2_max,
             const std::vector<std::string>                     & facies_names,
             const std::vector<Surface *>                       & faciesEstimInterval,
             Crava                                              * cravaResult,
             const std::vector<Grid2D *>                        & noiseScale,
             const ModelSettings                                * modelSettings,
             SpatialWellFilter                                  * filteredlogs,
             std::vector<WellData *>                              wells,
             int                                                  nWells = 0,
             const double                                         dz = 0.0,
             bool                                                 useFilter = false,
             bool                                                 relative = false);


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

  FFTGrid *              createLHCube(FFTGrid                  * likelihood,
                                      int                        fac,
                                      const std::vector<float> & priorFacies,
                                      std::vector<FFTGrid *>     priorFaciesCubes);

  void writeBWFaciesProb(std::vector<WellData *> wells,
                         int                     nWells);

private:
  int             nFacies_;
  int             syntWellLength_; // default is 1000 in the constructor, mean facies length is set to 10
  int             nBinsTrend_; //default is 20 in the constructor
  double          trend1BinSize_;
  double          trend2BinSize_;
  FFTGrid      ** faciesProb_;
  FFTGrid      *  faciesProbUndef_;

  void                   MakePosteriorElasticPDFRockPhysics(std::vector<std::vector<PosteriorElasticPDF *> >       & posteriorPdf,
                                                            SpatialWellFilter                                      * filteredlogs,
                                                            std::vector<FFTGrid *>                                 & priorFaciesCubes,
                                                            const ModelSettings                                    * modelSettings,
                                                            const std::map<std::string, DistributionsRock *>       & rock_distributions,
                                                            const std::vector<std::string>                         & facies_names,
                                                            const double                                           & dz,
                                                            const double                                           & trend1_min,
                                                            const double                                           & trend1_max,
                                                            const double                                           & trend2_min,
                                                            const double                                           & trend2_max,
                                                            bool                                                     useFilter);

  void                   MakePosteriorElasticPDF3D(std::vector<std::vector<PosteriorElasticPDF *> >       & posteriorPdf3d,
                                                 std::vector<Simbox*>                                     & volume,
                                                 const std::vector<double **>                             & sigmaEOrig,
                                                 bool                                                       useFilter,
                                                 std::vector<WellData *>                                    wells,
                                                 int                                                        nWells,
                                                 const std::vector<Surface *>                             & faciesEstimInterval,
                                                 const double                                               dz,
                                                 bool                                                       relative,
                                                 bool                                                       noVs,
                                                 std::vector<FFTGrid *>                                   & priorFaciesCubes,
                                                 Crava                                                    * cravaResult,
                                                 const std::vector<Grid2D *>                              & noiseScale,
                                                 const ModelSettings                                      * modelSettings,
                                                 const std::vector<std::string>                             facies_names);

  void                   CalculateFaciesProbFromPosteriorElasticPDF(FFTGrid                                           * alphagrid,
                                                            FFTGrid                                                   * betagrid,
                                                            FFTGrid                                                   * rhogrid,
                                                            const std::vector<std::vector<PosteriorElasticPDF *> >    & posteriorPdf3d,
                                                            const std::vector<Simbox *>                                 volume,
                                                            float                                                       p_undefined,
                                                            const std::vector<float>                                  & priorFacies,
                                                            std::vector<FFTGrid *>                                    & priorFaciesCubes,
                                                            const std::vector<Grid2D *>                               & noiseScale,
                                                            FFTGrid                                                   * seismicLH);



  void     makeFaciesProb(int                          nFac,
                       FFTGrid                       * alpha,
                       FFTGrid                       * beta,
                       FFTGrid                       * rho,
                       const std::vector<double **>  & sigmaEOrig,
                       bool                            useFilter,
                       std::vector<WellData *>         wells,
                       int                             nWells,
                       const std::vector<Surface *>  & faciesEstimInterval,
                       const double                    dz,
                       bool                            relative,
                       bool                            noVs,
                       float                           p_undef,
                       const std::vector<float>      & priorFacies,
                       std::vector<FFTGrid *>          priorFaciesCubes,
                       Crava                         * cravaResult,
                       const std::vector<Grid2D *>   & noiseScale,
                       const ModelSettings           * modelSettings,
                       FFTGrid                       * seismicLH,
                       const std::vector<std::string>  facies_names);

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

  void                   CalculateVariances(const std::vector<float> & alpha,
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

  float                  FindDensityFromPosteriorPDF(const double                                          & alpha,
                                                     const double                                          & beta,
                                                     const double                                          & rho,
                                                     const double                                          & s1,
                                                     const double                                          & s2,
                                                     const std::vector<std::vector<PosteriorElasticPDF*> > & posteriorPDF,
                                                     const int                                               facies,
                                                     const std::vector<Simbox *>                           & volume,
                                                     const std::vector<float>                              & t,
                                                     int                                                     nAng);

  void                   resampleAndWriteDensity(const FFTGrid     * const density,
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

  // This function draws random facies length following a geometric distribution with mean 10.0 ms, i.e. dz/10.0 bins
  void                   GenerateSyntWellData(std::vector<std::vector<float> >                        & syntWellData,
                                              std::vector<int>                                         & facies,
                                              const std::map<std::string, DistributionsRock *>         & rock_distributions,
                                              const std::vector<std::string>                           & facies_names,
                                              const std::vector<double>                                & trend1,
                                              const std::vector<double>                                & trend2,
                                              double                                                     dz);

  // shared routine for the calculateFaciesProb functions

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
