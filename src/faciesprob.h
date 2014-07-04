/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"
#include "libs/nrlib/flens/nrlib_flens.hpp"

#include <rplib/distributionsrock.h>
#include <src/posteriorelasticpdf.h>
#include <src/posteriorelasticpdf3d.h>
#include <rplib/syntwelldata.h>

#include <map>

class FFTGrid;
class AVOInversion;
class Simbox;
class RandomGen;
class CKrigingAdmin;
class KrigingData;
class ModelSettings;
class FilterWellLogs;
class SpatialWellFilter;
class ModelGeneral;
class CravaTrend;
class SeismicParametersHolder;
class BlockedLogsCommon;

class FaciesProb
{
public:
  FaciesProb(FFTGrid                                  * vp,
             FFTGrid                                  * vs,
             FFTGrid                                  * rho,
             int                                        nFac,
             float                                      p_undef,
             const std::vector<float>                 & priorFacies,
             std::vector<FFTGrid *>                     priorFaciesCubes,
             const std::vector<NRLib::Matrix>         & sigmaEOrig,
             bool                                       useFilter,
             std::map<std::string, BlockedLogsCommon *> blocked_logs,
             const std::vector<Surface *>             & faciesEstimInterval,
             const double                               dz,
             bool                                       relative,
             bool                                       noVs,
             AVOInversion                             * avoInversionResult,
             const std::vector<Grid2D *>              & noiseScale,
             const ModelSettings                      * modelSettings,
             FFTGrid                                  * seismicLH,
             const std::vector<std::string>             facies_names);

  FaciesProb(FFTGrid                                            * vp,
             FFTGrid                                            * vs,
             FFTGrid                                            * rho,
             int                                                  nFac,
             float                                                p_undef,
             const std::vector<float>                           & priorFacies,
             std::vector<FFTGrid *>                               priorFaciesCubes,
             FFTGrid                                            * seismicLH,
             const std::map<std::string, DistributionsRock *>   & rock_distributions,
             const std::vector<std::string>                     & facies_names,
             const std::vector<Surface *>                       & faciesEstimInterval,
             AVOInversion                                       * avoInversionResult,
             SeismicParametersHolder                            & seismicParameters,
             const std::vector<Grid2D *>                        & noiseScale,
             const ModelSettings                                * modelSettings,
             SpatialWellFilter                                  * filteredlogs,
             std::map<std::string, BlockedLogsCommon *>           blocked_wells,
             CravaTrend                                         & trend_cubes,
             //int                                                  nWells = 0,
             const double                                         dz = 0.0,
             bool                                                 useFilter = false,
             bool                                                 relative = false,
             const double                                         trend1_min = 0.0,
             const double                                         trend1_max = 0.0,
             const double                                         trend2_min = 0.0,
             const double                                         trend2_max = 0.0);


  ~FaciesProb();

  FFTGrid              * getFaciesProb(int i){return faciesProb_[i];};

  FFTGrid              * getFaciesProbUndef(){return faciesProbUndef_;};

  void                   calculateConditionalFaciesProb(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                                        const std::vector<Surface *>             & faciesEstimInterval,
                                                        const std::vector<std::string>           & faciesNames,
                                                        const double                               dz);

  void                   calculateFaciesProbGeomodel(const std::vector<float>           & priorFacies,
                                                     std::vector<FFTGrid *>               priorFaciesCubes);

  std::vector<double>    calculateChiSquareTest(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                                const std::vector<Surface *>             & faciesEstimInterval);

  FFTGrid *              createLHCube(FFTGrid                  * likelihood,
                                      int                        fac,
                                      const std::vector<float> & priorFacies,
                                      std::vector<FFTGrid *>     priorFaciesCubes);

  void                   writeBWFaciesProb(std::map<std::string, BlockedLogsCommon *> blocked_wells);

private:

  int                    MakePosteriorElasticPDFRockPhysics(std::vector<std::vector<PosteriorElasticPDF *> >       & posteriorPdf,
                                                            std::vector<Simbox*>                                   & volume,
                                                            AVOInversion                                           * avoInversionResult,
                                                            SeismicParametersHolder                                & seismicParameters,
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

  void                   checkProbabilities(const std::vector<std::string>  & faciesNames,
                                            FFTGrid                        ** faciesProb,
                                            int                               nFacies) const;

  int                    MakePosteriorElasticPDF3D(std::vector<std::vector<PosteriorElasticPDF *> >       & posteriorPdf3d,
                                                    std::vector<Simbox*>                                  & volume,
                                                    const std::vector<NRLib::Matrix>                      & sigmaEOrig,
                                                    bool                                                    useFilter,
                                                    std::map<std::string, BlockedLogsCommon *>              blocked_wells,
                                                    const std::vector<Surface *>                          & faciesEstimInterval,
                                                    const double                                            dz,
                                                    bool                                                    relative,
                                                    bool                                                    noVs,
                                                    std::vector<FFTGrid *>                                & priorFaciesCubes,
                                                    AVOInversion                                          * avoInversionResult,
                                                    const std::vector<Grid2D *>                           & noiseScale,
                                                    const ModelSettings                                   * modelSettings,
                                                    const std::vector<std::string>                          facies_names);

  void                   CalculateFaciesProbFromPosteriorElasticPDF(FFTGrid                                                   * vpgrid,
                                                                    FFTGrid                                                   * vsgrid,
                                                                    FFTGrid                                                   * rhogrid,
                                                                    const std::vector<std::vector<PosteriorElasticPDF *> >    & posteriorPdf,
                                                                    const std::vector<Simbox *>                               & volume,
                                                                    int                                                         nDimensions,
                                                                    float                                                       p_undefined,
                                                                    const std::vector<float>                                  & priorFacies,
                                                                    std::vector<FFTGrid *>                                    & priorFaciesCubes,
                                                                    const std::vector<Grid2D *>                               & noiseScale,
                                                                    FFTGrid                                                   * seismicLH,
                                                                    bool                                                        faciesProbFromRockPhysics,
                                                                    CravaTrend                                                & trend_cubes);


  void                   makeFaciesProb(int                                 nFac,
                                        FFTGrid                           * vp,
                                        FFTGrid                           * vs,
                                        FFTGrid                           * rho,
                                        const std::vector<NRLib::Matrix>  & sigmaEOrig,
                                        bool                                useFilter,
                                        std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                        const std::vector<Surface *>      & faciesEstimInterval,
                                        const double                        dz,
                                        bool                                relative,
                                        bool                                noVs,
                                        float                               p_undef,
                                        const std::vector<float>          & priorFacies,
                                        std::vector<FFTGrid *>              priorFaciesCubes,
                                        AVOInversion                      * avoInversionResult,
                                        const std::vector<Grid2D *>       & noiseScale,
                                        const ModelSettings               * modelSettings,
                                        FFTGrid                           * seismicLH,
                                        const std::vector<std::string>      facies_names);

  std::vector<FFTGrid *> makeFaciesHistAndSetPriorProb(const std::vector<double> & vp,
                                                       const std::vector<double> & vs,
                                                       const std::vector<double> & rho,
                                                       const std::vector<int>   & facies,
                                                       const Simbox             * volume);

  void                   makeFaciesDens(int                                nfac,
                                        const std::vector<NRLib::Matrix> & sigmaEOrig,
                                        bool                               useFilter,
                                        bool                               noVs,
                                        const std::vector<double>         & vpFiltered,
                                        const std::vector<double>         & vsFiltered,
                                        const std::vector<double>         & rhoFiltered,
                                        const std::vector<int>           & faciesLog,
                                        std::vector<FFTGrid *>           & density,
                                        Simbox                          ** volume,
                                        int                                index,
                                        NRLib::Matrix                    & G,
                                        AVOInversion                     * avoInversionResult,
                                        const std::vector<Grid2D *>      & noiseScale);

  void                   setNeededLogsSpatial(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                              const std::vector<Surface *> & faciesEstimInterval,
                                              const double                   dz,
                                              bool                           relative,
                                              bool                           noVs,
                                              bool                           useFilter,
                                              std::vector<double>           & vpFiltered,
                                              std::vector<double>           & vsFiltered,
                                              std::vector<double>           & rhoFiltered,
                                              std::vector<int>              & faciesLog);

  void                   CalculateVariances(const std::vector<double> & vp,
                                            const std::vector<double> & vs,
                                            const std::vector<double> & rho,
                                            const std::vector<int>    & facies,
                                            double                    & varVp,
                                            double                    & varVs,
                                            double                    & varRho);

  float                  findDensity(float                                       vp,
                                     float                                       vs,
                                     float                                       rho,
                                     const std::vector<std::vector<FFTGrid*> > & density,
                                     int                                         facies,
                                     const std::vector<Simbox *>               & volume,
                                     const std::vector<float>                  & t,
                                     int                                         nAng);

  float                  FindDensityFromPosteriorPDF(const double                                          & vp,
                                                     const double                                          & vs,
                                                     const double                                          & rho,
                                                     const double                                          & s1,
                                                     const double                                          & s2,
                                                     const std::vector<std::vector<PosteriorElasticPDF*> > & posteriorPDF,
                                                     const int                                               facies,
                                                     const std::vector<Simbox *>                           & volume,
                                                     const std::vector<float>                              & t,
                                                     int                                                     nAng,
                                                     bool                                                    faciesProbFromRockPhysics);

  void                   resampleAndWriteDensity(const FFTGrid     * const density,
                                                 const std::string & fileName,
                                                 const Simbox      * origVol,
                                                 Simbox            * volume,
                                                 int                 gridNo,
                                                 bool                writeSurface);
  Simbox *               createExpVol(const Simbox * volume);


  void                   calculateFaciesProb(FFTGrid                      * vpgrid,
                                             FFTGrid                      * vsgrid,
                                             FFTGrid                      * rhogrid,
                                             const std::vector<std::vector<FFTGrid *> > & density,
                                             const std::vector<Simbox *>    volume,
                                             float                          p_undef,
                                             const std::vector<float>     & priorFacies,
                                             std::vector<FFTGrid *>       & priorFaciesCubes,
                                             const std::vector<Grid2D *>   & noiseScale,
                                             FFTGrid                       * seismicLH);

  // This function draws random facies length following a geometric distribution with mean 10.0 ms, i.e. dz/10.0 bins
  void                   GenerateSyntWellData(std::vector<SyntWellData *>                              & syntWellData,
                                              const std::map<std::string, DistributionsRock *>         & rock_distributions,
                                              const std::vector<std::string>                           & facies_names,
                                              const std::vector<double>                                & trend1,
                                              const std::vector<double>                                & trend2,
                                              double                                                     dz,
                                              int                                                        nWellsPerCombinationOfTrendParams);

  // shared routine for the calculateFaciesProb functions

  void                   normalizeCubes(std::vector<FFTGrid *> & priorFaciesCubes);

  void                   checkConditionalProbabilities(float                         ** condFaciesProb,
                                                       const std::vector<std::string> & faciesNames,
                                                       const int                        nFacies,
                                                       const std::string              & identifier,
                                                       const bool                       accumulative,
                                                       bool                           & lowProbs,
                                                       int                            * faciesCount = NULL);

  void                   FillSigmaPriAndSigmaPost(AVOInversion * avoInversionResult,
                                                  double      ** sigma_pri_temp,
                                                  double      ** sigma_post_temp);

  void    SolveGEVProblem(double                           ** sigma_prior,
                          double                           ** sigma_posterior,
                          std::vector<std::vector<double> > & V);

  void CalculateTransform2D(const std::vector<double>  & d1,     //data vector 1
                            const std::vector<double>  & d2,     //data vector 2
                            const std::vector<double>  & d3,     //data vector 3
                            std::vector<double>        & x,      //output dimension 1
                            std::vector<double>        & y,      //output dimension 2
                            const std::vector<double>  & v1,     //linear transformation 1
                            const std::vector<double>  & v2);    //linear transformation 2

  int             nFacies_;
  FFTGrid      ** faciesProb_;
  FFTGrid       * faciesProbUndef_;
  int             syntWellLength_; // default is 1000 in the constructor, mean facies length is set to 10
  int             nBinsTrend_; //default is 20 in the constructor
  double          trend1BinSize_;
  double          trend2BinSize_;
};
#endif
