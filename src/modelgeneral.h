/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELGENERAL_H
#define MODELGENERAL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/cravatrend.h"
#include "src/seismicparametersholder.h"
#include "src/state4d.h"
#include "src/timeevolution.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionssolid.h"

struct irapgrid;
class Corr;
class Wavelet;
class Vario;
class Simbox;
class WellData;
class FFTGrid;
class RandomGen;
class GridMapping;
class InputFiles;
class TimeLine;
class WellData;
class SeismicParameters;

class ModelGeneral
{
public:
  ModelGeneral(ModelSettings           *& modelSettings,
               const InputFiles         * inputFiles,
               SeismicParametersHolder  & seismicParameters,
               Simbox                  *& timeBGSimbox);
  ~ModelGeneral();

  Simbox                   * getTimeSimbox()            const { return timeSimbox_             ;}
  Simbox                   * getTimeSimboxConstThick()  const { return timeSimboxConstThick_   ;}
  RandomGen                * getRandomGen()             const { return randomGen_              ;}
  GridMapping              * getTimeDepthMapping()      const { return timeDepthMapping_       ;}
  GridMapping              * getTimeCutMapping()        const { return timeCutMapping_         ;}
  CravaTrend               & getTrendCubes()                  { return trend_cubes_            ;}
  CravaTrend                 getTrendCubes()            const { return trend_cubes_            ;}
  Surface                  * getPriorCorrXY()           const { return priorCorrXY_            ;}

  bool                       getVelocityFromInversion() const { return velocityFromInversion_  ;}
  bool                       getFailed()                const { return failed_                 ;}
  std::vector<bool>          getFailedDetails()         const { return failed_details_         ;}

  void                       getCorrGradIJ(float & corrGradI, float &corrGradJ) const;
  Surface                  * getCorrelationDirection()  const { return correlationDirection_   ;}
  const State4D            & getState4D()               const { return state4d_                ;}

  TimeLine                 * getTimeLine()              const { return timeLine_               ;}
  std::vector<WellData *>  & getWells()             /*const*/ { return wells_                  ;}

  static void                readSegyFile(const std::string       & fileName,
                                          FFTGrid                *& target,
                                          Simbox                  * timeSimbox,
                                          Simbox                  * timeCutSimbox,
                                          ModelSettings          *& modelSettings,
                                          const SegyGeometry     *& geometry,
                                          int                       gridType,
                                          const std::string       & parName,
                                          float                     offset,
                                          const TraceHeaderFormat * format,
                                          std::string             & errText,
                                          bool                      nopadding = false);
  static void                checkThatDataCoverGrid(SegY        * segy,
                                                    float         offset,
                                                    Simbox      * timeCutSimbox,
                                                    float         guard_zone,
                                                    std::string & errText);
  static void                readStormFile(const std::string  & fileName,
                                           FFTGrid           *& target,
                                           const int            gridType,
                                           const std::string  & parName,
                                           Simbox             * timeSimbox,
                                           ModelSettings     *& modelSettings,
                                           std::string        & errText,
                                           bool                 isStorm  = true,
                                           bool                 nopadding = true);
  void                      mergeState4D(SeismicParametersHolder &  seismicParameters);
  void                      updateState4D(SeismicParametersHolder &  seismicParameters);

  void                      updateState4DWithSingleParameter(FFTGrid * EPost,
                                                             FFTGrid * CovPost,
                                                             int       parameterNumber);

  void                      updateState4DMu(FFTGrid * mu_vp_static,
                                            FFTGrid * mu_vs_static,
                                            FFTGrid * mu_rho_static,
                                            FFTGrid * mu_vp_dynamic,
                                            FFTGrid * mu_vs_dynamic,
                                            FFTGrid * mu_rho_dynamic);

  std::map<std::string, DistributionsRock *> getRockDistributionTime0() const;

  const std::vector<float>       & getPriorFacies()           /*const*/ { return priorFacies_          ;}
  const std::vector<FFTGrid *>   & getPriorFaciesCubes()      /*const*/ { return priorFaciesProbCubes_ ;}
  const std::vector<std::string> & getFaciesNames(void)           const { return faciesNames_          ;}
  std::vector<int>                 getFaciesLabel()               const { return faciesLabels_         ;}

  void addFaciesLabel(int faciesLabel)                    { faciesLabels_.push_back(faciesLabel) ;}
  void addFaciesName(const std::string & faciesName)      { faciesNames_.push_back(faciesName)   ;}

  bool getIs4DActive() const {return(do4DInversion_);}

  static FFTGrid  * createFFTGrid(int nx,
                                  int ny,
                                  int nz,
                                  int nxp,
                                  int nyp,
                                  int nzp,
                                  bool fileGrid);

  static void       readGridFromFile(const std::string       & fileName,
                                     const std::string       & parName,
                                     const float               offset,
                                     FFTGrid                *& grid,
                                     const SegyGeometry     *& geometry,
                                     const TraceHeaderFormat * format,
                                     int                       gridType,
                                     const Simbox            * timeSimbox,
                                     const Simbox            * timeCutSimbox,
                                     const ModelSettings     * modelSettings,
                                     std::string             & errorText,
                                     bool                      nopadding = false);

  static void       readSegyFile(const std::string       & fileName,
                                 FFTGrid                *& target,
                                 const Simbox            * timeSimbox,
                                 const Simbox            * timeCutSimbox,
                                 const ModelSettings     * modelSettings,
                                 const SegyGeometry     *& geometry,
                                 int                       gridType,
                                 const std::string       & parName,
                                 float                     offset,
                                 const TraceHeaderFormat * format,
                                 std::string             & errText,
                                 bool                      nopadding = false);

  static void       checkThatDataCoverGrid(const SegY   * segy,
                                           float         offset,
                                           const Simbox * timeCutSimbox,
                                           float         guard_zone,
                                           std::string & errText);
  static void       readStormFile(const std::string  & fileName,
                                  FFTGrid           *& target,
                                  const int            gridType,
                                  const std::string  & parName,
                                  const Simbox       * timeSimbox,
                                  const ModelSettings * modelSettings,
                                  std::string        & errText,
                                  bool                 isStorm  = true,
                                  bool                 nopadding = true);
  static void       loadVelocity(FFTGrid           *& velocity,
                                 const Simbox       * timeSimbox,
                                 const Simbox       * timeCutSimbox,
                                 const ModelSettings * modelSettings,
                                 const std::string  & velocityField,
                                 bool               & velocityFromInversion,
                                 std::string        & errText,
                                 bool               & failed);


  void              processWellLocation(FFTGrid                     ** seisCube,
                                        float                       ** reflectionMatrix,
                                        ModelSettings                * modelSettings,
                                        const std::vector<Surface *> & interval);              // Changes wells

  void              processPriorCorrelations(Background                     * background,
                                             std::vector<WellData *>          wells,
                                             const Simbox                   * timeSimbox,
                                             const ModelSettings            * modelSettings,
                                             const std::vector<float>       & priorFacies,
                                             FFTGrid                       ** seisCube,
                                             const InputFiles               * inputFiles,
                                             SeismicParametersHolder        & seismicParameters,
                                             std::string                    & errText,
                                             bool                           & failed);

   void             processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
                                          std::vector<WellData *>        wells,
                                          Simbox                       * timeSimbox,
                                          Simbox                       * timeCutSimbox,
                                          ModelSettings                * modelSettings,
                                          bool                         & failed,
                                          std::string                  & errTxt,
                                          const InputFiles             * inputFiles);

  void              generateRockPhysics3DBackground(const std::vector<DistributionsRock *>           & rock_distribution,
                                                    const std::vector<float>                         & probability,
                                                    FFTGrid                                          & vp,
                                                    FFTGrid                                          & vs,
                                                    FFTGrid                                          & rho);

  void              calculateCovariancesFromRockPhysics(const std::vector<DistributionsRock *>           & rock,
                                                        const std::vector<float>                         & probability,
                                                        NRLib::Grid2D<double>                            & param_corr,
                                                        std::string                                      & errTxt);

  void              complete4DBackground(const int nx,const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad,NRLib::Vector &initial_mean,NRLib::Matrix &initial_cov);

  //void            getInitial3DPriorFrom4D(SeismicParametersHolder & seismicParameters);
  bool              do4DRockPhysicsInversion(ModelSettings* modelSettings);

  void              mergeCovariance(std::vector<FFTGrid *> & sigma) {state4d_.mergeCov(sigma);}

  void              advanceTime(const int               & previous_vintage,
                                const double            & time_change,
                                SeismicParametersHolder & seismicParameters,
                                ModelSettings           * modelSettings);

  void              setTimeSimbox(Simbox * new_timeSimbox);

  void              setTimeDepthMapping(GridMapping * new_timeDepthMapping);
  void              dump4Dparameters(ModelSettings* modelSettings, std::string identifyer, int timestep,bool printPadding=false);
  void              dumpSeismicParameters(ModelSettings* modelSettings, std::string identifyer, int timestep,SeismicParametersHolder &  current_state);

private:
  void              processWells(std::vector<WellData *> & wells,
                                 Simbox                  * timeSimbox,
                                 ModelSettings          *& modelSettings,
                                 const InputFiles        * inputFiles,
                                 std::string             & errText,
                                 bool                    & failed);

  void              setFaciesNamesFromWells(std::vector<WellData *>        wells,
                                            ModelSettings               *& modelSettings,
                                            std::string                  & tmpErrText,
                                            int                          & error);


  void              setFaciesNamesFromRockPhysics();

  void              setUp3DPartOf4DBackground(const std::vector<DistributionsRock *>           & rock,
                                              const std::vector<float>                         & probability,
                                              const Simbox                                     & timeSimbox,
                                              const ModelSettings                              & modelSettings,
                                              SeismicParametersHolder                          & seismicParameters,
                                              State4D                                          & state4d,
                                              std::string                                      & errTxt);

  void              copyCorrelationsTo4DState(SeismicParametersHolder                    & seismicParameters,
                                              State4D                                    & state4d);

  bool              process4DBackground(ModelSettings           *& modelSettings,
                                        const InputFiles         * inputFiles,
                                        SeismicParametersHolder  & seismicParameters,
                                        std::string              & errText,
                                        bool                     & failed,
                                        NRLib::Vector            & initialMean,
                                        NRLib::Matrix            & initialCov);

  void              calculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
                                                       const std::vector<float>               & probability,
                                                       const std::vector<double>              & trend_position,
                                                       NRLib::Grid2D<double>                  & sigma_sum) const;

  void              makeTimeSimboxes(Simbox          *& timeSimbox,
                                     Simbox          *& timeCutSimbox,
                                     Simbox          *& timeBGSimbox,
                                     Simbox          *& timeSimboxConstThick,
                                     Surface         *& correlationDirection,
                                     ModelSettings   *& modelSettings,
                                     const InputFiles * inputFiles,
                                     std::string      & errText,
                                     bool             & failed);
  void              logIntervalInformation(const Simbox      * simbox,
                                           const std::string & header_text1,
                                           const std::string & header_text2);
  void              setupExtendedTimeSimbox(Simbox  * timeSimbox,
                                            Surface * corrSurf,
                                            Simbox *& timeCutSimbox,
                                            int       outputFormat,
                                            int       outputDomain,
                                            int       otherOutput);
  void              setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
                                                  Surface  * corrSurf,
                                                  Simbox  *& timeBGSimbox,
                                                  int        outputFormat,
                                                  int        outputDomain,
                                                  int        otherOutput);
  void              processDepthConversion(Simbox           * timeCutSimbox,
                                           Simbox           * timeSimbox,
                                           ModelSettings    * modelSettings,
                                           const InputFiles * inputFiles,
                                           std::string      & errText,
                                           bool             & failedVelocity);

  void              processRockPhysics(Simbox                        * timeSimbox,
                                       Simbox                        * timeCutSimbox,
                                       ModelSettings                 * modelSettings,
                                       bool                          & failed,
                                       std::string                   & errTxt,
                                       const std::vector<WellData *> & wells,
                                       const InputFiles              * inputFiles);



  void              printExpectationAndCovariance(const std::vector<double>   & expectation,
                                                  const NRLib::Grid2D<double> & covariance,
                                                  const bool                  & has_trend) const;

  void              setSimboxSurfaces(Simbox                        *& simbox,
                                      const std::vector<std::string> & surfFile,
                                      ModelSettings                  * modelSettings,
                                      std::string                    & errText,
                                      bool                           & failed);
  void              estimateXYPaddingSizes(Simbox         * timeSimbox,
                                           ModelSettings *& modelSettings);
  void              estimateZPaddingSize(Simbox         * timeSimbox,
                                         ModelSettings *& modelSettings);
  int               setPaddingSize(int    nx,
                                   double px);

  void              printSettings(ModelSettings    * modelSettings,
                                  const InputFiles * inputFiles);
  //Compute correlation gradient in terms of i,j and k in grid.
  NRLib::Vector      findPlane(Surface * surf); //Finds plane l2-closest to surface.
  //Create planar surface with same extent as template, p[0]+p[1]*x+p[2]*y
  Surface          * createPlaneSurface(const NRLib::Vector & planeParams,
                                        Surface             * templateSurf);
  void               writeAreas(const SegyGeometry * areaParams,
                                Simbox             * timeSimbox,
                                std::string        & text);
  void               findSmallestSurfaceGeometry(const double   x0,
                                                 const double   y0,
                                                 const double   lx,
                                                 const double   ly,
                                                 const double   rot,
                                                 double       & xMin,
                                                 double       & yMin,
                                                 double       & xMax,
                                                 double       & yMax);
  void              getGeometryFromGridOnFile(const std::string          seismicFile,
                                              const TraceHeaderFormat  * thf,
                                              SegyGeometry            *& geometry,
                                              int                      & fileType,
                                              std::string              & errText);
  SegyGeometry    * geometryFromCravaFile(const std::string & fileName);
  SegyGeometry    * geometryFromStormFile(const std::string & fileName, std::string & errText, bool scale = false);
  //void            processStructureParameters();

  void              estimateCorrXYFromSeismic(Surface *& CorrXY,
                                              FFTGrid ** seisCube,
                                              int numberOfAngles);
  Surface         * findCorrXYGrid(const Simbox * timeSimbox, const ModelSettings * modelSettings);

  int               computeTime(int year, int month, int day) const;

  void              checkFaciesNamesConsistency(ModelSettings     *& modelSettings,
                                                const InputFiles   * inputFiles,
                                                std::string        & tmpErrText) const;

  void              readPriorFaciesProbCubes(const InputFiles        * inputFiles,
                                             ModelSettings           * modelSettings,
                                             std::vector<FFTGrid *>  & priorFaciesProbCubes,
                                             Simbox                  * timeSimbox,
                                             Simbox                  * timeCutSimbox,
                                             std::string             & errTxt,
                                             bool                    & failed);

  void              validateCorrelationMatrix(float              ** C,
                                              const ModelSettings *  modelSettings,
                                              std::string         &  errTxt);
  void              makeCorr2DPositiveDefinite(Surface         * corrXY);


  Simbox                  * timeSimbox_;                 ///< Information about simulation area.
  Simbox                  * timeSimboxConstThick_;       ///< Simbox with constant thickness

  Surface                 * correlationDirection_;       ///< Grid giving the correlation direction.
  RandomGen               * randomGen_;                  ///< Random generator.

  double                    gradX_;                      ///< X-gradient of correlation rotation.
  double                    gradY_;                      ///< Y-gradient of correlation rotation.
                                                         ///< These are only used with correlation surfaces.

  CravaTrend                                                    trend_cubes_;            ///< Trend cubes used in rock phyiscs prior model
  std::map<std::string, std::vector<DistributionsRock *> >      rock_distributions_;     ///< Rocks used in rock physics model
  std::map<std::string, std::vector<DistributionWithTrend *> >  reservoir_variables_;    ///< Reservoir variables used in the rock physics model

  TimeEvolution             timeEvolution_;

  GridMapping             * timeDepthMapping_;           ///< Contains both simbox and mapping used for depth conversion
  GridMapping             * timeCutMapping_;             ///< Simbox and mapping for timeCut

  bool                      velocityFromInversion_;

  bool                      failed_;                     ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;             ///< Detailed failed information.

  TimeLine                * timeLine_;
  std::vector<WellData *>   wells_;                      ///< Well data

  bool                      forwardModeling_;
  int                       numberOfWells_;

  std::vector<float>        priorFacies_;                ///< Prior facies probabilities
  std::vector<FFTGrid *>    priorFaciesProbCubes_;       ///< Cubes for prior facies probabilities

  std::vector<int>          faciesLabels_;               ///< Facies labels, flyttes til blockedlogs
  std::vector<std::string>  faciesNames_;                ///< Facies names   (nFacies = faciesNames.size()). Use for ordering of facies

  bool                      do4DInversion_;
  bool                      do4DRockPhysicsInversion_;
  State4D                   state4d_;                    ///< State4D holds the 27 grdis needed for 4D inversion.

  Surface                 * priorCorrXY_;                ///< Lateral correlation

};

#endif
