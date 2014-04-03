/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELGENERAL_H
#define MODELGENERAL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "src/definitions.h"
//#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/cravatrend.h"
#include "src/seismicparametersholder.h"
#include "src/state4d.h"
#include "src/timeevolution.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionssolid.h"

#include "nrlib/grid/grid.hpp"
#include "src/blockedlogscommon.h"
#include "src/commondata.h"

struct irapgrid;
class Corr;
class Wavelet;
class Vario;
class Simbox;
//class WellData;
class FFTGrid;
class RandomGen;
class GridMapping;
class InputFiles;
class TimeLine;
//class WellData;
class SeismicParameters;
class CravaTrend;

//class MultiIntervalGrid;

class ModelGeneral
{
public:
  //ModelGeneral(ModelSettings           *& modelSettings,
  //             const InputFiles         * inputFiles,
  //             SeismicParametersHolder  & seismicParameters,
  //             Simbox                  *& timeBGSimbox);

  ModelGeneral(ModelSettings           *& modelSettings, //Multiple intervals
               const InputFiles         * inputFiles,
               SeismicParametersHolder  & seismicParameters,
               CommonData               * common_data,
               int                        i_interval);

  ~ModelGeneral();

  const Simbox                             * GetTimeSimbox()            const { return simbox_                  ;}
  RandomGen                                * GetRandomGen()             const { return random_gen_              ;}
  GridMapping                              * GetTimeDepthMapping()      const { return time_depth_mapping_      ;}
  CravaTrend                               & GetTrendCubes()                  { return trend_cubes_             ;}
  CravaTrend                                 GetTrendCubes()            const { return trend_cubes_             ;}
  bool                                       GetVelocityFromInversion() const { return velocity_from_inversion_ ;}
  bool                                       GetMultiInterval()         const { return multi_interval_          ;}
  State4D                                    GetState4D()               const { return state4d_                 ;}
  TimeLine                                 * GetTimeLine()              const { return time_line_               ;}
  std::map<std::string, BlockedLogsCommon *> GetBlockedWells()                { return blocked_logs_            ;}
  const std::vector<float>                 & GetPriorFacies()       /*const*/ { return prior_facies_            ;}
  const std::vector<FFTGrid *>             & GetPriorFaciesCubes()  /*const*/ { return prior_facies_prob_cubes_ ;}
  const std::vector<std::string>           & GetFaciesNames(void)       const { return facies_names_            ;}
  std::vector<int>                           GetFaciesLabel()           const { return facies_labels_           ;}
  bool                                       GetIs4DActive()            const { return(do_4D_inversion_)        ;}

  void AddFaciesLabel(int faciesLabel)                                        { facies_labels_.push_back(faciesLabel) ;}
  void AddFaciesName(const std::string & faciesName)                          { facies_names_.push_back(faciesName)   ;}

  //Simbox                   * getTimeSimboxConstThick()  const { return timeSimboxConstThick_   ;}
  //GridMapping              * getTimeCutMapping()        const { return timeCutMapping_         ;}
  //Surface                  * getPriorCorrXY()           const { return priorCorrXY_            ;}
  //bool                       getFailed()                const { return failed_                 ;}
  //std::vector<bool>          getFailedDetails()         const { return failed_details_         ;}
  //void                       getCorrGradIJ(float & corrGradI, float &corrGradJ) const;
  //Surface                  * getCorrelationDirection()  const { return correlationDirection_   ;}
  //std::vector<WellData *>  & getWells()             /*const*/ { return wells_                  ;}

  //static void                readSegyFile(const std::string       & fileName,
  //                                        FFTGrid                *& target,
  //                                        Simbox                  * timeSimbox,
  //                                        Simbox                  * timeCutSimbox,
  //                                        ModelSettings          *& modelSettings,
  //                                        const SegyGeometry     *& geometry,
  //                                        int                       gridType,
  //                                        const std::string       & parName,
  //                                        float                     offset,
  //                                        const TraceHeaderFormat * format,
  //                                        std::string             & errText,
  //                                        bool                      nopadding = false);
  //static void                checkThatDataCoverGrid(SegY        * segy,
  //                                                  float         offset,
  //                                                  Simbox      * timeCutSimbox,
  //                                                  float         guard_zone,
  //                                                  std::string & errText);
  //static void                readStormFile(const std::string  & fileName,
  //                                         FFTGrid           *& target,
  //                                         const int            gridType,
  //                                         const std::string  & parName,
  //                                         Simbox             * timeSimbox,
  //                                         ModelSettings     *& modelSettings,
  //                                         std::string        & errText,
  //                                         bool                 isStorm  = true,
  //                                         bool                 nopadding = true);

  std::map<std::string, DistributionsRock *> GetRockDistributionTime0() const;

  static FFTGrid  * CreateFFTGrid(int nx,
                                  int ny,
                                  int nz,
                                  int nxp,
                                  int nyp,
                                  int nzp,
                                  bool fileGrid);

  //static void       readGridFromFile(const std::string       & fileName,
  //                                   const std::string       & parName,
  //                                   const float               offset,
  //                                   FFTGrid                *& grid,
  //                                   const SegyGeometry     *& geometry,
  //                                   const TraceHeaderFormat * format,
  //                                   int                       gridType,
  //                                   const Simbox            * timeSimbox,
  //                                   const Simbox            * timeCutSimbox,
  //                                   const ModelSettings     * modelSettings,
  //                                   std::string             & errorText,
  //                                   bool                      nopadding = false);

  //static void       readSegyFile(const std::string       & fileName,
  //                               FFTGrid                *& target,
  //                               const Simbox            * timeSimbox,
  //                               const Simbox            * timeCutSimbox,
  //                               const ModelSettings     * modelSettings,
  //                               const SegyGeometry     *& geometry,
  //                               int                       gridType,
  //                               const std::string       & parName,
  //                               float                     offset,
  //                               const TraceHeaderFormat * format,
  //                               std::string             & errText,
  //                               bool                      nopadding = false);

  //static void       checkThatDataCoverGrid(const SegY   * segy,
  //                                         float         offset,
  //                                         const Simbox * timeCutSimbox,
  //                                         float         guard_zone,
  //                                         std::string & errText);
  //static void       readStormFile(const std::string  & fileName,
  //                                FFTGrid           *& target,
  //                                const int            gridType,
  //                                const std::string  & parName,
  //                                const Simbox       * timeSimbox,
  //                                const ModelSettings * modelSettings,
  //                                std::string        & errText,
  //                                bool                 isStorm  = true,
  //                                bool                 nopadding = true);
  //static void       loadVelocity(FFTGrid           *& velocity,
  //                               const Simbox       * timeSimbox,
  //                               const Simbox       * timeCutSimbox,
  //                               const ModelSettings * modelSettings,
  //                               const std::string  & velocityField,
  //                               bool               & velocityFromInversion,
  //                               std::string        & errText,
  //                               bool               & failed);

  //void              processWellLocation(FFTGrid                     ** seisCube,
  //                                      float                       ** reflectionMatrix,
  //                                      ModelSettings                * modelSettings,
  //                                      const std::vector<Surface *> & interval);              // Changes wells

  //void              processPriorCorrelations(Background                     * background,
  //                                           std::vector<WellData *>          wells,
  //                                           const Simbox                   * timeSimbox,
  //                                           const ModelSettings            * modelSettings,
  //                                           const std::vector<float>       & priorFacies,
  //                                           FFTGrid                       ** seisCube,
  //                                           const InputFiles               * inputFiles,
  //                                           SeismicParametersHolder        & seismicParameters,
  //                                           std::string                    & errText,
  //                                           bool                           & failed);

   //void             processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
   //                                       std::vector<WellData *>        wells,
   //                                       Simbox                       * timeSimbox,
   //                                       Simbox                       * timeCutSimbox,
   //                                       ModelSettings                * modelSettings,
   //                                       bool                         & failed,
   //                                       std::string                  & errTxt,
   //                                       const InputFiles             * inputFiles);

  //void              generateRockPhysics3DBackground(const std::vector<DistributionsRock *>           & rock_distribution,
  //                                                  const std::vector<float>                         & probability,
  //                                                  FFTGrid                                          & vp,
  //                                                  FFTGrid                                          & vs,
  //                                                  FFTGrid                                          & rho);

  void              CalculateCovariancesFromRockPhysics(const std::vector<DistributionsRock *>           & rock,
                                                        const std::vector<float>                         & probability,
                                                        NRLib::Grid2D<double>                            & param_corr,
                                                        std::string                                      & errTxt);

  void              Complete4DBackground(const int nx,const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad,NRLib::Vector &initial_mean,NRLib::Matrix &initial_cov);

  //void              getInitial3DPriorFrom4D(SeismicParametersHolder & seismicParameters);
  bool              Do4DRockPhysicsInversion(ModelSettings* modelSettings);

  void              MergeCovariance(std::vector<FFTGrid *> & sigma) {state4d_.mergeCov(sigma);}

  void              AdvanceTime(int time_step, SeismicParametersHolder & seismicParameters,ModelSettings* modelSettings);
  void              LastUpdateOfStaticAndDynamicParts(SeismicParametersHolder &  seismicParameters,ModelSettings* modelSettings);
  void              Dump4Dparameters(ModelSettings* modelSettings, std::string identifyer, int timestep);
  void              DumpSeismicParameters(ModelSettings* modelSettings, std::string identifyer, int timestep,SeismicParametersHolder &  current_state);

private:
  //void              processWells(std::vector<WellData *> & wells,
  //                               Simbox                  * timeSimbox,
  //                               ModelSettings          *& modelSettings,
  //                               const InputFiles        * inputFiles,
  //                               std::string             & errText,
  //                               bool                    & failed);

  //void              setFaciesNamesFromWells(std::vector<WellData *>        wells,
  //                                          ModelSettings               *& modelSettings,
  //                                          std::string                  & tmpErrText,
  //                                          int                          & error);

  //void              setFaciesNamesFromRockPhysics();

  //void              setUp3DPartOf4DBackground(const std::vector<DistributionsRock *>           & rock,
  //                                            const std::vector<float>                         & probability,
  //                                            const Simbox                                     & timeSimbox,
  //                                            const ModelSettings                              & modelSettings,
  //                                            SeismicParametersHolder                          & seismicParameters,
  //                                            State4D                                          & state4d,
  //                                            std::string                                      & errTxt);

  void              CopyCorrelationsTo4DState(SeismicParametersHolder                    & seismicParameters,
                                              State4D                                    & state4d);

  //bool              process4DBackground(ModelSettings           *& modelSettings,
  //                                      const InputFiles         * inputFiles,
  //                                      SeismicParametersHolder  & seismicParameters,
  //                                      std::string              & errText,
  //                                      bool                     & failed,
  //                                      NRLib::Vector            & initialMean,
  //                                      NRLib::Matrix            & initialCov);

  void              SetupState4D(SeismicParametersHolder & seismicParameters,
                                 const Simbox            * simbox,
                                 State4D                 & state4d,
                                 NRLib::Vector           & initialMean,
                                 NRLib::Matrix           & initialCov);

  void              CalculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
                                                       const std::vector<float>               & probability,
                                                       const std::vector<double>              & trend_position,
                                                       NRLib::Grid2D<double>                  & sigma_sum) const;

  //void              makeTimeSimboxes(Simbox          *& timeSimbox,
  //                                   Simbox          *& timeCutSimbox,
  //                                   Simbox          *& timeBGSimbox,
  //                                   Simbox          *& timeSimboxConstThick,
  //                                   Surface         *& correlationDirection,
  //                                   ModelSettings   *& modelSettings,
  //                                   const InputFiles * inputFiles,
  //                                   std::string      & errText,
  //                                   bool             & failed);
  //void              logIntervalInformation(const Simbox      * simbox,
  //                                         const std::string & header_text1,
  //                                         const std::string & header_text2);
  //void              setupExtendedTimeSimbox(Simbox  * timeSimbox,
  //                                          Surface * corrSurf,
  //                                          Simbox *& timeCutSimbox,
  //                                          int       outputFormat,
  //                                          int       outputDomain,
  //                                          int       otherOutput);
  //void              setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
  //                                                Surface  * corrSurf,
  //                                                Simbox  *& timeBGSimbox,
  //                                                int        outputFormat,
  //                                                int        outputDomain,
  //                                                int        otherOutput);
  //void              processDepthConversion(Simbox           * timeCutSimbox,
  //                                         Simbox           * timeSimbox,
  //                                         ModelSettings    * modelSettings,
  //                                         const InputFiles * inputFiles,
  //                                         std::string      & errText,
  //                                         bool             & failedVelocity);

  //void              processRockPhysics(Simbox                       * timeSimbox,
  //                                     Simbox                       * timeCutSimbox,
  //                                     ModelSettings                * modelSettings,
  //                                     bool                         & failed,
  //                                     std::string                  & errTxt,
  //                                     const InputFiles             * inputFiles);

  //void              printExpectationAndCovariance(const std::vector<double>   & expectation,
  //                                                const NRLib::Grid2D<double> & covariance,
  //                                                const bool                  & has_trend) const;

  //void              SetSimboxSurfaces(Simbox                        *& simbox,
  //                                    const std::vector<std::string> & surfFile,
  //                                    ModelSettings                  * modelSettings,
  //                                    std::string                    & errText,
  //                                    bool                           & failed);
  //void              EstimateXYPaddingSizes(Simbox         * timeSimbox,
  //                                         ModelSettings *& modelSettings);
  //void              EstimateZPaddingSize(Simbox         * timeSimbox,
  //                                       ModelSettings *& modelSettings);
  //int               SetPaddingSize(int    nx,
  //                                 double px);

  //void              PrintSettings(ModelSettings    * modelSettings,
  //                                const InputFiles * inputFiles);
  //Compute correlation gradient in terms of i,j and k in grid.
  //NRLib::Vector      FindPlane(Surface * surf); //Finds plane l2-closest to surface.
  //Create planar surface with same extent as template, p[0]+p[1]*x+p[2]*y
  //Surface          * CreatePlaneSurface(const NRLib::Vector & planeParams,
  //                                      Surface             * templateSurf);
  //void               WriteAreas(const SegyGeometry * areaParams,
  //                              Simbox             * timeSimbox,
  //                              std::string        & text);
  //void               FindSmallestSurfaceGeometry(const double   x0,
  //                                               const double   y0,
  //                                               const double   lx,
  //                                               const double   ly,
  //                                               const double   rot,
  //                                               double       & xMin,
  //                                               double       & yMin,
  //                                               double       & xMax,
  //                                               double       & yMax);
  //void              GetGeometryFromGridOnFile(const std::string         seismicFile,
  //                                            const TraceHeaderFormat * thf,
  //                                            SegyGeometry           *& geometry,
  //                                            std::string             & errText);
  //SegyGeometry    * GeometryFromCravaFile(const std::string & fileName);
  //SegyGeometry    * GeometryFromStormFile(const std::string & fileName, std::string & errText, bool scale = false);
  //void            processStructureParameters();

  void              EstimateCorrXYFromSeismic(Surface *& CorrXY,
                                              FFTGrid ** seisCube,
                                              int numberOfAngles);
  //Surface         * findCorrXYGrid(const Simbox * timeSimbox, const ModelSettings * modelSettings);

  //int               ComputeTime(int year, int month, int day) const;

  //void              CheckFaciesNamesConsistency(ModelSettings     *& modelSettings,
  //                                              const InputFiles   * inputFiles,
  //                                              std::string        & tmpErrText) const;

  //void              readPriorFaciesProbCubes(const InputFiles        * inputFiles,
  //                                           ModelSettings           * modelSettings,
  //                                           std::vector<FFTGrid *>  & priorFaciesProbCubes,
  //                                           Simbox                  * timeSimbox,
  //                                           Simbox                  * timeCutSimbox,
  //                                           std::string             & errTxt,
  //                                           bool                    & failed);

  void              ValidateCorrelationMatrix(float              ** C,
                                              const ModelSettings *  modelSettings,
                                              std::string         &  errTxt);
  void              MakeCorr2DPositiveDefinite(Surface         * corrXY);


  const Simbox                                                * simbox_;                 ///< Information about simulation area.

  RandomGen                                                   * random_gen_;             ///< Random generator.

  CravaTrend                                                    trend_cubes_;            ///< Trend cubes used in rock phyiscs prior model
  std::map<std::string, std::vector<DistributionsRock *> >      rock_distributions_;     ///< Rocks used in rock physics model
  std::map<std::string, std::vector<DistributionWithTrend *> >  reservoir_variables_;    ///< Reservoir variables used in the rock physics model

  TimeEvolution                                                 time_evolution_;

  GridMapping                                                 * time_depth_mapping_;           ///< Contains both simbox and mapping used for depth conversion

  bool                                                          velocity_from_inversion_;

  TimeLine                                                    * time_line_;

  std::map<std::string, BlockedLogsCommon *>                    blocked_logs_;

  std::vector<float>                                            prior_facies_;                ///< Prior facies probabilities
  std::vector<FFTGrid *>                                        prior_facies_prob_cubes_;     ///< Cubes for prior facies probabilities

  std::vector<int>                                              facies_labels_;               ///< Facies labels, flyttes til blockedlogs
  std::vector<std::string>                                      facies_names_;                ///< Facies names   (nFacies = faciesNames.size()). Use for ordering of facies

  bool                                                          do_4D_inversion_;
  bool                                                          do_4D_rock_physics_vnversion_;
  State4D                                                       state4d_;                     ///< State4D holds the 27 grdis needed for 4D inversion.

  bool                                                          multi_interval_;              ///< True if there is multiple intervals


  //Fjernes:
  //int                       numberOfWells_;

  //Simbox                  * timeSimboxConstThick_;       ///< Simbox with constant thickness  //Fjernes

  //Surface                 * correlationDirection_;       ///< Grid giving the correlation direction. //Fjernes

  //double                    gradX_;                      ///< X-gradient of correlation rotation. //Fjernes
  //double                    gradY_;                      ///< Y-gradient of correlation rotation. //Fjernes
  //                                                       ///< These are only used with correlation surfaces.

  //GridMapping             * timeCutMapping_;             ///< Simbox and mapping for timeCut  //Fjernes

  //bool                      failed_;                     ///< Indicates whether errors occured during construction. //Fjernes
  //std::vector<bool>         failed_details_;             ///< Detailed failed information. //Fjernes

  //std::vector<WellData *>   wells_;                      ///< Well data //Fjernes

  //bool                      forwardModeling_; //Flyttes til ModelAvoStatic

  //Surface                 * priorCorrXY_;                ///< Lateral correlation //Fjernes

};

#endif
