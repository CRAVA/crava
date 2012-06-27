#ifndef MODELGENERAL_H
#define MODELGENERAL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/cravatrend.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionsfluid.h"

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
class State4D;
class WellData;

class ModelGeneral
{
public:
  ModelGeneral(ModelSettings    *& modelSettings,
               const InputFiles  * inputFiles,
               Simbox           *& timeBGSimbox);
  ~ModelGeneral();

  Simbox            * getTimeSimbox()            const { return timeSimbox_             ;}
  Simbox            * getTimeSimboxConstThick()  const { return timeSimboxConstThick_   ;}
  RandomGen         * getRandomGen()             const { return randomGen_              ;}
  GridMapping       * getTimeDepthMapping()      const { return timeDepthMapping_       ;}
  GridMapping       * getTimeCutMapping()        const { return timeCutMapping_         ;}
  CravaTrend          getTrendCubes()            const { return trend_cubes_            ;}

  bool                getVelocityFromInversion() const { return velocityFromInversion_  ;}
  bool                getFailed()                const { return failed_                 ;}
  std::vector<bool>   getFailedDetails()         const { return failed_details_         ;}

  void                getCorrGradIJ(float & corrGradI, float &corrGradJ) const;
  Surface           * getCorrelationDirection()  const {return correlationDirection_    ;}

  TimeLine                    * getTimeLine()                  const {return(timeLine_) ;}
  std::vector<WellData *>     & getWells()                 /*const*/ { return wells_    ;}

  std::vector<DistributionsRock *> getRockDistributions(void)  const { return rock_distributions_;}

  const std::vector<float>       & getPriorFacies()           /*const*/ { return priorFacies_            ;}
  const std::vector<FFTGrid *>   & getPriorFaciesCubes()      /*const*/ { return priorFaciesProbCubes_   ;}

  static FFTGrid    * createFFTGrid(int nx,
                                    int ny,
                                    int nz,
                                    int nxp,
                                    int nyp,
                                    int nzp,
                                    bool fileGrid);
  static void         readGridFromFile(const std::string       & fileName,
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
  static void         readSegyFile(const std::string       & fileName,
                                   FFTGrid                *& target,
                                   const Simbox            * timeSimbox,
                                   const Simbox            * timeCutSimbox,
                                   const ModelSettings     * modelSettings,
                                   const SegyGeometry     *& geometry,
                                   int                       gridType,
                                   float                     offset,
                                   const TraceHeaderFormat * format,
                                   std::string             & errText,
                                   bool                      nopadding = false);
  static void         checkThatDataCoverGrid(const SegY   * segy,
                                             float         offset,
                                             const Simbox * timeCutSimbox,
                                             float         guard_zone,
                                             std::string & errText);
  static void         readStormFile(const std::string  & fileName,
                                    FFTGrid           *& target,
                                    const int            gridType,
                                    const std::string  & parName,
                                    const Simbox       * timeSimbox,
                                    const ModelSettings * modelSettings,
                                    std::string        & errText,
                                    bool                 isStorm  = true,
                                    bool                 nopadding = true);
  static void         loadVelocity(FFTGrid           *& velocity,
                                   const Simbox       * timeSimbox,
                                   const Simbox       * timeCutSimbox,
                                   const ModelSettings * modelSettings,
                                   const std::string  & velocityField,
                                   bool               & velocityFromInversion,
                                   std::string        & errText,
                                   bool               & failed);

  void                generateRockPhysics3DBackground(const std::vector<DistributionsRock *> & rock,
                                                      const std::vector<double>             & probability,
                                                      FFTGrid                               & vp,
                                                      FFTGrid                               & vs,
                                                      FFTGrid                               & rho,
                                                      double                                & varVp,
                                                      double                                & varVs,
                                                      double                                & varRho,
                                                      double                                & crossVpVs,
                                                      double                                & crossVpRho,
                                                      double                                & crossVsRho);

  void                generateRockPhysics4DBackground(const std::vector<DistributionsRock *> & rock,
                                                      const std::vector<double>             & probability,
                                                      int                                     lowCut,
                                                      Corr                                  & correlations, //The grids here get/set correctly.
                                                      const Simbox                          & timeSimbox,
                                                      const ModelSettings                   & modelSettings,
                                                      State4D                               & state4d);

  void             processWellLocation(FFTGrid                     ** seisCube,
                                       float                       ** reflectionMatrix,
                                       ModelSettings                * modelSettings,
                                       const std::vector<Surface *> & interval);              // Changes wells

  void             processPriorCorrelations(Corr                 *& correlations,
                                            Background            * background,
                                            std::vector<WellData *> wells,
                                            const Simbox          * timeSimbox,
                                            const ModelSettings   * modelSettings,
                                            FFTGrid              ** seisCube,
                                            const InputFiles      * inputFiles,
                                            std::string           & errText,
                                            bool                  & failed);

   void             processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
                                          std::vector<WellData *>        wells,
                                          Simbox                       * timeSimbox,
                                          Simbox                       * timeCutSimbox,
                                          ModelSettings                * modelSettings,
                                          bool                         & failed,
                                          std::string                  & errTxt,
                                          const InputFiles             * inputFiles);

private:
  void             processWells(std::vector<WellData *> & wells,
                                Simbox                  * timeSimbox,
                                ModelSettings          *& modelSettings,
                                const InputFiles        * inputFiles,
                                std::string             & errText,
                                bool                    & failed);

  void             setFaciesNames(std::vector<WellData *>     wells,
                                  ModelSettings               *& modelSettings,
                                  std::string                  & tmpErrText,
                                  int                          & error);

  void                makeTimeSimboxes(Simbox          *& timeSimbox,
                                       Simbox          *& timeCutSimbox,
                                       Simbox          *& timeBGSimbox,
                                       Simbox          *& timeSimboxConstThick,
                                       Surface         *& correlationDirection,
                                       ModelSettings   *& modelSettings,
                                       const InputFiles * inputFiles,
                                       std::string      & errText,
                                       bool             & failed);
  void                logIntervalInformation(const Simbox      * simbox,
                                             const std::string & header_text1,
                                             const std::string & header_text2);
  void                setupExtendedTimeSimbox(Simbox  * timeSimbox,
                                              Surface * corrSurf,
                                              Simbox *& timeCutSimbox,
                                              int       outputFormat,
                                              int       outputDomain,
                                              int       otherOutput);
  void                setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
                                                    Surface  * corrSurf,
                                                    Simbox  *& timeBGSimbox,
                                                    int        outputFormat,
                                                    int        outputDomain,
                                                    int        otherOutput);
  void                processDepthConversion(Simbox           * timeCutSimbox,
                                             Simbox           * timeSimbox,
                                             ModelSettings    * modelSettings,
                                             const InputFiles * inputFiles,
                                             std::string      & errText,
                                             bool             & failedVelocity);

  void             processRockPhysics(Simbox                       * timeSimbox,
                                      Simbox                       * timeCutSimbox,
                                      ModelSettings                * modelSettings,
                                      bool                         & failed,
                                      std::string                  & errTxt,
                                      const InputFiles             * inputFiles);
  void                setSimboxSurfaces(Simbox                        *& simbox,
                                        const std::vector<std::string> & surfFile,
                                        ModelSettings                  * modelSettings,
                                        std::string                    & errText,
                                        bool                           & failed);
  void                estimateXYPaddingSizes(Simbox         * timeSimbox,
                                             ModelSettings *& modelSettings);
  void                estimateZPaddingSize(Simbox         * timeSimbox,
                                           ModelSettings *& modelSettings);
  int                 setPaddingSize(int    nx,
                                     double px);

  void                printSettings(ModelSettings    * modelSettings,
                                    const InputFiles * inputFiles);
  //Compute correlation gradient in terms of i,j and k in grid.
  double           * findPlane(Surface * surf); //Finds plane l2-closest to surface.
  //Create planar surface with same extent as template, p[0]+p[1]*x+p[2]*y
  Surface          * createPlaneSurface(double  * planeParams,
                                        Surface * templateSurf);
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
  void              getGeometryFromGridOnFile(const std::string         seismicFile,
                                              const TraceHeaderFormat * thf,
                                              SegyGeometry           *& geometry,
                                              std::string             & errText);
  SegyGeometry    * geometryFromCravaFile(const std::string & fileName);
  SegyGeometry    * geometryFromStormFile(const std::string & fileName, std::string & errText, bool scale = false);
  //void             processStructureParameters();

  void             estimateCorrXYFromSeismic(Surface *& CorrXY,
                                             FFTGrid ** seisCube,
                                             int numberOfAngles);
  Surface        * findCorrXYGrid(const Simbox * timeSimbox, const ModelSettings * modelSettings);

  int              computeTime(int year, int month, int day) const;

  void             processFaciesInformation(ModelSettings     *& modelSettings,
                                            const InputFiles   * inputFiles,
                                            std::string        & tmpErrText,
                                            int                & error) const;

  void             readPriorFaciesProbCubes(const InputFiles        * inputFiles,
                                            ModelSettings           * modelSettings,
                                            std::vector<FFTGrid *>  & priorFaciesProbCubes,
                                            Simbox                  * timeSimbox,
                                            Simbox                  * timeCutSimbox,
                                            std::string             & errTxt,
                                            bool                    & failed);

  Simbox            * timeSimbox_;            ///< Information about simulation area.
  Simbox            * timeSimboxConstThick_;  ///< Simbox with constant thickness

  Surface           * correlationDirection_;  ///< Grid giving the correlation direction.
  RandomGen         * randomGen_;             ///< Random generator.

  double              gradX_;                 ///< X-gradient of correlation rotation.
  double              gradY_;                 ///< Y-gradient of correlation rotation.
                                              ///< These are only used with correlation surfaces.


  CravaTrend                       trend_cubes_;        ///< Trend cubes used in rock phyiscs prior model
  std::vector<DistributionsRock *> rock_distributions_; ///< Rock physics prior models
  std::map<std::string, DistributionsFluid *> fluid_distributions_; ///< Fluids used in rock physics model

  GridMapping       * timeDepthMapping_;      ///< Contains both simbox and mapping used for depth conversion
  GridMapping       * timeCutMapping_;        ///< Simbox and mapping for timeCut

  bool                velocityFromInversion_;

  bool                failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>   failed_details_;        ///< Detailed failed information.

  TimeLine                * timeLine_;
  std::vector<WellData *>   wells_;           ///< Well data

  bool                      forwardModeling_;
  int                       numberOfWells_;

  std::vector<float>                  priorFacies_;           ///< Prior facies probabilities
  std::vector<FFTGrid *>              priorFaciesProbCubes_;  ///< Cubes for prior facies probabilities

};

#endif
