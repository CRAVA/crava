#ifndef MODEL_H
#define MODEL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "lib/global_def.h"
#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"

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

class Model
{
public:
  Model(char * fileName);
  ~Model();

  ModelSettings  * getModelSettings()         const { return modelSettings_          ;}
  Simbox         * getTimeSimbox()            const { return timeSimbox_             ;}
  Simbox         * getTimeSimboxConstThick()  const { return timeSimboxConstThick_   ;}
  WellData      ** getWells()                 const { return wells_                  ;}
  FFTGrid        * getBackAlpha()             const { return background_->getAlpha() ;}
  FFTGrid        * getBackBeta()              const { return background_->getBeta()  ;}
  FFTGrid        * getBackRho()               const { return background_->getRho()   ;}
  Corr           * getCorrelations()          const { return correlations_           ;}
  float          * getPriorFacies()           const { return priorFacies_            ;}  
  FFTGrid       ** getSeisCubes()             const { return seisCube_               ;}
  Wavelet       ** getWavelets()              const { return wavelet_                ;}
  float         ** getAMatrix()               const { return reflectionMatrix_       ;}
  RandomGen      * getRandomGen()             const { return randomGen_              ;} 
  GridMapping    * getTimeDepthMapping()      const { return timeDepthMapping_       ;}
  GridMapping    * getTimeCutMapping()        const { return timeCutMapping_         ;}
  bool             getVelocityFromInversion() const { return velocityFromInversion_  ;}
  bool             getFailed()                const { return failed_                 ;}
  void             releaseGrids();                                        // Cuts connection to SeisCube_ and  backModel_
  void             getCorrGradIJ(float & corrGradI, float &corrGradJ) const;
  static void      writeSurfaceToFile(Surface           * surface,
                                      const std::string & name,
                                      int                 format);

  enum             backFileTypes{STORMFILE = -2, SEGYFILE = -1};

private:
  void             makeTimeSimboxes(Simbox        *& timeSimbox,
                                    Simbox        *& timeCutSimbox,
                                    Simbox        *& timeBGSimbox,
                                    Simbox        *& timeSimboxConstThick,
                                    Surface       *& correlationDirection,
                                    ModelSettings *& modelSettings, 
                                    InputFiles     * inputFiles,
                                    char           * errText,
                                    bool           & failed);
  void             setupExtendedTimeSimbox(Simbox  * timeSimbox, 
                                           Surface * corrSurf, 
                                           Simbox *& timeCutSimbox,
                                           int       outputFormat,
                                           int       outputFlag);
  void             setupExtendedBackgroundSimbox(Simbox   * timeSimbox, 
                                                 Surface  * corrSurf, 
                                                 Simbox  *& timeBGSimbox,
                                                 int        outputFormat,
                                                 int        outputFlag);
  void             processSeismic(FFTGrid      **& seisCube,
                                  Simbox        *& timeSimbox,
                                  ModelSettings *& modelSettings, 
                                  InputFiles     * inputFiles,
                                  char           * errText,
                                  bool           & failed);
  void             processWells(WellData     **& wells,
                                Simbox         * timeSimbox,
                                Simbox         * timeBGSimbox,
                                Simbox         * timeSimboxConstThick,
                                FFTGrid       ** seisCube,
                                RandomGen      * randomGen,
                                ModelSettings *& modelSettings, 
                                InputFiles     * inputFiles,
                                char           * errText,
                                bool           & failed);
  void             processBackground(Background   *& background, 
                                     WellData     ** wells,
                                     Simbox        * timeSimbox,
                                     Simbox        * timeBGSimbox,
                                     ModelSettings * modelSettings, 
                                     InputFiles    * inputFile,
                                     char          * errText,
                                     bool          & failed);
  void             processPriorCorrelations(Corr         *& correlations,
                                            Background    * background,
                                            WellData     ** wells,
                                            Simbox        * timeSimbox,
                                            ModelSettings * modelSettings, 
                                            InputFiles    * inputFiles,
                                            char          * errText,
                                            bool          & failed);
  void             processReflectionMatrix(float       **& reflectionMatrix,
                                           Background    * background,
                                           ModelSettings * modelSettings, 
                                           InputFiles    * inputfiles,
                                           char          * errText,
                                           bool          & failed);
  void             processWavelets(Wavelet     **& wavelet,
                                   FFTGrid      ** seisCube,
                                   WellData     ** wells,
                                   float        ** reflectionMatrix,
                                   Simbox        * timeSimbox,
                                   Surface      ** waveletEstimInterval,
                                   Surface      ** shiftGrids,
                                   Surface      ** gainGrids,
                                   ModelSettings * modelSettings, 
                                   InputFiles    * inputFiles,
                                   char          * errText,
                                   bool          & failed);
  void             processPriorFaciesProb(float        *& priorFacies,
                                          WellData     ** wells,
                                          RandomGen     * randomGen,
                                          int             nz,
                                          ModelSettings * modelSettings);
  void             processDepthConversion(Simbox        * timeCutSimbox, 
                                          Simbox        * timeSimbox_,
                                          ModelSettings * modelSettings_, 
                                          InputFiles    * inputFiles,
                                          char          * errText, 
                                          bool          & failedVelocity);
  void             loadVelocity(FFTGrid           *& velocity,
                                Simbox             * timeSimbox,
                                ModelSettings      * modelSettings, 
                                const std::string  & velocityField, 
                                char               * errText,
                                bool               & failed);
  void             setSimboxSurfaces(Simbox                        *& simbox, 
                                     const std::vector<std::string> & surfFile, 
                                     bool                             parallelSurfaces, 
                                     double                           dTop,
                                     double                           lz, 
                                     double                           dz, 
                                     int                              nz,
                                     int                              outputFormat,
                                     int                              outputFlag,
                                     char                           * errText,
                                     int                            & error);
  void             estimateXYPaddingSizes(Simbox         * timeSimbox,
                                          ModelSettings *& modelSettings);
  void             estimateZPaddingSize(Simbox         * timeSimbox,
                                        ModelSettings *& modelSettings);
  int              readSegyFile(const std::string   & fileName, 
                                FFTGrid            *& target, 
                                Simbox             *& timeSimbox, 
                                ModelSettings      *& modelSettings,
                                char                * errText,
                                const SegyGeometry *& geometry,
                                int                   gridType,
                                int                   i = IMISSING); 
  int              readStormFile(const std::string  & fName, 
                                 FFTGrid           *& target, 
                                 const char         * parName,
                                 Simbox             * timeSimbox,
                                 ModelSettings     *& modelSettings, 
                                 char               * errText);
  void             estimateCorrXYFromSeismic(Surface *& CorrXY,
                                             FFTGrid ** seisCube,
                                             int        nAngles);
  Surface        * findCorrXYGrid(ModelSettings * modelSettings);
  int              setPaddingSize(int    nx, 
                                  double px);
  void             loadExtraSurfaces(Surface  **& waveletEstimInterval,
                                     Surface  **& faciesEstimInterval,
                                     Simbox     * timeSimbox,
                                     InputFiles * inputFiles,
                                     char       * errText,
                                     bool       & failed);
  float         ** readMatrix(const std::string & fileName, 
                              int                 n1, 
                              int                 n2, 
                              const std::string & readReason, 
                              char              * errText);
  void             setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                                Background    * background,
                                                ModelSettings * modelSettings);
  int              findFileType(const std::string & fileName);
  void             checkAvailableMemory(Simbox            * timeSimbox,
                                        ModelSettings     * modelSettings,
                                        const std::string & seismicFile);
  void             checkFaciesNames(WellData      ** wells,
                                    ModelSettings *& modelSettings,
                                    char           * tmpErrText,
                                    int            & error);
  void             printSettings(ModelSettings * modelSettings,
                                 InputFiles    * inputFiles);
  int              getWaveletFileFormat(const std::string & fileName, 
                                        char              * errText);
  //Compute correlation gradient in terms of i,j and k in grid.
  double *         findPlane(Surface * surf); //Finds plane l2-closest to surface.             
  //Create planar surface with same extent as template, p[0]+p[1]*x+p[2]*y
  Surface *        createPlaneSurface(double  * planeParams, 
                                      Surface * templateSurf);
  void             writeAreas(const SegyGeometry * areaParams,
                              Simbox             * timeSimbox,
                              std::string        & text);
  void             findSmallestSurfaceGeometry(const double   x0, 
                                               const double   y0, 
                                               const double   lx, 
                                               const double   ly, 
                                               const double   rot,
                                               double       & xMin,
                                               double       & yMin,
                                               double       & xMax,
                                               double       & yMax); 
  ModelSettings  * modelSettings_;
  Simbox         * timeSimbox_;            ///< Information about simulation area.
  Simbox         * timeSimboxConstThick_;  ///< Simbox with constant thickness   
  
  WellData      ** wells_;                 ///< Well data
  Background     * background_;            ///< Holds the background model.
  Corr           * correlations_;          ///<
  FFTGrid       ** seisCube_;              ///< Seismic data cubes
  Wavelet       ** wavelet_;               ///< Wavelet for angle
  Surface       ** shiftGrids_;            ///< Grids containing shift data for wavelets
  Surface       ** gainGrids_;             ///< Grids containing gain data for wavelets.
  Surface       ** waveletEstimInterval_;  ///< Grids giving the wavelet estimation interval.
  Surface       ** faciesEstimInterval_;   ///< Grids giving the facies estimation intervals.
  Surface        * correlationDirection_;  ///< Grid giving the correlation direction.
  RandomGen      * randomGen_;             ///< Random generator.
  float          * priorFacies_;
  float         ** reflectionMatrix_;      ///< May specify own Zoeppritz-approximation. Default NULL,
                                           ///< indicating that standard approximation will be used.

  double           gradX_;                 ///< X-gradient of correlation rotation. 
  double           gradY_;                 ///< Y-gradient of correlation rotation.
                                           ///< These are only used with correaltion surfaces.

  bool             failed_;                ///< Indicates whether errors occured during construction. 
  GridMapping    * timeDepthMapping_;      ///< Contains both simbox and mapping used for depth conversion
  GridMapping    * timeCutMapping_;        ///< Simbox and mapping for timeCut

  bool             velocityFromInversion_;
};

#endif

