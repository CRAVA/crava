#ifndef MODELGENERAL_H
#define MODELGENERAL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"

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

class ModelGeneral
{
public:
  ModelGeneral(ModelSettings *& modelSettings, const InputFiles * inputFiles, Simbox *& timeBGSimbox);
  ~ModelGeneral();

  Simbox                      * getTimeSimbox()            const { return timeSimbox_             ;}
  Simbox                      * getTimeSimboxConstThick()  const { return timeSimboxConstThick_   ;}
  RandomGen                   * getRandomGen()             const { return randomGen_              ;}
  GridMapping                 * getTimeDepthMapping()      const { return timeDepthMapping_       ;}
  GridMapping                 * getTimeCutMapping()        const { return timeCutMapping_         ;}

  bool                          getVelocityFromInversion() const { return velocityFromInversion_  ;}
  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  void                          getCorrGradIJ(float & corrGradI, float &corrGradJ) const;

  static FFTGrid*               createFFTGrid(int nx,
                                              int ny,
                                              int nz,
                                              int nxp,
                                              int nyp,
                                              int nzp,
                                              bool fileGrid);
  static void      readGridFromFile(const std::string       & fileName,
                                    const std::string       & parName,
                                    const float               offset,
                                    FFTGrid                *& grid,
                                    const SegyGeometry     *& geometry,
                                    const TraceHeaderFormat * format,
                                    int                       gridType,
                                    Simbox                  * timeSimbox,
                                    ModelSettings           * modelSettings,
                                    std::string             & errorText,
                                    bool                      nopadding = false);
  static void      readSegyFile(const std::string       & fileName,
                                FFTGrid                *& target,
                                Simbox                 *& timeSimbox,
                                ModelSettings          *& modelSettings,
                                const SegyGeometry     *& geometry,
                                int                       gridType,
                                float                     offset,
                                const TraceHeaderFormat * format,
                                std::string             & errText,
                                bool                      nopadding = false);
  static void      readStormFile(const std::string  & fileName,
                                 FFTGrid           *& target,
                                 const int            gridType,
                                 const std::string  & parName,
                                 Simbox             * timeSimbox,
                                 ModelSettings     *& modelSettings,
                                 std::string        & errText,
                                 bool                 isStorm  = true,
                                 bool                 nopadding = true);
  static void      loadVelocity(FFTGrid           *& velocity,
                                Simbox             * timeSimbox,
                                ModelSettings      * modelSettings,
                                const std::string  & velocityField,
                                bool               & velocityFromInversion,
                                std::string        & errText,
                                bool               & failed);

private:
  void             makeTimeSimboxes(Simbox          *& timeSimbox,
                                    Simbox          *& timeCutSimbox,
                                    Simbox          *& timeBGSimbox,
                                    Simbox          *& timeSimboxConstThick,
                                    Surface         *& correlationDirection,
                                    ModelSettings   *& modelSettings,
                                    const InputFiles * inputFiles,
                                    std::string      & errText,
                                    bool             & failed);
  void             logIntervalInformation(const Simbox      * simbox,
                                          const std::string & header_text1,
                                          const std::string & header_text2);
  void             setupExtendedTimeSimbox(Simbox  * timeSimbox,
                                           Surface * corrSurf,
                                           Simbox *& timeCutSimbox,
                                           int       outputFormat,
                                           int       outputDomain,
                                           int       otherOutput);
  void             setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
                                                 Surface  * corrSurf,
                                                 Simbox  *& timeBGSimbox,
                                                 int        outputFormat,
                                                 int        outputDomain,
                                                 int        otherOutput);
  void             processDepthConversion(Simbox           * timeCutSimbox,
                                          Simbox           * timeSimbox,
                                          ModelSettings    * modelSettings,
                                          const InputFiles * inputFiles,
                                          std::string      & errText,
                                          bool             & failedVelocity);
  void             setSimboxSurfaces(Simbox                        *& simbox,
                                     const std::vector<std::string> & surfFile,
                                     ModelSettings                  * modelSettings,
                                     std::string                    & errText,
                                     bool                           & failed);
  void             estimateXYPaddingSizes(Simbox         * timeSimbox,
                                          ModelSettings *& modelSettings);
  void             estimateZPaddingSize(Simbox         * timeSimbox,
                                        ModelSettings *& modelSettings);
  int              setPaddingSize(int    nx,
                                  double px);

  void             printSettings(ModelSettings       * modelSettings,
                                 const InputFiles    * inputFiles);
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
  void             getGeometryFromGridOnFile(const std::string         seismicFile,
                                             const TraceHeaderFormat * thf,
                                             SegyGeometry           *& geometry,
                                             std::string             & errText);
  SegyGeometry   * geometryFromCravaFile(const std::string & fileName);
  SegyGeometry   * geometryFromStormFile(const std::string & fileName, std::string & errText, bool scale = false);

  Simbox                  * timeSimbox_;            ///< Information about simulation area.
  Simbox                  * timeSimboxConstThick_;  ///< Simbox with constant thickness

  Surface                 * correlationDirection_;  ///< Grid giving the correlation direction.
  RandomGen               * randomGen_;             ///< Random generator.

  double                    gradX_;                 ///< X-gradient of correlation rotation.
  double                    gradY_;                 ///< Y-gradient of correlation rotation.
                                                    ///< These are only used with correlation surfaces.

  GridMapping             * timeDepthMapping_;      ///< Contains both simbox and mapping used for depth conversion
  GridMapping             * timeCutMapping_;        ///< Simbox and mapping for timeCut

  bool                      velocityFromInversion_;

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.
};

#endif
