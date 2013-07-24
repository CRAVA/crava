/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELAVODYNAMIC_H
#define MODELAVODYNAMIC_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/vario.h"
#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"

struct irapgrid;
class Wavelet;
class Simbox;
class WellData;
class FFTGrid;
class RandomGen;
class GridMapping;
class InputFiles;
class ModelAVOStatic;
class ModelGeneral;
class SeismicParametersHolder;

class ModelAVODynamic
{
public:
  ModelAVODynamic(ModelSettings               *& modelSettings,
                  const InputFiles             * inputFiles,
                  std::vector<bool>              failedGeneralDetails,
                  std::vector<bool>              failedStaticDetails,
                  const Simbox                 * timeSimbox,
                  const Simbox                 * timeBGSimbox,
                  const Surface                * correlationDirection,
                  RandomGen                    * /*randomGen*/,
                  GridMapping                  * timeDepthMapping,
                  const GridMapping            * timeCutMapping,
                  const std::vector<Surface *> & waveletEstimInterval,
                  const std::vector<Surface *> & /*wellMoveInterval*/,
                  const std::vector<Surface *> & /*faciesEstimInterval*/,
                  ModelAVOStatic               * modelAVOstatic,
                  ModelGeneral                 * modelGeneral,
                  int                            t,
                  SeismicParametersHolder      & seismicParameters);    // modelAVOstatic::wells_ are altered. modelAVOstatic is deliberately sent in as un-const.

  ModelAVODynamic(ModelSettings          *& modelSettings,
                  const InputFiles        * inputFiles,
                  ModelAVOStatic          * modelAVOstatic,
                  ModelGeneral            * modelGeneral,
                  SeismicParametersHolder & seismicParameters,
                  const Simbox            * timeSimbox,
                  const Surface           * correlationDirection,
                  const GridMapping       * timeDepthMapping,
                  const GridMapping       * timeCutMapping,
                  int                       t);

  ~ModelAVODynamic();

  FFTGrid                    ** getSeisCubes()             const { return seisCube_               ;}
  Wavelet                    ** getWavelets()              const { return wavelet_                ;}
  float                      ** getAMatrix()               const { return reflectionMatrix_       ;}
  Grid2D                      * getLocalNoiseScale(int i)  const { return localNoiseScale_[i]     ;}
  const std::vector<Grid2D *> & getLocalNoiseScales()      const { return localNoiseScale_        ;}

  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  Vario                       * getAngularCorr()           const { return angularCorr_                    ;}
  float                         getSNRatio(int i)          const { return SNRatio_[i]                     ;}
  bool                          getUseLocalNoise()         const { return useLocalNoise_                  ;}
  float                         getAngle(int i)            const { return angle_[i]                       ;}
  bool                          getEstimateWavelet(int i)  const { return estimateWavelet_[i]             ;}
  bool                          getMatchEnergies(int i)    const { return matchEnergies_[i]               ;}
  int                           getNumberOfAngles()        const { return static_cast<int>(angle_.size()) ;}


  void                          releaseGrids();                        // Cuts connection to SeisCube_

  static float  ** readMatrix(const std::string & fileName,
                              int                 n1,
                              int                 n2,
                              const std::string & readReason,
                              std::string       & errText);

  static void             processBackground(Background                    *& background,
                                            const std::vector<WellData *>  & wells,
                                            const Simbox                   * timeSimbox,
                                            const Simbox                   * timeBGSimbox,
                                            GridMapping                   *& timeDepthMapping,
                                            const GridMapping              * timeCutMapping,
                                            const ModelSettings            * modelSettings,
                                            ModelGeneral                   * modelGeneral,
                                            const InputFiles               * inputFile,
                                            const int                      & thisTimeLapse,
                                            std::string                    & errText,
                                            bool                           & failed);
private:
  void             processSeismic(FFTGrid         **& seisCube,
                                  const Simbox      * timeSimbox,
                                  const GridMapping * timeDepthMapping,
                                  const GridMapping * timeCutMapping,
                                  const ModelSettings * modelSettings,
                                  const InputFiles  * inputFiles,
                                  std::string       & errText,
                                  bool              & failed);


  void             processReflectionMatrix(float               **& reflectionMatrix,
                                           const Background      * background,
                                           const std::vector<WellData *> & wells,
                                           const ModelSettings   * modelSettings,
                                           const InputFiles      * inputFiles,
                                           std::string           & errText,
                                           bool                  & failed);

  void             processWavelets(Wavelet                    **& wavelet,
                                   FFTGrid                     ** seisCube,
                                   std::vector<WellData *>        wells,
                                   const float          * const * reflectionMatrix,
                                   const Simbox                 * timeSimbox,
                                   const Surface                * correlationDirection,
                                   const std::vector<Surface *> & waveletEstimInterval,
                                   ModelSettings                * modelSettings,
                                   const InputFiles             * inputFiles,
                                   std::string                  & errText,
                                   bool                         & failed);

  int              process1DWavelet(const ModelSettings          * modelSettings,
                                    const InputFiles             * inputFiles,
                                    const Simbox                 * timeSimbox,
                                    //const FFTGrid        * const * seisCube,
                                    std::vector<WellData *>        wells,
                                    //const std::vector<Surface *> & waveletEstimInterval,
                                    const float                  * reflectionMatrix,
                                    std::string                  & errText,
                                    Wavelet                     *& wavelet,
                                    unsigned int                   i,
                                    bool                           useRickerWavelet);

 int               process3DWavelet(const ModelSettings                     * modelSettings,
                                    const InputFiles                        * inputFiles,
                                    const Simbox                            * timeSimbox,
                                    //const FFTGrid                   * const * seisCube,
                                    //const std::vector<WellData *>           & wells,
                                    //const std::vector<Surface *>            & waveletEstimInterval,
                                    const float                             * reflectionMatrix,
                                    std::string                             & errText,
                                    Wavelet                                *& wavelet,
                                    unsigned int                              i);

  void             setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                                double          vsvp,
                                                const ModelSettings * modelSettings);
  int              getWaveletFileFormat(const std::string & fileName,
                                        std::string & errText);
  double           vsvpFromWells(const std::vector<WellData *> & wells,
                                 int                     nWells);
  void             readAndWriteLocalGridsToFile(const std::string   & fileName,
                                                const std::string   & type,
                                                const float           scaleFactor,
                                                const ModelSettings * modelSettings,
                                                const unsigned int    i,
                                                const Simbox        * timeSimbox,
                                                const Grid2D        * grid);
  void             resampleSurfaceToGrid2D(const Simbox  * simbox,
                                           const Surface * surface,
                                           Grid2D        * outgrid);
  void             resampleGrid2DToSurface(const Simbox  * simbox,
                                           const Grid2D  * grid,
                                           Surface      *& surface);
  bool             findTimeGradientSurface(const std::string     & refTimeFile,
                                           const Simbox          * simbox,
                                           NRLib::Grid2D<float>  & refTimeGradX,
                                           NRLib::Grid2D<float>  & refTimeGradY);
  void             computeStructureDepthGradient(double                 v0,
                                                 double                 radius,
                                                 const Simbox         * timeSimbox,
                                                 const Surface        * t0Surf,
                                                 const Surface        * correlationDirection,
                                                 NRLib::Grid2D<float> & structureDepthGradX,
                                                 NRLib::Grid2D<float> & structureDepthGradY);
  void            computeReferenceTimeGradient( const Simbox         * timeSimbox,
                                                const Surface        * t0Surf,
                                                NRLib::Grid2D<float> & refTimeGradX,
                                                NRLib::Grid2D<float> & refTimeGradY);

  void              calculateSmoothGrad(const Surface * surf, double x, double y, double radius, double ds,  double& gx, double& gy);

  int                       numberOfAngles_;
  FFTGrid                ** seisCube_;              ///< Seismic data cubes
  Wavelet                ** wavelet_;               ///< Wavelet for angle

  float                  ** reflectionMatrix_;      ///< May specify own Zoeppritz-approximation. Default NULL,
                                                    ///< indicating that standard approximation will be used.

  GridMapping             * timeDepthMapping_;      ///< Contains both simbox and mapping used for depth conversion
  GridMapping             * timeCutMapping_;        ///< Simbox and mapping for timeCut*/

  std::vector<Grid2D *>     localNoiseScale_;       ///< Scale factors for local noise

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  Vario                   * angularCorr_;
  std::vector<float>        SNRatio_;
  std::vector<float>        angle_;
  std::vector<bool>         matchEnergies_;
  std::vector<bool>         estimateWavelet_;
  bool                      useLocalNoise_;
  int                       thisTimeLapse_;

};

#endif
