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
class Corr;
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
  ModelAVODynamic(ModelSettings       *& modelSettings,
                  const InputFiles     * inputFiles,
                  std::vector<bool>      failedGeneralDetails,
                  std::vector<bool>      failedStaticDetails,
                  Simbox               * timeSimbox,
                  Simbox              *& timeBGSimbox,
                  Surface              * correlationDirection,
                  RandomGen            * /*randomGen*/,
                  GridMapping          * timeDepthMapping,
                  GridMapping          * timeCutMapping,
                  std::vector<Surface *> waveletEstimInterval,
                  std::vector<Surface *> wellMoveInterval,
                  std::vector<Surface *> faciesEstimInterval,
                  ModelAVOStatic       * modelAVOstatic,
                  ModelGeneral         * modelGeneral,
                  int                    t);    // modelAVOstatic::wells_ are altered. modelAVOstatic is deliberately sent in as un-const.

  ModelAVODynamic(ModelSettings          *& modelSettings,
                  const InputFiles        * inputFiles,
                  ModelAVOStatic          * modelAVOstatic,
                  ModelGeneral            * modelGeneral,
                  SeismicParametersHolder & seismicParameters,
                  Simbox                  * timeSimbox,
                  Surface                 * correlationDirection,
                  GridMapping             * timeDepthMapping,
                  GridMapping             * timeCutMapping,
                  int                       t);

  ~ModelAVODynamic();

  FFTGrid                     * getBackAlpha()             const { return background_->getAlpha() ;}
  FFTGrid                     * getBackBeta()              const { return background_->getBeta()  ;}
  FFTGrid                     * getBackRho()               const { return background_->getRho()   ;}
  Corr                        * getCorrelations()          const { return correlations_           ;}
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


  void                          releaseGrids();                        // Cuts connection to SeisCube_ and  backModel_

  static float  ** readMatrix(const std::string & fileName,
                              int                 n1,
                              int                 n2,
                              const std::string & readReason,
                              std::string       & errText);

private:
  void             processSeismic(FFTGrid         **& seisCube,
                                  Simbox           *& timeSimbox,
                                  GridMapping      *& timeDepthMapping,
                                  GridMapping      *& timeCutMapping,
                                  ModelSettings    *& modelSettings,
                                  const InputFiles  * inputFiles,
                                  std::string       & errText,
                                  bool              & failed);

  void             processBackground(Background           *& background,
                                     std::vector<WellData *> wells,
                                     Simbox                * timeSimbox,
                                     Simbox                * timeBGSimbox,
                                     GridMapping          *& timeDepthMapping,
                                     GridMapping          *& timeCutMapping,
                                     ModelSettings         * modelSettings,
                                     const InputFiles      * inputFile,
                                     std::string           & errText,
                                     bool                  & failed);

  void             processReflectionMatrix(float               **& reflectionMatrix,
                                           Background            * background,
                                           std::vector<WellData *> wells,
                                           ModelSettings         * modelSettings,
                                           const InputFiles      * inputFiles,
                                           std::string           & errText,
                                           bool                  & failed);

  void             processWavelets(Wavelet                    **& wavelet,
                                   FFTGrid                     ** seisCube,
                                   std::vector<WellData *>        wells,
                                   float                       ** reflectionMatrix,
                                   Simbox                       * timeSimbox,
                                   const Surface                * correlationDirection,
                                   const std::vector<Surface *> & waveletEstimInterval,
                                   ModelSettings                * modelSettings,
                                   const InputFiles             * inputFiles,
                                   std::string                  & errText,
                                   bool                         & failed);

  int              process1DWavelet(ModelSettings                * modelSettings,
                                    const InputFiles             * inputFiles,
                                    Simbox                       * timeSimbox,
                                    FFTGrid                     ** seisCube,
                                    std::vector<WellData *>        wells,
                                    const std::vector<Surface *> & waveletEstimInterval,
                                    float                        * reflectionMatrix,
                                    std::string                  & errText,
                                    Wavelet                     *& wavelet,
                                    unsigned int                   i,
                                    bool                           useRickerWavelet);

 int               process3DWavelet(ModelSettings                           * modelSettings,
                                    const InputFiles                        * inputFiles,
                                    Simbox                                  * timeSimbox,
                                    FFTGrid                                ** seisCube,
                                    std::vector<WellData *>                   wells,
                                    const std::vector<Surface *>            & waveletEstimInterval,
                                    float                                   * reflectionMatrix,
                                    std::string                             & errText,
                                    Wavelet                                *& wavelet,
                                    unsigned int                              i,
                                    const NRLib::Grid2D<float>              & refTimeGradX,
                                    const NRLib::Grid2D<float>              & refTimeGradY,
                                    const std::vector<std::vector<double> > & tGradX,
                                    const std::vector<std::vector<double> > & tGradY);

  void             setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                                double          vsvp,
                                                ModelSettings * modelSettings);
  int              getWaveletFileFormat(const std::string & fileName,
                                        std::string & errText);
  double           vsvpFromWells(std::vector<WellData *> wells,
                                 int                     nWells);
  void             readAndWriteLocalGridsToFile(const std::string   & fileName,
                                                const std::string   & type,
                                                const float           scaleFactor,
                                                const ModelSettings * modelSettings,
                                                const unsigned int    i,
                                                Simbox              * timeSimbox,
                                                Grid2D             *& grid);
  void             resampleSurfaceToGrid2D(const Simbox  * simbox,
                                           const Surface * surface,
                                           Grid2D        * outgrid);
  void             resampleGrid2DToSurface(const Simbox  * simbox,
                                           const Grid2D  * grid,
                                           Surface      *& surface);
  bool             findTimeGradientSurface(const std::string     & refTimeFile,
                                           Simbox                * simbox,
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
  Background              * background_;            ///< Holds the background model.
  Corr                    * correlations_;          ///<
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
