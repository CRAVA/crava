#ifndef MODELAVODYNAMIC_H
#define MODELAVODYNAMIC_H

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
class ModelAVOStatic;

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
                  RandomGen            * randomGen,
                  GridMapping          * timeDepthMapping,
                  GridMapping          * timeCutMapping,
                  std::vector<Surface *> waveletEstimInterval,
                  std::vector<Surface *> wellMoveInterval,
                  std::vector<Surface *> faciesEstimInterval,
                  ModelAVOStatic       * modelAVOstatic);    // modelAVOstatic::wells_ are altered. modelAVOstatic is deliberately sent in as un-const.
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

  void                          releaseGrids();                        // Cuts connection to SeisCube_ and  backModel_
private:
  void             processSeismic(FFTGrid         **& seisCube,
                                  Simbox           *& timeSimbox,
                                  GridMapping      *& timeDepthMapping,
                                  GridMapping      *& timeCutMapping,
                                  ModelSettings    *& modelSettings,
                                  const InputFiles  * inputFiles,
                                  std::string       & errText,
                                  bool              & failed);

  void             processBackground(Background      *& background,
                                     WellData        ** wells,
                                     Simbox           * timeSimbox,
                                     Simbox           * timeBGSimbox,
                                     GridMapping     *& timeDepthMapping,
                                     GridMapping     *& timeCutMapping,
                                     ModelSettings    * modelSettings,
                                     const InputFiles * inputFile,
                                     std::string      & errText,
                                     bool             & failed);
  void             processPriorCorrelations(Corr            *& correlations,
                                            Background       * background,
                                            WellData        ** wells,
                                            Simbox           * timeSimbox,
                                            ModelSettings    * modelSettings,
                                            const InputFiles * inputFiles,
                                            std::string      & errText,
                                            bool             & failed);

  void             processReflectionMatrix(float           **& reflectionMatrix,
                                           Background        * background,
                                           WellData         ** wells,
                                           ModelSettings     * modelSettings,
                                           const InputFiles  * inputFiles,
                                           std::string       & errText,
                                           bool              & failed);

  void             processWavelets(Wavelet                    **& wavelet,
                                   FFTGrid                     ** seisCube,
                                   WellData                    ** wells,
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
                                    WellData                    ** wells,
                                    const std::vector<Surface *> & waveletEstimInterval,
                                    float                        * reflectionMatrix,
                                    std::string                  & errText,
                                    Wavelet                     *& wavelet,
                                    unsigned int                   i);

 int               process3DWavelet(ModelSettings                           * modelSettings,
                                    const InputFiles                        * inputFiles,
                                    Simbox                                  * timeSimbox,
                                    FFTGrid                                ** seisCube,
                                    WellData                               ** wells,
                                    const std::vector<Surface *>            & waveletEstimInterval,
                                    float                                   * reflectionMatrix,
                                    std::string                             & errText,
                                    Wavelet                                *& wavelet,
                                    unsigned int                              i,
                                    const NRLib::Grid2D<float>              & refTimeGradX,
                                    const NRLib::Grid2D<float>              & refTimeGradY,
                                    const std::vector<std::vector<double> > & tGradX,
                                    const std::vector<std::vector<double> > & tGradY);
  void             estimateCorrXYFromSeismic(Surface *& CorrXY,
                                             FFTGrid ** seisCube,
                                             int        nAngles);
  Surface        * findCorrXYGrid(Simbox * timeSimbox, ModelSettings * modelSettings);
  float         ** readMatrix(const std::string & fileName,
                              int                 n1,
                              int                 n2,
                              const std::string & readReason,
                              std::string       & errText);
  void             setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                                double          vsvp,
                                                ModelSettings * modelSettings);
  int              getWaveletFileFormat(const std::string & fileName,
                                        std::string & errText);
  double           vsvpFromWells(WellData ** wells,
                                 int         nWells);
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
};

#endif
