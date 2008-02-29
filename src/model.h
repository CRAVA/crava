#ifndef MODEL_H
#define MODEL_H

#include <stdio.h>

#include "lib/global_def.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"

struct irapgrid;
class Corr;
class ModelFile;
class Wavelet;
class Vario;
class Simbox;
class WellData;
class FFTGrid;
class RandomGen;

class Model{
public:
  Model(char * fileName);
  ~Model();

  ModelSettings  * getModelSettings()         const { return modelSettings_          ;}
  Simbox         * getTimeSimbox()            const { return timeSimbox_             ;}
  Simbox         * getDepthSimbox()           const { return depthSimbox_            ;}
  WellData      ** getWells()                 const { return wells_                  ;}
  FFTGrid        * getBackAlpha()             const { return background_->getAlpha() ;}
  FFTGrid        * getBackBeta()              const { return background_->getBeta()  ;}
  FFTGrid        * getBackRho()               const { return background_->getRho()   ;}
  Corr           * getPriorCorrelations()     const { return priorCorrelations_      ;}
  FFTGrid       ** getSeisCubes()             const { return seisCube_               ;}
  Wavelet       ** getWavelets()              const { return wavelet_                ;}
  float         ** getAMatrix()               const { return reflectionMatrix_       ;}
  RandomGen      * getRandomGen()             const { return randomGen_              ;} 
  irapgrid       * getCorrXYGrid();

  bool             hasSignalToNoiseRatio()    const { return hasSignalToNoiseRatio_  ;}
  bool             getFailed()                const { return failed_                 ;}
  void             releaseWells();                                        // Deallocates well data.
  void             releaseGrids();                                        // Cuts connection to SeisCube_ and  backModel_

private:
  void             makeTimeSimbox(char * errText);
  void             makeDepthSimbox(char * errText);
  void             processSeismic(char * errText);
  void             processWells(char * errText);
  void             processBackground(char * errText);
  void             processPriorCorrelations(char * errText);
  void             processReflectionMatrix(char * errText);
  void             processWavelets(void);
  void             setSimboxSurfaces(Simbox * simbox, char ** surfFile, bool parallelSurfaces, 
                                     double dTop, double lz, double dz, int nz, int & error);
  void             estimateXYPaddingSizes(void);
  void             estimateZPaddingSize(void);
  int              readSegyFiles(char ** fNames, int nFiles, FFTGrid ** target, char * errText,
                                 int gridType, int start = 0);
  int              readStormFile(char *fName, FFTGrid * & target, const char * parName, char * errText);

  void             estimateCorrXYFromSeismic(irapgrid * CorrXY);

  int              setPaddingSize(int nx, float px);
  float         ** readMatrix(char * fileName, int n1, int n2, const char * readReason, char * errText);
  float         ** getClassicAMatrix(void);
  void             setupDefaultReflectionMatrix(void);
  int              checkFileOpen(char ** fNames, int nFiles, const char * command, char * errText, 
                                 int start = 0, bool details = true);
  void             checkAvailableMemory(void);
  void             checkFaciesNames(void);
  void             printSettings(void);

  ModelFile      * modelFile_;
  ModelSettings  * modelSettings_;
  Simbox         * timeSimbox_;            // Information about simulation area.
  Simbox         * depthSimbox_;           // Simbox with depth surfaces
  WellData      ** wells_;                 // Well data
  Background     * background_;            // Holds the background model.
  Corr           * priorCorrelations_;     //
  FFTGrid       ** seisCube_;              // Seismic data cubes
  Wavelet       ** wavelet_;               // Wavelet for angle
  irapgrid      ** shiftGrids_;            // Grids containing shift data for wavelets
  irapgrid      ** gainGrids_;             // Grids containing gain data for wavelets.
  RandomGen      * randomGen_;             // Random generator.

  float         ** reflectionMatrix_;      // May specify own Zoeppritz-approximation. Default NULL,
                                           // indicating that standard approximation will be used.

  int              nxPad_;                 // total grid size size in x direction (padding included)
  int              nyPad_;
  int              nzPad_;

  bool             hasSignalToNoiseRatio_; // Use SN ratio instead of error variance in model file. 
  bool             failed_;                // Indicates whether errors ocuured during construction. 
};

#endif

