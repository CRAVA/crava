#ifndef MODELFILE_H
#define MODELFILE_H

#include <stdio.h>

class Vario;
class ModelSettings;

class ModelFile
{
public:
  ModelFile(char * fileName);
  ~ModelFile(void);
  ModelSettings  * getModelSettings(void)         const { return modelSettings_         ;}

  char           * getSeedFile(void)              const { return seedFile_              ;}
  char          ** getBackFile(void)              const { return backFile_              ;}
  char          ** getWellFile(void)              const { return wellFile_              ;}
  char          ** getHeaderList(void)            const { return headerList_            ;}
  char          ** getTimeSurfFile(void)          const { return timeSurfFile_          ;}
  char          ** getDepthSurfFile(void)         const { return depthSurfFile_         ;}
  char           * getVelocityField(void)         const { return velocityField_         ;}
  char          ** getSeismicFile(void)           const { return seismicFile_           ;}
  char          ** getWaveletFile(void)           const { return waveletFile_           ;}
  char          ** getWaveletEstIntFile(void)     const { return waveletEstIntFile_     ;}
  char          ** getFaciesEstIntFile(void)      const { return faciesEstIntFile_      ;}     
  char           * getReflMatrFile(void)          const { return reflMatrFile_          ;}
  char           * getParamCorrFile(void)         const { return paramCorrFile_         ;}
  char           * getCorrDirFile(void)           const { return corrDirFile_           ;}
  int              getSeisType(int i)             const { return seisType_[i]           ;}
  float          * getWaveletScale(void)          const { return waveScale_             ;}
  int              getNWaveletTransfArgs(void)    const { return nWaveletTransfArgs_    ;}
  int              getSeed(void)                  const { return seed_                  ;}
  float          * getConstBack(void)             const { return constBack_             ;}
  double           getTimeDTop(void)              const { return time_dTop_             ;}
  double           getTimeLz(void)                const { return time_lz_               ;}
  double           getTimeDz(void)                const { return time_dz_               ;}
  int              getTimeNz(void)                const { return time_nz_               ;}
  bool             getHasSignalToNoiseRatio(void) const { return hasSignalToNoiseRatio_ ;}
  bool             getFaciesLogGiven(void)        const { return faciesLogGiven_        ;}
  bool             getDoDepthConversion(void)     const { return doDepthConversion_     ;}
  bool             getGenerateBackground(void)    const { return generateBackground_    ;}
  bool             getParallelTimeSurfaces(void)  const { return parallelTimeSurfaces_  ;}
  bool             getParsingFailed(void)         const { return failed_                ;}
  
  enum             backFileTypes{STORMFILE = -2, SEGYFILE = -1};

private:
  
  int              readCommandWells(char ** params, int & pos, char * errText);
  int              readCommandBackground(char ** params, int & pos, char * errText);
  int              readCommandArea(char ** params, int & pos, char * errText);
  int              readCommandTimeSurfaces(char ** params, int & pos, char * errText);
  int              readCommandDepthConversion(char ** params, int & pos, char * errText);
  int              readCommandSeismic(char ** params, int & pos, char * errText, int seisType = 0);
  int              readCommandAngularCorr(char ** params, int & pos, char * errText);
  int              readCommandLateralCorr(char ** params, int & pos, char * errText);
  int              readCommandSimulate(char ** params, int & pos, char * errText);
  int              readCommandSeed(char ** params, int & pos, char * errText);
  int              readCommandPadding(char ** params, int & pos, char * errText);
  int              readCommandPrefix(char ** params, int & pos, char * errText);
  int              readCommandOutput(char ** params, int & pos, char * errText);
  int              readCommandWhiteNoise(char ** params, int & pos, char * errText);
  int              readCommandSegYOffset(char ** params, int & pos, char * errText);
  int              readCommandForceFile(char ** params, int & pos, char * errText);
  int              readCommandKriging(char ** params, int & pos, char * errText);
  int              readCommandDebug(char ** params, int & pos, char * errText);
  int              readCommandLocalWavelet(char ** params, int & pos, char * errText);
  int              readCommandEnergyTreshold(char ** params, int & pos, char * errText);
  int              readCommandParameterCorr(char ** params, int & pos, char * errText);
  int              readCommandReflectionMatrix(char ** params, int & pos, char * errText);
  int              readCommandFrequencyBand(char ** params, int & pos, char * errText);
  int              readCommandBackgroundControl(char ** params, int & pos, char * errText);
  int              readCommandGiveSignalToNoiseRatios(char ** params, int & pos, char * errText);
  int              readCommandSeismicResolution(char ** params, int & pos, char * errText);
  int              readCommandWaveletTaperingL(char ** params, int & pos, char * errText);
  int              readCommandPUndef(char ** params, int & pos, char * errText);
  int              readCommandMaxDeviationAngle(char ** params, int & pos, char * errText);
  int              readCommandAllowedParameterValues(char ** params, int & pos, char * errText);
  int              readCommandAllowedResidualVariances(char ** params, int & pos, char * errText);
  int              readCommandCorrelationDirection(char ** params, int & pos, char * errText);
  int              readCommandWaveletEstimationInterval(char ** params, int & pos, char * errText);
  int              readCommandFaciesEstimationInterval(char ** params, int & pos, char * errText);
  int              readCommandLogLevel(char ** params, int & pos, char * errText);
  int              readCommandTraceHeaderFormat(char ** params, int & pos, char * errText);
  int              getParNum(char ** params, int pos, int & error, char * errText,
                             const char * command, int min, int max = 0);
  int              checkFileOpen(char ** fNames, int nFiles, const char * command, char * errText, int start = 0,
                                 bool details = true);
  int              checkFileType(char * fileName);
  Vario          * createVario(char ** param, int nPar, const char * command, char * errText);

  ModelSettings  * modelSettings_;

  char          ** headerList_;            // The keywords to look for for time, sonic, shear sonic and density
  char          ** wellFile_;              // File names: wells
  char          ** backFile_;              // File names (temporarily stored).
  char          ** timeSurfFile_;          // File names: top and base time surfaces
  char          ** depthSurfFile_;         // File names: top and base depth surfaces
  char           * velocityField_;         // File names: velocity field, or command
  char          ** seismicFile_;           // File names: seismic data
  char          ** waveletFile_;           // File names: wavelets
  char          ** waveletEstIntFile_;     // File names: Wavelet estimation interval
  char          ** faciesEstIntFile_;      // File names: Facies estimation interval
  char           * seedFile_;              // File specifying the seed
  char           * reflMatrFile_;          // File name for reflection matrix file.
  char           * paramCorrFile_;         // File name for parameter correlations.
  char           * corrDirFile_;           // File name for correlation direction

  int            * seisType_;              // Type of seismic (STANDARDSEIS or PSSEIS)
  int              nSeisData_;             // Number of seismic cubes to condition on.

  float          * angle_;                 // Angle of cube.
  float          * noiseEnergy_;           // Noise Variance.
  float          * waveScale_;             // wavelet scaling/gain
  int              nWaveletTransfArgs_;    // Number of parameters in wavelet transform. Should be equal to nSeisData_

  int              seed_;

  char          ** fName_;                 // Facies names
  int            * fnr_;                   // Facies labels

  float          * constBack_;             // Constant background model (temporarily stored)
                                           // Note: Negative value ==> read from file (actual value gives format).

  double           time_dTop_;             // Used when top and base surfaces are parallel
  double           time_lz_;
  double           time_dz_;
  int              time_nz_;

  bool             faciesLogGiven_;
  bool             doDepthConversion_;     // Flag taking value 1 if depth surfaces are defined
  bool             parallelTimeSurfaces_;
  bool             generateBackground_;    // Generate background model from filtered wells.
  bool             hasSignalToNoiseRatio_; // Use SN ratio instead of error variance in model file. 

  bool             failed_;                // Indicates whether errors ocuured during construction. 
};

#endif

