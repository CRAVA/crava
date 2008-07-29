
#ifndef MODELSETTINGS_H
#define MODELSETTINGS_H

#include <stdio.h>

#include "lib/global_def.h"

class Vario;

class ModelSettings
{
public:
  ModelSettings(void);
  ~ModelSettings(void);

  Vario        * getAngularCorr(void)          const { return angularCorr_      ;} 
  Vario        * getLateralCorr(void)          const { return lateralCorr_      ;}
  Vario        * getBackgroundVario(void)      const { return backgroundVario_  ;} 
  float        * getKrigingParameters(void)    const { return krigingParams_    ;}
  float        * getAngle(void)                const { return angle_            ;}
  int            getNumberOfAngles(void)       const { return nAngles_          ;} 
  float        * getNoiseEnergy(void)          const { return noiseEnergy_      ;} 
  bool         * getMatchEnergies(void)        const { return matchEnergies_    ;} 
  char         * getFaciesName(int i)          const { return faciesNames_[i]   ;}
  int            getIndicatorBGTrend(int i)    const { return indBGTrend_[i]    ;}
  int            getIndicatorWavelet(int i)    const { return indWavelet_[i]    ;}
  int            getIndicatorFacies(int i)     const { return indFacies_[i]     ;}
  int            getNumberOfFacies(void)       const { return nFacies_          ;}
  int            getNumberOfWells(void)        const { return nWells_           ;}
  int            getNumberOfSimulations(void)  const { return nSimulations_     ;}
  float          getAlphaMin(void)             const { return alpha_min_        ;}
  float          getAlphaMax(void)             const { return alpha_max_        ;}
  float          getBetaMin(void)              const { return beta_min_         ;}
  float          getBetaMax(void)              const { return beta_max_         ;}
  float          getRhoMin(void)               const { return rho_min_          ;}
  float          getRhoMax(void)               const { return rho_max_          ;}
  float          getVarAlphaMin(void)          const { return var_alpha_min_    ;}
  float          getVarAlphaMax(void)          const { return var_alpha_max_    ;}
  float          getVarBetaMin(void)           const { return var_beta_min_     ;}
  float          getVarBetaMax(void)           const { return var_beta_max_     ;}
  float          getVarRhoMin(void)            const { return var_rho_min_      ;}
  float          getVarRhoMax(void)            const { return var_rho_max_      ;}
  float          getMaxHzBackground(void)      const { return maxHz_background_ ;}
  float          getMaxHzSeismic(void)         const { return maxHz_seismic_    ;}
  float          getMaxRankCorr(void)          const { return maxRankCorr_      ;}
  float          getMaxMergeDist(void)         const { return maxMergeDist_     ;}
  float          getMaxDevAngle(void)          const { return maxDevAngle_      ;}
  float          getLowCut(void)               const { return lowCut_           ;}
  float          getHighCut(void)              const { return highCut_          ;}
  float          getWNC(void)                  const { return wnc_              ;}
  float          getEnergyThreshold(void)      const { return energyThreshold_  ;}
  float          getMinRelWaveletAmp(void)     const { return minRelWaveletAmp_ ;}
  float          getMaxWaveletShift(void)      const { return maxWaveletShift_  ;}
  float          getWaveletTaperingL(void)     const { return waveletTaperingL_ ;}
  float          getXpad(void)                 const { return xPad_             ;}
  float          getYpad(void)                 const { return yPad_             ;}
  float          getZpad(void)                 const { return zPad_             ;}
  int            getNXpad(void)                const { return nxPad_            ;}
  int            getNYpad(void)                const { return nyPad_            ;}
  int            getNZpad(void)                const { return nzPad_            ;}
  float          getSegyOffset(void)           const { return segyOffset_       ;}
  float          getPundef(void)               const { return p_undef_          ;}
  double         getLzLimit(void)              const { return lzLimit_          ;}
  int            getOutputFlag(void)           const { return outputFlag_       ;}
  int            getFormatFlag(void)           const { return formatFlag_       ;}
  int            getDebugFlag(void)            const { return debugFlag_        ;}
  static int     getDebugLevel(void)                 { return debugFlag_        ;}
  int            getFileGrid(void)             const { return fileGrid_         ;}
  bool           getGenerateSeismic(void)      const { return generateSeismic_  ;}
  bool           getDoInversion(void);


  void           setAngularCorr(Vario * vario);    
  void           setLateralCorr(Vario * vario);    
  void           setBackgroundVario(Vario * vario);
  void           setKrigingParameters(float * krigingParams, int nParams);
  void           setAngle(float * angle, int nAngles);
  void           setNumberOfAngles(int nAngles)              { nAngles_          = nAngles          ;} 
  void           setNoiseEnergy(float * noiseEnergy, int nAngles);
  void           setMatchEnergies(float * waveletScale, int nAngles);
  void           setAllIndicatorsTrue(int nWells);
  void           setIndicatorBGTrend(int * indBGTrend, int nWells);
  void           setIndicatorWavelet(int * indWavelet, int nWells);
  void           setIndicatorFacies(int * indFacies, int nWells);
  void           setFaciesNames(char ** faciesNames, int nFacies);
  void           setNumberOfFacies(int nFacies)              { nFacies_          = nFacies          ;}
  void           setNumberOfWells(int nWells)                { nWells_           = nWells           ;} 
  void           setNumberOfSimulations(int nSimulations)    { nSimulations_     = nSimulations     ;} 
  void           setAlphaMin(float alpha_min)                { alpha_min_        = alpha_min        ;}
  void           setAlphaMax(float alpha_max)                { alpha_max_        = alpha_max        ;}
  void           setBetaMin(float beta_min)                  { beta_min_         = beta_min         ;}
  void           setBetaMax(float beta_max)                  { beta_max_         = beta_max         ;}
  void           setRhoMin(float rho_min)                    { rho_min_          = rho_min          ;}
  void           setRhoMax(float rho_max)                    { rho_max_          = rho_max          ;}
  void           setVarAlphaMin(float var_alpha_min)         { var_alpha_min_    = var_alpha_min    ;}
  void           setVarAlphaMax(float var_alpha_max)         { var_alpha_max_    = var_alpha_max    ;}
  void           setVarBetaMin(float var_beta_min)           { var_beta_min_     = var_beta_min     ;}
  void           setVarBetaMax(float var_beta_max)           { var_beta_max_     = var_beta_max     ;}
  void           setVarRhoMin(float var_rho_min)             { var_rho_min_      = var_rho_min      ;}
  void           setVarRhoMax(float var_rho_max)             { var_rho_max_      = var_rho_max      ;}
  void           setMaxHzBackground(float maxHz_background)  { maxHz_background_ = maxHz_background ;}
  void           setMaxHzSeismic(float maxHz_seismic)        { maxHz_seismic_    = maxHz_seismic    ;}
  void           setMaxRankCorr(float maxRankCorr)           { maxRankCorr_      = maxRankCorr      ;}
  void           setMaxMergeDist(float maxMergeDist)         { maxMergeDist_     = maxMergeDist     ;}
  void           setMaxDevAngle(float maxDevAngle)           { maxDevAngle_      = maxDevAngle      ;}
  void           setLowCut(float lowCut)                     { lowCut_           = lowCut           ;}
  void           setHighCut(float highCut)                   { highCut_          = highCut          ;}
  void           setWNC(float wnc)                           { wnc_              = wnc              ;}
  void           setEnergyThreshold(float energyThreshold)   { energyThreshold_  = energyThreshold  ;}
  void           setinRelWaveletAmp(float minRelWaveletAmp)  { minRelWaveletAmp_ = minRelWaveletAmp ;}
  void           setMaxWaveletShift(float maxWaveletShift)   { maxWaveletShift_  = maxWaveletShift  ;}
  void           setWaveletTaperingL(float waveletTaperingL) { waveletTaperingL_ = waveletTaperingL ;}
  void           setXpad(float xPad)                         { xPad_             = xPad             ;}
  void           setYpad(float yPad)                         { yPad_             = yPad             ;}
  void           setZpad(float zPad)                         { zPad_             = zPad             ;}
  void           setNXpad(int nxPad)                         { nxPad_            = nxPad            ;}
  void           setNYpad(int nyPad)                         { nyPad_            = nyPad            ;}
  void           setNZpad(int nzPad)                         { nzPad_            = nzPad            ;}
  void           setSegyOffset(float segyOffset)             { segyOffset_       = segyOffset       ;}
  void           setPundef(float p_undef)                    { p_undef_          = p_undef          ;}
  void           setLzLimit(double lzLimit)                  { lzLimit_          = lzLimit          ;}
  void           setOutputFlag(int outputFlag);
  void           setFormatFlag(int formatFlag)               { formatFlag_       = formatFlag       ;}
  void           setDebugFlag(int debugFlag)                 { debugFlag_        = debugFlag        ;}
  void           setFileGrid(int fileGrid)                   { fileGrid_         = fileGrid         ;}
  void           setGenerateSeismic(bool generateSeismic)    { generateSeismic_  = generateSeismic  ;}

  enum           outputGrids{PREDICTION        = 1, 
                             CORRELATION       = 2, 
                             RESIDUAL          = 4, 
                             VP                = 8, 
                             VS                = 16,
                             RHO               = 32, 
                             LAMELAMBDA        = 64, 
                             LAMEMU            = 128, 
                             POISSONRATIO      = 256,
                             AI                = 512, 
                             SI                = 1024, 
                             VPVSRATIO         = 2048, 
                             MURHO             = 4096, 
                             LAMBDARHO         = 8192, 
                             PRIORCORRELATIONS = 16384, 
                             BACKGROUND        = 32768, 
                             WELLS             = 65536, 
                             WAVELETS          = 131072, 
                             NOTIME            = 262144,
                             FACIESPROB        = 524288,
                             FACIESPROBRELATIVE = 1048576};
  enum             sseismicTypes{STANDARDSEIS = 0, PSSEIS = 1};

  static void    setFilePrefix(char * filePrefix);
  static char  * makeFullFileName(const char * name, const char * postfix = NULL);

private:

  Vario        * angularCorr_;           // Variogram for lateral error correlation
  Vario        * lateralCorr_;           // Variogram for lateral parameter correlation 
  Vario        * backgroundVario_;       // Used for lateral background correlation.

  float        * krigingParams_;   

  float        * angle_;                 // Angles
  int            nAngles_;               //
  float        * noiseEnergy_;           // Noise Variance .
  bool         * matchEnergies_;         // Let dataVariance_ = signalVariance_

  int          * indBGTrend_;            // Use well to estimate background trend? (1=yes,0=no)
  int          * indWavelet_;            // Use well to estimate wavelet? (1=yes,0=no)
  int          * indFacies_;             // Use well to estimate facies? (1=yes,0=no)

  char        ** faciesNames_;           // Facies names
  int            nFacies_;

  int            nWells_;
  int            nSimulations_;

  float          alpha_min_;             // Vp - smallest allowed value
  float          alpha_max_;             // Vp - largest allowed value
  float          beta_min_;              // Vs - smallest allowed value
  float          beta_max_;              // Vs - largest allowed value
  float          rho_min_;               // Rho - smallest allowed value
  float          rho_max_;               // Rho - largest allowed value
 
  float          var_alpha_min_;         //| These min and max values are used for consistency check. If  
  float          var_alpha_max_;         //| variances are outside these ranges there is probably a
  float          var_beta_min_;          //| problem with the logs.
  float          var_beta_max_;          //| 
  float          var_rho_min_;           //| The limits are for point variances. The minimum allowed variance 
  float          var_rho_max_;           //| for parameters will be scaled with 1/dt*dt

  float          maxHz_background_;      // Background resolution (high cut frequency)
  float          maxHz_seismic_;         // Seismic resolution (high cut frequency)

  float          maxRankCorr_;           // Vp-Vs correlation threshold for regarding Vs log synthetic
  float          maxMergeDist_;          // log entries closer than this will be merged
  float          maxDevAngle_;           // Wells with a local deviation larger than this is treated as deviated

  float          lowCut_;                // lower limit for frequency to be inverted
  float          highCut_;               // upper limit for frecuency to be inverted

  float          wnc_;                   // White noise component, see crava.h  

  float          energyThreshold_;       // If energy in reflection trace divided by mean energy
                                         // in reflection traces is lower than this, the reflections
                                         // will be interpolated. Default 0.

  float          minRelWaveletAmp_;      // Minimum relative wavelet amplitude. Smaller amplitudes are disregarded.
  float          maxWaveletShift_;       // Largest allowed shift when estimating wavelet
  float          waveletTaperingL_;      // Til Odds waveletestimering

  float          xPad_;                  // Padding factor in x direction
  float          yPad_;
  float          zPad_; 

  int            nxPad_;                 // Number of cells to pad in x direction
  int            nyPad_;
  int            nzPad_; 

  float          segyOffset_;            // Starttime for SegY cubes.

  float          p_undef_;               // NBNB-PAL: Hva gjør denne?

  double         lzLimit_;               // Minimum allowed value for (min interval thickness)/(max interval thickness)

  int            outputFlag_;            // Decides which grids to write (except simulation)
  int            formatFlag_;            // Decides output format, see fftgird.h
  int            fileGrid_;              // Indicator telling if grids are to be kept on file

  bool           generateSeismic_;        // Forward modelling

  static int     debugFlag_;
  static std::string  filePrefix_;            // Prefix (including path) for all output files
};

#endif
