#ifndef INPUTFILES_H
#define INPUTFILES_H

#include <string>
#include <vector>
#include <map>

class InputFiles
{
public:
  InputFiles(void);
  ~InputFiles(void);

  const std::string              & getSeedFile(void)             const { return seedFile_             ;}
  const std::string              & getWellFile(int i)            const { return wellFiles_[i]         ;}
  const std::string              & getSeismicFile(int i)         const { return seismicFiles_[i]      ;}
  const std::string              & getWaveletFile(int i)         const { return waveletFiles_[i]      ;}
  const std::string              & getShiftFile(int i)           const { return waveletShiftFiles_[i] ;}
  const std::string              & getScaleFile(int i)           const { return waveletScaleFiles_[i] ;}
  const std::string              & getWaveletEstIntFile(int i)   const { return waveletEstIntFile_[i] ;}
  const std::string              & getWaveletFilterFile(int i)   const { return waveletFilterFiles_[i];}
  const std::string              & getWaveletCorrFile(int i)     const { return waveletCorrFiles_[i]  ;}
  const std::string              & getWellMoveIntFile(int i)     const { return wellMoveIntFile_[i]   ;}
  const std::string              & getFaciesEstIntFile(int i)    const { return faciesEstIntFile_[i]  ;}
  const std::vector<std::string> & getTimeSurfFiles(void)        const { return timeSurfFiles_        ;}
  const std::string              & getTimeSurfFile(int i)        const { return timeSurfFiles_[i]     ;}
  const std::vector<std::string> & getDepthSurfFiles(void)       const { return depthSurfFiles_       ;}
  const std::string              & getDepthSurfFile(int i)       const { return depthSurfFiles_[i]    ;}
  const std::string              & getVelocityField(void)        const { return velocityField_        ;}
  const std::string              & getBackFile(int i)            const { return backFile_[i]          ;}
  const std::string              & getBackVelFile(void)          const { return backVelFile_          ;}
  const std::string              & getCorrDirFile(void)          const { return corrDirFile_          ;}
  const std::string              & getReflMatrFile(void)         const { return reflMatrFile_         ;}
  const std::string              & getParamCorrFile(void)        const { return paramCorrFile_        ;}
  const std::string              & getLocalNoiseFile(int i)      const { return localNoiseFiles_[i]   ;}
  const std::string              & getTempCorrFile(void)         const { return tempCorrFile_         ;}
  const std::string              & getRefSurfaceFile(void)       const { return refSurfaceFile_       ;}
  const std::string              & getInputDirectory(void)       const { return inputDirectory_       ;}
  const std::map<std::string,std::string> & getPriorFaciesProbFile(void)   const {return priorFaciesProb_ ;}
  const std::string              & getAreaSurfaceFile(void)      const { return areaSurfaceFile_      ;}

  int                              getNumberOfSeismicFiles(void) const { return static_cast<int>(seismicFiles_.size());}

  void setSeedFile(const std::string & seedFile)                          { seedFile_             = seedFile          ;}
  void addWellFile(const std::string & wellFile)                          { wellFiles_.push_back(wellFile)            ;}
  void addSeismicFile(const std::string & seismicFile)                    { seismicFiles_.push_back(seismicFile)      ;}
  void addWaveletFile(const std::string & waveletFile)                    { waveletFiles_.push_back(waveletFile)      ;}
  void addShiftFile(const std::string & shiftFile)                        { waveletShiftFiles_.push_back(shiftFile)   ;}
  void addScaleFile(const std::string & scaleFile)                        { waveletScaleFiles_.push_back(scaleFile)   ;}
  void addNoiseFile(const std::string &noiseFile)                         { localNoiseFiles_.push_back(noiseFile)     ;}
  void setAreaSurfaceFile(const std::string &areaFile)                    { areaSurfaceFile_      = areaFile          ;}
  void setWaveletEstIntFile(int i, const std::string & waveletEstIntFile) { waveletEstIntFile_[i] = waveletEstIntFile ;}
  void addWaveletFilterFile(const std::string & filterFile)               { waveletFilterFiles_.push_back(filterFile) ;}
  void addWaveletCorrFile(const std::string & corrFile)                   { waveletCorrFiles_.push_back(corrFile)     ;}
  void setWellMoveIntFile(int i, const std::string & wellMoveIntFile)     { wellMoveIntFile_[i]   = wellMoveIntFile   ;}
  void setFaciesEstIntFile(int i, const std::string & faciesEstIntFile)   { faciesEstIntFile_[i]  = faciesEstIntFile  ;}
  void addTimeSurfFile(const std::string & timeSurfFile)                  { timeSurfFiles_.push_back(timeSurfFile)    ;}
  void setDepthSurfFile(int i, const std::string & depthSurfFile)         { depthSurfFiles_[i]    = depthSurfFile     ;}
  void addMultizoneSurfaceFile(const std::string & fileName)              { multizoneSurfaceFiles_.push_back(fileName);}
  void setVelocityField(const std::string & velocityField)                { velocityField_        = velocityField     ;}
  void setBackFile(int i, const std::string & backFile)                   { backFile_[i]          = backFile          ;}
  void setBackVelFile(const std::string & backVelFile)                    { backVelFile_          = backVelFile       ;}
  void setReflMatrFile(const std::string & reflMatrFile)                  { reflMatrFile_         = reflMatrFile      ;}
  void setCorrDirFile(const std::string & corrDirFile)                    { corrDirFile_          = corrDirFile       ;}
  void setParamCorrFile(const std::string & paramCorrFile)                { paramCorrFile_        = paramCorrFile     ;}
  void setTempCorrFile(const std::string & tempCorrFile)                  { tempCorrFile_         = tempCorrFile      ;}
  void setRefSurfaceFile(const std::string & refSurfaceFile)              { refSurfaceFile_       = refSurfaceFile    ;}
  void setInputDirectory(std::string inputDirectory)                      { inputDirectory_       = inputDirectory    ;}
  void setPriorFaciesProb(std::string name,std::string file)              { priorFaciesProb_[name] = file                ;}
  std::string addInputPathAndCheckFiles();

private:
  std::string addPathAndCheck(std::string & fileName, const bool possiblyNumber = false);

  std::string                seedFile_;              ///< File specifying the seed
  std::vector<std::string>   wellFiles_;             ///< File names: wells
  std::vector<std::string>   seismicFiles_;          ///< File names: seismic data
  std::vector<std::string>   waveletFiles_;          ///< File names: wavelets
  std::vector<std::string>   waveletShiftFiles_;     ///< File names: wavelets
  std::vector<std::string>   waveletScaleFiles_;     ///< File names: wavelets
  std::vector<std::string>   waveletEstIntFile_;     ///< File names: Wavelet estimation interval
  std::vector<std::string>   waveletFilterFiles_;    ///< File names: Filter files for 3D wavelets
  std::vector<std::string>   waveletCorrFiles_;      ///< File names: Filter with correction factors for 3D wavelets
  std::vector<std::string>   wellMoveIntFile_;       ///< File names: Well move interval
  std::vector<std::string>   faciesEstIntFile_;      ///< File names: Facies estimation interval
  std::vector<std::string>   timeSurfFiles_;         ///< File names: top and base time surfaces
  std::vector<std::string>   depthSurfFiles_;        ///< File names: top and base depth surfaces
  std::vector<std::string>   multizoneSurfaceFiles_; ///< File names: top and base surfaces for multizone background model
  std::string                velocityField_;         ///< File names: velocity field, or command
  std::vector<std::string>   backFile_;              ///< File names (temporarily stored).
  std::string                backVelFile_;           ///< File names Vp-velocity for background modelling
  std::string                reflMatrFile_;          ///< File name for reflection matrix file.
  std::string                corrDirFile_;           ///< File name for correlation direction
  std::string                paramCorrFile_;         ///< File name for correlation between parameters.
  std::string                tempCorrFile_;          ///< File name for temporal parameter correlations.
  std::string                refSurfaceFile_;        ///< File name for reference time surface corresponding to z0.

  std::string                inputDirectory_;        ///< Base directory for input files.
  std::vector<std::string>   localNoiseFiles_;       ///< File names: local noise
  std::map<std::string, std::string> priorFaciesProb_; ///< File names for locally varying prior facies probability
  std::string                areaSurfaceFile_;
};

#endif
