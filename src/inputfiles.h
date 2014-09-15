/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef INPUTFILES_H
#define INPUTFILES_H

#include <vector>
#include <map>
#include <string>

class InputFiles
{
public:
  InputFiles(void);
  ~InputFiles(void);

  const std::string              & getSeedFile(void)             const { return seedFile_             ;}
  const std::string              & getWellFile(int i)            const { return wellFiles_[i]         ;}
  const std::vector<std::string> & getWellFiles(void)            const { return wellFiles_            ;}
  const std::string              & getSeismicFile(int i, int j)  const { return timeLapseSeismicFiles_[i][j]      ;}
  const std::string              & getWaveletFile(int i, int j)  const { return timeLapseWaveletFiles_[i][j]      ;}
  const std::string              & getShiftFile(int i, int j)    const { return timeLapseWaveletShiftFiles_[i][j] ;}
  const std::string              & getScaleFile(int i, int j)    const { return timeLapseWaveletScaleFiles_[i][j] ;}
  const std::string              & getLocalNoiseFile(int i, int j)const{ return timeLapseLocalNoiseFiles_[i][j]   ;}
  const std::vector<std::string> & getTravelTimeHorizons(int i)  const { return timeLapseTravelTimeHorizons_[i]   ;}
  const std::string              & getRmsVelocities(int i)       const { return rmsVelocities_[i]                 ;}
  const std::string              & getWaveletEstIntFileTop(int i)const { return waveletEstIntFileTop_[i]          ;}
  const std::string              & getWaveletEstIntFileBase(int i)const{ return waveletEstIntFileBase_[i]         ;}
  const std::string              & getWaveletFilterFile(int i)   const { return waveletFilterFiles_[i];}
  const std::string              & getWaveletCorrFile(int i)     const { return waveletCorrFiles_[i]  ;}
  const std::string              & getWellMoveIntFile(int i)     const { return wellMoveIntFile_[i]   ;}
  const std::string              & getFaciesEstIntFile(int i)    const { return faciesEstIntFile_[i]  ;}
  const std::string              & getTimeSurfTopFile()          const { return timeSurfTopFile_      ;}
  const std::string              & getDepthSurfTopFile()         const { return depthSurfTopFile_     ;}
  const std::string              & getVelocityField(void)        const { return velocityField_        ;}
  const std::string              & getBackFile(int i)            const { return backFile_[i]          ;}
  const std::string              & getBackVelFile(void)          const { return backVelFile_          ;}
  const std::string              & getReflMatrFile(void)         const { return reflMatrFile_         ;}
  const std::string              & getParamCovFile(void)         const { return paramCovFile_         ;}
  const std::string              & getTempCorrFile(void)         const { return tempCorrFile_         ;}
  const std::string              & getParamAutoCovFile(void)     const { return paramAutoCovFile_     ;}
  const std::string              & getRefSurfaceFile(void)       const { return refSurfaceFile_       ;}
  const std::string              & getInputDirectory(void)       const { return inputDirectory_       ;}
  const std::map<std::string,std::string> & getPriorFaciesProbFile(void) const {return priorFaciesProb_ ;}
  const std::string              & getAreaSurfaceFile(void)      const { return areaSurfaceFile_      ;}
  //const std::vector<std::string> & getMultizoneSurfaceFiles()    const { return multizoneSurfaceFiles_;}
  const std::string              & getTrendCube(int i)           const { return trendCubes_[i]        ;}
  const std::string              & getGravimetricData(int i)     const { return gravimetricData_[i]   ;}
  const std::vector<std::string> & getSeismicFiles(void)         const { return seismicFiles_         ;}
  const std::map<std::string, std::string> & getCorrDirFiles(void)                        const { return corrDirFiles_                        ;}
  const std::string                        & getCorrDirFile(std::string name)             const { return corrDirFiles_.find(name)->second     ;}
  const std::map<std::string, std::string> & getCorrDirTopSurfaceFiles(void)              const { return corrDirTopFiles_                     ;}
  const std::string                        & getCorrDirTopSurfaceFile(std::string name)   const { return corrDirTopFiles_.find(name)->second  ;}
  const std::map<std::string, std::string> & getCorrDirBaseSurfaceFiles(void)             const { return corrDirBaseFiles_                    ;}
  const std::string                        & getCorrDirBaseSurfaceFile(std::string name)  const { return corrDirBaseFiles_.find(name)->second ;}
  const std::vector<std::vector<std::string> > & getTimeLapseSeismicFiles(void)    const { return timeLapseSeismicFiles_                 ;}
  const std::map<std::string, std::string> & getBaseTimeSurfaces(void)             const { return timeSurfBaseFiles_                     ;}
  const std::string                        & getBaseTimeSurface(std::string name)  const { return timeSurfBaseFiles_.find(name)->second  ;}
  const std::map<std::string, std::string> & getBaseDepthSurfaces(void)            const { return depthSurfBaseFiles_                    ;}
  const std::string                        & getBaseDepthSurface(std::string name) const { return depthSurfBaseFiles_.find(name)->second ;}

  int                              getNumberOfSeismicFiles(int i)const { return static_cast<int>(timeLapseSeismicFiles_[i].size());}

  void setSeedFile(const std::string & seedFile)                          { seedFile_             = seedFile          ;}
  void addWellFile(const std::string & wellFile)                          { wellFiles_.push_back(wellFile)            ;}
  void addSeismicFile(const std::string & seismicFile)                    { seismicFiles_.push_back(seismicFile)      ;}
  void addRmsVelocity(const std::string & rmsVelocityFile)                { rmsVelocities_.push_back(rmsVelocityFile) ;}
  void addTravelTimeHorizon(const std::string & horizon)                  { travelTimeHorizons_.push_back(horizon)    ;}
  void addGravimetricData(const std::string & dataFile)                   { gravimetricData_.push_back(dataFile)      ;}
  void addWaveletFile(const std::string & waveletFile)                    { waveletFiles_.push_back(waveletFile)      ;}
  void addShiftFile(const std::string & shiftFile)                        { waveletShiftFiles_.push_back(shiftFile)   ;}
  void addScaleFile(const std::string & scaleFile)                        { waveletScaleFiles_.push_back(scaleFile)   ;}
  void addNoiseFile(const std::string &noiseFile)                         { localNoiseFiles_.push_back(noiseFile)     ;}
  void setAreaSurfaceFile(const std::string &areaFile)                    { areaSurfaceFile_      = areaFile          ;}
  void addWaveletEstIntFileTop(const std::string & fileTop)               { waveletEstIntFileTop_.push_back(fileTop)  ;}
  void addWaveletEstIntFileBase(const std::string & fileBase)             { waveletEstIntFileBase_.push_back(fileBase);}
  void addWaveletFilterFile(const std::string & filterFile)               { waveletFilterFiles_.push_back(filterFile) ;}
  void addWaveletCorrFile(const std::string & corrFile)                   { waveletCorrFiles_.push_back(corrFile)     ;}
  void setWellMoveIntFile(int i, const std::string & wellMoveIntFile)     { wellMoveIntFile_[i]   = wellMoveIntFile   ;}
  void setFaciesEstIntFile(int i, const std::string & faciesEstIntFile)   { faciesEstIntFile_[i]  = faciesEstIntFile  ;}
  void setTimeSurfTopFile(const std::string & timeSurfTopFile)            { timeSurfTopFile_      = timeSurfTopFile   ;}
  void setDepthSurfTopFile(const std::string & depthSurfTopFile)          { depthSurfTopFile_     = depthSurfTopFile  ;}
  //void addMultizoneSurfaceFile(const std::string & fileName)              { multizoneSurfaceFiles_.push_back(fileName);}
  void setVelocityField(const std::string & velocityField)                { velocityField_        = velocityField     ;}
  void setBackFile(int i, const std::string & backFile)                   { backFile_[i]          = backFile          ;}
  void setBackVelFile(const std::string & backVelFile)                    { backVelFile_          = backVelFile       ;}
  void setReflMatrFile(const std::string & reflMatrFile)                  { reflMatrFile_         = reflMatrFile      ;}
  void setParamCovFile(const std::string & paramCovFile)                  { paramCovFile_         = paramCovFile      ;}
  void setTempCorrFile(const std::string & tempCorrFile)                  { tempCorrFile_         = tempCorrFile      ;}
  void setParamAutoCovFile(const std::string & paramAutoCovFile)          { paramAutoCovFile_     = paramAutoCovFile  ;}
  void setRefSurfaceFile(const std::string & refSurfaceFile)              { refSurfaceFile_       = refSurfaceFile    ;}
  void setInputDirectory(std::string inputDirectory)                      { inputDirectory_       = inputDirectory    ;}
  void setPriorFaciesProb(std::string name,std::string file)              { priorFaciesProb_[name] = file             ;}
  void addTrendCubes(std::string trendParameterFile)                      { trendCubes_.push_back(trendParameterFile) ;}
  void addDefaultWaveletEstIntFileTop(void)                               { waveletEstIntFileTop_.push_back("")       ;}
  void addDefaultWaveletEstIntFileBase(void)                              { waveletEstIntFileBase_.push_back("")      ;}
  void setCorrDirFile(const std::string & interval_name, const std::string & file_name)            { corrDirFiles_[interval_name]       = file_name ;}
  void setCorrDirTopSurfaceFile(const std::string & interval_name, const std::string & file_name)  { corrDirTopFiles_[interval_name]    = file_name ;}
  void setCorrDirBaseSurfaceFile(const std::string & interval_name, const std::string & file_name) { corrDirBaseFiles_[interval_name]   = file_name ;}
  void setBaseTimeSurface(const std::string & interval_name, const std::string & file_name)        { timeSurfBaseFiles_[interval_name]  = file_name ;}
  void setBaseDepthSurface(const std::string & interval_name, const std::string & file_name)       { depthSurfBaseFiles_[interval_name] = file_name ;}
  std::string addInputPathAndCheckFiles();

  void clearTimeLapse(void)                                               { seismicFiles_.clear();
                                                                            waveletShiftFiles_.clear();
                                                                            waveletScaleFiles_.clear();
                                                                            waveletFiles_.clear();
                                                                            localNoiseFiles_.clear()                  ;}

  void addTimeLapse(void)                                                { timeLapseSeismicFiles_.push_back(seismicFiles_);
                                                                           timeLapseWaveletShiftFiles_.push_back(waveletShiftFiles_);
                                                                           timeLapseWaveletScaleFiles_.push_back(waveletScaleFiles_);
                                                                           timeLapseWaveletFiles_.push_back(waveletFiles_);
                                                                           timeLapseLocalNoiseFiles_.push_back(localNoiseFiles_)   ;}

  void clearTimeLapseTravelTime()                                        { travelTimeHorizons_.clear()                             ;}

  void addTimeLapseTravelTime()                                          { timeLapseTravelTimeHorizons_.push_back(travelTimeHorizons_);}

private:
  std::string addPathAndCheck(std::string & fileName, const bool possiblyNumber = false);

  std::string                seedFile_;              ///< File specifying the seed
  std::vector<std::string>   wellFiles_;             ///< File names: wells
  std::vector<std::string>   seismicFiles_;          ///< File names: seismic data
  std::vector<std::string>   rmsVelocities_;         ///< File names: RMS velocities U^2, vector over time lapses
  std::vector<std::string>   travelTimeHorizons_;    ///< File names: Horizons used in travel time inversion
  std::vector<std::string>   gravimetricData_;       ///< File names: Gravimetric data
  std::vector<std::string>   waveletFiles_;          ///< File names: wavelets
  std::vector<std::string>   waveletShiftFiles_;     ///< File names: wavelets
  std::vector<std::string>   waveletScaleFiles_;     ///< File names: wavelets
  std::vector<std::string>   waveletEstIntFileTop_;  ///< File names: Top of wavelet estimation interval for each time lapse
  std::vector<std::string>   waveletEstIntFileBase_; ///< File names: Base of wavelet estimation interval for each time lapse
  std::vector<std::string>   waveletFilterFiles_;    ///< File names: Filter files for 3D wavelets
  std::vector<std::string>   waveletCorrFiles_;      ///< File names: Filter with correction factors for 3D wavelets
  std::vector<std::string>   wellMoveIntFile_;       ///< File names: Well move interval
  std::vector<std::string>   faciesEstIntFile_;      ///< File names: Facies estimation interval
  //std::vector<std::string>   multizoneSurfaceFiles_; ///< File names: top and base surfaces for multizone background model
  std::string                velocityField_;         ///< File names: velocity field, or command
  std::vector<std::string>   backFile_;              ///< File names (temporarily stored).
  std::string                backVelFile_;           ///< File names Vp-velocity for background modelling
  std::string                reflMatrFile_;          ///< File name for reflection matrix file.
  std::string                paramCovFile_;          ///< File name for correlation between parameters.
  std::string                tempCorrFile_;          ///< File name for temporal parameter correlations.
  std::string                paramAutoCovFile_;      ///< File name for parameter autocovariance.
  std::string                refSurfaceFile_;        ///< File name for reference time surface corresponding to z0.

  std::string                        timeSurfTopFile_;    ///< File names: top time surface
  std::map<std::string, std::string> timeSurfBaseFiles_;  ///< File names: base time surfaces per interval
  std::string                        depthSurfTopFile_;   ///< File names: top depth surface
  std::map<std::string, std::string> depthSurfBaseFiles_; ///< File names: base depth surfaces per interval

  std::vector<std::string>           trendCubes_;         ///< File name for the trend cubes in the rock physics model
  std::string                        inputDirectory_;     ///< Base directory for input files.
  std::vector<std::string>           localNoiseFiles_;    ///< File names: local noise
  std::map<std::string, std::string> priorFaciesProb_;    ///< File names for locally varying prior facies probability
  std::string                        areaSurfaceFile_;

  std::map<std::string, std::string>     corrDirFiles_;                ///< File names: Map between interval name and correlation direction
  std::map<std::string, std::string>     corrDirTopFiles_;             ///< File names: Map between interval name and correlation direction for top-surface
  std::map<std::string, std::string>     corrDirBaseFiles_;            ///< File names: Map between interval name and correlation direction for base-surface

  std::vector<std::vector<std::string> > timeLapseWaveletShiftFiles_;  ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseWaveletScaleFiles_;  ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseWaveletFiles_;       ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseSeismicFiles_;       ///< File names: seismic data for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseLocalNoiseFiles_;    ///< File names: local noise for each time lapse and angle gather

  std::vector<std::vector<std::string> > timeLapseTravelTimeHorizons_; ///< File names: Horizons for each travel time time lapse
};

#endif
