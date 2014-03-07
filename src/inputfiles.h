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
  const std::vector<std::string> & getTimeSurfFiles(void)        const { return timeSurfFiles_        ;}
  const std::string              & getTimeSurfFile(int i)        const { return timeSurfFiles_[i]     ;}
  const std::vector<std::string> & getDepthSurfFiles(void)       const { return depthSurfFiles_       ;}
  const std::string              & getDepthSurfFile(int i)       const { return depthSurfFiles_[i]    ;}
  const std::string              & getVelocityField(void)        const { return velocityField_        ;}
  const std::string              & getBackFile(int i)            const { return backFile_[i]          ;}
  const std::string              & getBackVelFile(void)          const { return backVelFile_          ;}
  const std::string              & getCorrDirFile(void)          const { return corrDirFile_          ;}
  const std::string              & getCorrDirTopFile(void)       const { return corrDirTopFile_       ;}
  const std::string              & getCorrDirBaseFile(void)      const { return corrDirBaseFile_      ;}
  const std::map<std::string, std::string> & getCorrDirIntervalFiles(void)                        const { return interval_corrDirFiles_     ;}
  const std::string                        & getCorrDirIntervalFile(std::string name)             const { return interval_corrDirFiles_.find(name)->second     ;}
  const std::map<std::string, std::string> & getCorrDirIntervalTopSurfaceFiles(void)              const { return interval_corrDirTopFiles_  ;}
  const std::string                        & getCorrDirIntervalTopSurfaceFile(std::string name)   const { return interval_corrDirTopFiles_.find(name)->second  ;}
  const std::map<std::string, std::string> & getCorrDirIntervalBaseSurfaceFiles(void)             const { return interval_corrDirBaseFiles_ ;}
  const std::string                        & getCorrDirIntervalBaseSurfaceFile(std::string name)  const { return interval_corrDirBaseFiles_.find(name)->second ;}
  const std::string              & getReflMatrFile(void)         const { return reflMatrFile_         ;}
  const std::string              & getParamCovFile(void)         const { return paramCovFile_        ;}
  const std::string              & getTempCorrFile(void)         const { return tempCorrFile_         ;}
  const std::string              & getRefSurfaceFile(void)       const { return refSurfaceFile_       ;}
  const std::string              & getInputDirectory(void)       const { return inputDirectory_       ;}
  const std::map<std::string,std::string> & getPriorFaciesProbFile(void)   const {return priorFaciesProb_ ;}
  const std::string              & getAreaSurfaceFile(void)      const { return areaSurfaceFile_      ;}
  const std::vector<std::string> & getMultizoneSurfaceFiles()    const { return multizoneSurfaceFiles_;}
  const std::string              & getTrendCube(int i)           const { return trendCubes_[i]        ;}
  const std::string              & getGravimetricData(int i)     const { return gravimetricData_[i]   ;}
  const std::vector<std::string> & getSeismicFiles(void)         const { return seismicFiles_         ;}
  const std::vector<std::vector<std::string> > & getTimeLapseSeismicFiles(void) const { return timeLapseSeismicFiles_ ;}
  const std::map<std::string, std::string> & getIntervalBaseTimeSurfaces(void)             const { return interval_base_time_surface_                     ;}
  const std::string                        & getIntervalBaseTimeSurface(std::string name)  const { return interval_base_time_surface_.find(name)->second  ;}
  const std::map<std::string, std::string> & getIntervalBaseDepthSurfaces(void)            const { return interval_base_depth_surface_                    ;}
  const std::string                        & getIntervalBaseDepthSurface(std::string name) const { return interval_base_depth_surface_.find(name)->second ;}

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
  void addTimeSurfFile(const std::string & timeSurfFile)                  { timeSurfFiles_.push_back(timeSurfFile)    ;}
  void setDepthSurfFile(int i, const std::string & depthSurfFile)         { depthSurfFiles_[i]    = depthSurfFile     ;}
  void addMultizoneSurfaceFile(const std::string & fileName)              { multizoneSurfaceFiles_.push_back(fileName);}
  void setVelocityField(const std::string & velocityField)                { velocityField_        = velocityField     ;}
  void setBackFile(int i, const std::string & backFile)                   { backFile_[i]          = backFile          ;}
  void setBackVelFile(const std::string & backVelFile)                    { backVelFile_          = backVelFile       ;}
  void setReflMatrFile(const std::string & reflMatrFile)                  { reflMatrFile_         = reflMatrFile      ;}
  void setCorrDirFile(const std::string & corrDirFile)                    { corrDirFile_          = corrDirFile       ;}
  void setCorrDirTopSurfaceFile(const std::string & corrDirFile)          { corrDirTopFile_       = corrDirFile       ;}
  void setCorrDirBaseSurfaceFile(const std::string & corrDirFile)         { corrDirBaseFile_      = corrDirFile       ;}
  void setCorrDirIntervalFile(const std::string & interval_name, const std::string & file_name)            { interval_corrDirFiles_[interval_name]     = file_name ;}
  void setCorrDirIntervalTopSurfaceFile(const std::string & interval_name, const std::string & file_name)  { interval_corrDirTopFiles_[interval_name]  = file_name ;}
  void setCorrDirIntervalBaseSurfaceFile(const std::string & interval_name, const std::string & file_name) { interval_corrDirBaseFiles_[interval_name] = file_name ;}
  void setParamCovFile(const std::string & paramCovFile)                  { paramCovFile_        = paramCovFile     ;}
  void setTempCorrFile(const std::string & tempCorrFile)                  { tempCorrFile_         = tempCorrFile      ;}
  void setRefSurfaceFile(const std::string & refSurfaceFile)              { refSurfaceFile_       = refSurfaceFile    ;}
  void setInputDirectory(std::string inputDirectory)                      { inputDirectory_       = inputDirectory    ;}
  void setPriorFaciesProb(std::string name,std::string file)              { priorFaciesProb_[name] = file             ;}
  void addTrendCubes(std::string trendParameterFile)                      { trendCubes_.push_back(trendParameterFile) ;}
  std::string addInputPathAndCheckFiles();
  void setIntervalBaseTimeSurface(const std::string & interval_name, const std::string & file_name)  { interval_base_time_surface_[interval_name] = file_name   ;}
  void setIntervalBaseDepthSurface(const std::string & interval_name, const std::string & file_name) { interval_base_depth_surface_[interval_name] = file_name  ;}
  void addDefaultWaveletEstIntFileTop(void)                               { waveletEstIntFileTop_.push_back("")       ;}
  void addDefaultWaveletEstIntFileBase(void)                              { waveletEstIntFileBase_.push_back("")      ;}

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
  std::vector<std::string>   timeSurfFiles_;         ///< File names: top and base time surfaces
  std::vector<std::string>   depthSurfFiles_;        ///< File names: top and base depth surfaces
  std::vector<std::string>   multizoneSurfaceFiles_; ///< File names: top and base surfaces for multizone background model
  std::string                velocityField_;         ///< File names: velocity field, or command
  std::vector<std::string>   backFile_;              ///< File names (temporarily stored).
  std::string                backVelFile_;           ///< File names Vp-velocity for background modelling
  std::string                reflMatrFile_;          ///< File name for reflection matrix file.
  std::string                corrDirFile_;           ///< File name for correlation direction
  std::string                corrDirTopFile_;        ///< File name for correlation direction for top-surface
  std::string                corrDirBaseFile_;       ///< File name for correlation direction for base-surface
  std::string                paramCovFile_;          ///< File name for correlation between parameters.
  std::string                tempCorrFile_;          ///< File name for temporal parameter correlations.
  std::string                refSurfaceFile_;        ///< File name for reference time surface corresponding to z0.

  std::vector<std::string>                         trendCubes_;               ///< File name for the trend cubes in the rock physics model

  std::string                inputDirectory_;        ///< Base directory for input files.
  std::vector<std::string>   localNoiseFiles_;       ///< File names: local noise
  std::map<std::string, std::string> priorFaciesProb_; ///< File names for locally varying prior facies probability
  std::string                areaSurfaceFile_;

  std::map<std::string, std::string>    interval_base_time_surface_;     ///< File names: Map between interval name and base time surface file name/value.
  std::map<std::string, std::string>    interval_base_depth_surface_;    ///< File names: Map between interval name and base depth surface file name.
  std::map<std::string, std::string>    interval_corrDirFiles_;           ///< File names: Map between interval name and correlation direction
  std::map<std::string, std::string>    interval_corrDirTopFiles_;        ///< File names: Map between interval name and correlation direction for top-surface
  std::map<std::string, std::string>    interval_corrDirBaseFiles_;       ///< File names: Map between interval name and correlation direction for base-surface

  std::vector<std::vector<std::string> > timeLapseWaveletShiftFiles_; ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseWaveletScaleFiles_; ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseWaveletFiles_;      ///< File names: wavelets for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseSeismicFiles_;      ///< File names: seismic data for each time lapse and angle gather
  std::vector<std::vector<std::string> > timeLapseLocalNoiseFiles_;   ///< File names: local noise for each time lapse and angle gather

  std::vector<std::vector<std::string> > timeLapseTravelTimeHorizons_; ///< File names: Horizons for each travel time time lapse
};

#endif
