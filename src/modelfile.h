#ifndef MODELFILE_H
#define MODELFILE_H

#include <stdio.h>

class Vario;
class ModelSettings;
class InputFiles;

class ModelFile
{
public:
  ModelFile(char * fileName);
  ~ModelFile(void);

  ModelSettings  * getModelSettings(void)         const { return modelSettings_ ;}
  InputFiles     * getInputFiles(void)            const { return inputFiles_    ;}

  bool             getParsingFailed(void)         const { return failed_        ;}

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
  int              readCommandBackgroundVelocity(char ** params, int & pos, char * errText);
  void             parse_old(char * fileName, char **& params, int & nParam);
  void             parse(char * fileName, char **& params, int & nParam);
  int              getParNum(char ** params, int pos, int & error, char * errText,
                             const char * command, int min, int max = 0);
  int              checkFileOpen(char ** fNames, int nFiles, const char * command, char * errText, int start = 0,
                                 bool details = true);
  int              checkFileOpen(const std::string & fName, const char * command, char * errText);


  Vario          * createVario(char ** param, int nPar, const char * command, char * errText);
  void             trimString(std::string & str);
  void             findApostropheEnclosedParts(std::string & str);
  void             splitString(std::vector<std::string> & entries, const std::string & str);
  void             reintroduceBlanks(std::vector<std::string> & entries, unsigned int iStart);

  ModelSettings  * modelSettings_;
  InputFiles     * inputFiles_;

  bool             failed_;                // Indicates whether errors ocuured during construction. 
};

#endif
