/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef PARAMETEROUTPUT_H
#define PARAMETEROUTPUT_H

#include <string>

class FFTFileGrid;
class FFTGrid;
class Simbox;
class ModelGeneral;
class ModelSettings;

class ParameterOutput
{
public:
  //Conventions for writeParameters:
  // simNum = -1 indicates prediction, otherwise filename ends with n+1.
  // All grids are in normal domain, and on log scale.
  static void      writeParameters(const Simbox  * simbox,
                                   ModelGeneral  * modelGeneral,
                                   const ModelSettings * modelSettings,
                                   FFTGrid       * postAlpha,
                                   FFTGrid       * postBeta,
                                   FFTGrid       * postRho,
                                   int             outputFlag,
                                   bool            fileGrid,
                                   int             simNum,
                                   bool            kriged);
  static void      writeToFile(const Simbox      * simbox,
                               ModelGeneral      * modelGeneral,
                               const ModelSettings     * modelSettings,
                               FFTGrid           * grid,
                               const std::string & fileName,
                               const std::string & sgriLabel,
                               bool padding=false);


private:
  static void      computeAcousticImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                            FFTGrid * Alpha, FFTGrid * Rho,
                                            bool fileGrid, const std::string & fileName);
  static void      computeShearImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                         FFTGrid * Beta, FFTGrid * Rho,
                                         bool fileGrid, const std::string & fileName);
  static void      computeVpVsRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                    FFTGrid * Alpha, FFTGrid * Beta,
                                    bool fileGrid, const std::string & fileName);
  static void      computePoissonRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                       FFTGrid * Alpha, FFTGrid * Beta,
                                       bool fileGrid, const std::string & fileName);
  static void      computeLameMu(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                 FFTGrid * Beta, FFTGrid * Rho,
                                 bool fileGrid, const std::string & FileName);
  static void      computeLameLambda(const Simbox * simbox,ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                     FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                                     bool fileGrid, const std::string & fileName);
  static void      computeMuRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                                bool fileGrid, const std::string & fileName);
  static void      computeLambdaRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
                                    FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                                    bool fileGrid, const std::string & fileName);

  static FFTGrid * createFFTGrid(FFTGrid * referenceGrid, bool fileGrid);
};
#endif
