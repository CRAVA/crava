#ifndef PARAMETEROUTPUT_H
#define PARAMETEROUTPUT_H

class FFTFileGrid;
class FFTGrid;
class Simbox;
class Model;

class ParameterOutput
{
public:
  //Conventions for writeParameters:
  // simNum = -1 indicates prediction, otherwise filename ends with n+1.
  // All grids are in normal domain, and on log scale.
  static void      writeParameters(const Simbox  * simbox,
                                   Model         * model,
                                   FFTGrid       * postAlpha, 
                                   FFTGrid       * postBeta, 
                                   FFTGrid       * postRho, 
                                   int             outputFlag,
                                   int             fileGrid,
                                   int             simNum); 
  static void      writeToFile(const Simbox      * simbox, 
                               Model             * model, 
                               FFTGrid           * grid, 
                               const std::string & fileName, 
                               const std::string & sgriLabel);
  
private: 
  static void      computeAcousticImpedance(const Simbox * simbox, Model * model, 
                                            FFTGrid * Alpha, FFTGrid * Rho, 
                                            int fileGrid, const std::string & fileName);
  static void      computeShearImpedance(const Simbox * simbox, Model * model,
                                         FFTGrid * Beta, FFTGrid * Rho, 
                                         int fileGrid, const std::string & fileName);
  static void      computeVpVsRatio(const Simbox * simbox, Model * model,
                                    FFTGrid * Alpha, FFTGrid * Beta, 
                                    int fileGrid, const std::string & fileName);
  static void      computePoissonRatio(const Simbox * simbox, Model * model,
                                       FFTGrid * Alpha, FFTGrid * Beta, 
                                       int fileGrid, const std::string & fileName);
  static void      computeLameMu(const Simbox * simbox, Model * model,
                                 FFTGrid * Beta, FFTGrid * Rho, 
                                 int fileGrid, const std::string & FileName);
  static void      computeLameLambda(const Simbox * simbox,Model * model,
                                     FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, 
                                     int fileGrid, const std::string & fileName);
  static void      computeMuRho(const Simbox * simbox, Model * model,
                                FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, 
                                int fileGrid, const std::string & fileName);
  static void      computeLambdaRho(const Simbox * simbox, Model * model, 
                                    FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, 
                                    int fileGrid, const std::string & fileName);
  
  static FFTGrid * createFFTGrid(FFTGrid * referenceGrid, int fileGrid);
};
#endif
