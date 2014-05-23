/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef PARAMETEROUTPUT_H
#define PARAMETEROUTPUT_H

#include <string>

#include "src/definitions.h"

class FFTFileGrid;
class FFTGrid;
class Simbox;
//class ModelGeneral;
class ModelSettings;
class GridMapping;

class ParameterOutput
{
public:
  //Conventions for writeParameters:
  // simNum = -1 indicates prediction, otherwise filename ends with n+1.
  // All grids are in normal domain, and on log scale.
  static void      WriteParameters(const Simbox        * simbox,
                                   GridMapping         * time_depth_mapping,
                                   const ModelSettings * modelSettings,
                                   StormContGrid       * postAlpha,
                                   StormContGrid       * postBeta,
                                   StormContGrid       * postRho,
                                   int                   outputFlag,
                                   //bool                  fileGrid,
                                   int                   simNum,
                                   bool                  kriged);

  static void      WriteToFile(const Simbox        * simbox,
                               GridMapping         * time_depth_mapping,
                               const ModelSettings * modelSettings,
                               StormContGrid       * grid,
                               const std::string   & fileName,
                               const std::string   & sgriLabel,
                               bool                  padding=false);

  //static void      WriteToFile(const Simbox        * simbox,
  //                             GridMapping         * time_depth_mapping,
  //                             const ModelSettings * model_settings,
  //                             FFTGrid             * grid,
  //                             const std::string   & file_name,
  //                             const std::string   & sgri_label,
  //                             bool                  padding = false);

  static void     WriteFile(const ModelSettings     * model_settings,
                            StormContGrid           * storm_grid,
                            const std::string       & fName,
                            const std::string       & subDir,
                            const Simbox            * simbox,
                            const std::string         label = "NO_LABEL",
                            const float               z0 = 0.0,
                            const GridMapping       * depthMap = NULL,
                            const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS),
                            bool                      padding = false);

private:
  //static void      computeAcousticImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                          FFTGrid * Alpha, FFTGrid * Rho,
  //                                          bool fileGrid, const std::string & fileName);
  static void      ComputeAcousticImpedance(const Simbox        * simbox,
                                            GridMapping         * time_depth_mapping,
                                            const ModelSettings * model_settings,
                                            StormContGrid       * vp,
                                            StormContGrid       * rho,
                                            const std::string   & file_name);

  //static void      computeShearImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                       FFTGrid * Beta, FFTGrid * Rho,
  //                                       bool fileGrid, const std::string & fileName);
  static void      ComputeShearImpedance(const Simbox        * simbox,
                                         GridMapping         * time_depth_mapping,
                                         const ModelSettings * model_settings,
                                         StormContGrid       * vs,
                                         StormContGrid       * rho,
                                         const std::string   & file_name);

  //static void      computeVpVsRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                  FFTGrid * Alpha, FFTGrid * Beta,
  //                                  bool fileGrid, const std::string & fileName);
  static void     ComputeVpVsRatio(const Simbox        * simbox,
                                   GridMapping         * time_depth_mapping,
                                   const ModelSettings * model_settings,
                                   StormContGrid       * vp,
                                   StormContGrid       * vs,
                                   const std::string   & file_name);

  //static void      computePoissonRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                     FFTGrid * Alpha, FFTGrid * Beta,
  //                                     bool fileGrid, const std::string & fileName);
  static void      ComputePoissonRatio(const Simbox        * simbox,
                                       GridMapping         * time_depth_mapping,
                                       const ModelSettings * model_settings,
                                       StormContGrid       * vp,
                                       StormContGrid       * vs,
                                       const std::string   & file_name);

  //static void      computeLameMu(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                               FFTGrid * Beta, FFTGrid * Rho,
  //                               bool fileGrid, const std::string & FileName);
  static void      ComputeLameMu(const Simbox        * simbox,
                                 GridMapping         * time_depth_mapping,
                                 const ModelSettings * model_settings,
                                 StormContGrid       * vs,
                                 StormContGrid       * rho,
                                 const std::string   & file_name);

  //static void      computeLameLambda(const Simbox * simbox,ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                   FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
  //                                   bool fileGrid, const std::string & fileName);
  static void      ComputeLameLambda(const Simbox        * simbox,
                                     GridMapping         * time_depth_mapping,
                                     const ModelSettings * model_settings,
                                     StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     StormContGrid       * rho,
                                     const std::string   & file_name);

  //static void      computeMuRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                              FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
  //                              bool fileGrid, const std::string & fileName);
  static void      ComputeMuRho(const Simbox        * simbox,
                                GridMapping         * time_depth_mapping,
                                const ModelSettings * model_settings,
                                StormContGrid       * vp,
                                StormContGrid       * vs,
                                StormContGrid       * rho,
                                const std::string   & file_name);

  //static void      computeLambdaRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * modelSettings,
  //                                  FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
  //                                  bool fileGrid, const std::string & fileName);
  static void      ComputeLambdaRho(const Simbox        * simbox,
                                    GridMapping         * time_depth_mapping,
                                    const ModelSettings * model_settings,
                                    StormContGrid       * vp,
                                    StormContGrid       * vs,
                                    StormContGrid       * rho,
                                    const std::string   & file_name);

  //static FFTGrid * createFFTGrid(FFTGrid * referenceGrid, bool fileGrid);

  static void      ExpTransf(StormContGrid * grid);

  static void      WriteResampledStormCube(const StormContGrid * storm_grid,
                                           const GridMapping   * gridmapping,
                                           const std::string   & file_name,
                                           const Simbox        * simbox,
                                           const int             format,
                                           float                 z0);
};
#endif
