/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef PARAMETEROUTPUT_H
#define PARAMETEROUTPUT_H

#include <string>

#include "src/definitions.h"
#include "libs/fft/include/fftw.h"

//class FFTFileGrid;
//class FFTGrid;
class Simbox;
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
                                   const ModelSettings * model_settings,
                                   StormContGrid       * vp,
                                   StormContGrid       * vs,
                                   StormContGrid       * rho,
                                   int                   output_flag,
                                   int                   sim_num,
                                   bool                  kriged);

  static void      WriteToFile(const Simbox        * simbox,
                               GridMapping         * time_depth_mapping,
                               const ModelSettings * model_settings,
                               StormContGrid       * grid,
                               const std::string   & file_name,
                               const std::string   & sgri_label,
                               bool                  padding = false);

  //static void      WriteToFile(const Simbox        * simbox,
  //                             GridMapping         * time_depth_mapping,
  //                             const ModelSettings * model_settings,
  //                             FFTGrid             * grid,
  //                             const std::string   & file_name,
  //                             const std::string   & sgri_label,
  //                             bool                  padding = false);

  static void     WriteFile(const ModelSettings     * model_settings,
                            StormContGrid           * storm_grid,
                            const std::string       & f_name,
                            const std::string       & sub_dir,
                            const Simbox            * simbox,
                            //Note: Used to mark synthetic seismic, where value in cell is actually value at top of cell.
                            //Real seismic data when resampled are values at base of cell, and must be shifted an index before this is used.
                            bool                      is_seismic = false,
                            const std::string         label = "NO_LABEL",
                            const GridMapping       * depth_map = NULL,
                            bool                      padding = false);

private:

  static void      ComputeAcousticImpedance(const Simbox        * simbox,
                                            GridMapping         * time_depth_mapping,
                                            const ModelSettings * model_settings,
                                            StormContGrid       * vp,
                                            StormContGrid       * rho,
                                            const std::string   & file_name);

  static void      ComputeShearImpedance(const Simbox        * simbox,
                                         GridMapping         * time_depth_mapping,
                                         const ModelSettings * model_settings,
                                         StormContGrid       * vs,
                                         StormContGrid       * rho,
                                         const std::string   & file_name);

  static void     ComputeVpVsRatio(const Simbox        * simbox,
                                   GridMapping         * time_depth_mapping,
                                   const ModelSettings * model_settings,
                                   StormContGrid       * vp,
                                   StormContGrid       * vs,
                                   const std::string   & file_name);

  static void      ComputePoissonRatio(const Simbox        * simbox,
                                       GridMapping         * time_depth_mapping,
                                       const ModelSettings * model_settings,
                                       StormContGrid       * vp,
                                       StormContGrid       * vs,
                                       const std::string   & file_name);

  static void      ComputeLameMu(const Simbox        * simbox,
                                 GridMapping         * time_depth_mapping,
                                 const ModelSettings * model_settings,
                                 StormContGrid       * vs,
                                 StormContGrid       * rho,
                                 const std::string   & file_name);

  static void      ComputeLameLambda(const Simbox        * simbox,
                                     GridMapping         * time_depth_mapping,
                                     const ModelSettings * model_settings,
                                     StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     StormContGrid       * rho,
                                     const std::string   & file_name);

  static void      ComputeMuRho(const Simbox        * simbox,
                                GridMapping         * time_depth_mapping,
                                const ModelSettings * model_settings,
                                StormContGrid       * vp,
                                StormContGrid       * vs,
                                StormContGrid       * rho,
                                const std::string   & file_name);

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
                                           float                 z0,
                                           bool                  is_depth);

  static void     SeismicShift(NRLib::Grid<float> * grid);
};
#endif
