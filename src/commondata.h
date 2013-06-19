/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef COMMONDATA_H
#define COMMONDATA_H

#include <math.h>
#include <string>
#include "src/simbox.h"
#include "src/welldata.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "src/blockedlogscommon.h"
#include "src/tasklist.h"

#include "src/seismicstorage.h"

class ModelSettings;
class ModelGeneral;
class InputFiles;
//class FFTGrid;
//class Surface;

class CommonData{
public:
  CommonData(ModelSettings    * model_settings,
             InputFiles       * input_files);

  ~ CommonData();

  //SET FUNCTIONS

  void        setEstimationSimbox();
  void        setFullInversionVoulume();

  //GET FUNCTIONS

  const Simbox          & getEstimationSimbox()     const { return estimation_simbox_     ;}
  const NRLib::Volume   & getFullInversionVolume()  const { return full_inversion_volume_ ;}
  


private:

  void getGeometryFromGridOnFile(const std::string                grid_file,
                                        const TraceHeaderFormat * thf,
                                        SegyGeometry           *& geometry,
                                        std::string             & err_text);

  SegyGeometry * geometryFromCravaFile(const std::string & file_name);

  SegyGeometry * geometryFromStormFile(const std::string    & file_name,
                                       std::string          & err_text,
                                       bool                   scale = false);

  bool createOuterTemporarySimbox(ModelSettings           * model_settings,
                                  InputFiles              * input_files,
                                  Simbox                  & estimation_simbox,
                                  NRLib::Volume           & full_inversion_interval,
                                  std::string             & err_text);

  void writeAreas(const SegyGeometry  * area_params,
                  Simbox              * time_simbox,
                  std::string         & err_text);

  void findSmallestSurfaceGeometry(const double   x0,
                                   const double   y0,
                                   const double   lx,
                                   const double   ly,
                                   const double   rot,
                                   double       & x_min,
                                   double       & y_min,
                                   double       & x_max,
                                   double       & y_max);

  void setSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
                                    NRLib::Volume                  & full_inversion_volume,
                                    const InputFiles               * input_files,
                                    const ModelSettings            * model_settings,
                                    std::string                    & err_text,
                                    bool                           & failed);

  void setSurfacesSingleInterval(Simbox                           & estimation_simbox,
                                 NRLib::Volume                    & full_inversion_volume,
                                 const std::vector<std::string>   & surf_file,
                                 ModelSettings                    * model_settings,
                                 std::string                      & err_text,
                                 bool                             & failed);

  bool readSeismicData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles);

  bool readWellData(ModelSettings     * modelSettings,
                    InputFiles        * inputFiles);

  bool blockWellsForEstimation(const ModelSettings                            * const model_settings, 
                               const InputFiles                               * const input_files, 
                               const Simbox                                   & estimation_simbox, 
                               const NRLib::Volume                            & full_inversion_volume, 
                               const std::vector<NRLib::Well>                 & wells,
                               std::vector<BlockedLogsCommon *>               & blocked_logs_common,
                               std::string                                    & err_text);

  bool setupReflectionMatrixAndTempWavelet();
  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();

  int computeTime(int year, int month, int day) const;

  // Bool variables indicating whether the corresponding data processing
  // succeeded
  bool outer_temp_simbox_;
  bool read_seismic_;
  bool read_wells_;
  bool block_wells_;
  bool setup_reflection_matrix_;
  bool optimize_well_location_;
  bool wavelet_estimation_shape_;
  bool prior_corr_estimation_;
  bool setup_estimation_rock_physics_;

  Simbox          estimation_simbox_;
  NRLib::Volume   full_inversion_volume_;

  //std::vector<NRLib::SegY *>   segy_files_;
  //std::vector<StormContGrid *> stormgrids_;

  std::vector<SeismicStorage> seismic_data_;

  //FFTGrid                   ** seisCube_;                 ///< Seismic data cubes
  std::vector<WellData *>           welldata_;              ///< Well data
  std::vector<BlockedLogsCommon *>  blocked_logs_common_;   ///< Blocked wells for estimation
  std::vector<NRLib::Well>          wells_;
  //Simbox                     * estimation_simbox_;
  //NRLib::Volume              * full_inversion_volume_;

};

#endif
