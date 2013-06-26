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
#include "lib/utils.h"
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

  void        SetEstimationSimbox();
  void        SetFullInversionVoulume();

  //GET FUNCTIONS

  const Simbox          & GetEstimationSimbox()     const { return estimation_simbox_     ;}
  const NRLib::Volume   & GetFullInversionVolume()  const { return full_inversion_volume_ ;}



private:

  void GetGeometryFromGridOnFile(const std::string         grid_file,
                                 const TraceHeaderFormat * thf,
                                 SegyGeometry           *& geometry,
                                 std::string             & err_text);

  SegyGeometry * GetGeometryFromCravaFile(const std::string & file_name);

  SegyGeometry * GetGeometryFromStormFile(const std::string    & file_name,
                                       std::string          & err_text,
                                       bool                   scale = false);

  bool CreateOuterTemporarySimbox(ModelSettings           * model_settings,
                                  InputFiles              * input_files,
                                  Simbox                  & estimation_simbox,
                                  NRLib::Volume           & full_inversion_interval,
                                  std::string             & err_text);

  void WriteAreas(const SegyGeometry  * area_params,
                  Simbox              * time_simbox,
                  std::string         & err_text);

  void FindSmallestSurfaceGeometry(const double   x0,
                                   const double   y0,
                                   const double   lx,
                                   const double   ly,
                                   const double   rot,
                                   double       & x_min,
                                   double       & y_min,
                                   double       & x_max,
                                   double       & y_max);

  void SetSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
                                    NRLib::Volume                  & full_inversion_volume,
                                    const InputFiles               * input_files,
                                    const ModelSettings            * model_settings,
                                    std::string                    & err_text,
                                    bool                           & failed);

  void SetSurfacesSingleInterval(Simbox                           & estimation_simbox,
                                 NRLib::Volume                    & full_inversion_volume,
                                 const std::vector<std::string>   & surf_file,
                                 ModelSettings                    * model_settings,
                                 std::string                      & err_text,
                                 bool                             & failed);

  bool ReadSeismicData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles);

  bool ReadWellData(ModelSettings     * modelSettings,
                    InputFiles        * inputFiles,
                    std::string       & err_text);

  bool BlockWellsForEstimation(const ModelSettings                            * const model_settings,
                               //const InputFiles                               * const input_files,
                               const Simbox                                   & estimation_simbox,
                               const std::vector<NRLib::Well>                 & wells,
                               std::vector<BlockedLogsCommon *>               & blocked_logs_common,
                               std::string                                    & err_text);

  bool       CheckThatDataCoverGrid(const SegY   * segy,
                                    float         offset,
                                    const Simbox * timeCutSimbox,
                                    float         guard_zone);

  void ProcessLogsNorsarWell(NRLib::Well    & new_well,
                             std::string    & error_text,
                             bool           & failed);

  void ProcessLogsRMSWell(NRLib::Well       & new_well,
                          std::string       & error_text,
                          bool              & failed);


  //void        readNorsarWell(const std::string              & wellFileName,
  //                           NRLib::Well                    & new_well,
  //                           const std::vector<std::string> & logNames,
  //                           const std::vector<bool>        & inverseVelocity,
  //                           bool                             faciesLogGiven,
  //                           std::string                    & error);

  //void          readRMSWell(const std::string              & wellFileName,
  //                          NRLib::Well                    & new_well,
  //                          const std::vector<std::string> & logNames,
  //                          const std::vector<bool>        & inverseVelocity,
  //                          bool                             faciesLogGiven,
  //                          std::string                    & error);

  bool SetupReflectionMatrixAndTempWavelet(ModelSettings * model_settings,
                                           InputFiles *    input_files);

  float  ** ReadMatrix(const std::string & fileName,
                       int                 n1,
                       int                 n2,
                       const std::string & readReason,
                       std::string       & errText);

  void  SetupDefaultReflectionMatrix(float             **& reflectionMatrix,
                                     double                vsvp,
                                     const ModelSettings * modelSettings,
                                     int                   numberOfAngles,
                                     int                   thisTimeLapse);

  bool waveletHandling(ModelSettings * model_settings,
                       InputFiles * input_files);

  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();

  int ComputeTime(int year, int month, int day) const;

  // Bool variables indicating whether the corresponding data processing
  // succeeded
  bool outer_temp_simbox_;
  bool read_seismic_;
  bool read_wells_;
  bool block_wells_;
  bool setup_reflection_matrix_;
  bool temporary_wavelet_;
  bool optimize_well_location_;
  bool wavelet_estimation_shape_;
  bool prior_corr_estimation_;
  bool setup_estimation_rock_physics_;

  Simbox          estimation_simbox_;
  NRLib::Volume   full_inversion_volume_;

  //std::vector<SeismicStorage> seismic_data_;
  std::map<int, std::vector<SeismicStorage> >   seismic_data_;
  std::vector<NRLib::Well>                      wells_;
  std::vector<BlockedLogsCommon *>              blocked_logs_common_;   ///< Blocked wells for estimation

  //Well variables not contained in NRlib::Well
  //std::map<std::string, int>                        timemissing_;
  //std::map<std::string, double>                     xpos0_;
  //std::map<std::string, double>                     ypos0_;
  //std::map<std::string, std::vector<std::string> >  faciesnames_;
  //std::map<std::string, std::vector<int> >          faciesNr_;
  //std::map<std::string, int>                        faciesok_; //Bool?
  //std::map<std::string, int>                        nFacies_;


  std::map<int, float **> reflectionMatrix_;
  bool        reflection_matrix_from_file_; //False: created from global vp/vs

  std::vector<Wavelet*> temporary_wavelets_; //One wavelet per angle

  //float                  ** reflectionMatrix_;


  //Simbox                     * estimation_simbox_;
  //NRLib::Volume              * full_inversion_volume_;

};

#endif
