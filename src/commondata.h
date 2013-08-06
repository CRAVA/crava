/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef COMMONDATA_H
#define COMMONDATA_H

#include "src/simbox.h"
#include "src/welldata.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "lib/utils.h"
#include "src/blockedlogscommon.h"
#include "src/tasklist.h"
#include "src/seismicstorage.h"
#include "src/multiintervalgrid.h"

class CravaTrend;
class InputFiles;
class ModelSettings;
class ModelGeneral;


class CommonData{
public:
  CommonData(ModelSettings    * model_settings,
             InputFiles       * input_files);

  ~ CommonData();

  //GET FUNCTIONS

  const Simbox                    & GetEstimationSimbox()     const { return estimation_simbox_     ;}
  const NRLib::Volume             & GetFullInversionVolume()  const { return full_inversion_volume_ ;}
  const std::vector<NRLib::Well>  & GetWells()                const { return wells_                 ;}


private:

  void LoadWellMoveInterval(const InputFiles             * input_files,
                            const Simbox                 * estimation_simbox,
                            std::vector<Surface *>       & well_move_interval,
                            std::string                  & err_text,
                            bool                         & failed);

  void  OptimizeWellLocations(ModelSettings                                 * model_settings,
                              InputFiles                                    * input_files,
                              const Simbox                                  * estimation_simbox,
                              //const NRLib::Volume                           & volume,
                              std::vector<NRLib::Well>                      & wells,
                              std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                              std::map<int, std::vector<SeismicStorage> >   & seismic_data,
                              std::map<int, float **>                       & reflection_matrix,
                              std::string                                   & err_text,
                              bool                                          & failed);

  void MoveWell(const NRLib::Well & well,
                const Simbox      * simbox,
                double              delta_X,
                double              delta_Y,
                double              k_move);

  void  CalculateDeviation(NRLib::Well            & new_well,
                           const ModelSettings    * const model_settings,
                           float                  & dev_angle,
                           Simbox                 * simbox,
                           int                      use_for_wavelet_estimation);

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

  void SetSurfacesMultipleIntervals(const ModelSettings            * const model_settings,
                                    NRLib::Volume                  & full_inversion_volume,
                                    Simbox                         & estimation_simbox,
                                    const InputFiles               * input_files,
                                    std::string                    & err_text,
                                    bool                           & failed);

  void SetSurfacesSingleInterval(const ModelSettings              * const model_settings,
                                 NRLib::Volume                    & full_inversion_volume,
                                 Simbox                           & estimation_simbox,
                                 const std::vector<std::string>   & surf_file,
                                 std::string                      & err_text,
                                 bool                             & failed);

  bool ReadSeismicData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles);

  bool ReadWellData(ModelSettings     * modelSettings,
                    Simbox            * estimation_simbox,
                    InputFiles        * inputFiles,
                    std::string       & err_text);

  bool BlockWellsForEstimation(const ModelSettings                            * const model_settings,
                               const Simbox                                   & estimation_simbox,
                               std::vector<NRLib::Well>                       & wells,
                               std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_common,
                               std::string                                    & err_text);

  bool       CheckThatDataCoverGrid(const SegY            * segy,
                                    float                   offset,
                                    //const Simbox          * timeCutSimbox,
                                    float                   guard_zone);

  void ProcessLogsNorsarWell(NRLib::Well                  & new_well,
                             std::string                  & error_text,
                             bool                         & failed);

  void ProcessLogsRMSWell(NRLib::Well                     & new_well,
                          std::string                     & error_text,
                          bool                            & failed);


  bool  SetupReflectionMatrixAndTempWavelet(ModelSettings  * model_settings,
                                            InputFiles     * input_files);

  float  ** ReadMatrix(const std::string                  & fileName,
                       int                                  n1,
                       int                                  n2,
                       const std::string                  & readReason,
                       std::string                        & errText);

  void  SetupDefaultReflectionMatrix(float              **& reflectionMatrix,
                                     double                 vsvp,
                                     const ModelSettings  * modelSettings,
                                     int                    numberOfAngles,
                                     int                    thisTimeLapse);

  bool  WaveletHandling(ModelSettings                     * model_settings,
                        InputFiles                        * input_files);

  void  SetupTrendCubes(ModelSettings                     * model_settings,
                        InputFiles                        * input_files,
                        MultiIntervalGrid                 * multiple_interval_grid,
                        std::string                       & error_text,
                        bool                              & failed);

  void SetupRockPhysics(const ModelSettings                                 * model_settings,
                        const InputFiles                                    * input_files,
                        const MultiIntervalGrid                             * multiple_interval_grid,
                        const std::vector<CravaTrend>                       & trend_cubes,
                        const std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                        std::string                                         & error_text,
                        bool                                                & failed);

  int Process1DWavelet(const ModelSettings                * modelSettings,
                       const InputFiles             * inputFiles,
                       const SeismicStorage           * seismic_data,
                       std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                       //std::vector<BlockedLogsCommon *> blocked_logs_common,
                       const std::vector<Surface *> & waveletEstimInterval,
                       std::string                  & err_text,
                       Wavelet                     *& wavelet,
                       Grid2D                      *& local_noise_scale,
                       Grid2D                      *& local_noise_shift,
                       Grid2D                      *& local_noise_estimate,
                       unsigned int                   i_timelapse,
                       unsigned int                   j_angle,
                       const float                    angle,
                       float                          sn_ratio,
                       bool                           estimate_wavlet,
                       bool                           use_ricker_wavelet,
                       bool                           use_local_noise);

  int Process3DWavelet(const ModelSettings                     * model_settings,
                       const InputFiles                        * input_files,
                       const SeismicStorage           * seismic_data,
                       //std::vector<BlockedLogsCommon *> blocked_logs_common,
                       std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                       const std::vector<Surface *>            & wavelet_estim_interval,
                       std::string                             & err_text,
                       Wavelet                                *& wavelet,
                       unsigned int                              i_timelapse,
                       unsigned int                              j_angle,
                       float                                     angle,
                       float                                     sn_ratio,
                       const NRLib::Grid2D<float>              & ref_time_grad_x,
                       const NRLib::Grid2D<float>              & ref_time_grad_y,
                       const std::vector<std::vector<double> > & t_grad_x,
                       const std::vector<std::vector<double> > & t_grad_y,
                       bool                                      estimate_wavelet);

  void FindWaveletEstimationInterval(InputFiles             * input_files,
                                     std::vector<Surface *> & wavelet_estim_interval,
                                     std::string            & err_text);

  void ComputeStructureDepthGradient(double                 v0,
                                     double                 radius,
                                     const Surface        * t0_surf,
                                     const Surface        * correlation_direction,
                                     NRLib::Grid2D<float> & structure_depth_grad_x,
                                     NRLib::Grid2D<float> & structure_depth_grad_y);

  void ComputeReferenceTimeGradient(const Surface       * t0_surf,
                                    NRLib::Grid2D<float> &ref_time_grad_x,
                                    NRLib::Grid2D<float> &ref_time_grad_y);

  void CalculateSmoothGrad(const Surface * surf, double x, double y, double radius, double ds,  double& gx, double& gy);


  void ResampleSurfaceToGrid2D(const Surface * surface,
                               Grid2D        * outgrid);

  int  GetWaveletFileFormat(const std::string & fileName, std::string & errText);

  void ReadAndWriteLocalGridsToFile(const std::string   & fileName,
                                    const std::string   & type,
                                    const float           scaleFactor,
                                    const ModelSettings * modelSettings,
                                    const Grid2D        * grid,
                                    const float           angle);

  void ResampleGrid2DToSurface(const Simbox   * simbox,
                               const Grid2D   * grid,
                               Surface       *& surface);

  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();

  int ComputeTime(int year, int month, int day) const;

  // CLASS VARIABLES ---------------------------------------------------

  // Bool variables indicating whether the corresponding data processing
  // succeeded
  bool outer_temp_simbox_;
  bool read_seismic_;
  bool read_wells_;
  bool block_wells_;
  bool setup_reflection_matrix_;
  bool temporary_wavelet_;
  bool optimize_well_location_;
  bool wavelet_handling_;
  bool wavelet_estimation_shape_;
  bool prior_corr_estimation_;
  bool setup_estimation_rock_physics_;
  bool setup_multigrid_;
  bool setup_trend_cubes_;

  MultiIntervalGrid       * multiple_interval_grid_;
  Simbox                    estimation_simbox_;
  NRLib::Volume             full_inversion_volume_;

  std::map<int, std::vector<SeismicStorage> >   seismic_data_;

  // Well logs
  std::vector<NRLib::Well>                      wells_;

  // Blocked well logs
  std::map<std::string, BlockedLogsCommon *>    mapped_blocked_logs_; //
  std::vector<std::string>                      continuous_logs_to_be_blocked_; // Continuous logs that should be blocked
  std::vector<std::string>                      discrete_logs_to_be_blocked_;   // Discrete logs that should be blocked

  // trend cubes
  int                                                           n_trend_cubes_;
  std::vector<CravaTrend>                                       trend_cubes_;
  std::map<std::string, std::vector<DistributionsRock *> >      rock_distributions_;     ///< Rocks used in rock physics model
  std::map<std::string, std::vector<DistributionWithTrend *> >  reservoir_variables_;    ///< Reservoir variables used in the rock physics model

  //Well variables not contained in NRlib::Well
  //std::map<std::string, int>                        timemissing_;
  //std::map<std::string, double>                     xpos0_;
  //std::map<std::string, double>                     ypos0_;
  //std::map<std::string, std::vector<std::string> >  faciesnames_;
  //std::map<std::string, std::vector<int> >          faciesNr_;
  //std::map<std::string, int>                        faciesok_; //Bool?
  //std::map<std::string, int>                        nFacies_;


  std::map<int, float **>   reflection_matrix_;
  bool                      reflection_matrix_from_file_; //False: created from global vp/vs

  std::vector<Wavelet*>     temporary_wavelets_; //One wavelet per angle

  std::map<int, Wavelet**>   wavelets_;
  std::map<int, std::vector<Grid2D *> > local_noise_scale_;
  std::map<int, std::vector<Grid2D *> > local_shift_;
  std::map<int, std::vector<Grid2D *> > local_scale_;
  std::map<int, std::vector<float> > global_noise_estimate_;
  std::map<int, std::vector<float> > sn_ratio_;

  //Simbox                     * estimation_simbox_;
  //NRLib::Volume              * full_inversion_volume_;

};
#endif
