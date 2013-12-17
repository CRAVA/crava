/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef COMMONDATA_H
#define COMMONDATA_H

#include "src/simbox.h"
//#include "src/welldata.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "lib/utils.h"
#include "src/blockedlogscommon.h"
#include "src/tasklist.h"
#include "src/seismicstorage.h"
#include "src/multiintervalgrid.h"
#include "src/background.h"
//#include "src/timeevolution.h"
#include "src/gridmapping.h"

class CravaTrend;
class InputFiles;
class ModelSettings;
//class ModelGeneral;
class TimeLine;
//class ModelGravityStatic;
class MultiIntervalGrid;
class GridMapping;
class BlockedLogsCommon;

class CommonData{
public:
  CommonData(ModelSettings                            * model_settings,
             InputFiles                               * input_files);

  ~ CommonData();

  enum                 gridTypes{CTMISSING, DATA, PARAMETER};

  //GET FUNCTIONS

  const Simbox                    & GetEstimationSimbox()     const { return estimation_simbox_      ;}
  const NRLib::Volume             & GetFullInversionVolume()  const { return full_inversion_volume_  ;}
  const std::vector<NRLib::Well>  & GetWells()                const { return wells_                  ;}
  const MultiIntervalGrid         * GetMultipleIntervalGrid() const { return multiple_interval_grid_ ;}
  MultiIntervalGrid               * GetMultipleIntervalGrid(void)   { return multiple_interval_grid_ ;}

  const std::map<std::string, std::vector<DistributionsRock *> >     & GetDistributionsRock()  const { return rock_distributions_  ;}
  const std::map<std::string, std::vector<DistributionWithTrend *> > & GetReservoirVariables() const { return reservoir_variables_ ;}

  TimeLine                               * GetTimeLine()                        { return time_line_       ;}
  const std::vector<std::string>         & GetFaciesNames()               const { return facies_names_    ;}
  const std::vector<int>                 & GetFaciesLabels()              const { return facies_labels_   ;}
  const std::vector<std::vector<float> > & GetPriorFacies()               const { return prior_facies_    ;}
  const std::vector<float>               & GetPriorFaciesInterval(int i)  const { return prior_facies_[i] ;}

  const std::map<std::string, BlockedLogsCommon *> GetBlockedLogs()       const { return mapped_blocked_logs_ ;}

  std::vector<Surface *>                 & GetFaciesEstimInterval()             { return facies_estim_interval_ ;}

  std::vector<float>       & GetGravityObservationUtmxTimeLapse(int time_lapse)  { return observation_location_utmx_[time_lapse]  ;}
  std::vector<float>       & GetGravityObservationUtmyTimeLapse(int time_lapse)  { return observation_location_utmy_[time_lapse]  ;}
  std::vector<float>       & GetGravityObservationDepthTimeLapse(int time_lapse) { return observation_location_depth_[time_lapse] ;}
  std::vector<float>       & GetGravityResponseTimeLapse(int time_lapse)         { return gravity_response_[time_lapse]           ;}
  std::vector<float>       & GetGravityStdDevTimeLapse(int time_lapse)           { return gravity_std_dev_[time_lapse]            ;}

  GridMapping              * GetTimeDepthMapping()                             { return time_depth_mapping_            ;}
  bool                       GetVelocityFromInversion()                        { return velocity_from_inversion_       ;}

  bool                                    GetUseLocalNoise()                            { return use_local_noise_                            ;}
  std::map<int, std::vector<Grid2D *> > & GetLocalNoiseScale()                          { return local_noise_scale_                          ;}
  std::vector<Grid2D *>                 & GetLocalNoiseScaleTimeLapse(int time_lapse)   { return local_noise_scale_.find(time_lapse)->second ;}
  std::vector<SeismicStorage>           & GetSeismicDataTimeLapse(int time_lapse)       { return seismic_data_.find(time_lapse)->second      ;}
  std::map<int, std::vector<float> >    & GetSNRatio()                                  { return sn_ratio_                                   ;}
  std::vector<float>                    & GetSNRatioTimeLapse(int time_lapse)           { return sn_ratio_.find(time_lapse)->second          ;}

  bool                                    GetRefMatFromFileGlobalVpVs()                 { return refmat_from_file_global_vpvs_               ;}
  float **                                GetReflectionMatrixTimeLapse(int time_lapse)  { return reflection_matrix_.find(time_lapse)->second ;}

  std::vector<std::vector<float> >      & GetAngularCorrelation(int time_lapse)         { return angular_correlations_[time_lapse]           ;}
  //Vario                                 * GetAngularCorrelation(int time_lapse)         { return angular_correlations_[time_lapse]           ;}

  std::vector<Wavelet *>                  & GetWavelet(int time_lapse)                  { return wavelets_.find(time_lapse)->second          ;}
  const std::vector<std::vector<double> > & GetSyntSeis(int time_lapse)                 { return synt_seis_.find(time_lapse)->second         ;}

  const std::vector<std::vector<double> > & GetTGradX()                           const { return t_grad_x_                                   ;}
  const std::vector<std::vector<double> > & GetTGradY()                           const { return t_grad_y_                                   ;}
  const NRLib::Grid2D<float>              & GetRefTimeGradX()                     const { return ref_time_grad_x_                            ;}
  const NRLib::Grid2D<float>              & GetRefTimeGradY()                     const { return ref_time_grad_y_                            ;}

  const Surface                           * GetPriorCorrXY(int i_interval)              { return prior_corr_XY_[i_interval]                  ;}



  void  SetupDefaultReflectionMatrix(float              **& reflection_matrix,
                                     double                 vsvp,
                                     const ModelSettings  * model_settings,
                                     int                    number_of_angles,
                                     int                    this_time_lapse);

  void FillInData(NRLib::Grid<double> & grid,
                  FFTGrid             * fft_grid_new,
                  const Simbox        & simbox,
                  StormContGrid       * storm_grid,
                  const SegY          * segy,
                  FFTGrid             * fft_grid_old,
                  float                 smooth_length,
                  int                 & missing_traces_simbox,
                  int                 & missing_traces_padding,
                  int                 & dead_traces_simbox,
                  int                   grid_type,
                  bool                  scale    = false,
                  bool                  is_segy  = true,
                  bool                  is_storm = false);

private:

  void LoadWellMoveInterval(const InputFiles             * input_files,
                            const Simbox                 * estimation_simbox,
                            std::vector<Surface *>       & well_move_interval,
                            std::string                  & err_text);
                            //bool                         & failed);

  bool  OptimizeWellLocations(ModelSettings                                 * model_settings,
                              InputFiles                                    * input_files,
                              const Simbox                                  * estimation_simbox,
                              //const NRLib::Volume                           & volume,
                              std::vector<NRLib::Well>                      & wells,
                              std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                              std::map<int, std::vector<SeismicStorage> >   & seismic_data,
                              std::map<int, float **>                       & reflection_matrix,
                              std::string                                   & err_text);
                              //bool                                          & failed);

  void MoveWell(const NRLib::Well & well,
                const Simbox      * simbox,
                double              delta_X,
                double              delta_Y,
                double              k_move);

  void  CalculateDeviation(NRLib::Well            & new_well,
                           const ModelSettings    * const model_settings,
                           float                  & dev_angle,
                           Simbox                 * simbox);

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

  void  EstimateXYPaddingSizes(Simbox         * simbox,
                               ModelSettings  * model_settings) const;

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
                                    std::string                    & err_text);
                                    //bool                           & failed);

  void SetSurfacesSingleInterval(const ModelSettings              * const model_settings,
                                 NRLib::Volume                    & full_inversion_volume,
                                 Simbox                           & estimation_simbox,
                                 const std::vector<std::string>   & surf_file,
                                 std::string                      & err_text);
                                 //bool                             & failed);

  bool ReadSeismicData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles,
                       std::string    & err_text);

  bool ReadWellData(ModelSettings                   * model_settings,
                    Simbox                          * estimation_simbox,
                    InputFiles                      * input_files,
                    std::vector<std::string>          & log_names,
                    const std::vector<std::string>  & log_names_from_user,
                    const std::vector<bool>         & inverse_velocity,
                    bool                              facies_log_given,
                    std::string                     & err_text);

  bool BlockWellsForEstimation(const ModelSettings                            * const model_settings,
                               const Simbox                                   & estimation_simbox,
                               const MultiIntervalGrid                        & multiple_interval_grid,
                               std::vector<NRLib::Well>                       & wells,
                               std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_common,
                               std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_for_correlation,
                               std::string                                    & err_text);

  //bool       CheckThatDataCoverGrid(const SegY            * segy,
  //                                  float                   offset,
  //                                  const Simbox          * timeCutSimbox,
  //                                  float                   guard_zone) const;

  void CutWell(std::string           well_file_name,
               NRLib::Well         & well,
               const NRLib::Volume & full_inversion_volume);

  void ProcessLogsNorsarWell(NRLib::Well                      & new_well,
                             std::vector<std::string>         & log_names_from_user,
                             const std::vector<bool>          & inverse_velocity,
                             bool                               facies_log_given,
                             std::string                      & error_text);
                             //bool                             & failed);

  void ProcessLogsRMSWell(NRLib::Well                     & new_well,
                          std::vector<std::string>        & log_names_from_user,
                          const std::vector<bool>         & inverse_velocity,
                          bool                              facies_log_given,
                          std::string                     & error_text);
                          //bool                            & failed);

  bool  SetupReflectionMatrix(ModelSettings  * model_settings,
                              InputFiles     * input_files,
                              std::string    & err_text);

  bool  SetupTemporaryWavelet(ModelSettings * model_settings,
                              //InputFiles    * input_files,
                              std::string   & err_text);

  bool  WaveletHandling(ModelSettings * model_settings,
                        InputFiles    * input_files,
                        std::string   & err_text);

  //bool       CheckThatDataCoverGrid(const SegY            * segy,
  //                                  float                   offset,
  //                                  //const Simbox          * timeCutSimbox,
  //                                  float                   guard_zone);

  bool       CheckThatDataCoverGrid(const SegY   * segy,
                                    float          offset,
                                    const Simbox & time_cut_simbox,
                                    float          guard_zone,
                                    std::string &  err_text) const;

  bool       CheckThatDataCoverGrid(StormContGrid * stormgrid,
                                    const Simbox  & simbox,
                                    float           guard_zone,
                                    std::string   & err_text) const;

  void ProcessLogsNorsarWell(NRLib::Well                  & new_well,
                             std::string                  & error_text,
                             bool                         & failed);

  void ProcessLogsRMSWell(NRLib::Well                     & new_well,
                          std::string                     & error_text,
                          bool                            & failed);

  //void ReadFaciesNamesFromWellFile(ModelSettings            * model_settings,
  //                                 std::string                well_file_name,
  //                                 std::vector<int>         & facies_nr,
  //                                 std::vector<std::string> & facies_names,
  //                                 std::string              & err_txt);

  void ReadFaciesNamesFromWellLogs(NRLib::Well              & well,
                                   std::vector<int>         & facies_nr,
                                   std::vector<std::string> & facies_names);
                                   //std::string              & err_txt);

  void SetFaciesNamesFromWells(ModelSettings            *& model_settings,
                               std::string               & err_text);
                               //bool                      & failed);

  void GetMinMaxFnr(int            & min,
                    int            & max,
                    const int        n_facies,
                    std::vector<int> facies_nr);


  float  ** ReadMatrix(const std::string                  & file_name,
                       int                                  n1,
                       int                                  n2,
                       const std::string                  & read_reason,
                       std::string                        & err_text);

  //void  SetupDefaultReflectionMatrix(float              **& reflection_matrix,
  //                                   double                 vsvp,
  //                                   const ModelSettings  * model_settings,
  //                                   int                    number_of_angles,
  //                                   int                    this_time_lapse);

  int Process1DWavelet(const ModelSettings                      * modelSettings,
                       const InputFiles                         * inputFiles,
                       const SeismicStorage                     * seismic_data,
                       std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                       const std::vector<Surface *>             & waveletEstimInterval,
                       const float                              * reflection_matrix,
                       std::vector<double>                      & synt_seic,
                       std::string                              & err_text,
                       Wavelet                                 *& wavelet,
                       Grid2D                                  *& local_noise_scale,
                       Grid2D                                  *& local_noise_shift,
                       Grid2D                                  *& local_noise_estimate,
                       unsigned int                               i_timelapse,
                       unsigned int                               j_angle,
                       const float                                angle,
                       float                                    & sn_ratio,
                       bool                                       estimate_wavlet,
                       bool                                       use_ricker_wavelet,
                       bool                                       use_local_noise);

  int Process3DWavelet(const ModelSettings                      * model_settings,
                       const InputFiles                         * input_files,
                       const SeismicStorage                     * seismic_data,
                       std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                       const std::vector<Surface *>             & wavelet_estim_interval,
                       const float                              * reflection_matrix,
                       std::string                              & err_text,
                       Wavelet                                 *& wavelet,
                       unsigned int                               i_timelapse,
                       unsigned int                               j_angle,
                       float                                      angle,
                       float                                      sn_ratio,
                       const NRLib::Grid2D<float>               & ref_time_grad_x,
                       const NRLib::Grid2D<float>               & ref_time_grad_y,
                       const std::vector<std::vector<double> >  & t_grad_x,
                       const std::vector<std::vector<double> >  & t_grad_y,
                       bool                                       estimate_wavelet);

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

  bool SetupTrendCubes(ModelSettings                     * model_settings,
                       InputFiles                        * input_files,
                       MultiIntervalGrid                 * multiple_interval_grid,
                       std::string                       & error_text);
                       //bool                              & failed);

  bool SetupRockPhysics(const ModelSettings                                 * model_settings,
                        const InputFiles                                    * input_files,
                        const MultiIntervalGrid                             * multiple_interval_grid,
                        const std::vector<CravaTrend>                       & trend_cubes,
                        const std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                        int                                                   n_trend_cubes,
                        std::string                                         & error_text);
                        //bool                                                & failed);

  void PrintExpectationAndCovariance(const std::vector<double>   & expectation,
                                     const NRLib::Grid2D<double> & covariance,
                                     const bool                  & has_trend) const;

  bool EstimateWaveletShape();

  bool SetupPriorFaciesProb(ModelSettings                               * model_settings,
                            InputFiles                                  * input_files,
                            std::string                                 & err_text);

  void FindFaciesEstimationInterval(InputFiles             * input_files,
                                    std::vector<Surface *> & facies_estim_interval,
                                    std::string            & err_text);

  void CheckFaciesNamesConsistency(ModelSettings     *& model_settings,
                                   const InputFiles   * input_files,
                                   std::string        & tmp_err_text) const;
                                   //int                  i_interval) const;

  void SetFaciesNamesFromRockPhysics();

  void ReadPriorFaciesProbCubes(const InputFiles                                * input_files,
                                ModelSettings                                   * model_settings,
                                std::vector<std::vector<NRLib::Grid<double> > > & prior_facies_prob_cubes,
                                const std::vector<Simbox>                       & interval_simboxes,
                                std::string                                     & err_text);

  //static FFTGrid  * CreateFFTGrid(int nx,
  //                                int ny,
  //                                int nz,
  //                                int nxp,
  //                                int nyp,
  //                                int nzp,
  //                                bool fileGrid);

  void ReadGridFromFile(const std::string                 & file_name,
                        const std::string                 & par_name,
                        const float                         offset,
                        std::vector<NRLib::Grid<double> > & grid,
                        const SegyGeometry               *& geometry,
                        const TraceHeaderFormat           * format,
                        int                                 grid_type,
                        const std::vector<Simbox>         & interval_simboxes,
                        const Simbox                      & estimation_simbox,
                        const ModelSettings               * model_settings,
                        std::string                       & err_text,
                        bool                                nopadding = true);

  void ReadGridFromFile(const std::string                 & file_name,
                        const std::string                 & par_name,
                        const float                         offset,
                        NRLib::Grid<double>               & interval_grids,
                        const SegyGeometry               *& geometry,
                        const TraceHeaderFormat           * format,
                        int                                 grid_type,
                        const Simbox                      & simboxes,
                        const ModelSettings               * model_settings,
                        std::string                       & err_text,
                        bool                                nopadding = true);

  void ReadCravaFile(NRLib::Grid<double> & grid,
                     const std::string   & file_name,
                     std::string         & err_text);

  void GetZPaddingFromCravaFile(const std::string & file_name,
                                std::string       & err_text,
                                int               & nz_pad);

  void ReadSegyFile(const std::string                 & file_name,
                    std::vector<NRLib::Grid<double> > & interval_grids,
                    const std::vector<Simbox>         & interval_simboxes,
                    const Simbox                      & simbox,
                    const ModelSettings               * model_settings,
                    const SegyGeometry               *& geometry,
                    int                                 grid_type,
                    const std::string                 & par_name,
                    float                               offset,
                    const TraceHeaderFormat           * format,
                    std::string                       & err_text,
                    bool                                nopadding = true);

  int GetFillNumber(int i, int n, int np);

  int FindClosestFactorableNumber(int leastint);

  void SmoothTraceInGuardZone(std::vector<float> & data_trace,
                              float                dz_data,
                              float                smooth_length);

  void ResampleTrace(const std::vector<float> & data_trace,
                     const rfftwnd_plan       & fftplan1,
                     const rfftwnd_plan       & fftplan2,
                     fftw_real                * rAmpData,
                     fftw_real                * rAmpFine,
                     int                        cnt,
                     int                        rnt,
                     int                        cmt,
                     int                        rmt);

  void InterpolateGridValues(std::vector<float> & grid_trace,
                             float                z0_grid,
                             float                dz_grid,
                             fftw_real          * rAmpFine,
                             float                z0_data,
                             float                dz_fine,
                             int                  n_fine,
                             int                  nz,
                             int                  nzp);
                             //int                  grid_nk);

  void InterpolateAndShiftTrend(std::vector<float>       & interpolated_trend,
                                float                      z0_grid,
                                float                      dz_grid,
                                const std::vector<float> & trend_long,
                                float                      z0_data,
                                float                      dz_fine,
                                int                        n_fine,
                                int                        nz,
                                int                        nzp);
                                //int                        grid_nk);

  int GetZSimboxIndex(int k,
                      int nz,
                      int nzp);
                      //int grid_nk);

  void SetTrace(const std::vector<float> & trace,
                NRLib::Grid<double>      & grid,
                size_t                     i,
                size_t                     j);

  void SetTrace(float                 value,
                NRLib::Grid<double> & grid,
                size_t                i,
                size_t                j);

  void SetTrace(const std::vector<float> & trace,
                FFTGrid                  * grid,
                size_t                     i,
                size_t                     j);

  void SetTrace(float     value,
                FFTGrid * grid,
                size_t    i,
                size_t    j);

  void ReadStormFile(const std::string                 & file_name,
                     std::vector<NRLib::Grid<double> > & interval_grids,
                     const int                           grid_type,
                     const std::string                 & par_name,
                     const std::vector<Simbox>         & interval_simboxes,
                     const ModelSettings               * model_settings,
                     std::string                       & err_text,
                     bool                                is_storm = true,
                     bool                                nopadding = true);

  bool SetupDepthConversion(ModelSettings  * model_settings,
                            InputFiles     * input_files,
                            std::string    & err_text);

  bool SetupBackgroundModel(ModelSettings  * model_settings,
                            InputFiles     * input_files,
                            std::string    & err_text);

  double FindMeanVsVp(NRLib::Grid<double> & vp,
                      NRLib::Grid<double> & vs);

  //void GetAvgMinMaxGrid(const NRLib::Grid<double> & grid,
  //                      double                    & avg,
  //                      double                    & min,
  //                      double                    & max);

  void SetUndefinedCellsToGlobalAverageGrid(NRLib::Grid<double> & grid,
                                            const double          avg);

  //void LogTransformGrid(NRLib::Grid<double> & grid);

  void SubtractGrid(NRLib::Grid<double>       & to_grid,
                    const NRLib::Grid<double> & from_grid);

  void ChangeSignGrid(NRLib::Grid<double> & grid);

  //void LoadVelocity(FFTGrid             *& velocity,
  //                  const Simbox         * interval_simbox, //timeSimbox,
  //                  const Simbox         * simbox, //timeCutSimbox,
  void LoadVelocity(NRLib::Grid<double>  & velocity,
                    const Simbox         & interval_simbox,
                    const ModelSettings  * model_settings,
                    const std::string    & velocity_field,
                    bool                 & velocity_from_inversion,
                    std::string          & err_text);

  //std::vector<std::map<std::string, DistributionsRock *> > GetRockDistributionTime0() const;
  std::map<std::string, DistributionsRock *> GetRockDistributionTime0() const;

  void GenerateRockPhysics3DBackground(const std::vector<DistributionsRock *> & rock_distribution,
                                       const std::vector<float>               & probability,
                                       //std::vector<NRLib::Grid<double> >      & parameters,
                                       NRLib::Grid<double>                    & vp,
                                       NRLib::Grid<double>                    & vs,
                                       NRLib::Grid<double>                    & rho,
                                       Simbox                                 & simbox,
                                       const CravaTrend                       & trend_cube);
                                       //int                                      i_interval);

  void SetupExtendedBackgroundSimbox(Simbox   * simbox,
                                     Surface  * corr_surf,
                                     Simbox  *& bg_simbox,
                                     int        output_format,
                                     int        output_domain,
                                     int        other_output);

  void SetupExtendedBackgroundSimbox(Simbox   * simbox,
                                     Surface  * top_corr_surf,
                                     Surface  * base_corr_surf,
                                     Simbox  *& bg_simbox,
                                     int        output_format,
                                     int        output_domain,
                                     int        other_output);


  bool SetupPriorCorrelation(ModelSettings                                                * model_settings,
                             const InputFiles                                             * input_files,
                             const std::vector<NRLib::Well>                               & wells,
                             std::map<std::string, BlockedLogsCommon *>                   & mapped_blocked_logs_for_correlation,
                             const std::vector<Simbox>                                    & interval_simboxes,
                             const std::vector<std::vector<float> >                       & prior_facies_prob,
                             const std::vector<std::string>                               & facies_names,
                             const std::vector<CravaTrend>                                & trend_cubes,
                             const std::map<int, std::vector<SeismicStorage> >            & seismic_data,
                             std::string                                                  & err_text);

  void            GetCorrGradIJ(float         & corr_grad_I, 
                                float         & corr_grad_J,
                                const Simbox  * simbox) const;

  void  CalculateCorrelationsFromRockPhysics(const std::vector<DistributionsRock *>           & rock_distribution,
                                             const std::map<std::string, float>               & probability,
                                             const std::vector<std::string>                   & facies_names,
                                             const CravaTrend                                 & trend_cubes,
                                             NRLib::Grid2D<double>                            & param_corr,
                                             std::string                                      & err_txt);

  void  CalculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
                                           const std::vector<float>               & probability,
                                           const std::vector<double>              & trend_position,
                                           NRLib::Grid2D<double>                  & sigma_sum) const;

  void ValidateCorrelationMatrix(float               ** C,
                                 const ModelSettings *  model_settings,
                                 std::string         &  err_txt);

  Surface * FindCorrXYGrid(const Simbox           * time_simbox,
                           const ModelSettings    * model_settings) const;

  bool  SetupTimeLine(ModelSettings * model_settings,
                      InputFiles    * input_files,
                      std::string   & err_text_common);

  bool  SetupGravityInversion(ModelSettings * model_settings,
                              InputFiles    * input_files,
                              std::string   & err_text_common);

  void  ReadGravityDataFile(const std::string   & file_name,
                            const std::string   & read_reason,
                            int                   n_obs,
                            int                   n_columns,
                            std::vector <float> & obs_loc_utmx,
                            std::vector <float> & obs_loc_utmy,
                            std::vector <float> & obs_loc_depth,
                            std::vector <float> & gravity_response,
                            std::vector <float> & gravity_std_dev,
                            std::string         & err_text);

  //void  SetUpscaledPaddingSize(ModelSettings * model_settings);

  //int   SetPaddingSize(int original_nxp, int upscaling_factor);

  //std::vector<int> FindClosestFactorableNumber(int leastint);

  bool  SetupTravelTimeInversion(ModelSettings * model_settings,
                                 InputFiles    * input_files,
                                 std::string   & err_text_common);

  void  ProcessHorizons(std::vector<Surface>   & horizons,
                        const InputFiles       * input_files,
                        std::string            & err_text,
                        bool                   & failed,
                        int                      i_timelapse);

  void ReadAngularCorrelations(ModelSettings * model_settings,
                               std::string   & err_text);

  bool optimizeWellLocations();
  bool estimateWaveletShape();

  int ComputeTime(int year, int month, int day) const;

  FFTGrid * createFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool file_grid);

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
  //bool wavelet_estimation_shape_;
  bool setup_multigrid_;
  bool setup_trend_cubes_;
  bool setup_estimation_rock_physics_;
  bool setup_prior_facies_probabilities_;
  bool setup_depth_conversion_;
  bool setup_background_model_;
  bool setup_prior_correlation_;
  bool setup_timeline_;
  //bool prior_corr_estimation_;
  bool setup_gravity_inversion_;
  bool setup_traveltime_inversion_;

  MultiIntervalGrid                           * multiple_interval_grid_;
  Simbox                                        estimation_simbox_;
  NRLib::Volume                                 full_inversion_volume_;

  std::map<int, std::vector<SeismicStorage> >   seismic_data_; //Map timelapse

  // Well logs
  std::vector<std::string>                      log_names_;
  std::vector<NRLib::Well>                      wells_;

  // Blocked well logs
  std::map<std::string, BlockedLogsCommon *>    mapped_blocked_logs_;           ///< Blocked logs with estimation simbox
  std::map<std::string, BlockedLogsCommon *>    mapped_blocked_logs_for_correlation_; ///< Blocked logs for estimation of vertical corr
  std::vector<std::string>                      continuous_logs_to_be_blocked_; ///< Continuous logs that should be blocked
  std::vector<std::string>                      discrete_logs_to_be_blocked_;   ///< Discrete logs that should be blocked

  // Trend cubes and rock physics
  int                                           n_trend_cubes_;
  //std::vector<CravaTrend>                       trend_cubes_;  //Trend cubes per interval.

  std::map<std::string, std::vector<DistributionsRock *> >     rock_distributions_;     ///< Rocks used in rock physics model
  std::map<std::string, std::vector<DistributionWithTrend *> > reservoir_variables_;    ///< Reservoir variables used in the rock physics model

  // prior facies
  std::vector<std::vector<float> >               prior_facies_;                  ///< Prior facies probabilities
  std::vector<Surface *>                         facies_estim_interval_;

  //std::vector<std::vector<FFTGrid *> >          prior_facies_prob_cubes_;       ///< Cubes for prior facies probabilities //H Need to move this to multi_interval_grid_
  //std::vector<std::vector<NRLib::Grid<double> > > prior_facies_prob_cubes_;       ///< Cubes for prior facies probabilities. Vector(facies) vector(intervals).

  // Timeline
  TimeLine                                    * time_line_;

  bool                                          forward_modeling_;

  std::map<int, float **>                       reflection_matrix_; //Map timelapse
  bool                                          refmat_from_file_global_vpvs_;  //True if reflection matrix is from file or set up from global vp/vs value.

  std::map<int, std::vector<Wavelet *> >             wavelets_; //Map time_lapse, vector angles
  std::map<int, std::vector<Grid2D *> >              local_noise_scale_;
  std::map<int, std::vector<Grid2D *> >              local_shift_;
  std::map<int, std::vector<Grid2D *> >              local_scale_;
  std::map<int, std::vector<float> >                 global_noise_estimate_;
  std::map<int, std::vector<float> >                 sn_ratio_;
  bool                                               use_local_noise_;
  std::map<int, std::vector<std::vector<double> > >  synt_seis_; //Map time_lapse, vector angles

  std::vector<std::vector<double> >             t_grad_x_;
  std::vector<std::vector<double> >             t_grad_y_;
  NRLib::Grid2D<float>                          ref_time_grad_x_; ///< Time gradient in x-direction for reference time surface (t0)
  NRLib::Grid2D<float>                          ref_time_grad_y_; ///< Time gradient in x-direction for reference time surface (t0)

  std::vector<std::vector<int> >                facies_nr_wells_;               ///< Facies Numbers per well.
  std::vector<std::vector<std::string> >        facies_names_wells_;            ///< Facies Names per well
  std::vector<bool>                             facies_log_wells_;              ///< True if this well has a facies log

  std::vector<std::string>                      facies_names_;                  ///< Facies names combined for wells. (Intervals?)
  std::vector<int>                              facies_labels_;

  // Prior correlation
  std::vector<Surface *>                        prior_corr_XY_;
  //std::vector<NRLib::Matrix>                      prior_var_0_;
  //std::vector<std::vector<NRLib::Grid<double> > > prior_cov_; //Vp, vs, rho
  //std::vector<std::vector<NRLib::Grid<double> > > prior_corr_; //Vp-vs, Vp-Rho, Vs-Rho

  std::vector<Wavelet*>                         temporary_wavelets_;            ///< Wavelet per angle

  //Gravimetry parameters per timelapse
  std::vector<std::vector<float> >              observation_location_utmx_;     ///< Vectors to store observation location coordinates
  std::vector<std::vector<float> >              observation_location_utmy_;
  std::vector<std::vector<float> >              observation_location_depth_;
  std::vector<std::vector<float> >              gravity_response_;              ///< Vector to store base line gravity response
  std::vector<std::vector<float> >              gravity_std_dev_;               ///< Vector to store base line gravity standard deviation

  //Traveltime parameters per timelapse
  std::vector<std::vector<Surface> >            horizons_;                      ///< Horizons used for horizon inversion
  std::vector<NRLib::Grid<double> >             rms_data_;                      ///< RMS data U^2

  //Depth conversion
  GridMapping                                 * time_depth_mapping_;
  bool                                          velocity_from_inversion_;

  //Angular correlations
  std::vector<std::vector<std::vector<float> > > angular_correlations_;
  //std::vector<Vario *>                          angular_correlations_;


};
#endif
