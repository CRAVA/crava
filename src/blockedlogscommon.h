/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKEDLOGSCOMMON_H
#define BLOCKEDLOGSCOMMON_H

#include "src/seismicstorage.h"
#include "fftw.h"
#include <map>
#include <string>
#include "nrlib/well/well.hpp"
#include "src/definitions.h"
#include "src/fftgrid.h"

class BlockedLogsCommon{
public:
  BlockedLogsCommon();

  BlockedLogsCommon(NRLib::Well     * well_data,
                    const Simbox    * const estimation_simbox,
                    bool              interpolate,
                    //bool            & failed,
                    std::string     & err_text,
                    float             max_hz_background = 0.0,
                    float             max_hz_seismic = 0.0);

  BlockedLogsCommon(const NRLib::Well                * well_data,
                    const std::vector<std::string>   & cont_logs_to_be_blocked,
                    const std::vector<std::string>   & disc_logs_to_be_blocked,
                    const Simbox                     * const estimation_simbox,
                    bool                               interpolate,
                    //bool                             & failed,
                    std::string                      & err_text);

  BlockedLogsCommon(const NRLib::Well   * well_data,
                    const StormContGrid & stormgrid,
                    float                 max_hz_background = 0.0,
                    float                 max_hz_seismic = 0.0);

  ~BlockedLogsCommon();

  void  FindOptimalWellLocation(std::vector<SeismicStorage>   & seismic_data,
                              const Simbox                  * time_simbox,
                              float                        ** refl_coef,
                              int                             n_angles,
                              const std::vector<float>      & angle_weight,
                              float                           max_shift,
                              int                             i_max_offset,
                              int                             j_max_offset,
                              const std::vector<Surface *>    limits,
                              int                           & i_move,
                              int                           & j_move,
                              float                         & k_move);


  //GET FUNCTIONS --------------------------------


  std::string                            GetWellName()         const   { return well_name_                                                            ;}
  int                                    GetNumberOfBlocks()   const   { return n_blocks_                                                             ;}
  int                                    GetNFacies()          const   { return static_cast<int>(facies_map_.size())                                  ;}
  const std::vector<double>            & GetXpos(void)         const   { return x_pos_blocked_                                                        ;}
  const std::vector<double>            & GetYpos(void)         const   { return y_pos_blocked_                                                        ;}
  const std::vector<double>            & GetZpos(void)         const   { return z_pos_blocked_                                                        ;}

  const std::vector<double>            & GetTVD(void)          const   { return z_pos_blocked_                                                        ;}
  const std::vector<double>            & GetTWT(void)          const   { return twt_blocked_                                                          ;}

  const std::map<int, std::string>     & GetFaciesMap()       const    { return facies_map_                                                           ;}
  const std::vector<int>               & GetIposVector()       const   { return i_pos_                                                                ;}
  const std::vector<int>               & GetJposVector()       const   { return j_pos_                                                                ;}
  const std::vector<int>               & GetKposVector()       const   { return k_pos_                                                                ;}

  const std::vector<double>            & GetVpBlocked(void)    const   { return continuous_logs_blocked_.find("Vp")->second                           ;}
  const std::vector<double>            & GetVsBlocked(void)    const   { return continuous_logs_blocked_.find("Vs")->second                           ;}
  const std::vector<int>               & GetFaciesBlocked(void) const  { return facies_blocked_                                                       ;}
  const std::vector<double>            & GetRhoBlocked(void)   const   { return continuous_logs_blocked_.find("Rho")->second                          ;}
  const std::vector<double>            & GetVpUnblocked(void)  const   { return continuous_logs_unblocked_.find("Vp")->second                         ;}
  const std::vector<double>            & GetVsUnblocked(void)  const   { return continuous_logs_unblocked_.find("Vs")->second                         ;}
  const std::vector<double>            & GetRhoUnblocked(void) const   { return continuous_logs_unblocked_.find("Rho")->second                        ;}

  const std::vector<double>            & GetVpSeismicResolution(void)  const { return cont_logs_seismic_resolution_.find("Vp")->second                ;}
  const std::vector<double>            & GetVsSeismicResolution(void)  const { return cont_logs_seismic_resolution_.find("Vs")->second                ;}
  const std::vector<double>            & GetRhoSeismicResolution(void) const { return cont_logs_seismic_resolution_.find("Rho")->second               ;}

  const std::vector<double>            & GetVpBackgroundResolution(void)  const { return cont_logs_background_resolution_.find("Vp")->second          ;}
  const std::vector<double>            & GetVsBackgroundResolution(void)  const { return cont_logs_background_resolution_.find("Vs")->second          ;}
  const std::vector<double>            & GetRhoBackgroundResolution(void) const { return cont_logs_background_resolution_.find("Rho")->second         ;}

  const std::vector<double>            & GetVpHighCutBackground(void)  const { return cont_logs_highcut_background_.find("Vp")->second                ;}
  const std::vector<double>            & GetVsHighCutBackground(void)  const { return cont_logs_highcut_background_.find("Vs")->second                ;}
  const std::vector<double>            & GetRhoHighCutBackground(void) const { return cont_logs_highcut_background_.find("Rho")->second               ;}

  bool                                   HasContLog(std::string s)     { return (continuous_logs_blocked_.find(s) != continuous_logs_blocked_.end())  ;}
  bool                                   HasDiscLog(std::string s)     { return (discrete_logs_blocked_.find(s) != discrete_logs_blocked_.end())      ;}

  void                                   GetVerticalTrend(const std::vector<double>  & blocked_log,
                                                          std::vector<double>        & trend);

  void                                   GetVerticalTrend(const int        * blocked_log,
                                                          std::vector<int> & trend);

  void                                   GetVerticalTrend(const std::vector<double> & blocked_log,
                                                          float                     * trend) const;

  void                                   GetVerticalTrendLimited(const std::vector<double>       & blocked_log,
                                                                 std::vector<double>             & trend,
                                                                 const std::vector<Surface *>    & limits);

  void                                   GetBlockedGrid(const SeismicStorage   * grid,
                                                        const Simbox           * estimation_simbox,
                                                        std::vector<double>    & blocked_log,
                                                        int                      i_offset = 0,
                                                        int                      j_offset = 0);

  void                                   GetBlockedGrid(const FFTGrid       * grid,
                                                        std::vector<double> & blocked_log,
                                                        int                   i_offset = 0,
                                                        int                   j_offset = 0);


  // FUNCTIONS -----------------------------------

  void                                   EstimateCor(fftw_complex * var1_c,
                                                     fftw_complex * var2_c,
                                                     fftw_complex * ccor_1_2_c,
                                                     int            cnzp) const;

  void                                   FillInCpp(const float * coeff,
                                                   int           start,
                                                   int           length,
                                                   fftw_real   * cpp_r,
                                                   int           nzp);

  void                                   SetLogFromVerticalTrend(std::vector<double>   & vertical_trend,
                                                                 double                  z0,              // z-value of center in top layer
                                                                 double                  dz,              // dz in vertical trend
                                                                 int                     nz,              // layers in vertical trend
                                                                 std::string             type,
                                                                 int                     i_angle);

  void                                   FillInSeismic(std::vector<double> & seismic_data,
                                                       int                   start,
                                                       int                   length,
                                                       fftw_real           * seis_r,
                                                       int                   nzp) const;

  void                                   SetSeismicGradient(double                            v0,
                                                            const NRLib::Grid2D<float>   &    structure_depth_grad_x,
                                                            const NRLib::Grid2D<float>   &    structure_depth_grad_y,
                                                            const NRLib::Grid2D<float>   &    ref_time_grad_x,
                                                            const NRLib::Grid2D<float>   &    ref_time_grad_y,
                                                            std::vector<double>          &    x_gradient,
                                                            std::vector<double>          &    y_gradient);

  void                                   SetTimeGradientSettings(float distance, float sigma_m);

  void                                   FindSeismicGradient(const std::vector<SeismicStorage> & seismic_data,
                                                             const Simbox                      * const estimation_simbox,
                                                             int                                 n_angles,
                                                             std::vector<double>               & x_gradient,
                                                             std::vector<double>               & y_gradient,
                                                             std::vector<std::vector<double> > & sigma_gradient);

  void                                   FindContinuousPartOfData(const std::vector<bool> & hasData,
                                                                  int                       nz,
                                                                  int                     & start,
                                                                  int                     & length) const;

  int                                    FindMostProbable(const int * count,
                                                          int         n_facies,
                                                          int         block_index);

  void                                   FilterLogs(float max_hz_background,
                                                    float max_hz_seismic);

  void                                   CreateHighCutBackground(std::vector<int> b_ind);

private:

  // FUNCTIONS------------------------------------

  void                  SetLogFromVerticalTrend(std::vector<double>   & blocked_log,
                                                std::vector<double>   & z_pos,
                                                int                     n_blocks,
                                                std::vector<double>   & vertical_trend,
                                                double                  z0,
                                                double                  dzVt,
                                                int                     nz);

  void         InterpolateTrend(const double   * blocked_log,
                                double         * trend);

  void         InterpolateTrend(const std::vector<double> & blocked_log,
                                double                    * trend);

  void         InterpolateTrend(const std::vector<double>    & blocked_log,
                                std::vector<double>          & trend);

  void         InterpolateTrend(const std::vector<double>      & blocked_log,
                                 std::vector<double>           & trend,
                                 const std::vector<Surface *>  & limits);

  void         InterpolateTrend(const int        * blocked_log,
                                std::vector<int> & trend);

  double ComputeElasticImpedance(double         vp,
                                 float         vs,
                                 float         rho,
                                 const float * coeff) const;

  void    RemoveMissingLogValues(const NRLib::Well                            * well_data,
                                 std::vector<double>                          & x_pos_unblocked,
                                 std::vector<double>                          & y_pos_unblocked,
                                 std::vector<double>                          & z_pos_unblocked,
                                 std::vector<double>                          & twt_unblocked,
                                 std::vector<int>                             & facies_unblocked,
                                 std::map<std::string, std::vector<double> >  & continuous_logs_unblocked,
                                 std::map<std::string, std::vector<int> >     & discrete_logs_unblocked,
                                 const std::vector<std::string>               & cont_logs_to_be_blocked,
                                 const std::vector<std::string>               & disc_logs_to_be_blocked,
                                 unsigned int                                 & n_data,
                                 bool                                         & failed,
                                 std::string                                  & err_text);

  void    BlockWell(const Simbox                                        * const estimation_simbox,
                    const std::map<std::string, std::vector<double> >   & continuous_logs_unblocked,
                    const std::map<std::string, std::vector<int> >      & discrete_logs_unblocked,
                    std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                    std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                    unsigned int                                          n_data,
                    bool                                                  facies_log_defined,
                    const std::map<int, std::string>                    & facies_map,
                    bool                                                  interpolate,
                    bool                                                & failed,
                    std::string                                         & err_text);

  void    FindSizeAndBlockPointers(const Simbox         * const estimation_simbox,
                                   std::vector<int>     & b_ind);

  void    FindSizeAndBlockPointers(const StormContGrid  & stormgrid,
                                   std::vector<int>     & bInd);

  void    FindBlockIJK(const Simbox                     * const estimation_simbox,
                       const std::vector<int>           & b_ind);

  void    FindBlockIJK(const StormContGrid    & stormgrid,
                       const std::vector<int>   b_ind);

  void    BlockCoordinateLog(const std::vector<int>     &  b_ind,
                             const std::vector<double>  &  coord,
                             std::vector<double>        &  blocked_coord);

  void    BlockContinuousLog(const std::vector<int>     &  b_ind,
                             const std::vector<double>  &  well_log,
                             std::vector<double>        &  blocked_log);

  void    BlockFaciesLog(const std::vector<int>     & b_ind,
                           const std::vector<int>     & well_log,
                           const std::map<int,std::string>     & facies_map,
                           int                          n_facies,
                           std::vector<int>           &  blocked_log);

  void    InterpolateContinuousLog(std::vector<double>   & blocked_log,
                                   int                     start,
                                   int                     end,
                                   int                     index,
                                   float                   rel);

  int   FindMostProbable(const std::vector<int>  & count,
                         int                       n_facies,
                         int                       block_index);

  void    SmoothTrace(std::vector<float> &trace);

  void    FindPeakTrace(std::vector<float> &trace, std::vector<double> &z_peak, std::vector<double> &peak,
                        std::vector<double> &b, double dz, double z_top);

  void    PeakMatch(std::vector<double> &zPeak, std::vector<double> &peak, std::vector<double> &b,
                    std::vector<double> &zPeakW, std::vector<double> &peakW, std::vector<double> &bW);

  double  ComputeShift(std::vector<double> &z_peak, std::vector<double> &z_peak_w, double z0);

  void    ComputeGradient(std::vector<double> &q_epsilon, std::vector<double> &q_epsilon_data,
                          std::vector<double> &z_shift, int nx, int ny, double dx, double dy);

  void    SmoothGradient(std::vector<double>               & x_gradient,
                         std::vector<double>               & y_gradient,
                         std::vector<double>               & q_epsilon,
                         std::vector<double>               & q_epsilon_data,
                         std::vector<std::vector<double> > & sigma_gradient);

  void    ComputePrecisionMatrix(double &a, double &b, double &c);


  float   ComputeElasticImpedance(double         alpha,
                                  double         beta,
                                  double         rho,
                                  const float * coeff) const;

  void    SetLogFromVerticalTrend(float     *& blocked_log,
                                  const std::vector<double> & zpos,
                                  int          nBlocks,
                                  float      * vertical_trend,
                                  double       z0,
                                  double       dzVt,
                                  int          nz);

  bool    ResampleTime(std::vector<double> & time_resampled,
                       int                   nd,
                       double              & dt);

  void    ResampleLog(std::vector<double>       & log_resampled,
                      const std::vector<double> & log,
                      const std::vector<double> & time,
                      const std::vector<double> & time_resampled,
                      int                         nd,
                      double                      dt);

  void    InterpolateLog(std::vector<double>       & log_interpolated,
                         const std::vector<double> & log_resampled,
                         int                         nd);

  void    ApplyFilter(std::vector<double> & log_filtered,
                      std::vector<double> & log_interpolated,
                      int                   n_time_samples,
                      double                dt_milliseconds,
                      float                 max_hz);

  // CLASS VARIABLES -----------------------------

  unsigned int       n_blocks_;         // number of blocks
  unsigned int       n_data_; //Number of non-missing
  std::string        well_name_;

  std::map<int, std::string>  facies_map_;
  bool               facies_log_defined_;
  //std::vector<int>   facies_numbers_;
  //std::vector<int>   facies_;
  // Blocked values

  std::vector<double> x_pos_blocked_;  // Blocked x position
  std::vector<double> y_pos_blocked_;  // Blocked y position
  std::vector<double> z_pos_blocked_;  // Blocked z position
  std::vector<double> twt_blocked_;    // Blocked twt value
  std::vector<int>    facies_blocked_; // Blocked facies log

  std::vector<int> i_pos_;    // Simbox i position for block
  std::vector<int> j_pos_;    // Simbox j position for block
  std::vector<int> k_pos_;    // Simbox k position for block

  int n_continuous_logs_;     // Number of continuous logs
  int n_discrete_logs_;       // Number of discrete logs

  std::map<std::string, std::vector<double> > continuous_logs_blocked_;  // Map between variable name and blocked continuous log
  std::map<std::string, std::vector<int> > discrete_logs_blocked_;       // Map between variable name and blocked discrete log

  std::map<std::string, std::vector<double> > cont_logs_seismic_resolution_;  // Map between variable name and blocked continuous log
  std::map<std::string, std::vector<int> > disc_logs_seismic_resolution_;  // Map between variable name and blocked discrete log

  std::map<std::string, std::vector<double> > cont_logs_background_resolution_;  // Map between variable name and blocked continuous log
  std::map<std::string, std::vector<int> > disc_logs_background_resolution_;  // Map between variable name and blocked discrete log

  std::map<std::string, std::vector<double> > cont_logs_highcut_background_;  // Continuous logs high-cut filtered to background resolution (log-domain)
  std::map<std::string, std::vector<int> > disc_logs_highcut_background_; // Discrete logs high-cut filtered to background resolution (log-domain)

  std::vector<std::vector<double> > actual_synt_seismic_data_; ///< Forward modelled seismic data using local wavelet
  std::vector<std::vector<double> > well_synt_seismic_data_;   ///< Forward modelled seismic data using wavelet estimated in well

  float lateral_threshold_gradient_; //Minimum lateral distance where gradient lines must not cross
  float sigma_m_; //Smoothing factor for the gradients

  // Unblocked values

  std::vector<double> x_pos_unblocked_;
  std::vector<double> y_pos_unblocked_;
  std::vector<double> z_pos_unblocked_;
  std::vector<double> twt_unblocked_;
  std::vector<int>    facies_unblocked_;

  std::map<std::string, std::vector<double> > continuous_logs_unblocked_;  // Map between variable name and unblocked continuous log
  std::map<std::string, std::vector<int> > discrete_logs_unblocked_;       // Map between variable name and unblocked discrete log

  ///H Copied from blockedlogs.h, they are only set in SetLogFromVerticalTrend. Are they needed later? Incorporated into cont_logs_seismic_resolution_
  //float        * alpha_seismic_resolution_; ///<
  //float        * beta_seismic_resolution_;  ///< Logs filtered to resolution of inversion result
  //float        * rho_seismic_resolution_;   ///<
  //float       ** actual_synt_seismic_data_; ///< Forward modelled seismic data using local wavelet
  //float       ** well_synt_seismic_data_;   ///< Forward modelled seismic data using wavelet estimated in well
  int            n_angles_;                  ///< Number of AVA stacks

  bool                      interpolate_;              ///<

  //int                       n_angles_;                 ///< Number of angles
  int                       n_layers_;                 ///< Number of layers in estimation_simbox
  float                     dz_;                       ///< Simbox dz value for block

  //int                       n_facies_;

  int                       first_M_;                   ///< First well log entry contributing to blocked well
  int                       last_M_;                    ///< Last well log entry contributing to blocked well

  int                       first_B_;                   ///< First block with contribution from well log
  int                       last_B_;                    ///< Last block with contribution from well log



};
#endif
