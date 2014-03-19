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

class MultiIntervalGrid;

class BlockedLogsCommon{
public:

  // Constructor for wavelet estimation blocking
  //BlockedLogsCommon(const NRLib::Well * well_data,
  //                  const Simbox      * const estimation_simbox,
  //                  bool                interpolate,
  //                  std::string       & err_text,
  //                  float               max_hz_background = 0.0,
  //                  float               max_hz_seismic = 0.0);


  //Copy constructor
  BlockedLogsCommon(const BlockedLogsCommon & logs);

  // Constructor for correlation estimation blocking
  BlockedLogsCommon(NRLib::Well                      * well_data,
                    const std::vector<std::string>   & cont_logs_to_be_blocked,
                    const std::vector<std::string>   & disc_logs_to_be_blocked,
                    const MultiIntervalGrid          * multiple_interval_grid,
                    bool                               interpolate,
                    std::string                      & err_text);

  // Constructor for blocking in the surrounding estimation simbox
  BlockedLogsCommon(NRLib::Well                      * well_data,
                    const std::vector<std::string>   & cont_logs_to_be_blocked,
                    const std::vector<std::string>   & disc_logs_to_be_blocked,
                    const Simbox                     * const estimation_simbox,
                    bool                               interpolate,
                    std::string                      & err_text);

  // Constructor for wavelet estimation blocking
  BlockedLogsCommon(const NRLib::Well   * well_data,
                    const StormContGrid & stormgrid);

  ~BlockedLogsCommon();


  //GET FUNCTIONS --------------------------------

  int                                    GetNBlocksWithData(std::string s) const { return n_blocks_with_data_.find(s)->second                         ;}
  int                                    GetNBlocksWithDataTot() const { return n_blocks_with_data_tot_                                               ;}
  std::string                            GetWellName()         const   { return well_name_                                                            ;}
  int                                    GetNumberOfBlocks()   const   { return n_blocks_                                                             ;}
  int                                    GetNFacies()          const   { return static_cast<int>(facies_map_.size())                                  ;}
  const std::vector<double>            & GetXposBlocked(void)  const   { return x_pos_blocked_                                                        ;}
  const std::vector<double>            & GetYposBlocked(void)  const   { return y_pos_blocked_                                                        ;}
  const std::vector<double>            & GetZposBlocked(void)  const   { return z_pos_blocked_                                                        ;}
  const std::vector<double>            & GetXposRawLogs()      const   { return x_pos_raw_logs_                                                       ;}
  const std::vector<double>            & GetYposRawLogs()      const   { return y_pos_raw_logs_                                                       ;}
  const std::vector<double>            & GetZposRawLogs()      const   { return z_pos_raw_logs_                                                       ;}

  const std::map<int, std::string>     & GetFaciesMap()        const   { return facies_map_                                                           ;}
  const std::vector<int>               & GetIposVector()       const   { return i_pos_                                                                ;}
  const std::vector<int>               & GetJposVector()       const   { return j_pos_                                                                ;}
  const std::vector<int>               & GetKposVector()       const   { return k_pos_                                                                ;}

  const std::vector<double>            & GetVpBlocked(void)    const   { return continuous_logs_blocked_.find("Vp")->second                           ;}
        std::vector<double>            & GetVpBlocked(void)            { return continuous_logs_blocked_.find("Vp")->second                           ;}
  const std::vector<double>            & GetVsBlocked(void)    const   { return continuous_logs_blocked_.find("Vs")->second                           ;}
        std::vector<double>            & GetVsBlocked(void)            { return continuous_logs_blocked_.find("Vs")->second                           ;}
  const std::vector<double>            & GetRhoBlocked(void)   const   { return continuous_logs_blocked_.find("Rho")->second                          ;}
        std::vector<double>            & GetRhoBlocked(void)           { return continuous_logs_blocked_.find("Rho")->second                          ;}
  const std::vector<double>            & GetMDBlocked(void)    const   { return continuous_logs_blocked_.find("MD")->second                           ;}
  const std::vector<int>               & GetFaciesBlocked(void) const  { return facies_blocked_                                                       ;}
  const std::vector<double>            & GetVpRawLogs(void)  const     { return continuous_raw_logs_.find("Vp")->second                               ;}
  const std::vector<double>            & GetVsRawLogs(void)  const     { return continuous_raw_logs_.find("Vs")->second                               ;}
  const std::vector<double>            & GetRhoRawLogs(void) const     { return continuous_raw_logs_.find("Rho")->second                              ;}

  const std::vector<double>            & GetVpSeismicResolution(void)  const { return cont_logs_seismic_resolution_.find("Vp")->second                ;}
  const std::vector<double>            & GetVsSeismicResolution(void)  const { return cont_logs_seismic_resolution_.find("Vs")->second                ;}
  const std::vector<double>            & GetRhoSeismicResolution(void) const { return cont_logs_seismic_resolution_.find("Rho")->second               ;}

  const std::vector<double>            & GetVpHighCutBackground(void)  const { return cont_logs_highcut_background_.find("Vp")->second                ;}
  const std::vector<double>            & GetVsHighCutBackground(void)  const { return cont_logs_highcut_background_.find("Vs")->second                ;}
  const std::vector<double>            & GetRhoHighCutBackground(void) const { return cont_logs_highcut_background_.find("Rho")->second               ;}

  const std::vector<double>            & GetVpHighCutSeismic(void)  const { return cont_logs_highcut_seismic_.find("Vp")->second                      ;}
  const std::vector<double>            & GetVsHighCutSeismic(void)  const { return cont_logs_highcut_seismic_.find("Vs")->second                      ;}
  const std::vector<double>            & GetRhoHighCutSeismic(void) const { return cont_logs_highcut_seismic_.find("Rho")->second                     ;}

  const std::vector<double>            & GetVpPredicted(void)  const { return continuous_logs_predicted_.find("Vp")->second                           ;}
  const std::vector<double>            & GetVsPredicted(void)  const { return continuous_logs_predicted_.find("Vs")->second                           ;}
  const std::vector<double>            & GetRhoPredicted(void) const { return continuous_logs_predicted_.find("Rho")->second                          ;}

  const std::vector<double>            & GetVpForFacies(void)  const { return vp_for_facies_                                                          ;}
  const std::vector<double>            & GetRhoForFacies(void) const { return rho_for_facies_                                                         ;}

  const std::vector<double>            & GetCpp(int angle)             const { return cpp_.find(angle)->second                                        ;}
  const std::vector<double>            & GetFaciesProb(int angle)      const { return facies_prob_.find(angle)->second                                ;}
  const std::vector<double>            & GetRealSeismicData(int angle) const { return real_seismic_data_.find(angle)->second                          ;}
  //const std::vector<double>            & GetActualSyntSeismicData(int angle) const { return actual_synt_seismic_data_.find(angle)->second ;}
  //const std::vector<double>            & GetWellSyntSeismicData(int angle) const { return well_synt_seismic_data_.find(angle)->second ;}

  bool                                   GetIsDeviated(void)                const { return is_deviated_                                               ;}

  bool                                   GetUseForFaciesProbabilities(void) const { return (use_for_facies_probabilities_ > 0)                        ;}
  bool                                   GetUseForWaveletEstimation(void)   const { return (use_for_wavelet_estimation_ > 0)                          ;}
  bool                                   GetUseForFiltering(void)           const { return (use_for_filtering_ > 0)                                   ;}
  bool                                   GetUseForBackgroundTrend(void)     const { return (use_for_background_trend_ > 0)                            ;}

  float                                  GetDz(void)                        const { return dz_                                                        ;}

  bool                                   HasSyntheticVsLog(void)            const { return (real_vs_log_ == 0)                                        ;}
  bool                                   HasContLog(std::string s)     { return (continuous_logs_blocked_.find(s) != continuous_logs_blocked_.end())  ;}
  bool                                   HasDiscLog(std::string s)     { return (discrete_logs_blocked_.find(s) != discrete_logs_blocked_.end())      ;}

  void                                   GetVerticalTrend(const std::vector<double>  & blocked_log,
                                                          std::vector<double>        & trend);

  void                                   GetVerticalTrend(const int        * blocked_log,
                                                          std::vector<int> & trend);

  //void                                   GetVerticalTrend(const std::vector<double> & blocked_log,
  //                                                        float                     * trend) const;

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

  void                                   GetBlockedGrid(const NRLib::Grid<double> * grid,
                                                        std::vector<double>       & blocked_log,
                                                        int                         i_offset = 0,
                                                        int                         j_offset = 0);


  //SET FUNCTIONS --------------------------------

  void                                   SetDeviated(bool deviated)                                     { is_deviated_ = deviated                                      ;}
  void                                   SetRealVsLog(int real_vs_log)                                  { real_vs_log_ = real_vs_log                                   ;}
  void                                   SetUseForFaciesProbabilities(int use_for_facies_probabilities) { use_for_facies_probabilities_ = use_for_facies_probabilities ;}
  void                                   SetUseForBackgroundTrend(int use_for_background_trend)         { use_for_background_trend_ = use_for_background_trend         ;}
  void                                   SetUseForFiltering(int use_for_filtering)                      { use_for_filtering_ = use_for_filtering                       ;}
  void                                   SetUseForWaveletEstimation(int use_for_wavelet_estimation)     { use_for_wavelet_estimation_ = use_for_wavelet_estimation     ;}

  void                                   SetNAngles(int n_angles)                                       { n_angles_  = n_angles                                        ;}


  // FUNCTIONS -----------------------------------

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

  void                                   EstimateCor(fftw_complex * var1_c,
                                                     fftw_complex * var2_c,
                                                     fftw_complex * ccor_1_2_c,
                                                     int            cnzp) const;

  void                                   FillInCpp(const float * coeff,
                                                   int           start,
                                                   int           length,
                                                   fftw_real   * cpp_r,
                                                   int           nzp);

  void                                   SetLogFromVerticalTrend(const std::vector<double> & vertical_trend,
                                                                 double                      z0,              // z-value of center in top layer
                                                                 double                      dz,              // dz in vertical trend
                                                                 int                         nz,              // layers in vertical trend
                                                                 std::string                 type,
                                                                 int                         i_angle,
                                                                 int                         n_angles);

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

  void                                   WriteWell(const int                        formats,
                                                   const float                      max_hz_background,
                                                   const float                      max_hz_seismic,
                                                   const std::vector<std::string> & facies_name,
                                                   const std::vector<int>         & facies_label);

  void                                   GenerateSyntheticSeismic(const float   * const * refl_coef,
                                                                  int                     n_angles,
                                                                  std::vector<Wavelet *> & wavelet,
                                                                  int                     nz,
                                                                  int                     nzp,
                                                                  const Simbox          * timeSimbox);

  void                                   SetLogFromGrid(FFTGrid    * grid,
                                                        int          i_angle,
                                                        int          n_angles,
                                                        std::string  type);

  void                                   SetSpatialFilteredLogs(std::vector<double>       & filtered_log,
                                                                int                         n_data,
                                                                std::string                 type,
                                                                const std::vector<double> & bg);

  void                                   FindMeanVsVp(const NRLib::Surface<double> & top,
                                                      const NRLib::Surface<double> & bot,
                                                      double                         mean_vs_vp,
                                                      int                            n_vs_vp);

  bool                                   VolumeFocus(const NRLib::Volume     & volume); //Sets all observations outside volume to missing.
                                                                                        //Returns false if none left - in that case, object is not modified.

private:

  // FUNCTIONS------------------------------------

  void                  SetLogFromVerticalTrend(std::vector<double>       & blocked_log,
                                                std::vector<double>       & z_pos,
                                                int                         n_blocks,
                                                const std::vector<double> & vertical_trend,
                                                double                      z0,
                                                double                      dzVt,
                                                int                         nz);

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
                                 std::vector<double>                          & x_pos_raw_logs,
                                 std::vector<double>                          & y_pos_raw_logs,
                                 std::vector<double>                          & z_pos_raw_logs,
                                 std::vector<int>                             & facies_raw_logs,
                                 std::map<std::string, std::vector<double> >  & continuous_logs_raw_logs,
                                 std::map<std::string, std::vector<int> >     & discrete_logs_raw_logs,
                                 const std::vector<std::string>               & cont_logs_to_be_blocked,
                                 const std::vector<std::string>               & disc_logs_to_be_blocked,
                                 unsigned int                                 & n_data,
                                 bool                                         & failed,
                                 std::string                                  & err_text);

  void    BlockWell(const Simbox                                        * const estimation_simbox,
                    const NRLib::Well                                   * well,
                    const std::map<std::string, std::vector<double> >   & continuous_logs_raw_logs,
                    const std::map<std::string, std::vector<int> >      & discrete_logs_raw_logs,
                    std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                    std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                    unsigned int                                          n_data,
                    bool                                                  facies_log_defined,
                    const std::map<int, std::string>                    & facies_map,
                    bool                                                  interpolate,
                    bool                                                & failed,
                    std::string                                         & err_text);

  void    BlockWellForCorrelationEstimation(const MultiIntervalGrid                             * multiple_interval_grid,
                                            const NRLib::Well                                   * well,
                                            const std::map<std::string, std::vector<double> >   & continuous_logs_raw_logs,
                                            const std::map<std::string, std::vector<int> >      & discrete_logs_raw_logs,
                                            std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                                            std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                                            unsigned int                                          n_data,
                                            bool                                                  facies_log_defined,
                                            const std::map<int, std::string>                    & facies_map,
                                            bool                                                  interpolate,
                                            int                                                 & n_layers,
                                            bool                                                & failed,
                                            std::string                                         & err_text);

  void    FindSizeAndBlockPointers(const MultiIntervalGrid       * multiple_interval_grid,
                                   std::vector<int>              & b_ind,
                                   int                             n_data,
                                   int                           & n_layers,
                                   unsigned int                  & n_blocks,
                                   std::map<std::string, int>    & n_layers_adjusted_per_interval);

  void    FindSizeAndBlockPointers(const Simbox         * const estimation_simbox,
                                   std::vector<int>     & b_ind,
                                   unsigned int         & n_blocks);

  void    FindSizeAndBlockPointers(const StormContGrid  & stormgrid,
                                   std::vector<int>     & bInd,
                                   unsigned int         & n_blocks);

  void    FindBlockIJK(const Simbox                     * const estimation_simbox,
                       const std::vector<int>           & b_ind);

  void    FindBlockIJK(const MultiIntervalGrid          * multiple_interval_grid,
                       const std::vector<double>        & x_pos,
                       const std::vector<double>        & y_pos,
                       const std::vector<double>        & z_pos,
                       const std::vector<int>           & b_ind);

  void    FindBlockIJK(const StormContGrid              & stormgrid,
                       const std::vector<int>             b_ind);

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

  void    CountBlocksWithData(const std::vector<double>                          & x_pos_blocked,
                              const std::vector<double>                          & y_pos_blocked,
                              const std::vector<double>                          & z_pos_blocked,
                              const std::map<std::string, std::vector<double> >  & continuous_logs_blocked,
                              unsigned int                                         n_blocks,
                              std::map<std::string, int>                         & n_blocks_with_data,
                              int                                                & n_blocks_with_data_tot);

  void    CountBlocksWithDataPerInterval(const MultiIntervalGrid                            * multiple_interval_grid,
                                         const std::vector<double>                          & x_pos_blocked,
                                         const std::vector<double>                          & y_pos_blocked,
                                         const std::vector<double>                          & z_pos_blocked,
                                         const std::map<std::string, std::vector<double> >  & continuous_logs_blocked,
                                         unsigned int                                         n_blocks,
                                         std::map<std::string, int>                         & n_blocks_with_data,
                                         int                                                & n_blocks_with_data_tot,
                                         int                                                  n_intervals) const;

  void    WriteRMSWell(const float                      max_hz_background,
                       const float                      max_hz_seismic,
                       const std::vector<std::string> & facies_name,
                       const std::vector<int>         & facies_label);

  void    WriteNorsarWell(const float max_hz_background,
                          const float max_hz_seismic);

  void    FindXYZForVirtualPart(const Simbox * simbox);

  void    UpdateLog(std::vector<double>                           & data,
                    const std::vector<std::pair<size_t, size_t> > & intervals);

  void    UpdateLog(std::vector<int>                              & data,
                    const std::vector<std::pair<size_t, size_t> > & intervals);

  // CLASS VARIABLES -----------------------------

  std::map<std::string, int>  n_layers_adjusted_per_interval_;   // Number of layers per interval using dz from the first interval everywhere
  std::map<std::string, int>  n_blocks_with_data_;               // Number of blocks with data per interval
  int                         n_blocks_with_data_tot_;           // Total number of blocks with data
  unsigned int                n_blocks_;                         // Number of blocks
  unsigned int                n_data_;                           // Number of non-missing
  std::string                 well_name_;

  std::map<int, std::string>  facies_map_;
  bool                        facies_log_defined_;
  //std::vector<int>   facies_numbers_;
  //std::vector<int>   facies_;

  // Blocked values
  std::vector<double> x_pos_blocked_;  // Blocked x position
  std::vector<double> y_pos_blocked_;  // Blocked y position
  std::vector<double> z_pos_blocked_;  // Blocked z position
  std::vector<int>    facies_blocked_; // Blocked facies log

  std::vector<int>    s_pos_;    // Simbox number for block
  std::vector<int>    i_pos_;    // Simbox i position for block
  std::vector<int>    j_pos_;    // Simbox j position for block
  std::vector<int>    k_pos_;    // Simbox k position for block

  int                 n_continuous_logs_;     // Number of continuous logs
  int                 n_discrete_logs_;       // Number of discrete logs

  std::map<std::string, std::vector<double> > continuous_logs_blocked_;         // Map between variable name and blocked continuous log
  std::map<std::string, std::vector<int> >    discrete_logs_blocked_;           // Map between variable name and blocked discrete log

  std::map<std::string, std::vector<double> > cont_logs_seismic_resolution_;    // Continuous logs filtered to resolution of inversion result
  std::map<std::string, std::vector<double> > cont_logs_highcut_background_;    // Continuous logs high-cut filtered to background resolution (log-domain)
  std::map<std::string, std::vector<double> > cont_logs_highcut_seismic_;       // Continuous logs high-cut filtered to approx. seismic resolution (log-domain)

  std::vector<std::vector<double> >           actual_synt_seismic_data_;        ///< Forward modelled seismic data using local wavelet
  std::vector<std::vector<double> >           well_synt_seismic_data_;          ///< Forward modelled seismic data using wavelet estimated in well

  float lateral_threshold_gradient_; //Minimum lateral distance where gradient lines must not cross
  float sigma_m_; //Smoothing factor for the gradients

  // Logs from well, not blocked
  std::vector<double> x_pos_raw_logs_;
  std::vector<double> y_pos_raw_logs_;
  std::vector<double> z_pos_raw_logs_;
  std::vector<int>    facies_raw_logs_;

  std::map<std::string, std::vector<double> > continuous_raw_logs_;  // Map between variable name and raw continuous logs
  std::map<std::string, std::vector<int> >    discrete_raw_logs_;    // Map between variable name and raw discrete logs

  //Variables needed in SetLogFromGrid and later used in WriteWell
  std::map<int, std::vector<double> >         cpp_;                        ///< Mapped reflection coefficients for angles
  std::map<std::string, std::vector<double> > continuous_logs_predicted_;  ///< Map between variable name and predicted continuous log
  std::map<int, std::vector<double> >         real_seismic_data_;          ///< Map between angle and real seismic data
  std::map<int, std::vector<double> >         facies_prob_;                ///< Map between angle and facies prob

  std::vector<double> vp_for_facies_;
  std::vector<double> rho_for_facies_;

  int                       n_angles_;                  ///< Number of AVA stacks

  bool                      interpolate_;              ///<

  int                       n_layers_;                 ///< Number of layers in estimation_simbox
  float                     dz_;                       ///< Simbox dz value for block

  int                       first_S_;                   ///< First simbox that the well log appears in
  int                       last_S_;                    ///< Last simbox that the well log appears in

  int                       first_M_;                   ///< First well log entry contributing to blocked well
  int                       last_M_;                    ///< Last well log entry contributing to blocked well

  int                       first_B_;                   ///< First block with contribution from well log
  int                       last_B_;                    ///< Last block with contribution from well log

  std::vector<int>          n_well_log_obs_in_interval_;      ///< Number of well log observations in each interval

  bool                      is_deviated_;

  //Used in logging in crava.cpp. Copied from well
  int                       real_vs_log_;                    //Uses the indicator enum from Modelsettings
  int                       use_for_facies_probabilities_;   //Uses the indicator enum from Modelsettings

  int                       use_for_background_trend_;       //Uses the indicator enum from Modelsettings
  int                       use_for_filtering_;              //Uses the indicator enum from Modelsettings
  int                       use_for_wavelet_estimation_;     //Uses the indicator enum from Modelsettings

};
#endif
