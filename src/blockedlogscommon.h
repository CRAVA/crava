/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef BLOCKEDLOGSCOMMON_H
#define BLOCKEDLOGSCOMMON_H

#include "src/seismicstorage.h"
#include "nrlib/well/well.hpp"

class BlockedLogsCommon
{

public:
  BlockedLogsCommon(const NRLib::Well   * const well_data,
                    const Simbox        * const estimation_simbox,
                    bool                  interpolate);

  ~ BlockedLogsCommon(void);


  //GET FUNCTIONS --------------------------------

  int                                    GetNumberOfBlocks()   const   { return n_blocks_                                                     ;}
  const std::vector<double>            & GetXpos(void)         const   { return continuous_logs_[continuous_log_names_.find("X_pos")->second] ;}
  const std::vector<double>            & GetYpos(void)         const   { return continuous_logs_[continuous_log_names_.find("Y_pos")->second] ;}
  const std::vector<double>            & GetZpos(void)         const   { return continuous_logs_[continuous_log_names_.find("TVD")->second]   ;}
  const std::vector<double>            & GetTVD(void)          const   { return continuous_logs_[continuous_log_names_.find("TVD")->second]   ;}
  const std::vector<double>            & GetTWT(void)          const   { return continuous_logs_[continuous_log_names_.find("TWT")->second]   ;}
  const std::vector<double>            & GetVp(void)          const    { return continuous_logs_[continuous_log_names_.find("Vp")->second]   ;}
  const std::vector<double>            & GetVs(void)          const    { return continuous_logs_[continuous_log_names_.find("Vs")->second]   ;}
  const std::vector<double>            & GetRho(void)          const   { return continuous_logs_[continuous_log_names_.find("Rho")->second]   ;}
  bool                                   HasContLog(std::string s)     { return continuous_log_names_.find(s) != continuous_log_names_.end()  ;}
  bool                                   HasDiscLog(std::string s)     { return discrete_log_names_.find(s) != discrete_log_names_.end()      ;}
  const std::vector<int>                 GetIposVector(void)   const   { return i_pos_                                                        ;}
  const std::vector<int>                 GetJposVector(void)   const   { return j_pos_                                                        ;}
  const std::string                      GetWellName(void)     const   { return well_name_                                                    ;}


  // FUNCTIONS -----------------------------------

  void    SetSeismicGradient(double                            v0,
                             const NRLib::Grid2D<float>   &    structure_depth_grad_x,
                             const NRLib::Grid2D<float>   &    structure_depth_grad_y,
                             const NRLib::Grid2D<float>   &    ref_time_grad_x,
                             const NRLib::Grid2D<float>   &    ref_time_grad_y,
                             std::vector<double>          &    x_gradient,
                             std::vector<double>          &    y_gradient);

  void    SetTimeGradientSettings(float distance, float sigma_m);

  void    FindSeismicGradient(const std::vector<SeismicStorage> & seismic_data,
                              const Simbox                      * const estimation_simbox,
                              int                                 n_angles,
                              std::vector<double>               & x_gradient,
                              std::vector<double>               & y_gradient,
                              std::vector<std::vector<double> > & sigma_gradient);

  void    GetBlockedGrid(const Simbox         * estimation_simbox,
                         const SeismicStorage * seismic_data,
                         double               * blockedLog,
                         int                    iOffset = 0,
                         int                    jOffset = 0);

  void    GetVerticalTrend(const std::vector<double> & blocked_log,
                           double * trend);

  void    FindContiniousPartOfData(const std::vector<bool> & has_data,
                                   int                       nz,
                                   int                     & start,
                                   int                     & length);

  void    FillInCpp(const float * coeff,
                    int           start,
                    int           length,
                    fftw_real   * cpp_r,
                    int           nzp);

  void    FillInSeismic(double    * seismicData,
                        int         start,
                        int         length,
                        fftw_real * seis_r,
                        int         nzp) const;

  void    EstimateCor(fftw_complex * var1_c,
                      fftw_complex * var2_c,
                      fftw_complex * ccor_1_2_c,
                      int            cnzp) const;

  void    SetLogFromVerticalTrend(float      * vertical_trend,
                                  double       z0,              // z-value of center in top layer
                                  double       dz,              // dz in vertical trend
                                  int          nz,              // layers in vertical trend
                                  std::string  type,
                                  int          i_angle);

private:

  // FUNCTIONS------------------------------------

  void    BlockWell(const NRLib::Well                  * const well_data,
                    const Simbox                       * const estimation_simbox,
                    std::map<std::string, int>         & continuous_log_names,
                    std::map<std::string, int>         & discrete_log_names,
                    std::vector<std::vector<double> >  & continuous_logs,
                    std::vector<std::vector<int> >     & discrete_logs,
                    bool                                 interpolate);

  void    FindSizeAndBlockPointers(const Simbox         * const estimation_simbox,
                                   std::vector<int>     & bInd);

  void    FindBlockIJK(const Simbox                     * const estimation_simbox,
                       const std::vector<int>           & bInd);

  void    BlockCoordinateLog(const std::vector<int>     &  b_ind,
                             const std::vector<double>  &  coord,
                             std::vector<double>        &  blocked_coord);

  void    BlockContinuousLog(const int                 *  b_ind,
                             const std::vector<double> &  well_log,
                             std::vector<double>       &  blocked_log);

  void    InterpolateContinuousLog(std::vector<double>   & blocked_log,
                                   int                     start,
                                   int                     end,
                                   int                     index,
                                   float                   rel);

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

  void    InterpolateTrend(const std::vector<double> & blocked_log,
                           double * trend);

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

  // CLASS VARIABLES -----------------------------

  std::string        well_name_;        ///< Name of well

  int                n_blocks_;         // number of blocks in well log

  std::vector<double> x_pos_;  // Simbox x position
  std::vector<double> y_pos_;  // Simbox y position
  std::vector<double> z_pos_;  // Simbox z position

  //std::vector<float> twt_;    // Two-way travel time (ms)
  //std::vector<float> tvd_;    // True vertical depth

  std::vector<int> i_pos_;    // Simbox i position
  std::vector<int> j_pos_;    // Simbox j position
  std::vector<int> k_pos_;    // Simbox k position

  std::map<std::string, int> continuous_log_names_; //
  std::map<std::string, int> discrete_log_names_;   //

  std::map<std::string, int> continuous_log_names_blocked_; //
  std::map<std::string, int> discrete_log_names_blocked_;   //

  int n_continuous_logs_;
  int n_discrete_logs_;

  std::vector<std::vector<double> > continuous_logs_;
  std::vector<std::vector<int> >    discrete_logs_;

  std::vector<std::vector<double> > continuous_logs_blocked_;
  std::vector<std::vector<int> >    discrete_logs_blocked_;

  ///H Copied from blockedlogs.h, they are only set in SetLogFromVerticalTrend. Are they needed later?
  float        * alpha_seismic_resolution_; ///<
  float        * beta_seismic_resolution_;  ///< Logs filtered to resolution of inversion result
  float        * rho_seismic_resolution_;   ///<
  float       ** actual_synt_seismic_data_; ///< Forward modelled seismic data using local wavelet
  float       ** well_synt_seismic_data_;   ///< Forward modelled seismic data using wavelet estimated in well
  int            n_angles_;                  ///< Number of AVA stacks

  int                       n_layers_;                 ///< Number of layers in estimation_simbox

  float                     dz_;                       ///< Simbox dz value for block

  int                       first_M_;                   ///< First well log entry contributing to blocked well
  int                       last_M_;                    ///< Last well log entry contributing to blocked well


  int                       first_B_;                   ///< First block with contribution from well log
  int                       last_B_;                    ///< Last block with contribution from well log

  float                     lateral_threshold_gradient_;  //Minimum lateral distance where gradient lines must not cross
  float                     sigma_m_;                     //Smoothing factor for the gradients

  bool                      interpolate_;                 //Should vertical interpolation be done? (for Roxar)
};

#endif
