/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/blockedlogscommon.h"
#include "src/simbox.h"
#include "lib/utils.h"
#include "nrlib/flens/nrlib_flens.hpp"
#include "src/multiintervalgrid.h"
#include "src/wavelet1D.h"
#include "src/cravatrend.h"

BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well                * well_data,
                                     const std::vector<std::string>   & cont_logs_to_be_blocked,
                                     const std::vector<std::string>   & disc_logs_to_be_blocked,
                                     const Simbox                     * const estimation_simbox,
                                     bool                               interpolate,
                                     bool                               restrict_to_visible,
                                     bool                             & is_inside,
                                     std::string                      & err_text) {

  n_angles_                     = 0;
  well_name_                    = well_data->GetWellName();
  n_layers_                     = estimation_simbox->getnz();
  n_blocks_                     = 0;
  n_blocks_with_data_.insert(std::pair<std::string, int>(estimation_simbox->GetIntervalName(),0));
  n_blocks_with_data_tot_       = 0;
  interpolate_                  = interpolate;
  bool failed                   = false;
  is_deviated_                  = well_data->IsDeviated();
  use_for_facies_probabilities_ = well_data->GetUseForFaciesProbabilities();
  real_vs_log_                  = well_data->GetRealVsLog();
  use_for_wavelet_estimation_   = well_data->GetUseForWaveletEstimation();
  use_for_background_trend_     = well_data->GetUseForBackgroundTrend();
  use_for_filtering_            = well_data->GetUseForFiltering();
  use_for_rock_physics_         = well_data->GetUseForRockPhysics();
  well_name_                    = well_data->GetWellName();
  facies_log_defined_           = false;
  first_S_                      = IMISSING;
  last_S_                       = IMISSING;
  first_M_                      = IMISSING;
  last_M_                       = IMISSING;
  first_B_                      = IMISSING;
  last_B_                       = IMISSING;
  lateral_threshold_gradient_   = RMISSING;
  sigma_m_                      = RMISSING;
  //n_blocks_with_data_.resize(1,0);

  std::string err_text_tmp = "";

  // FACIES
  if (well_data->HasDiscLog("Facies")) {
    facies_log_defined_ = true;
    facies_map_ = well_data->GetFaciesMap();
  }
  lateral_threshold_gradient_ = 0; //Minimum lateral distance where gradient lines must not cross
  sigma_m_ = 0; //Smoothing factor for the gradients

  // 20130627 EN: Missing data are removed upon construction of a well_data object, whereas
  // NRLib::Well objects, which are used here, keep the logs as they are in the input files.
  RemoveMissingLogValues(well_data, x_pos_raw_logs_, y_pos_raw_logs_, z_pos_raw_logs_,
                         facies_raw_logs_, continuous_raw_logs_, discrete_raw_logs_, cont_logs_to_be_blocked,
                         disc_logs_to_be_blocked, n_data_, failed, err_text_tmp);

  if (failed)
    err_text_tmp += "Logs were not successfully read from well " + well_name_ +".\n";

  is_inside = true; //If well is insde simbox

  if (!failed)
    BlockWell(estimation_simbox, well_data, continuous_raw_logs_, discrete_raw_logs_,
              facies_map_, continuous_logs_blocked_, cont_logs_highcut_seismic_,
              cont_logs_highcut_background_, discrete_logs_blocked_,
              x_pos_blocked_, y_pos_blocked_, z_pos_blocked_, facies_blocked_,
              n_data_, i_pos_, j_pos_, k_pos_, first_M_, last_M_, first_B_, last_B_, n_blocks_, n_blocks_with_data_,
              n_blocks_with_data_tot_, facies_log_defined_, interpolate, restrict_to_visible, dz_, failed, is_inside, err_text_tmp);

  if (err_text_tmp != "") {
    LogKit::LogFormatted(LogKit::Low,"\nBlocking of well " + well_name_ + " in simbox failed:\n");
    LogKit::LogFormatted(LogKit::Low, err_text_tmp + "\n");

    err_text += err_text_tmp;
  }

  n_continuous_logs_ = static_cast<int>(continuous_logs_blocked_.size());
  n_discrete_logs_   = static_cast<int>(discrete_logs_blocked_.size());

}

// For correlation estimation
BlockedLogsCommon::BlockedLogsCommon(NRLib::Well                      * well_data,
                                     const std::vector<std::string>   & cont_logs_to_be_blocked,
                                     const std::vector<std::string>   & disc_logs_to_be_blocked,
                                     const MultiIntervalGrid          * multiple_interval_grid,
                                     bool                               interpolate,
                                     bool                             & is_inside,
                                     std::string                      & err_text) {
  n_angles_                     = 0;
  well_name_                    = well_data->GetWellName();
  n_layers_                     = 0;
  n_blocks_                     = 0;
  n_blocks_with_data_tot_       = 0;
  real_vs_log_                  = well_data->GetRealVsLog();
  use_for_facies_probabilities_ = well_data->GetUseForFaciesProbabilities();
  interpolate_                  = interpolate;
  is_deviated_                  = well_data->IsDeviated();
  use_for_wavelet_estimation_   = well_data->GetUseForWaveletEstimation();
  use_for_background_trend_     = well_data->GetUseForBackgroundTrend();
  use_for_filtering_            = well_data->GetUseForFiltering();
  use_for_rock_physics_         = well_data->GetUseForRockPhysics();
  facies_log_defined_           = false;
  lateral_threshold_gradient_   = RMISSING;
  sigma_m_                      = RMISSING;
  first_S_                      = IMISSING;
  last_S_                       = IMISSING;
  first_M_                      = IMISSING;
  last_M_                       = IMISSING;
  first_B_                      = IMISSING;
  last_B_                       = IMISSING;

  std::string err_text_tmp = "";

  //n_blocks_with_data_.resize(multiple_interval_grid->GetNIntervals(),0);
  const std::vector<Simbox *> interval_simboxes = multiple_interval_grid->GetIntervalSimboxes();
  facies_log_defined_ = false;
  for (int i = 0; i < multiple_interval_grid->GetNIntervals(); i++) {
    n_blocks_with_data_.insert(std::pair<std::string, int>(multiple_interval_grid->GetIntervalName(i), 0));
  }

  bool failed = false;

  // FACIES
  if (well_data->HasDiscLog("Facies")) {
    facies_log_defined_ = true;
    facies_map_ = well_data->GetFaciesMap();
  }

  RemoveMissingLogValues(well_data, x_pos_raw_logs_, y_pos_raw_logs_, z_pos_raw_logs_,
                         facies_raw_logs_, continuous_raw_logs_, discrete_raw_logs_, cont_logs_to_be_blocked,
                         disc_logs_to_be_blocked, n_data_, failed, err_text_tmp);

  if (failed)
    err_text_tmp += "Logs were not successfully read from well " + well_name_ +".\n";

  if (!failed) {
    if (interval_simboxes.size() == 1) {
      //
      // If there is only one interval
      //
      n_layers_ = interval_simboxes[0]->getnz();
      BlockWell(interval_simboxes[0],
                well_data,
                continuous_raw_logs_,
                discrete_raw_logs_,
                facies_map_,
                continuous_logs_blocked_,
                cont_logs_highcut_seismic_,
                cont_logs_highcut_background_,
                discrete_logs_blocked_,
                x_pos_blocked_,
                y_pos_blocked_,
                z_pos_blocked_,
                facies_blocked_,
                n_data_,
                i_pos_,
                j_pos_,
                k_pos_,
                first_M_,
                last_M_,
                first_B_,
                last_B_,
                n_blocks_,
                n_blocks_with_data_,
                n_blocks_with_data_tot_,
                facies_log_defined_,
                interpolate,
                false,
                dz_,
                failed,
                is_inside,
                err_text_tmp);
    }
    else {
      //
      // If there are multiple intervals: Neighbouring cells vertically
      // must have the same correlation within each trace for correlation estimation.
      //
      BlockWellForCorrelationEstimation(multiple_interval_grid,
                                        well_data,
                                        continuous_raw_logs_,
                                        discrete_raw_logs_,
                                        facies_map_,
                                        x_pos_blocked_,
                                        y_pos_blocked_,
                                        z_pos_blocked_,
                                        facies_blocked_,
                                        continuous_logs_blocked_,
                                        cont_logs_highcut_seismic_,
                                        cont_logs_highcut_background_,
                                        discrete_logs_blocked_,
                                        n_data_,
                                        n_blocks_,
                                        n_blocks_with_data_,
                                        n_blocks_with_data_tot_,
                                        n_well_log_obs_in_interval_,
                                        n_layers_adjusted_per_interval_,
                                        i_pos_,
                                        j_pos_,
                                        k_pos_,
                                        s_pos_,
                                        first_M_,
                                        last_M_,
                                        first_S_,
                                        last_S_,
                                        first_B_,
                                        last_B_,
                                        facies_log_defined_,
                                        interpolate,
                                        n_layers_,
                                        dz_,
                                        failed,
                                        err_text_tmp);
    }
  }

  if (err_text_tmp != "") {
    LogKit::LogFormatted(LogKit::Low,"\nBlocking of well " + well_name_ + " in the outer estimation simbox failed:\n");
    LogKit::LogFormatted(LogKit::Low, err_text_tmp + "\n");

    err_text += err_text_tmp;
  }

  n_continuous_logs_ = static_cast<int>(continuous_logs_blocked_.size());
  n_discrete_logs_   = static_cast<int>(discrete_logs_blocked_.size());
}

//Blocked logs for RockPhysics
BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well              * well,
                                     const Simbox                   * simbox,
                                     const CravaTrend               & trend_cubes,
                                     const std::vector<std::string> & cont_logs_to_be_blocked,
                                     const std::vector<std::string> & disc_logs_to_be_blocked,
                                     std::string                    & err_text)
: lateral_threshold_gradient_(RMISSING),
  sigma_m_(RMISSING),
  first_S_(IMISSING),
  last_S_ (IMISSING),
  first_M_(IMISSING),
  last_M_ (IMISSING),
  first_B_(IMISSING),
  last_B_ (IMISSING)

{
  use_for_facies_probabilities_ = well->GetUseForFaciesProbabilities();
  use_for_wavelet_estimation_   = well->GetUseForWaveletEstimation();
  use_for_background_trend_     = well->GetUseForBackgroundTrend();
  use_for_filtering_            = well->GetUseForFiltering();
  use_for_rock_physics_         = well->GetUseForRockPhysics();
  real_vs_log_                  = well->GetRealVsLog();
  is_deviated_                  = well->IsDeviated();

  bool failed          = false;
  std::string tmp_err_text = "";

  // FACIES
  if (well->HasDiscLog("Facies")) {
    facies_log_defined_ = true;
    facies_map_ = well->GetFaciesMap();
  }
  int n_facies = well->GetNFacies();
  facies_names_.resize(n_facies);

  std::vector<int> facies_nr(n_facies);

  std::map<int, std::string>  facies_map = well->GetFaciesMap();
  int i = 0;
  for (std::map<int,std::string>::const_iterator it = facies_map.begin(); it != facies_map.end(); it++) {
    facies_names_[i] = it->second;
    facies_nr[i] = i;
    i++;
  }

  //First run RemoveMissingLogValues since x_pos etc. are needed in FindSizeAndBlockPointers
  RemoveMissingLogValues(well, x_pos_raw_logs_, y_pos_raw_logs_, z_pos_raw_logs_,
                         facies_raw_logs_, continuous_raw_logs_, discrete_raw_logs_, cont_logs_to_be_blocked,
                         disc_logs_to_be_blocked, n_data_, failed, tmp_err_text);

  std::vector<int> b_ind(n_data_); // Gives which block each well log entry contributes to
  n_layers_ = simbox->getnz();

  bool is_inside = true;
  FindSizeAndBlockPointers(simbox, true, b_ind, n_layers_, first_M_, last_M_, n_blocks_, is_inside);

  if (is_inside == false)
    tmp_err_text += "Error when blocking logs for RockPhysics: Well "+well->GetWellName()+" was not found within the estimation simbox surrounding the inversion intervals.\n";

  if (tmp_err_text != "") {
    double dz;
    FindBlockIJK(simbox, b_ind, first_M_, last_M_, first_B_, last_B_, i_pos_, j_pos_, k_pos_, dz);

    FindTrendPositions(i_pos_, j_pos_, k_pos_, n_blocks_, trend_cubes, s1_, s2_);

    std::vector<double> blocked_vp;
    std::vector<double> blocked_vs;
    std::vector<double> blocked_rho;
    std::vector<double> blocked_porosity;
    std::vector<int>    blocked_facies;

    if (well->HasContLog("Vp"))
      BlockContinuousLog(b_ind, well->GetContLog("Vp"),       blocked_vp);

    if (well->HasContLog("Vs"))
      BlockContinuousLog(b_ind, well->GetContLog("Vs"),       blocked_vs);

    if (well->HasContLog("Rho"))
      BlockContinuousLog(b_ind, well->GetContLog("Rho"),      blocked_rho);

    if (well->HasContLog("Porosity"))
      BlockPorosityLog(b_ind, well->GetContLog("Porosity"), blocked_porosity);

    BlockFaciesLog(b_ind, facies_raw_logs_, facies_map, static_cast<int>(facies_map.size()), blocked_facies);

    vp_for_facies_.resize(n_facies,       std::vector<double>(n_blocks_, RMISSING));
    vs_for_facies_.resize(n_facies,       std::vector<double>(n_blocks_, RMISSING));
    rho_for_facies_.resize(n_facies,      std::vector<double>(n_blocks_, RMISSING));
    porosity_for_facies_.resize(n_facies, std::vector<double>(n_blocks_, RMISSING));

    AssignToFacies(blocked_vp,       blocked_facies, facies_nr, vp_for_facies_);
    AssignToFacies(blocked_vs,       blocked_facies, facies_nr, vs_for_facies_);
    AssignToFacies(blocked_rho,      blocked_facies, facies_nr, rho_for_facies_);
    AssignToFacies(blocked_porosity, blocked_facies, facies_nr, porosity_for_facies_);

    CalculateBulkShear(n_blocks_, n_facies, bulk_modulus_, shear_modulus_);
  }
  else
    err_text += tmp_err_text;


}

BlockedLogsCommon::BlockedLogsCommon(const BlockedLogsCommon & logs)
{
  n_layers_adjusted_per_interval_ = logs.n_layers_adjusted_per_interval_;
  n_blocks_with_data_             = logs.n_blocks_with_data_;
  n_blocks_with_data_tot_         = logs.n_blocks_with_data_tot_;
  n_blocks_                       = logs.n_blocks_;
  n_data_                         = logs.n_data_;
  well_name_                      = logs.well_name_;

  facies_map_         = logs.facies_map_;
  facies_log_defined_ = logs.facies_log_defined_;

  // Blocked values
  x_pos_blocked_  = logs.x_pos_blocked_;
  y_pos_blocked_  = logs.y_pos_blocked_;
  z_pos_blocked_  = logs.z_pos_blocked_;
  facies_blocked_ = logs.facies_blocked_;

  s_pos_ = logs.s_pos_;
  i_pos_ = logs.i_pos_;
  j_pos_ = logs.j_pos_;
  k_pos_ = logs.k_pos_;

  n_continuous_logs_ = logs.n_continuous_logs_;
  n_discrete_logs_   = logs.n_discrete_logs_;

  continuous_logs_blocked_ = logs.continuous_logs_blocked_;
  discrete_logs_blocked_   = logs.discrete_logs_blocked_;

  cont_logs_seismic_resolution_ = logs.cont_logs_seismic_resolution_;
  cont_logs_highcut_background_ = logs.cont_logs_highcut_background_;
  cont_logs_highcut_seismic_    = logs.cont_logs_highcut_seismic_;

  actual_synt_seismic_data_ = logs.actual_synt_seismic_data_;
  well_synt_seismic_data_   = logs.well_synt_seismic_data_;

  lateral_threshold_gradient_ = logs.lateral_threshold_gradient_;
  sigma_m_                    = logs.sigma_m_;

  // Logs from well, not blocked
  x_pos_raw_logs_  = logs.x_pos_raw_logs_;
  y_pos_raw_logs_  = logs.y_pos_raw_logs_;
  z_pos_raw_logs_  = logs.z_pos_raw_logs_;
  facies_raw_logs_ = logs.facies_raw_logs_;

  continuous_raw_logs_ = logs.continuous_raw_logs_;
  discrete_raw_logs_   = logs.discrete_raw_logs_;

  //Variables needed in SetLogFromGrid and later used in WriteWell
  continuous_logs_predicted_ = logs.continuous_logs_predicted_;
  real_seismic_data_         = logs.real_seismic_data_;
  facies_prob_               = logs.facies_prob_;

  vp_facies_filtered_         = logs.vp_facies_filtered_;
  rho_facies_filtered_        = logs.rho_facies_filtered_;

  n_angles_                   = logs.n_angles_;

  interpolate_                = logs.interpolate_;

  n_layers_                   = logs.n_layers_;
  dz_                         = logs.dz_;

  first_S_                    = logs.first_S_;
  last_S_                     = logs.last_S_;

  first_M_                    = logs.first_M_;
  last_M_                     = logs.last_M_;

  first_B_                    = logs.first_B_;
  last_B_                     = logs.last_B_;

  n_well_log_obs_in_interval_ = logs.n_well_log_obs_in_interval_;

  is_deviated_                = logs.is_deviated_;

  //Used in logging in crava.cpp. Copied from well
  real_vs_log_                  = logs.real_vs_log_;
  use_for_facies_probabilities_ = logs.use_for_facies_probabilities_;
  use_for_background_trend_     = logs.use_for_background_trend_;
  use_for_filtering_            = logs.use_for_filtering_;
  use_for_wavelet_estimation_   = logs.use_for_wavelet_estimation_;
  use_for_rock_physics_         = logs.use_for_rock_physics_;

  //Variables used for rockphysics
  facies_names_                 = logs.facies_names_;
  vp_for_facies_                = logs.vp_for_facies_;
  vs_for_facies_                = logs.vs_for_facies_;
  rho_for_facies_               = logs.rho_for_facies_;
  porosity_for_facies_          = logs.porosity_for_facies_;
  bulk_modulus_                 = logs.bulk_modulus_;
  shear_modulus_                = logs.shear_modulus_;
  s1_                           = logs.s1_;
  s2_                           = logs.s2_;
}

BlockedLogsCommon::~BlockedLogsCommon() {

}

void BlockedLogsCommon::BlockWellForCorrelationEstimation(const MultiIntervalGrid                             * multiple_interval_grid,
                                                          const NRLib::Well                                   * well,
                                                          const std::map<std::string, std::vector<double> >   & continuous_logs_raw_logs,
                                                          const std::map<std::string, std::vector<int> >      & discrete_raw_logs,
                                                          const std::map<int, std::string>                    & facies_map,
                                                          std::vector<double>                                 & x_pos_blocked,
                                                          std::vector<double>                                 & y_pos_blocked,
                                                          std::vector<double>                                 & z_pos_blocked,
                                                          std::vector<int>                                    & facies_blocked,
                                                          std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                                                          std::map<std::string, std::vector<double> >         & cont_logs_highcut_seismic,
                                                          std::map<std::string, std::vector<double> >         & cont_logs_highcut_background,
                                                          std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                                                          unsigned int                                          n_data,
                                                          unsigned int                                        & n_blocks,
                                                          std::map<std::string, int>                          & n_blocks_with_data,
                                                          int                                                 & n_blocks_with_data_tot,
                                                          std::vector<int>                                    & n_well_log_obs_in_interval,
                                                          std::map<std::string, int>                          & n_layers_adjusted_per_interval,
                                                          std::vector<int>                                    & i_pos,
                                                          std::vector<int>                                    & j_pos,
                                                          std::vector<int>                                    & k_pos,
                                                          std::vector<int>                                    & s_pos,
                                                          int                                                 & first_M,
                                                          int                                                 & last_M,
                                                          int                                                 & first_S,
                                                          int                                                 & last_S,
                                                          int                                                 & first_B,
                                                          int                                                 & last_B,
                                                          bool                                                  facies_log_defined,
                                                          bool                                                  interpolate,
                                                          int                                                 & n_layers,
                                                          double                                              & dz,
                                                          bool                                                & failed,
                                                          std::string                                         & err_text) const{
  std::vector<int> b_ind(n_data); // Which block each well log entry contributes to

  try{

    FindSizeAndBlockPointers(multiple_interval_grid,
                             b_ind,
                             n_data,
                             n_layers,
                             first_M,
                             last_M,
                             first_S,
                             last_S,
                             n_blocks,
                             n_layers_adjusted_per_interval,
                             err_text);
    if (err_text == "") {
      FindBlockIJK(multiple_interval_grid, b_ind, GetXposRawLogs(), GetYposRawLogs(), GetZposRawLogs(),
                   first_M, last_M, first_S, last_S, first_B, last_B, n_well_log_obs_in_interval,
                   i_pos, j_pos, k_pos, s_pos, dz);

      // Coordinate logs and necessary logs

      BlockCoordinateLog(b_ind, x_pos_raw_logs_, x_pos_blocked);
      BlockCoordinateLog(b_ind, y_pos_raw_logs_, y_pos_blocked);
      BlockCoordinateLog(b_ind, z_pos_raw_logs_, z_pos_blocked);

      //Extrapolates if missing values in the beginning or end.Not needed for estimation of correlations.
      //FindXYZForVirtualPart(multiple_interval_grid, i_pos, j_pos, k_pos, n_blocks, first_B, last_B,
      //                      first_S, last_S, x_pos_blocked, y_pos_blocked, z_pos_blocked);

      //HighCutSeismic blocked logs
      const std::map<std::string, std::vector<double> > & seismic_resolution_logs = well->GetContLogSeismicResolution();
      for (std::map<std::string, std::vector<double> >::const_iterator it = seismic_resolution_logs.begin(); it!=seismic_resolution_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        cont_logs_highcut_seismic.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      //Create HighCutBackground blocked logs
      const std::map<std::string, std::vector<double> > & background_resolution_logs = well->GetContLogBackgroundResolution();
      for (std::map<std::string, std::vector<double> >::const_iterator it = background_resolution_logs.begin(); it!=background_resolution_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        cont_logs_highcut_background.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      // Continuous logs
      for (std::map<std::string, std::vector<double> >::const_iterator it = continuous_logs_raw_logs.begin(); it!=continuous_logs_raw_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        continuous_logs_blocked.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      CountBlocksWithDataPerInterval(multiple_interval_grid, x_pos_blocked_, y_pos_blocked_, z_pos_blocked_, continuous_logs_blocked,
                                      n_blocks_, n_blocks_with_data, n_blocks_with_data_tot);

      // Discrete logs
      if (facies_log_defined)
        BlockFaciesLog(b_ind, facies_raw_logs_, facies_map, static_cast<int>(facies_map.size()), facies_blocked);

      (void) discrete_raw_logs;
      (void) discrete_logs_blocked;

      if (interpolate) {
        for (unsigned int i=1;i<n_data;i++) {
          if (abs(b_ind[i]-b_ind[i-1]) > 1) {
            int start, end;
            if (b_ind[i] > b_ind[i-1]) {
              start = b_ind[i-1];
              end   = b_ind[i];
            }
            else {
              start = b_ind[i];
              end   = b_ind[i-1];
            }
            for (int j = start+1;j<end;j++) {
              float t = static_cast<float>(j-start)/static_cast<float>(end-start);
              // Coordinate logs
              InterpolateContinuousLog(x_pos_blocked, start, end, j, t);
              InterpolateContinuousLog(y_pos_blocked, start, end, j, t);
              InterpolateContinuousLog(z_pos_blocked, start, end, j, t);

              // interpolate all blocked continuous logs
              for (std::map<std::string,std::vector<double> >::iterator it=continuous_logs_blocked.begin(); it!=continuous_logs_blocked.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }
              // interpolate all continuous logs with upper frequency cut
              for (std::map<std::string,std::vector<double> >::iterator it=cont_logs_highcut_seismic.begin(); it!=cont_logs_highcut_seismic.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }
              // interpolate all continuous logs with upper frequency cut
              for (std::map<std::string,std::vector<double> >::iterator it=cont_logs_highcut_background.begin(); it!=cont_logs_highcut_background.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }

            }
          }
        }
      }
    }
  }catch(NRLib::Exception & e) {
    err_text += "Error blocking logs for correlation estimation:\n";
    err_text += std::string(e.what());
    failed = true;
  }
}

void BlockedLogsCommon::BlockWell(const Simbox                                        * estimation_simbox,
                                  const NRLib::Well                                   * well,
                                  const std::map<std::string, std::vector<double> >   & continuous_logs_raw_logs,
                                  const std::map<std::string, std::vector<int> >      & discrete_logs_raw_logs,
                                  const std::map<int, std::string>                    & facies_map,
                                  std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                                  std::map<std::string, std::vector<double> >         & cont_logs_highcut_seismic,
                                  std::map<std::string, std::vector<double> >         & cont_logs_highcut_background,
                                  std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                                  std::vector<double>                                 & x_pos_blocked,
                                  std::vector<double>                                 & y_pos_blocked,
                                  std::vector<double>                                 & z_pos_blocked,
                                  std::vector<int>                                    & facies_blocked,
                                  unsigned int                                          n_data,
                                  std::vector<int>                                    & i_pos,
                                  std::vector<int>                                    & j_pos,
                                  std::vector<int>                                    & k_pos,
                                  int                                                 & first_M,
                                  int                                                 & last_M,
                                  int                                                 & first_B,
                                  int                                                 & last_B,
                                  unsigned int                                        & n_blocks,
                                  std::map<std::string, int>                          & n_blocks_with_data,
                                  int                                                 & n_blocks_with_data_tot,
                                  bool                                                  facies_log_defined,
                                  bool                                                  interpolate,
                                  bool                                                  visible_only,
                                  double                                              & dz,
                                  bool                                                & failed,
                                  bool                                                & is_inside,
                                  std::string                                         & err_text) const{

  std::vector<int> b_ind(n_data); // Which block each well log entry contributes to

  try {

    FindSizeAndBlockPointers(estimation_simbox, visible_only, b_ind, n_layers_, first_M, last_M, n_blocks, is_inside);
    if (is_inside == true) {
      FindBlockIJK(estimation_simbox,
                   b_ind,
                   first_M,
                   last_M,
                   first_B,
                   last_B,
                   i_pos,
                   j_pos,
                   k_pos,
                   dz);

      // Coordinate logs and necessary logs

      BlockCoordinateLog(b_ind, x_pos_raw_logs_, x_pos_blocked);
      BlockCoordinateLog(b_ind, y_pos_raw_logs_, y_pos_blocked);
      BlockCoordinateLog(b_ind, z_pos_raw_logs_, z_pos_blocked);

      //Extrapolates if missing values in the beginning or end.
      FindXYZForVirtualPart(estimation_simbox, i_pos, j_pos, k_pos, n_blocks, first_B, last_B,
                                 x_pos_blocked, y_pos_blocked, z_pos_blocked);

      // Continuous logs
      for (std::map<std::string, std::vector<double> >::const_iterator it = continuous_logs_raw_logs.begin(); it!=continuous_logs_raw_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        continuous_logs_blocked.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      //HighCutSeismic blocked logs
      const std::map<std::string, std::vector<double> > & seismic_resolution_logs = well->GetContLogSeismicResolution();
      for (std::map<std::string, std::vector<double> >::const_iterator it = seismic_resolution_logs.begin(); it!=seismic_resolution_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        cont_logs_highcut_seismic.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      //Create HighCutBackground blocked logs
      const std::map<std::string, std::vector<double> > & background_resolution_logs = well->GetContLogBackgroundResolution();
      for (std::map<std::string, std::vector<double> >::const_iterator it = background_resolution_logs.begin(); it!=background_resolution_logs.end(); it++) {
        std::vector<double> temp_vector_blocked;
        BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
        cont_logs_highcut_background.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
      }

      CountBlocksWithData(x_pos_blocked, y_pos_blocked, z_pos_blocked, continuous_logs_blocked,
                          estimation_simbox->GetIntervalName(), n_blocks, n_blocks_with_data, n_blocks_with_data_tot);

      // Discrete logs
      if (facies_log_defined)
        BlockFaciesLog(b_ind, facies_raw_logs_, facies_map, static_cast<int>(facies_map.size()), facies_blocked);

      (void) discrete_logs_raw_logs;
      (void) discrete_logs_blocked;

      if (interpolate) {
        for (unsigned int i=1;i<n_data;i++) {
          if (abs(b_ind[i]-b_ind[i-1]) > 1) {
            int start, end;
            if (b_ind[i] > b_ind[i-1]) {
              start = b_ind[i-1];
              end   = b_ind[i];
            }
            else {
              start = b_ind[i];
              end   = b_ind[i-1];
            }
            for (int j = start+1;j<end;j++) {
              float t = static_cast<float>(j-start)/static_cast<float>(end-start);
              // Coordinate logs
              InterpolateContinuousLog(x_pos_blocked, start, end, j, t);
              InterpolateContinuousLog(y_pos_blocked, start, end, j, t);
              InterpolateContinuousLog(z_pos_blocked, start, end, j, t);

              // all blocked continuous logs
              for (std::map<std::string,std::vector<double> >::iterator it=continuous_logs_blocked.begin(); it!=continuous_logs_blocked.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }
              for (std::map<std::string,std::vector<double> >::iterator it=cont_logs_highcut_seismic.begin(); it!=cont_logs_highcut_seismic.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }
              for (std::map<std::string,std::vector<double> >::iterator it=cont_logs_highcut_background.begin(); it!=cont_logs_highcut_background.end(); ++it) {
                InterpolateContinuousLog(it->second, start, end, j, t);
              }

            }
          }
        }
      }
    }

  }catch(NRLib::Exception & e) {
    err_text += "Error blocking wells:\n";
    err_text += std::string(e.what());
    failed = true;
  }

}

// For correlation estimation, multiple intervals -------------------------------------------------------------------------------
void  BlockedLogsCommon::FindSizeAndBlockPointers(const MultiIntervalGrid       * multiple_interval_grid,
                                                  std::vector<int>              & b_ind,
                                                  int                             n_data,
                                                  int                           & n_layers,
                                                  int                           & first_M,
                                                  int                           & last_M,
                                                  int                           & first_S,
                                                  int                           & last_S,
                                                  unsigned int                  & n_blocks,
                                                  std::map<std::string, int>    & n_layers_adjusted_per_interval,
                                                  std::string                   & err_text) const{

  const std::vector<Simbox *> interval_simboxes = multiple_interval_grid->GetIntervalSimboxes();
  int   n_intervals                             = static_cast<int>(interval_simboxes.size());
  const std::vector<double> & x_pos             = GetXposRawLogs();
  const std::vector<double> & y_pos             = GetYposRawLogs();
  const std::vector<double> & z_pos             = GetZposRawLogs();

  std::vector<int> nz(n_intervals);           // nz is the number of layers in each simbox
  for (size_t s = 0; s < interval_simboxes.size(); s++) {
    nz[s] = interval_simboxes[s]->getnz();
    n_layers += nz[s];
  }
  //const std::vector<double> dz_rel = multiple_interval_grid->GetDzRel();  // dz_rel is the dz of simbox i relative to the smallest dz

  std::vector<int> k_offset(n_intervals, 0);
  for (int s = 0; s < n_intervals; s++) {
    std::string interval_name = multiple_interval_grid->GetIntervalName(s);
    int n_lay = static_cast<int>(nz[s]);
    n_layers_adjusted_per_interval.insert(std::pair<std::string, int>(interval_name, n_lay)) ; // number of blocks per interval. Casting double to int cuts the decimals
    if(s < n_intervals-1)
      k_offset[s+1] = k_offset[s]+n_lay;
  }
  //
  // Find the adjusted number of blocks per interval
  //

  //
  // Find first cell in the first simbox where the well log is observed
  //
  int first_I(IMISSING);
  int first_J(IMISSING);
  int first_K(IMISSING);
  for (int m = 0 ; m < n_data ; m++) {
    for (int n = 0 ; n < n_intervals; n++) {
      // The intervals are sometimes overlapping
      if (interval_simboxes[n]->IsPointBetweenVisibleSurfaces(x_pos[m], y_pos[m], z_pos[m])) {
        interval_simboxes[n]->getIndexes(x_pos[m], y_pos[m], z_pos[m], first_I, first_J, first_K);
        if (first_I != IMISSING && first_J != IMISSING && first_K != IMISSING) {
          first_K = static_cast<int>(first_K); // the vertical blocks must be equally spaced for corr estimation
          first_S = n;
          first_M = m;
          break;
        }
      }
    }
    if (first_I != IMISSING && first_J != IMISSING && first_K != IMISSING)
      break;
  }
  if (first_S_ == IMISSING) {
    err_text += "Well " + GetWellName() + " does not pass through any of the inversion interval(s). \n";
    return;
  }
  //
  // Find last cell in last simbox where the well log is observed
  //
  int last_I(IMISSING);
  int last_J(IMISSING);
  int last_K(IMISSING);
  for (int m = n_data - 1 ; m > 0 ; m--) {
    for (int n = n_intervals-1 ; n >= 0; n--) {
      // The intervals are sometimes overlapping
      if (interval_simboxes[n]->IsPointBetweenVisibleSurfaces(x_pos[m], y_pos[m], z_pos[m])) {
        interval_simboxes[n]->getIndexes(x_pos[m], y_pos[m], z_pos[m], last_I, last_J, last_K);
        if (last_I != IMISSING && last_J != IMISSING && last_K != IMISSING) {
          last_K = static_cast<int>(last_K);
          last_S = n;
          last_M = m;
          break;
        }
      }
    }
    if (last_I != IMISSING && last_J != IMISSING && last_K != IMISSING)
      break;
  }
  //
  // Count number of blocks needed for the defined part of the well.
  //
  for (int m = 0 ; m < n_data ; m++) {
    b_ind[m] = IMISSING;
  }
  int new_I, new_J, new_K;
  int old_I = first_I;
  int old_J = first_J;
  int old_K = first_K;
  //int old_simbox = 0;

  int n_defined_blocks = 0;
  b_ind[first_M] = first_K + k_offset[first_S]; // The first defined well log entry (first_M) contributes to
                                                // the first defined block (first_K)

  //
  // The well positions used to be given in float rather than double. Unfortunately, this
  // allowed a well to oscillate between two or more cells, leading to a breakdown of the
  // algorithm below. To remedy for this we introduced the array simbox_ind which records
  // the indices of the simbox cells that are already accounted for, so that these are not
  // enlisted more than one time.
  //
  // ASSUMPTION: all simboxes have same nx and ny
  int * simbox_ind = new int[n_data];                                     // help hack
  const int nx    = interval_simboxes[0]->getnx();                         // help hack
  const int ny    = interval_simboxes[0]->getny();                         // help hack

  simbox_ind[0] = nx*ny*old_K + nx*old_J + old_I;                         // help hack

  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    // Find which interval simbox we are in (this function simply iterates through the simboxes
    // to find out where the x,y,z-coordinates belong
    int simbox_number = multiple_interval_grid->WhichSimbox(x_pos[m], y_pos[m], z_pos[m]);

    interval_simboxes[simbox_number]->getIndexes(x_pos[m], y_pos[m], z_pos[m] , new_I , new_J , new_K);

    int tot_K = new_K+k_offset[simbox_number];

    if (new_I != old_I || new_J != old_J || tot_K != old_K) {

      int  this_ind = nx*ny*tot_K + nx*new_J + new_I;
      bool block_not_listed = true;
      for (int l = 0 ; l < n_defined_blocks ; l++) {
        if (this_ind == simbox_ind[l]) {
          block_not_listed = false;
          break;
        }
      }
      if (block_not_listed) {
        simbox_ind[n_defined_blocks+1] = this_ind;
        old_I = new_I;
        old_J = new_J;
        old_K = tot_K;
        n_defined_blocks++;
      }
    }
    b_ind[m] = first_K + n_defined_blocks;
  }
  n_defined_blocks++;
  //
  // Why we cannot use n_blocks_ = n_layers:
  //
  // When we calculate the background model for each parameter we first
  // estimate a vertical trend in the total volume, and then we interpolate
  // the blocked log intop this trend volume. To avoid sharp contrast we
  // ensure that the blocked log is defined from top to base of the volume.
  // In regions where the log is undefined we generate it by kriging from
  // the rest of the log. Likewise, in regions where there is no blocked
  // log at all because the well was too short, we have to make a virtual
  // well.
  //
  //n_blocks = first_K + n_defined_blocks + (n_layers_ - last_K - 1);
  n_blocks = 0;
  for (int i=0; i<first_S_; i++)
    n_blocks += n_layers_adjusted_per_interval.find(interval_simboxes[i]->GetIntervalName())->second;                   // 1. Add number of blocks from intervals above the first well obs
  for (int i=last_S_+1; i<n_intervals; i++)
    n_blocks += n_layers_adjusted_per_interval.find(interval_simboxes[i]->GetIntervalName())->second;                   // 2. Add number of blocks from intervals below the last well obs
  n_blocks += first_K;                                                                                                  // 3. Add number of layers above the first well observation in the simbox with the first well obs
  n_blocks += n_layers_adjusted_per_interval.find(interval_simboxes[last_S_]->GetIntervalName())->second - last_K;      // 4. Add remaining layers below the last well observation in the simbox with the last well obs
  n_blocks += n_defined_blocks;                                                                                         // 5. Add number of defined blocks between first_K and last_K


  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"first_S_, last_S_          = %d, %d    \n",first_S_,last_S_);
    LogKit::LogFormatted(LogKit::Low,"first_M_, last_M_          = %d, %d    \n",first_M_,last_M_);
    LogKit::LogFormatted(LogKit::Low,"n_layers_                  = %d        \n",n_layers_);
    LogKit::LogFormatted(LogKit::Low,"first_I,first_J,first_K     = %d, %d, %d\n",first_I,first_J,first_K);
    LogKit::LogFormatted(LogKit::Low,"last_I,last_J,last_K        = %d, %d, %d\n",last_I,last_J,last_K);
    LogKit::LogFormatted(LogKit::Low,"n_defined_blocks, n_blocks_ = %d, %d    \n",n_defined_blocks,n_blocks_);
  }
  delete [] simbox_ind;
}

void  BlockedLogsCommon::FindSizeAndBlockPointers(const Simbox                  * const estimation_simbox,
                                                  bool                            visible_only,
                                                  std::vector<int>              & b_ind,
                                                  const int                     & n_layers,
                                                  int                           & first_M,
                                                  int                           & last_M,
                                                  unsigned int                  & n_blocks,
                                                  bool                          & is_inside) const{
  int   nd = static_cast<int>(b_ind.size());
  const std::vector<double> & x_pos = GetXposRawLogs();
  const std::vector<double> & y_pos = GetYposRawLogs();
  const std::vector<double> & z_pos = GetZposRawLogs();

  //
  // Find first cell in Simbox that the well hits
  //
  int first_I(IMISSING);
  int first_J(IMISSING);
  int first_K(IMISSING);
  for (int m = 0 ; m < nd ; m++) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], z_pos[m], first_I, first_J, first_K, visible_only);
    if (first_I != IMISSING && first_J != IMISSING && first_K != IMISSING) {
      first_M = m;
      break;
    }
  }
  if (first_M == IMISSING) {
    //err_text += "Well "+GetWellName()+" was not found within the estimation simbox surrounding the inversion intervals.\n";
    is_inside = false;
    return;
  }
  //
  // Find last cell in Simbox that the well hits
  //
  int last_I(IMISSING);
  int last_J(IMISSING);
  int last_K(IMISSING);
  for (int m = nd - 1 ; m > 0 ; m--) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], z_pos[m], last_I, last_J, last_K, visible_only);
    if (last_I != IMISSING && last_J != IMISSING && last_K != IMISSING) {
      last_M = m;
      break;
    }
  }
  //
  // Count number of blocks needed for the defined part of well.
  //
  for (int m = 0 ; m < nd ; m++) {
    b_ind[m] = IMISSING;
  }
  int new_I, new_J, new_K;
  int old_I = first_I;
  int old_J = first_J;
  int old_K = first_K;

  int n_defined_blocks = 0;
  b_ind[first_M] = first_K; // The first defined well log entry contributes to this block.

  //
  // The well positions used to be given in float rather than double. Unfortunately, this
  // allowed a well to oscillate between two or more cells, leading to a breakdown of the
  // algorithm below. To remedy for this we introduced array simboxInd which records the
  // indices of the simbox cells that are already accounted for, so that these are not
  // enlisted more than one time.
  //
  int * simbox_ind = new int[nd];
  const int nx    = estimation_simbox->getnx();
  const int ny    = estimation_simbox->getny();
  simbox_ind[0] = nx*ny*old_K + nx*old_J + old_I;

  for (int m = first_M + 1 ; m < last_M + 1 ; m++) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], z_pos[m], new_I ,new_J, new_K);

    if (new_I != old_I || new_J != old_J || new_K != old_K) {

      int  this_ind = nx*ny*new_K + nx*new_J + new_I;
      bool block_not_listed = true;
      for (int l = 0 ; l < n_defined_blocks ; l++) {
        if (this_ind == simbox_ind[l]) {
          block_not_listed = false;
          break;
        }
      }
      if (block_not_listed) {
        simbox_ind[n_defined_blocks+1] = this_ind;
        old_I = new_I;
        old_J = new_J;
        old_K = new_K;
        n_defined_blocks++;
      }
    }
    b_ind[m] = first_K + n_defined_blocks;
  }
  n_defined_blocks++;
  //
  // Why we cannot use n_blocks_ = nDefined:
  //
  // When we calculate the background model for each parameter we first
  // estimate a vertical trend in the total volume, anf then we interpolate
  // the blocked log intop this trend volume. To avoid sharp contrast we
  // ensure that the blocked log is defined from top to base of the volume.
  // In regions where the log is undefined we generate it by kriging from
  // the rest of the log. Likewise, in regions where there is no blocked
  // log at all because the well was too short, we have to make a virtual
  // well.
  //
  n_blocks = first_K + n_defined_blocks + (n_layers - last_K - 1);

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"first_M_, last_M_          = %d, %d    \n",first_M,last_M);
    LogKit::LogFormatted(LogKit::Low,"n_layers_                  = %d        \n",n_layers);
    LogKit::LogFormatted(LogKit::Low,"first_I,first_J,first_K     = %d, %d, %d\n",first_I,first_J,first_K);
    LogKit::LogFormatted(LogKit::Low,"last_I,last_J,last_K        = %d, %d, %d\n",last_I,last_J,last_K);
    LogKit::LogFormatted(LogKit::Low,"n_defined_blocks, n_blocks_ = %d, %d    \n",n_defined_blocks,n_blocks);
  }
  delete [] simbox_ind;
}

void    BlockedLogsCommon::CountBlocksWithData(const std::vector<double>                          & x_pos_blocked,
                                               const std::vector<double>                          & y_pos_blocked,
                                               const std::vector<double>                          & z_pos_blocked,
                                               const std::map<std::string, std::vector<double> >  & continuous_logs_blocked,
                                               std::string                                          interval_name,
                                               unsigned int                                         n_blocks,
                                               std::map<std::string, int>                         & n_blocks_with_data,
                                               int                                                & n_blocks_with_data_tot) const
{

  const std::vector<double> vp_log_blocked = continuous_logs_blocked.find("Vp")->second;
  n_blocks_with_data_tot = 0;

  for (size_t i = 0; i < n_blocks; i++) {
    if (vp_log_blocked[i] != RMISSING && x_pos_blocked[i] != RMISSING && y_pos_blocked[i] != RMISSING && z_pos_blocked[i] != RMISSING) {
      n_blocks_with_data.find(interval_name)->second++;
      n_blocks_with_data_tot++;
    }
  }

}

// Counts the number of blocks with data per interval --------------------------
void    BlockedLogsCommon::CountBlocksWithDataPerInterval(const MultiIntervalGrid                            * multiple_interval_grid,
                                                          const std::vector<double>                          & x_pos_blocked,
                                                          const std::vector<double>                          & y_pos_blocked,
                                                          const std::vector<double>                          & z_pos_blocked,
                                                          const std::map<std::string, std::vector<double> >  & continuous_logs_blocked,
                                                          unsigned int                                         n_blocks,
                                                          std::map<std::string, int>                         & n_blocks_with_data,
                                                          int                                                & n_blocks_with_data_tot) const
{
  // Use Vp to test that data is not missing
  const std::vector<double> vp_log_blocked = continuous_logs_blocked.find("Vp")->second;
  n_blocks_with_data_tot = 0;

  for (size_t i = 0; i < n_blocks; i++) {
    if (vp_log_blocked[i] != RMISSING && x_pos_blocked[i] != RMISSING && y_pos_blocked[i] != RMISSING && z_pos_blocked[i] != RMISSING) {
      int s = multiple_interval_grid->WhichSimbox(x_pos_blocked[i], y_pos_blocked[i], z_pos_blocked[i]);
      const std::string interval_name = multiple_interval_grid->GetIntervalName(s);
      n_blocks_with_data.find(interval_name)->second++;
      n_blocks_with_data_tot++;
    }
  }
}

// Function used for blocked logs for correlation estimation in multiple interval setting
//------------------------------------------------------------------------------
void    BlockedLogsCommon::FindBlockIJK(const MultiIntervalGrid          * multiple_interval_grid,
                                        const std::vector<int>           & bInd,
                                        const std::vector<double>        & x_pos_raw_logs,
                                        const std::vector<double>        & y_pos_raw_logs,
                                        const std::vector<double>        & z_pos_raw_logs,
                                        const int                        & first_M,
                                        const int                        & last_M,
                                        const int                        & first_S,
                                        const int                        & last_S,
                                        int                              & first_B,
                                        int                              & last_B,
                                        std::vector<int>                 & n_well_log_obs_in_interval,
                                        std::vector<int>                 & i_pos,
                                        std::vector<int>                 & j_pos,
                                        std::vector<int>                 & k_pos,
                                        std::vector<int>                 & s_pos,
                                        double                           & dz) const
{
  const std::vector<Simbox *> interval_simboxes = multiple_interval_grid->GetIntervalSimboxes();
  //const std::vector<double> dz_rel              = multiple_interval_grid->GetDzRel();
  i_pos.resize(n_blocks_);
  j_pos.resize(n_blocks_);
  k_pos.resize(n_blocks_);
  s_pos.resize(n_blocks_);

  //
  // Preliminary: Count number of well observations in each
  //
  double x,y,z;
  n_well_log_obs_in_interval.resize(interval_simboxes.size(), 0);

  for (size_t m = 0; m < x_pos_raw_logs_.size(); m++) {
    x = x_pos_raw_logs[m];
    y = y_pos_raw_logs[m];
    z = z_pos_raw_logs[m];

    //H Raw logs may be outside simbox
    int simbox_number = multiple_interval_grid->WhichSimbox(x,y,z);
    if (simbox_number > -1)
      n_well_log_obs_in_interval[simbox_number]++;
  }

  std::vector<int> k_offset(interval_simboxes.size(), 0);
  for (size_t s = 0; s < interval_simboxes.size(); s++) {
    int n_lay = static_cast<int>(interval_simboxes[s]->getnz());
    if(s < interval_simboxes.size()-1)
      k_offset[s+1] = k_offset[s]+n_lay;
  }

  //
  // 1. Set IJK for virtual part of well in upper part of the first simbox where it is found
  //
  int b = -1;   // block counter;
  int wl = -1;  // well log counter
  int first_I, first_J, first_K;
  interval_simboxes[first_S]->getIndexes(x_pos_raw_logs[first_M], y_pos_raw_logs[first_M], z_pos_raw_logs[first_M], first_I, first_J, first_K);
  for (int s = 0;  s < first_S_; s++) {
    for (int k = 0; k < interval_simboxes[s]->getnz(); k++) {
      b++;
      s_pos[b] = s;
      i_pos[b] = first_I;
      j_pos[b] = first_J;
      k_pos[b] = k + k_offset[s];
    }
  }
  //
  // IJK in the simbox with the first well observation
  //
  //b = -1;
  for (int k = 0; k < first_K; k++) {
    b++;
    s_pos[b] = first_S_;
    i_pos[b] = first_I;
    j_pos[b] = first_J;
    k_pos[b] = k + k_offset[first_S];
  }

  //
  // 2. Set IJK for the defined part of the well
  //
  b++;
  wl++;
  first_B = b;           // first block that the well log contributes to
  s_pos[b] = first_S_;
  i_pos[b] = first_I;
  j_pos[b] = first_J;
  k_pos[b] = static_cast<int>(first_K);
  int i = 0;
  int j = 0;
  int k = 0;
  int max_m = 0;
  if (first_S_ == last_S_) {
    max_m = last_M_+1;                              // 1. If the last well obs is in interval number first_S_
  }
  else {
    for (int s = 0; s <= first_S_; s++)
      max_m += n_well_log_obs_in_interval[s];      // 2. If the last well obs is in another interval
  }
  // loop over the first simbox where the well is observed
  for (int m = first_M_ + 1; m < max_m; m++) {
    wl++;
    if (bInd[wl] != bInd[wl - 1]) {
      b++;
      interval_simboxes[first_S]->getIndexes(x_pos_raw_logs[m], y_pos_raw_logs[m], z_pos_raw_logs[m], i, j, k);
      s_pos[b] = first_S;
      i_pos[b] = i;
      j_pos[b] = j;
      k_pos[b] = k + k_offset[first_S];
    }
  }
  //first_B_ = static_cast<int>(first_K*dz_rel[first_S_]);
  //last_B_  = b;
  max_m = 0;
  int min_m = 0;
  // the remaining simboxes where the well is observed
  for (int s = first_S+1; s <= last_S ; s++) {
    min_m = 0;
    max_m = 0;
    if (s == last_S_)
      max_m = last_M_+1;                          // 1. If the last well obs is in interval number s
    else {
      for (int t = 0; t < s+1; t++)
        max_m += n_well_log_obs_in_interval[t];  // 2. If the last well obs is in another interval
    }
    for (int t = 0; t < s; t++) {
      min_m += n_well_log_obs_in_interval[t];
    }
    int start = first_M_ + min_m; //Start of this simbox.
    for (int m = start; m < max_m; m++) {
      wl++;
      if (bInd[wl] != bInd[wl - 1]) {
        b++;
        k = 0;
        interval_simboxes[s]->getIndexes(x_pos_raw_logs[m], y_pos_raw_logs[m], z_pos_raw_logs[m], i, j, k);
        s_pos[b] = s;
        i_pos[b] = i;
        j_pos[b] = j;
        k_pos[b] = k + k_offset[s];;
      }
    }
  }
  last_B = b;

  //
  // 3. Set IJK for the virtual part of the well in the lower simboxes
  //

  int last_I, last_J, last_K;

  interval_simboxes[last_S]->getIndexes(x_pos_raw_logs[last_M], y_pos_raw_logs[last_M], z_pos_raw_logs[last_M], last_I, last_J, last_K);
  min_m = last_K+1;
  max_m = interval_simboxes[last_S]->getnz();
  for (k = last_K + 1; k < max_m; k++) {
    b++;
    s_pos[b] = last_S;
    i_pos[b] = last_I;
    j_pos[b] = last_J;
    k_pos[b] = k + k_offset[last_S];;
  }

  for (int s = last_S+1; s < static_cast<int>(interval_simboxes.size()); s++) {
    min_m = 0;
    max_m = interval_simboxes[s]->getnz();
    for (k = min_m; k< max_m; k++) {
      b++;
      s_pos[b] = s;
      i_pos[b] = last_I;
      j_pos[b] = last_J;
      k_pos[b] = k + k_offset[s];
    }
  }

  dz = static_cast<float>(interval_simboxes[0]->getRelThick(i_pos_[0],j_pos_[0])*interval_simboxes[0]->getdz());
  //dz_ = static_cast<float>(interval_simboxes[last_S_].getRelThick(i_pos_[0],j_pos_[0])*interval_simboxes[last_S_].getdz());

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"firstB_, lastB_        = %d, %d    \n",first_B,last_B);
    LogKit::LogFormatted(LogKit::Low,"firstI, firstJ, firstK = %d, %d, %d\n",first_I, first_J, first_K);
    LogKit::LogFormatted(LogKit::Low,"lastI,  lastJ,  lastK  = %d, %d, %d\n",last_I, last_J, last_K);
    for (unsigned int bb = 0; bb < n_blocks_; bb++)
      LogKit::LogFormatted(LogKit::Low,"b=%d   i,j,k=%d,%d,%d\n",bb,i_pos[bb],j_pos[bb],k_pos[bb]);
  }
}

//------------------------------------------------------------------------------
void    BlockedLogsCommon::FindBlockIJK(const Simbox                     * estimation_simbox,
                                        const std::vector<int>           & bInd,
                                        const int                        & first_M,
                                        const int                        & last_M,
                                        int                              & first_B,
                                        int                              & last_B,
                                        std::vector<int>                 & i_pos,
                                        std::vector<int>                 & j_pos,
                                        std::vector<int>                 & k_pos,
                                        double                           & dz) const
{
  i_pos.resize(n_blocks_);
  j_pos.resize(n_blocks_);
  k_pos.resize(n_blocks_);

  const std::vector<double> & x_pos = this->GetXposRawLogs();
  const std::vector<double> & y_pos = this->GetYposRawLogs();
  const std::vector<double> & z_pos = this->GetZposRawLogs();

  //
  // Set IJK for virtual part of well in upper part of simbox
  //
  int b = -1; // block counter;
  int first_I, first_J, first_K;
  estimation_simbox->getIndexes(x_pos[first_M], y_pos[first_M], z_pos[first_M], first_I, first_J, first_K);
  for (int k = 0 ; k < first_K ; k++) {
    b++;
    i_pos[b] = first_I;
    j_pos[b] = first_J;
    k_pos[b] = k;
  }

  //
  // Set IJK for the defined part of the well
  //
  b = first_K;
  i_pos[b] = first_I;
  j_pos[b] = first_J;
  k_pos[b] = first_K;
  int i, j, k;
  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    if (bInd[m] != bInd[m - 1]) {
      b++;
      estimation_simbox->getIndexes(x_pos[m], y_pos[m], z_pos[m], i, j, k);
      i_pos[b] = i;
      j_pos[b] = j;
      k_pos[b] = k;
    }
  }
  first_B = first_K;
  last_B  = b;

  //
  // Set IJK for the virtual part of well in lower part of simbox
  //
  int last_I,  last_J,  last_K;
  estimation_simbox->getIndexes(x_pos[last_M], y_pos[last_M], z_pos[last_M], last_I, last_J, last_K);
  for (k = last_K + 1 ; k < n_layers_ ; k++) {
    b++;
    i_pos[b] = last_I;
    j_pos[b] = last_J;
    k_pos[b] = k;
  }

  dz = static_cast<float>(estimation_simbox->getRelThick(i_pos[0],j_pos[0])*estimation_simbox->getdz());

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"firstB_, lastB_        = %d, %d    \n",first_B,last_B);
    LogKit::LogFormatted(LogKit::Low,"firstI, firstJ, firstK = %d, %d, %d\n",first_I, first_J, first_K);
    LogKit::LogFormatted(LogKit::Low,"lastI,  lastJ,  lastK  = %d, %d, %d\n",last_I, last_J, last_K);
    for (unsigned int bb = 0 ; bb < n_blocks_ ; bb++)
      LogKit::LogFormatted(LogKit::Low,"b=%d   i,j,k=%d,%d,%d\n",bb,i_pos[bb],j_pos[bb],k_pos[bb]);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogsCommon::FindBlockIJK(const StormContGrid             & stormgrid,
                                const std::vector<int>          & b_ind,
                                const int                       & first_M,
                                const int                       & last_M,
                                int                             & first_B,
                                int                             & last_B,
                                std::vector<int>                & i_pos,
                                std::vector<int>                & j_pos,
                                std::vector<int>                & k_pos) const
{
  i_pos.resize(n_blocks_);
  j_pos.resize(n_blocks_);
  k_pos.resize(n_blocks_);


  //int   dummy;
  const std::vector<double> & x = GetXposRawLogs();
  const std::vector<double> & y = GetYposRawLogs();
  const std::vector<double> & z = GetZposRawLogs();

  //
  // Set IJK for virtual part of well in upper part of stormgrid
  //
  int b = -1; // block counter;
  size_t first_I;
  size_t first_J;
  size_t first_K;
  stormgrid.FindIndex(x[first_M], y[first_M], z[first_M], first_I, first_J, first_K);

  for (size_t k = 0; k < first_K; k++) {
    b++;
    i_pos[b] = static_cast<int>(first_I);
    j_pos[b] = static_cast<int>(first_J);
    k_pos[b] = static_cast<int>(k);
  }

  //
  // Set IJK for the defined part of the well
  //
  b = static_cast<int>(first_K);
  i_pos[b] = static_cast<int>(first_I);
  j_pos[b] = static_cast<int>(first_J);
  k_pos[b] = static_cast<int>(first_K);
  size_t i, j, k;
  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    if (b_ind[m] != b_ind[m - 1]) {
      b++;
      stormgrid.FindIndex(x[m], y[m], z[m], i, j, k);
      i_pos[b] = static_cast<int>(i);
      j_pos[b] = static_cast<int>(j);
      k_pos[b] = static_cast<int>(k);
    }
  }
  first_B = static_cast<int>(first_K);
  last_B  = b;

  //
  // Set IJK for the virtual part of well in lower part of simbox
  //
  size_t last_I,  last_J,  last_K;
  stormgrid.FindIndex(x[last_M], y[last_M], z[last_M], last_I, last_J, last_K);

  for (int kk = static_cast<int>(last_K) + 1 ; kk < n_layers_ ; kk++) {
    b++;
    i_pos[b] = static_cast<int>(last_I);
    j_pos[b] = static_cast<int>(last_J);
    k_pos[b] = kk;
  }
}

void BlockedLogsCommon::BlockCoordinateLog(const std::vector<int>    &  b_ind,
                                           const std::vector<double> &  coord,
                                           std::vector<double>       &  blocked_coord) const
{
  //
  // Initialise arrays
  //
  blocked_coord.resize(n_blocks_, 0.0f);
  std::vector<int> count(n_blocks_, 0);

  //
  // Block log
  //
  // 1. Block defined part of log
  for (int m = first_M_ ; m < last_M_ + 1 ; m++) {
    blocked_coord[b_ind[m]] += coord[m];
    count[b_ind[m]]++;
  }
  // 2. Set virtual part of well log to RMISSING
  for (unsigned int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0)
      blocked_coord[l] /= count[l];
    else
      blocked_coord[l]  = RMISSING;
  }

}

//------------------------------------------------------------------------------
void BlockedLogsCommon::BlockFaciesLog(const std::vector<int>          & b_ind,
                                       const std::vector<int>          & well_log,
                                       const std::map<int,std::string> & facies_map,
                                       int                               n_facies,
                                       std::vector<int>                & blocked_log) const
{
  if (well_log.size() > 0) {
    //
    // Set undefined
    //

    std::vector<int> facies_numbers;
    for (std::map<int,std::string>::const_iterator it = facies_map.begin(); it != facies_map.end(); it++) {
      facies_numbers.push_back(it->first);
    }

    blocked_log.resize(n_blocks_);
    for (unsigned int m = 0 ; m < n_blocks_ ; m++)
      blocked_log[m] = IMISSING;

    int   max_allowed_value = 10000;  // Largest allowed value (facies number).
    std::vector<int> count(n_facies);
    std::vector<int> table(max_allowed_value);

    int max_facies_number = -1;
    for (int i = 0 ; i < n_facies ; i++)
      if (facies_numbers[i] > max_facies_number)
        max_facies_number = facies_numbers[i];

    if (max_facies_number > max_allowed_value) {
      LogKit::LogFormatted(LogKit::Error,"ERROR the maximum facies number " + NRLib::ToString(max_facies_number) + " is higher then the allowed limit in Crava " + NRLib::ToString(max_allowed_value) +".\n");
      LogKit::LogFormatted(LogKit::Error,"  Crava has to stop\n");
      exit(1);
    }

    //
    // Set up facies-to-position table.
    //
    // Example: If log values range from 2 to 4 the table looks like
    //
    // table[0] = IMISSING
    // table[1] = IMISSING
    // table[2] =    0
    // table[3] =    1
    // table[4] =    2
    // table[5] = IMISSING
    //    .          .
    //    .          .
    //
    for (int i = 0 ; i < max_allowed_value ; i++)
      table[i] = IMISSING;
    for (int i = 0 ; i < n_facies ; i++)
      table[facies_numbers[i]] = i;

    //
    // Block log
    //
    for (int i = 0 ; i < n_facies ; i++)
      count[i] = 0;
    int value = well_log[first_M_];
    if (value!=IMISSING)
      count[table[value]]++;

    for (int m = first_M_+1 ; m < last_M_ + 1 ; m++) {
      if (b_ind[m] != b_ind[m - 1]) { // bInd[m] is the block number which sample 'm' lies in
        blocked_log[b_ind[m-1]] = FindMostProbable(count, n_facies, b_ind[m-1]);
        for (int i = 0 ; i < n_facies ; i++)
          count[i] = 0;
      }
    value = well_log[m];
    if (value!=IMISSING)
      count[table[value]]++;
    }
    blocked_log[b_ind[last_M_]] = FindMostProbable(count, n_facies, b_ind[last_M_]);

    //
    // NOTE: The blocked log contains internal numbers 0, 1, 2, ... and
    //       are NOT the facies labels.
    //
  }
}

//------------------------------------------------------------------------------
void BlockedLogsCommon::BlockContinuousLog(const std::vector<int>     & b_ind,
                                           const std::vector<double>  & well_log,
                                           std::vector<double>        & blocked_log) const
{
  //
  // Initialise arrays
  //

  blocked_log.resize(n_blocks_, 0.0f);
  std::vector<int> count(n_blocks_, 0);

  //
  // Block log
  //
  for (int m = first_M_ ; m < last_M_ + 1 ; m++) {
    if (well_log[m] != RMISSING && well_log[m] != WELLMISSING) {
      blocked_log[b_ind[m]] += log(well_log[m]); //NBNB-PAL: Flytt denne logaritmen nedover...
      //blocked_log[b_ind[m]] += well_log[m];
      count[b_ind[m]]++;
      //LogKit::LogFormatted(LogKit::Low,"m=%d bInd[m]  log(wellLog[m])  %d  %.5f \n",m,bInd[m],log(wellLog[m]));
    }
  }
  for (unsigned int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0) {
      blocked_log[l] /= count[l];
      //blocked_log[l] = log(blocked_log[l]); //tmp
      //LogKit::LogFormatted(LogKit::Low,"l=%d   count[l]=%d  sum=%.3f  blocked_log[l]=%.4f \n",l,count[l],sum, blocked_log[l]);
    }
    else
      blocked_log[l] = RMISSING;
  }

}

//------------------------------------------------------------------------------
void BlockedLogsCommon::BlockPorosityLog(const std::vector<int>     & b_ind,
                                         const std::vector<double>  & well_log,
                                         std::vector<double>        & blocked_log) const
{
  //
  // Initialise arrays
  //

  blocked_log.resize(n_blocks_, 0.0f);
  std::vector<int> count(n_blocks_, 0);

  //
  // Block log
  //
  for (int m = first_M_ ; m < last_M_ + 1 ; m++) {
    if (well_log[m] != RMISSING && well_log[m] != WELLMISSING) {
      blocked_log[b_ind[m]] += well_log[m];
      count[b_ind[m]]++;
    }
  }

  for (unsigned int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0)
      blocked_log[l] /= count[l];
    else
      blocked_log[l] = RMISSING;
  }
}

//------------------------------------------------------------------------------
void  BlockedLogsCommon::InterpolateContinuousLog(std::vector<double>   & blocked_log,
                                                  int                     start,
                                                  int                     end,
                                                  int                     index,
                                                  float                   rel) const
{
  if (blocked_log[start] != RMISSING && blocked_log[end] != RMISSING && blocked_log[index] == RMISSING)
    blocked_log[index] = rel*blocked_log[end]+(1-rel)*blocked_log[start];
}


void BlockedLogsCommon::SetSeismicGradient(double                            v0,
                                           const NRLib::Grid2D<float>   &    structure_depth_grad_x,
                                           const NRLib::Grid2D<float>   &    structure_depth_grad_y,
                                           const NRLib::Grid2D<float>   &    ref_time_grad_x,
                                           const NRLib::Grid2D<float>   &    ref_time_grad_y,
                                           std::vector<double>          &    x_gradient,
                                           std::vector<double>          &    y_gradient) const
{
  x_gradient.resize(n_blocks_);
  y_gradient.resize(n_blocks_);

  double mp= 2.0/(v0*0.001);
  for (unsigned int k = 0; k < n_blocks_; k++) {
    int i = i_pos_[k];
    int j = j_pos_[k];
    x_gradient[k]= structure_depth_grad_x(i,j)*mp+ref_time_grad_x(i,j);
    y_gradient[k]= structure_depth_grad_y(i,j)*mp+ref_time_grad_y(i,j);
  }
}

void BlockedLogsCommon::SetTimeGradientSettings(float distance, float sigma_m)
{
  lateral_threshold_gradient_ = distance;
  sigma_m_ = sigma_m;
}

void BlockedLogsCommon::FindSeismicGradient(std::vector<SeismicStorage *>     & seismic_data,
                                            const Simbox                * const estimation_simbox,
                                            int                                 n_angles,
                                            std::vector<double>               & x_gradient,
                                            std::vector<double>               & y_gradient,
                                            std::vector<std::vector<double> > & sigma_gradient)
{
  int i, j, l;
  unsigned int k;
  int xEx = 2;
  int yEx = 2;
  int nZx = (2 * xEx + 1);
  int nZy = (2 * yEx + 1);

  x_gradient.resize(n_blocks_);
  y_gradient.resize(n_blocks_);
  std::vector<double> tmp_x_grad(n_blocks_);
  std::vector<double> tmp_y_grad(n_blocks_);
  std::vector<double> q_epsilon(4*n_blocks_);
  std::vector<double> q_epsilon_data(2*n_blocks_);

  std::vector<float> seis_trace;
  std::vector<double> z_shift(int(nZx*nZy*n_blocks_));

  //seismic peak position and characteristics in well
  std::vector<double> z_peak_well;
  std::vector<double> peak_well;
  std::vector<double> b_well;

  //seismic peak position and characterisitics in trace
  std::vector<double> z_peak;
  std::vector<double> peak;
  std::vector<double> b;

  int i0 = i_pos_[0];
  int j0 = j_pos_[0];

  //Check if well needs to change position in order for the whole
  //shift region to be contained in the seismic cube
  //NBNB marita M� testes om det fungerer for forskjellige br�nner
  int nx = estimation_simbox->getnx();
  int ny = estimation_simbox->getny();
  int di_neg = 0; int di_pos = 0; int dj_neg = 0; int dj_pos = 0; //The max replacement in well in x and y direction.
  for (k = 1; k < n_blocks_; k++) {
    di_pos = std::max(i0 - i_pos_[k] , di_pos);
    di_neg = std::min(i0 - i_pos_[k] , di_neg);
    dj_pos = std::max(j0 - j_pos_[k] , dj_pos);
    dj_neg = std::min(j0 - j_pos_[k] , dj_neg);
  }
  i0 = std::max(xEx, i0 + di_neg) - di_neg;
  i0 = std::min(nx - xEx - 1, i0 + di_pos) - di_pos;
  j0 = std::max(yEx, j0 + dj_neg) - dj_neg;
  j0 = std::min(ny - yEx - 1, j0 + dj_pos) - dj_pos;


  int di = i0 - i_pos_[0];
  int dj = j0 - j_pos_[0];
  if (di != 0 || dj != 0) {
    //adjust the well location
    for (k = 0; k < n_blocks_; k++) {
      i_pos_[k] += di;
      j_pos_[k] += dj;
    }
  }

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/traces");
  std::ofstream out(buffer);
  delete [] buffer;

  double dz, ztop, dzW, ztopW;
  for (l = 0; l < n_angles; l++) {
    seismic_data[l]->SetRandomAccess();
    for (j = -yEx; j <= yEx; j++) {
      for (i = -xEx; i <= xEx; i++) {

        seis_trace = seismic_data[l]->GetRealTrace(estimation_simbox, i0, j0); ///H Correct values returned?
        //seis_trace = seisCube[l]->getRealTrace2(i0, j0);

        SmoothTrace(seis_trace);
        //if (j == 0 ) {
        //  for (int s = 0; s < seisTrace.size(); s++)
        //    out << seisTrace[s] << std::endl;
        //}

        dzW =  estimation_simbox->getdz(i0,j0);
        ztopW =  estimation_simbox->getTop(i0,j0);
        FindPeakTrace(seis_trace, z_peak_well, peak_well, b_well, dzW, ztopW);

        //seis_trace = seisCube[l]->getRealTrace2(i0+i, j0+j);
        seis_trace = seismic_data[l]->GetRealTrace(estimation_simbox, i0, j0);
        SmoothTrace(seis_trace);
        if (i==0) {
          for (size_t s = 0; s < seis_trace.size(); s++)
            out << seis_trace[s] << std::endl;
        }

        dz =  estimation_simbox->getdz(i0+i, j0+j);
        ztop =  estimation_simbox->getTop(i0+i, j0+j);
        FindPeakTrace(seis_trace, z_peak, peak, b, dz, ztop);

        PeakMatch(z_peak,peak,b,z_peak_well,peak_well,b_well);//Finds the matching peaks in two traces
        z_shift[(i+2) + (j+2)*nZx] = ComputeShift(z_peak,z_peak_well,z_pos_blocked_[0]);

        for (k = 1; k < n_blocks_; k++) {
          //Check if well changes lateral position
          if ((i_pos_[k]- i_pos_[k-1] == 0) && (j_pos_[k] - j_pos_[k-1] == 0))
            z_shift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = ComputeShift(z_peak,z_peak_well,z_pos_blocked_[k]);
          else {
            //well has changed lateral position and we adapt to the new well position
            //seis_trace = seisCube[l]->getRealTrace2(i_pos_[k],j_pos_[k]);
            seis_trace = seismic_data[l]->GetRealTrace(estimation_simbox, i_pos_[k], j_pos_[k]);
            SmoothTrace(seis_trace);
            dzW = estimation_simbox->getdz(i_pos_[k],j_pos_[k]);
            ztopW = estimation_simbox->getTop(i_pos_[k],j_pos_[k]);
            FindPeakTrace(seis_trace, z_peak_well, peak_well, b_well, dzW, ztopW);

            //seis_trace = seisCube[l]->getRealTrace2(i_pos_[k]+i, j_pos_[k]+j);
            seis_trace = seismic_data[l]->GetRealTrace(estimation_simbox, i_pos_[k]+i, j_pos_[k]+j);
            SmoothTrace(seis_trace);
            dz = estimation_simbox->getdz(i_pos_[k]+i, j_pos_[k]+j);
            ztop = estimation_simbox->getTop(i_pos_[k]+i, j_pos_[k]+j);
            FindPeakTrace(seis_trace, z_peak, peak, b, dz, ztop);

            PeakMatch(z_peak,peak,b,z_peak_well,peak_well,b_well);

            z_shift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = ComputeShift(z_peak, z_peak_well, z_pos_blocked_[k]);
          }
        }
      }
    }
    seismic_data[l]->EndAccess();
    double dx = estimation_simbox->getdx();
    double dy = estimation_simbox->getdy();

    ComputeGradient(q_epsilon, q_epsilon_data, z_shift, nZx, nZx, dx, dy);
  }

  SmoothGradient(x_gradient, y_gradient, q_epsilon, q_epsilon_data, sigma_gradient);
  // NBNB Odd sl�r av estimeringen for � teste om det gir bedre resultat
 /* for (k = 0; k < nBlocks_; k++) {
    xGradient[k]=0.0;
    yGradient[k]=0.0;
  }*/

}

void BlockedLogsCommon::SmoothTrace(std::vector<float> &trace) const
{
  float smoothing_distance = 40; //ms in each direction
  int  L = static_cast<int>(ceil(smoothing_distance/dz_)); //number of lags in the gauss kernel
  double sigma = 10 /dz_; // ms / (ms/cell)

  unsigned int n_trace = static_cast<unsigned int>(trace.size());
  std::vector<float> gk(2*L+1);
  std::vector<float> s_trace(n_trace);
  unsigned int i;
  int j;
  float tmp;
  for (j = -L; j <= L; j++) {
    tmp = static_cast<float>((j*j)/(2*sigma*sigma));
    gk[j+L] = exp(-tmp);
  }

  float N;
  for (i = 0; i < n_trace; i++) {
    N = 0;
    for (j = -L; j <= L; j++) {
      if (i+j >= 0 && i+j < n_trace) {
        s_trace[i] += gk[j+L]*trace[i+j];
        N += gk[j+L];
      }
    }
    s_trace[i] /= N;
  }

  for (i = 0; i < n_trace; i++)
    trace[i] = s_trace[i];
}

void BlockedLogsCommon::FindPeakTrace(std::vector<float>    & trace,
                                      std::vector<double>   & z_peak,
                                      std::vector<double>   & peak,
                                      std::vector<double>   & b,
                                      double                  dz,
                                      double                  z_top) const
{
  int k;
  double x1, x2, x3, y1, y2, y3, y11, y12, y21;
  int N = static_cast<int>(trace.size());
  z_peak.resize(N); peak.resize(N); b.resize(N);

  double a;
  double c;
  int counter = 0;
  for (k = 1; k < N-1; k++) {
    if ((trace[k] >= trace[k-1] && trace[k] > trace[k+1]) || (trace[k] <= trace[k-1] && trace[k] < trace[k+1])) {
      //Data point for interpolation
      x1 = -dz; x2 = 0; x3 = dz;
      y1 = static_cast<double>(trace[k-1]); y2 = static_cast<double>(trace[k]); y3 = static_cast<double>(trace[k+1]);

      //Newton interpolation method
      y11 = (y2 - y1) / (x2 - x1); y12 = (y3 - y2) / (x3 - x2);
      y21 = (y12 - y11) / (x3 - x1);

      //y = ax + bx^2 + c
      c = y1 - y11*x1 + y21*x1*x2;
      a = y11 - y21*x1 - y21*x2;
      b[counter] = y21;
      //zPeak = -a/2b

      z_peak[counter] = - a/(2.0*b[counter]);

      double tmp = a*z_peak[counter] + b[counter]*z_peak[counter]*z_peak[counter] + c;
      peak[counter] = tmp;
      //Transform back to original z-axis
      z_peak[counter] += z_top + k*dz;
      counter++;

    }
  }
  z_peak.resize(counter); b.resize(counter); peak.resize(counter);

}

void BlockedLogsCommon::PeakMatch(std::vector<double>   & z_peak,
                                  std::vector<double>   & peak,
                                  std::vector<double>   & b,
                                  std::vector<double>   & z_peak_w,
                                  std::vector<double>   & peak_w,
                                  std::vector<double>   & bW) const
{
  //This routine matches the peaks from two traces and returns the set of peak positions that matches.
  unsigned int i, j;
  std::vector<double> pW(z_peak.size());
  std::vector<double> p(z_peak.size());

  double diffz;

  double maxdiffz = 5 * dz_; //matching criteria: Peaks must be no longer that 5 cells apart. (marita: input parameter?)
  double diffp = 0.5; //matcing the size of the peaks NBNB-Frode: This should maybe be an input parameter!

  int counter = 0;
  unsigned int lim = 0;
  for (i = 0; i < z_peak_w.size(); i++) {
    for (j = lim; j < z_peak.size(); j++) {
      diffz = fabs(z_peak_w[i] - z_peak[j]);
      if (diffz < maxdiffz) {
        //Check if the peaks point in the same direction
        if ((bW[i] < 0 && b[j] < 0)||(bW[i] >= 0 && b[j] >= 0)) {
          // Check for difference in peak size
          if ((fabs(peak_w[i] - peak[j]))/(fabs(peak_w[i]) + fabs(peak[j])) < diffp) {
            pW[counter] = z_peak_w[i];
            p[counter] =  z_peak[j];
            counter++;
            lim = j + 1;
          }
        }
      }
    }
  }
  z_peak_w.resize(counter); z_peak.resize(counter);
  for (i = 0; i < static_cast<unsigned int>(counter); i++) {
    z_peak_w[i] = pW[i];
    z_peak[i] =  p[i];
  }

}

double BlockedLogsCommon::ComputeShift(std::vector<double>  & z_peak,
                                       std::vector<double>  & z_peak_w,
                                       double                 z0) const
{
  //This routine computes the position of z0 between two peaks in the well and finds the corresponding distance in
  //the other trace. Then zShift is the difference in z between the two.
  unsigned int N = static_cast<unsigned int>(z_peak.size());
  if (N == 0)
    return RMISSING; //The case of no match in the traces
  else {
    unsigned int i;
    int pos = 0;
    double zShift;
    if (z0 < z_peak_w[0])
      zShift = z_peak_w[0] - z_peak[0];
    else if (z0 >= z_peak_w[N-1])
      zShift = z_peak_w[N-1] - z_peak[N-1];
    else {
      for (i = 0; i < N-1; i++) {
        if (z0 >= z_peak_w[i] && z0 < z_peak_w[i+1]) {
          pos = i;
          i = N;
        }
      }
      zShift = z0 - (z_peak[pos] + (z0 - z_peak_w[pos])/(z_peak_w[pos+1]-z_peak_w[pos])*(z_peak[pos+1]-z_peak[pos]));
    }

    return zShift;
  }
}

void BlockedLogsCommon::ComputeGradient(std::vector<double>     & q_epsilon,
                                        std::vector<double>     & q_epsilon_data,
                                        std::vector<double>     & z_shift,
                                        int                       nx,
                                        int                       ny,
                                        double                    dx,
                                        double                    dy) const
{
  //This fit the model zshift(x,y)= beta0 + beta1*x + beta2*y  ==> beta1 is x-gradient and beta2 is y-gradient
  int i, j, k;
  size_t l;
  std::vector<double> Z(3*nx*ny);
  std::vector<double> Y(nx*ny);
  std::vector<double> cov(9);
  std::vector<double> invcov(9);
  std::vector<double> regM(3*nx*ny);

  static bool append = false;

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/gradNoSmooth");
  std::ofstream out;
  if (append) {
    out.open(buffer, std::ios::app|std::ios::out);
  }
  else {
    out.open(buffer);
    append = true;
  }
  delete [] buffer;

  int ndata;
  double data;

  int cy, cz;
  int counter1 = 0;
  int counter2 = 0;
  for (l = 0; l < n_blocks_; l++) {
    cy = 0; cz = 0;
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        data = z_shift[i + j*nx + l*nx*ny];
        if (data != RMISSING) {
          Y[cy] = data;
          Z[cz] = 1.0;
          Z[cz + 1] = (i-(nx-1)/2)*dx;
          Z[cz + 2] = (j-(ny-1)/2)*dy;
          cy++;
          cz += 3;
        }
      }
    }
    Y.resize(cy);
    Z.resize(cz);

    ndata = cy;

    //Compute inverse covariance (ZtZ)^-1 (diagonal matrix for our purpose)
    double tmp;
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) {
        tmp = 0;
        for (k = 0; k < ndata; k++)
          tmp += Z[i + 3*k] * Z[j + 3*k];
        cov[i + 3*j] = tmp;
      }
    }
    double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  +   cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

    if (det != 0) {
      invcov[0] = (cov[4]*cov[8] - cov[5]*cov[7]) / det;
      invcov[1] = (cov[2]*cov[7] - cov[1]*cov[8]) / det;
      invcov[2] = (cov[1]*cov[5] - cov[2]*cov[4]) / det;
      invcov[3] = (cov[5]*cov[6] - cov[3]*cov[8]) / det;
      invcov[4] = (cov[0]*cov[8] - cov[2]*cov[6]) / det;
      invcov[5] = (cov[2]*cov[3] - cov[0]*cov[5]) / det;
      invcov[6] = (cov[3]*cov[7] - cov[4]*cov[6]) / det;
      invcov[7] = (cov[1]*cov[6] - cov[0]*cov[7]) / det;
      invcov[8] = (cov[0]*cov[4] - cov[1]*cov[3]) / det;



      //Compute regression matrix (ZtZ)^-1Zt
      regM.resize(static_cast<unsigned int>(3*ndata));
      for (j = 0; j < 3; j++) {
        for (i = 0; i < ndata; i++) {
          tmp = 0;
          for (k = 0; k < 3; k++)
            tmp += invcov[k + 3*j]*Z[k + 3*i];
          regM[i + j*ndata] = tmp;
        }
      }


      //Compute beta_1(gradientx) og beta_2(gradienty), beta_0 not necessary
      double beta0 = 0;
      double beta1 = 0;
      double beta2 = 0;
      for (j = 0; j < ndata; j++) {
        beta0 += regM[j]*Y[j];
        beta1 += regM[j + ndata]*Y[j];
        beta2 += regM[j + 2*ndata]*Y[j];}
      double sigma2 = 0;
      double sigmatmp;
      for (j = 0; j < ndata; j++) {
        sigmatmp = Y[j] - beta0*Z[j*3] - beta1*Z[j*3 + 1] - beta2*Z[j*3 + 2];
        sigma2 += sigmatmp*sigmatmp;
      }

      double qa = sigma2*invcov[4];
      double qb = sigma2*invcov[8];

      out << beta1 << " " << beta2 << " " << qa << " " << qb << std::endl;

      //cov(beta) = sigma2*(ZtZ)^{-1}
      q_epsilon[counter1] += cov[4]/sigma2;
      q_epsilon[counter1+1] += cov[5]/sigma2;
      q_epsilon[counter1+2] += cov[8]/sigma2;
      q_epsilon[counter1+3] += 0;
      counter1 += 4;

      q_epsilon_data[counter2] += (cov[4]* beta1 + cov[5]*beta2) / sigma2;
      q_epsilon_data[counter2+1] += (cov[5]* beta1 + cov[8]*beta2) / sigma2;
      counter2 += 2;

    }
  }
}

void BlockedLogsCommon::SmoothGradient(std::vector<double>               & x_gradient,
                                       std::vector<double>               & y_gradient,
                                       std::vector<double>               & q_epsilon,
                                       std::vector<double>               & q_epsilon_data,
                                       std::vector<std::vector<double> > & sigma_gradient) const
{
  int i, j;
  int n_beta = n_blocks_* 2;

  NRLib::SymmetricMatrix q_beta_data = NRLib::SymmetricZeroMatrix(n_beta);

  //Set the prior precicion values
  double a, b, c;
  ComputePrecisionMatrix(a,b,c);

  //Initialize the Qm|d matrix
  for (i = 0; i < n_beta - 2; i++) {
    q_beta_data(i,i)   = c + q_epsilon[2*i];
    q_beta_data(i,i+2) = b;
    if (i % 2 == 0)
      q_beta_data(i,i+1) = q_epsilon[2*i+1];
    else
      q_beta_data(i,i+1) = 0;
  }
  //Edge effects
  q_beta_data(0,0             )  += a - c;
  q_beta_data(1,1)               += a - c;
  q_beta_data(n_beta-2,n_beta-2)  = c + q_epsilon[2*(n_beta-2)];
  q_beta_data(n_beta-1,n_beta-1)  = c + q_epsilon[2*(n_beta-1)];
  q_beta_data(n_beta-2,n_beta-1)  = q_epsilon[2*(n_beta-2)+1];

  // Inversion of precision matrix

  NRLib::Matrix I = NRLib::IdentityMatrix(n_beta);

  // Compute the inverse of Qbeta_data
  // First we do cholesky factorization LL^T = covD
  // Solve Q inverse
  NRLib::CholeskySolve(q_beta_data, I);

  //Compute the product (Qbeta_Data)^-1 Qepsilon_data
  NRLib::Vector res(n_beta);

  for (i = 0; i < n_beta; i++) {
    double tmp = 0;
    for (j = 0; j < n_beta; j++)
      tmp += I(i,j)*q_epsilon_data[j];
    res(i) = tmp;
  }

  // Return Sigma_gradient
  sigma_gradient.resize(n_beta);
  std::vector<double> tmp_vec(n_beta);
  for (i = 0; i < n_beta; i++) {
    for (j = 0; j < n_beta; j++)
      tmp_vec[j] = I(i,j);
    sigma_gradient[i] = tmp_vec;
  }


  int counter = 0;
  for (i = 0; i < static_cast<int>(n_blocks_); i++) {
    x_gradient[i] = res(counter);
    y_gradient[i] = res(counter+1);
    counter += 2;
  }

  /*
   char* buffer2 = new char[1000];
   sprintf(buffer2,"%s.txt", "C:/Outputfiles/gradients");
   std::ofstream out2(buffer2);
   for (i = 0; i < nBlocks_; i++) {
       out2 << xGradient[i] << " " << yGradient[i];
     out2 << std::endl;
   }
  delete [] buffer2;
  */
}

void BlockedLogsCommon::ComputePrecisionMatrix(double   & a,
                                               double   & b,
                                               double   & c) const
{
  double minDist = static_cast<double>(lateral_threshold_gradient_);
  double sigma_m =  static_cast<double>(sigma_m_);
  double zQ = 1.645; //default 95% confidence interval

  double K = (dz_ /(minDist*zQ))*(dz_ /(minDist*zQ));
  double alpha = 1 - K/(2*sigma_m);

  double alpha_2 = alpha*alpha;
  double gamma_2 = sigma_m*sigma_m*(1-alpha_2);

  a = 1/sigma_m + alpha_2/gamma_2;
  b = - sqrt(alpha_2)/gamma_2;
  c = (alpha_2 + 1)/gamma_2;
}

void
BlockedLogsCommon::InterpolateTrend(const std::vector<double>   & blocked_log,
                                    double                      * trend) const
{
  for (unsigned int m = 1 ; m < n_blocks_ ; m++) {
    if (abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if (k_pos_[m] < k_pos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(k_pos_[m]-k_pos_[m-1]);
      float t = step_mult;
      for (int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if (trend[j] == RMISSING)
          trend[j] = t*blocked_log[m]+(1-t)*blocked_log[m-1];
        t += step_mult;
      }
    }
  }
}

float BlockedLogsCommon::ComputeElasticImpedance(double         vp,
                                                 double         vs,
                                                 double         rho,
                                                 const float  * coeff) const
{
  // vp, vs, rho are logtransformed
  float ang_imp;

  ang_imp = float(coeff[0]*vp+coeff[1]*vs+coeff[2]*rho );

  return(ang_imp);
}

void BlockedLogsCommon::EstimateCor(fftw_complex * var1_c,
                                    fftw_complex * var2_c,
                                    fftw_complex * ccor_1_2_c,
                                    int            cnzp) const
{
  for (int i=0;i<cnzp;i++) {
    ccor_1_2_c[i].re =  var1_c[i].re*var2_c[i].re + var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }
}

//------------------------------------------------------------------------------
void    BlockedLogsCommon::RemoveMissingLogValues(const NRLib::Well                            * well_data,
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
                                                  std::string                                  & err_text) const
{

  // Get size of data vector including WELLMISSING data points
  unsigned int n_data_with_wellmissing = well_data->GetNData();
  // The number of legal data excluding WELLMISSING is initially set to 0
  n_data = 0;
  // Value from Well_data class
  double OPENWORKS_MISSING = -999.25;

  // Get the continuous and discrete logs from the Well object

  try{
    std::map<std::string,std::vector<double> >  continuous_logs_well = well_data->GetContLog();
    std::map<std::string,std::vector<int> >     discrete_logs_well   = well_data->GetDiscLog();

    std::map<std::string, std::vector<double> >   continuous_logs_raw_logs_temp;
    std::map<std::string, std::vector<int> >      discrete_logs_raw_logs_temp;

    // Pick only the variables that are requested in the constructor --------------------------------------------------------
    // Find the continuous vectors in the wells that are to be blocked
    for (unsigned int i=0; i<cont_logs_to_be_blocked.size(); i++) {
      std::map<std::string,std::vector<double> >::iterator it = continuous_logs_well.find(cont_logs_to_be_blocked[i]);
      // if the well log contains this continuous log
      if (it != continuous_logs_well.end()) {
        continuous_logs_raw_logs_temp.insert(std::pair<std::string, std::vector<double> >(it->first, it->second));
      }
    }
    // Find the discrete vectors in the wells that are to be blocked
    for (unsigned int i=0; i<disc_logs_to_be_blocked.size(); i++) {
      std::map<std::string,std::vector<int> >::iterator it = discrete_logs_well.find(disc_logs_to_be_blocked[i]);
      // if the well log contains this discrete log
      if (it != discrete_logs_well.end()) {
        discrete_logs_raw_logs_temp.insert(std::pair<std::string, std::vector<int> >(it->first, it->second));
      }
    }
    for (std::map<std::string,std::vector<double> >::iterator it = continuous_logs_raw_logs_temp.begin(); it!=continuous_logs_raw_logs_temp.end(); it++) {
      std::vector<double> temp_vector;
      continuous_logs_raw_logs.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector));
    }

    for (std::map<std::string,std::vector<int> >::iterator it = discrete_logs_raw_logs_temp.begin(); it!=discrete_logs_raw_logs_temp.end(); it++) {
      std::vector<int> temp_vector;
      discrete_logs_raw_logs.insert(std::pair<std::string, std::vector<int> >(it->first, temp_vector));
    }

    // Remove WELLMISSING data from wells
    if (continuous_logs_well.find("Z_pos") != continuous_logs_well.end()) {

      for (unsigned int i = 0; i < n_data_with_wellmissing; i++) {
        // TWT/Z_pos is the variable we are testing for WELLMISSING values
        double dummy = continuous_logs_well.find("Z_pos")->second[i];
        if (dummy != WELLMISSING && dummy != OPENWORKS_MISSING) {

          x_pos_raw_logs.push_back(continuous_logs_well.find("X_pos")->second[i]);
          y_pos_raw_logs.push_back(continuous_logs_well.find("Y_pos")->second[i]);
          z_pos_raw_logs.push_back(continuous_logs_well.find("Z_pos")->second[i]);

          if (facies_log_defined_)
            facies_raw_logs.push_back(discrete_logs_well.find("Facies")->second[i]);

          // Loop over continuous variables and push back this element
          for (std::map<std::string,std::vector<double> >::iterator it = continuous_logs_raw_logs_temp.begin(); it!=continuous_logs_raw_logs_temp.end(); it++) {
            continuous_logs_raw_logs.find(it->first)->second.push_back(it->second[i]);
          }

          // Loop over discrete variables and push back this element
          for (std::map<std::string,std::vector<int> >::iterator it = discrete_logs_raw_logs_temp.begin(); it!=discrete_logs_raw_logs_temp.end(); it++) {
            discrete_logs_raw_logs.find(it->first)->second.push_back(it->second[i]);
          }

          n_data++;
        }
      }
    }
    else {
      failed = true;
      err_text += "\nError: TWT is not in part of logs in well " + well_data->GetWellName();
    }
  }catch(NRLib::Exception & e) {
    err_text += "Error: "+std::string(e.what());
    failed = true;
  }

}

void BlockedLogsCommon::FindOptimalWellLocation(std::vector<SeismicStorage *> & seismic_data,
                                                const Simbox                  * estimation_simbox,
                                                const Simbox                  & inversion_simbox,
                                                const NRLib::Matrix           & refl_matrix,
                                                int                             n_angles,
                                                const std::vector<float>      & angle_weight,
                                                float                           max_shift,
                                                int                             i_max_offset,
                                                int                             j_max_offset,
                                                const std::vector<Surface *>  & limits,
                                                int                           & i_move,
                                                int                           & j_move,
                                                float                         & k_move) const
{
  int   polarity;
  int   i,j,k,l,m;
  int   start, length;
  float sum;
  float shift_F;
  float max_tot;
  float f1,f2,f3;

  int nx              = estimation_simbox->getnx();
  int ny              = estimation_simbox->getny();
  int nz              = estimation_simbox->getnz();
  int nzp             = 2*nz; //Should not hurt, only single traces here.
  int cnzp            = nzp/2+1;
  int rnzp            = 2*cnzp;
  int i_tot_offset    = 2*i_max_offset+1;
  int j_tot_offset    = 2*j_max_offset+1;
  int polarityMax     = 0;
  float shift         = 0.0f;
  float max_value_tot = 0;
  float total_weight  = 0;
  float dz            = static_cast<float>(estimation_simbox->getdz());

  std::vector<double> seis_log(n_blocks_);
  std::vector<double> seis_data(n_layers_);
  std::vector<double> vp_vert(n_layers_);
  std::vector<double> vs_vert(n_layers_);
  std::vector<double> rho_vert(n_layers_);

  std::vector<int>   i_offset(i_tot_offset);
  std::vector<int>   j_offset(j_tot_offset);
  std::vector<int>   shift_I(n_angles);
  std::vector<int>   shift_I_max(n_angles);
  std::vector<float> max_value(n_angles);
  std::vector<float> max_value_max(n_angles);

  fftw_real    ** cpp_r               = new fftw_real*[n_angles];
  fftw_complex ** cpp_c               = reinterpret_cast<fftw_complex**>(cpp_r);

  fftw_real    ** cor_cpp_r           = new fftw_real*[n_angles];
  fftw_complex ** cor_cpp_c           = reinterpret_cast<fftw_complex**>(cor_cpp_r);

  fftw_real    ** seis_r              = new fftw_real*[n_angles];
  fftw_complex ** seis_c              = reinterpret_cast<fftw_complex**>(seis_r);

  fftw_real    ** ccor_seis_cpp_r     = new fftw_real*[n_angles];
  fftw_complex ** ccor_seis_cpp_c     = reinterpret_cast<fftw_complex**>(ccor_seis_cpp_r);

  fftw_real    ** ccor_seis_cpp_Max_r = new fftw_real*[n_angles];

  for (i=0; i < n_angles; i++) {
    max_value_max[i] = 0.0f;
    shift_I_max[i]   = 0;
  }

  // make offset vectors
  for (i=-i_max_offset; i < i_max_offset+1; i++) {
    i_offset[i+i_max_offset]=i;
  }
  for (j=-j_max_offset; j < j_max_offset+1; j++) {
    j_offset[j+j_max_offset]=j;
  }

  GetVerticalTrendLimited(GetVpBlocked(), vp_vert, limits);
  GetVerticalTrendLimited(GetVsBlocked(), vs_vert, limits);
  GetVerticalTrendLimited(GetRhoBlocked(), rho_vert, limits);

  std::vector<bool> has_data(n_layers_);
  for (i = 0 ; i < n_layers_ ; i++) {
    has_data[i] = vp_vert[i] != RMISSING && vs_vert[i] != RMISSING && rho_vert[i] != RMISSING;
  }
  FindContinuousPartOfData(has_data, n_layers_, start, length);

  for ( j=0; j<n_angles; j++ ) {
    seis_r[j]              = new fftw_real[rnzp];
    cpp_r[j]               = new fftw_real[rnzp];
    cor_cpp_r[j]           = new fftw_real[rnzp];
    ccor_seis_cpp_r[j]     = new fftw_real[rnzp];
    ccor_seis_cpp_Max_r[j] = new fftw_real[rnzp];
  }

  // Calculate reflection coefficients
  for ( j=0; j<n_angles; j++ ) {
    for (i=0; i<rnzp; i++) {
      cpp_r[j][i] = 0;
    }
    float * refl_coefficients = new float[3];
    refl_coefficients[0] = static_cast<float>(refl_matrix(j,0));
    refl_coefficients[1] = static_cast<float>(refl_matrix(j,1));
    refl_coefficients[2] = static_cast<float>(refl_matrix(j,2));
    FillInCpp(refl_coefficients,start,length,cpp_r[j],nzp);
    Utils::fft(cpp_r[j],cpp_c[j],nzp);
    EstimateCor(cpp_c[j],cpp_c[j],cor_cpp_c[j],cnzp);
    Utils::fftInv(cor_cpp_c[j],cor_cpp_r[j],nzp);
    delete [] refl_coefficients;
  }

  std::vector<NRLib::Grid<float> > seis_cube_small(n_angles,NRLib::Grid<float> (i_tot_offset,j_tot_offset,n_blocks_));

  for (j = 0 ; j < n_angles ; j++)
  {
    //seismic_data[j]->setAccessMode(FFTGrid::RANDOMACCESS);
    for (k = 0; k < i_tot_offset; k++)
    {
      for (l = 0; l < j_tot_offset; l++)
      {
        GetBlockedGrid(seismic_data[j], estimation_simbox, seis_log, i_offset[k], j_offset[l]);
        for (m = 0; m < static_cast<int>(n_blocks_); m++)
        {
          seis_cube_small[j](k, l, m) = static_cast<float>(seis_log[m]);
        }
      }
    }
    //seis_cube[j]->endAccess();
  }

  // Loop through possible well locations
  for (k=0; k<i_tot_offset; k++) {
    int i_index = i_pos_[0]+i_offset[k];
    if (i_index<0 || i_index>nx-1) //Check if position is within seismic range
      continue;

    for (l=0; l<j_tot_offset; l++) {
      int j_index = j_pos_[0]+j_offset[l];
      if (j_index<0 || j_index>ny-1) //Check if position is within seismic range
        continue;
      else {  //Check if position is within inversion simbox.
        double xp, yp, zp;
        estimation_simbox->getCoord(i_index, j_index, 0, xp, yp, zp);
        if (inversion_simbox.isInside(xp, yp) == false)
          continue;
      }

      for (j = 0; j < n_angles; j++) {

        for (m=0; m<static_cast<int>(n_blocks_); m++)
          seis_log[m] = seis_cube_small[j](k,l,m);

        GetVerticalTrend(seis_log, seis_data);
        FillInSeismic(seis_data,start,length,seis_r[j],nzp);

        Utils::fft(seis_r[j],seis_c[j],nzp);
        EstimateCor(seis_c[j],cpp_c[j],ccor_seis_cpp_c[j],cnzp);
        Utils::fftInv(ccor_seis_cpp_c[j],ccor_seis_cpp_r[j],nzp);
      }

      // if the sum from -max_shift to max_shift ms is
      // positive then polarity is positive
      dz = static_cast<float>(estimation_simbox->getRelThick(i_index,j_index)*estimation_simbox->getdz());
      sum = 0;
      for ( j=0; j<n_angles; j++ ) {
        if (angle_weight[j] > 0) {
          for (i=0;i<ceil(max_shift/dz);i++)//zero included
            sum+=ccor_seis_cpp_r[j][i];
          for (i=0;i<floor(max_shift/dz);i++)
            sum+=ccor_seis_cpp_r[j][nzp-i-1];
        }
      }
      polarity=-1;
      if (sum > 0)
        polarity=1;

      // Find maximum correlation and corresponding shift for each angle
      max_tot = 0.0;
      for ( j=0; j<n_angles; j++ ) {
        if (angle_weight[j]>0) {
          max_value[j] = 0.0f;
          shift_I[j]=0;
          for (i=0;i<ceil(max_shift/dz);i++) {
            if (ccor_seis_cpp_r[j][i]*polarity > max_value[j]) {
              max_value[j] = ccor_seis_cpp_r[j][i]*polarity;
              shift_I[j] = i;
            }
          }
          for (i=0;i<floor(max_shift/dz);i++) {
            if (ccor_seis_cpp_r[j][nzp-1-i]*polarity > max_value[j]) {
              max_value[j] = ccor_seis_cpp_r[j][nzp-1-i]*polarity;
              shift_I[j] = -1-i;
            }
          }
          max_tot += angle_weight[j]*max_value[j]; //Find weighted total maximum correlation
        }
      }

      if (max_tot > max_value_tot) {
        max_value_tot = max_tot;
        polarityMax = polarity;
        i_move       = i_offset[k];
        j_move       = j_offset[l];
        for (m=0; m<n_angles; m++) {
          shift_I_max[m]   = shift_I[m];
          max_value_max[m] = max_value[m];
          for (i=0;i<rnzp;i++)
            ccor_seis_cpp_Max_r[m][i] = ccor_seis_cpp_r[m][i];
        }
      }
    }
  }

  for (i=0; i<n_angles; i++) {
    shift_I[i] = shift_I_max[i];
    max_value[i] = max_value_max[i];
    for (j=0;j<rnzp;j++)
      ccor_seis_cpp_r[i][j] = ccor_seis_cpp_Max_r[i][j];
  }
  polarity = polarityMax;

  // Find kMove in optimal location
  for (j=0; j<n_angles; j++) {
    if (angle_weight[j]>0) {
      if (shift_I[j] < 0) {
        if (ccor_seis_cpp_r[j][nzp+shift_I[j]-1]*polarity < max_value[j]) //then local max
        {
          f1 = ccor_seis_cpp_r[j][nzp+shift_I[j]-1];
          f2 = ccor_seis_cpp_r[j][nzp+shift_I[j]];
          int ind3;
          if (shift_I[j]==-1)
            ind3 = 0;
          else
            ind3=nzp+shift_I[j]+1;
          f3 = ccor_seis_cpp_r[j][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shift_F=shift_I[j]+x0;
        }
        else  // do as good as we can
          shift_F=float(shift_I[j]);
      }
      else //positive or zero shift
      {
        if (ccor_seis_cpp_r[j][shift_I[j]+1]*polarity < max_value[j]) //then local max
        {
          f3 = ccor_seis_cpp_r[j][shift_I[j]+1];
          f2 = ccor_seis_cpp_r[j][shift_I[j]];
          int ind1;
          if (shift_I[j]==0)
            ind1 = nzp-1;
          else
            ind1=shift_I[j]-1;
          f1 = ccor_seis_cpp_r[j][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shift_F=shift_I[j]+x0;
        }
        else  // do as good as we can
          shift_F=float(shift_I[j]);
      }
      shift += angle_weight[j]*shift_F*dz;//weigthing shift according to well_weight
      total_weight += angle_weight[j];
    }
  }

  shift/=total_weight;
  k_move = shift;

  for ( j=0; j<n_angles; j++ ) {
    delete [] ccor_seis_cpp_Max_r[j];
    delete [] ccor_seis_cpp_r[j];
    delete [] cor_cpp_r[j];
    delete [] seis_r[j];
    delete [] cpp_r[j];
  }


  delete [] ccor_seis_cpp_Max_r;
  delete [] ccor_seis_cpp_r;
  delete [] cor_cpp_r;
  delete [] seis_r;
  delete [] cpp_r;
}

void BlockedLogsCommon::GetVerticalTrendLimited(const std::vector<double>          & log,
                                                std::vector<double>                & trend,
                                                const std::vector<Surface *>       & limits) const{
  if (log.size() > 0 && trend.size() > 0) {
    std::vector<int> count(n_layers_);
    for (int k = 0 ; k < n_layers_ ; k++) {
      trend[k] = 0.0;
      count[k] = 0;
    }
    for (int m = 0 ; m < static_cast<int>(n_blocks_); m++) {
      if (log[m] != RMISSING) {
        if (limits.size() == 0 ||
           (limits[0]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) <= z_pos_blocked_[m] &&
            limits[1]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) >= z_pos_blocked_[m])) {
          trend[k_pos_[m]] += log[m];
          count[k_pos_[m]]++;
        }
      }
    }
    for (int k = 0 ; k < n_layers_; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
    if (interpolate_)
      InterpolateTrend(log,trend,limits);

  }
  else {
    if (log.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrendLimited(): Trying to use an undefined log\n");
    if (trend.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrendLimited(): Trying to use an undefined trend\n");
    exit(1);
  }
}

void
BlockedLogsCommon::InterpolateTrend(const double  * blocked_log,
                                    double        * trend) const
{
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if (abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if (k_pos_[m] < k_pos_[m-1])
        delta = -1;
      double step_mult = static_cast<double>(delta)/static_cast<double>(k_pos_[m]-k_pos_[m-1]);
      double t = step_mult;
      for (int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if (trend[j] == RMISSING)
          trend[j] = t*blocked_log[m]+(1-t)*blocked_log[m-1];
        t += step_mult;
      }
    }
  }
}

void BlockedLogsCommon::InterpolateTrend(const std::vector<double>    & blocked_log,
                                         std::vector<double>          & trend) const
{
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if (abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if (k_pos_[m] < k_pos_[m-1])
        delta = -1;
      for (int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if (trend[j] == RMISSING)
          trend[j] = blocked_log[m-1];
      }
    }
  }
}

void  BlockedLogsCommon::InterpolateTrend(const std::vector<double>      & blocked_log,
                                          std::vector<double>            & trend,
                                          const std::vector<Surface *>   & limits) const
{
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if (abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if (k_pos_[m] < k_pos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(k_pos_[m]-k_pos_[m-1]);
      float t = step_mult;
      for (int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if (trend[j] == RMISSING) {
          if (limits.size() == 0 ||
             (limits[0]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) <= z_pos_blocked_[m] &&
              limits[1]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) >= z_pos_blocked_[m])) {
            trend[j] = t*blocked_log[m]+(1-t)*blocked_log[m-1];
          }
        }
        t += step_mult;
      }
    }
  }
}

void BlockedLogsCommon::FindContinuousPartOfData(const std::vector<bool> & hasData,
                                                 int                       nz,
                                                 int                     & start,
                                                 int                     & length) const
{
  int  i;
  int  l_pice            =  0;
  int  length_max_pice    = -1;
  int  start_longest_pice =  0;
  bool previous_had_data  = false;

  for (i = 0; i < nz; i++) {
    if (hasData[i]) {
      if (! previous_had_data)
        l_pice=1;
      else
        l_pice++;
      previous_had_data = true;
    }
    else {
      if (previous_had_data) {
        if (length_max_pice < l_pice) {
          length_max_pice  = l_pice;
          start_longest_pice = i-l_pice;
        }
      }
      previous_had_data=false;
    }
  }

  if (previous_had_data) {
    if (length_max_pice < l_pice) {
      length_max_pice  = l_pice;
      start_longest_pice = i-l_pice;
    }
  }

  start  = start_longest_pice;
  length = length_max_pice;

  if (length == -1)
    length = 0;
}

void BlockedLogsCommon::FillInCpp(const float * coeff,
                                  int           start,
                                  int           length,
                                  fftw_real   * cpp_r,
                                  int           nzp) const
{
  int i;

  for (i=0; i < nzp; i++)
    cpp_r[i]=0;

  std::vector<double> vp_vert(n_layers_);
  std::vector<double> vs_vert(n_layers_);
  std::vector<double> rho_vert(n_layers_);

  GetVerticalTrend(GetVpBlocked(), vp_vert);
  GetVerticalTrend(GetVsBlocked(), vs_vert);
  GetVerticalTrend(GetRhoBlocked(), rho_vert);

  //Make reflection coefficients consistent with seismic data indexes. Seismic data at index i is the response at top of cell i, so
  //the reflection coefficients in a cell must also be associated with the top of the cell.
  for (i = start+1; i < start+length; i++) {
    double ei1 = ComputeElasticImpedance(vp_vert[i-1],   static_cast<float>(vs_vert[i-1]),  static_cast<float>(rho_vert[i-1]),   coeff);
    double ei2 = ComputeElasticImpedance(vp_vert[i], static_cast<float>(vs_vert[i]),static_cast<float>(rho_vert[i]), coeff);
    cpp_r[i] =  static_cast<fftw_real>(ei2-ei1);
  }

}

double BlockedLogsCommon::ComputeElasticImpedance(double        vp,
                                                  float         vs,
                                                  float         rho,
                                                  const float * coeff) const
{
  // vp, vs, rho are logtransformed
  double ang_imp;

  ang_imp = coeff[0]*vp+coeff[1]*vs+coeff[2]*rho;

  return ang_imp;
}

void BlockedLogsCommon::GetVerticalTrend(const std::vector<double>  & log,
                                         std::vector<double>        & trend) const
{
  if (log.size() > 0 && trend.size() > 0) {
    std::vector<double> count(n_layers_);
    for (int k = 0 ; k < n_layers_ ; k++) {
      trend[k] = 0.0;
      count[k] = 0;
    }
    for (int m = 0 ; m < static_cast<int>(n_blocks_) ; m++) {
      if (log[m] != RMISSING) {
        trend[k_pos_[m]] += log[m];
        count[k_pos_[m]]++;
      }
    }
    for (int k = 0 ; k < n_layers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
    if (interpolate_ == true)
      InterpolateTrend(log, trend);

  }
  else {
    if (log.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrend(): Trying to use an undefined log\n");
    if (trend.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrend(): Trying to use an undefined trend\n");
    exit(1);
  }
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::GetVerticalTrend(const int         * blocked_log,
                                         std::vector<int>  & trend) const
{
  if (blocked_log != NULL && trend.size() > 0) {

    int * count = new int[GetNFacies()];
    for (int k = 0 ; k < GetNFacies() ; k++) {
      count[k] = 0;
    }

    for (int k = 0 ; k < n_layers_ ; k++) {
      for (int i = 0 ; i < GetNFacies() ; i++)
        count[i] = 0;
      for (size_t m = 0 ; m < n_blocks_ ; m++) {
        if (k_pos_[m] == k) {
          if (blocked_log[m] != IMISSING) {
            count[blocked_log[m]]++;        // Count the number of times a facies occurs in layer 'k'
          }
        }
      }
      trend[k] = FindMostProbable(count, GetNFacies(), k);
    }
    if (interpolate_ == true)
      InterpolateTrend(blocked_log, trend);

    delete [] count;
  }
  else {
    if (blocked_log == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }
}

//--------------------------------------------------------------------------------------
int BlockedLogsCommon::FindMostProbable(const int * count,
                                        int         n_facies,
                                        int         block_index) const
{
  int  max_index     = IMISSING;
  int  max_count     = 0;
  bool inconclusive = false;

  for (int i=0 ; i < n_facies ; i++ ) {
    if (count[i] > 0 && count[i] > max_count) {
      max_count     = count[i];
      max_index     = i;
      inconclusive = false;
    }
    else if (count[i] > 0 && count[i] == max_count) {
      inconclusive = true;
    }
  }

  if (inconclusive) {
    std::vector<int> equal;
    for (int i=0 ; i < n_facies ; i++ ) {
      if (count[i] == max_count) {
        equal.push_back(i);
      }
    }
    int j = (block_index + 1) % equal.size();
    max_index = equal[j];
  }

  return (max_index);
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::InterpolateTrend(const int          * blocked_log,
                                         std::vector<int>   & trend) const
{
  for (size_t m = 1 ; m < n_blocks_ ; m++) {
    if (abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if (k_pos_[m] < k_pos_[m-1])
        delta = -1;
      for (int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if (trend[j] == RMISSING)
          trend[j] = blocked_log[m-1];
      }
    }
  }
}


//--------------------------------------------------------------------------------------
void BlockedLogsCommon::GetBlockedGrid(SeismicStorage         * grid,
                                       const Simbox           * estimation_simbox,
                                       std::vector<double>    & blocked_log,
                                       int                      i_offset,
                                       int                      j_offset) const
{
  grid->SetRandomAccess(); //Needed if grid is FFTGrid.
  for (int m = 0; m < static_cast<float>(n_blocks_); m++) {
    //LogKit::LogFormatted(LogKit::Low,"m=%d  ipos_[m], jpos_[m], kpos_[m] = %d %d %d\n",m,ipos_[m], jpos_[m], kpos_[m]);
    blocked_log[m] = grid->GetRealTraceValue(estimation_simbox, i_pos_[m]+i_offset, j_pos_[m]+j_offset, k_pos_[m]);
  }
  grid->EndAccess();      //Needed if grid is FFTGrid.
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::GetBlockedGrid(const NRLib::Grid<float> * grid,
                                       std::vector<double>      & blocked_log,
                                       int                        i_offset,
                                       int                        j_offset) const
{
  for (size_t m = 0 ; m < n_blocks_ ; m++) {
    //LogKit::LogFormatted(LogKit::Low,"m=%d  ipos_[m], jpos_[m], kpos_[m] = %d %d %d\n",m,ipos_[m], jpos_[m], kpos_[m]);
    blocked_log[m] = static_cast<float>(grid->GetValue(i_pos_[m]+i_offset, j_pos_[m]+j_offset, k_pos_[m]));
  }
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::FillInSeismic(std::vector<double>   & seismic_data,
                                      int                     start,
                                      int                     length,
                                      fftw_real             * seis_r,
                                      int                     nzp,
                                      bool                    top_value) const
{
  int i;
  for (i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for (i=start; i<start+length; i++)
    seis_r[i] = static_cast<fftw_real>(seismic_data[i]);

  if (top_value == true)
    Utils::ShiftTrace(&(seis_r[start]),length, false);
}

//--------------------------------------------------------------------------------------
void  BlockedLogsCommon::SetLogFromVerticalTrend(const std::vector<double>                    & vertical_trend,
                                                 std::map<std::string, std::vector<double> >  & cont_logs_seismic_resolution,
                                                 std::vector<std::vector<double> >            & actual_synt_seismic_data,
                                                 std::vector<std::vector<double> >            & well_synt_seismic_data,
                                                 std::string                                    type,
                                                 int                                            i_angle) const
{

  if (type != "WELL_SYNTHETIC_SEISMIC")
  {
    std::vector<double> blocked_log(n_blocks_);

    if (type == "ACTUAL_SYNTHETIC_SEISMIC") {
      std::vector<fftw_real> trace(vertical_trend.size());
      for (size_t i=0;i<vertical_trend.size();i++)
        trace[i] = static_cast<float>(vertical_trend[i]);
      Utils::ShiftTrace(&(trace[0]), trace.size());
      std::vector<double> v_trend(vertical_trend.size());
      for (size_t i=0;i<vertical_trend.size();i++)
        v_trend[i] = static_cast<double>(trace[i]);
      for (size_t i=0 ; i<n_blocks_ ; i++)
        blocked_log[i] = v_trend[k_pos_[i]];
      if (actual_synt_seismic_data_.size() == 0)
        actual_synt_seismic_data.resize(n_angles_); // nAngles is set along with real_seismic_data_
      actual_synt_seismic_data[i_angle] = blocked_log;
    }
    else {
      for (size_t i=0 ; i<n_blocks_ ; i++)
        blocked_log[i] = vertical_trend[k_pos_[i]];
      if (type == "VP_SEISMIC_RESOLUTION")
        cont_logs_seismic_resolution.insert(std::pair<std::string, std::vector<double> >("Vp", blocked_log));
      else if (type == "VS_SEISMIC_RESOLUTION")
        cont_logs_seismic_resolution.insert(std::pair<std::string, std::vector<double> >("Vs", blocked_log));
      else if (type == "RHO_SEISMIC_RESOLUTION")
        cont_logs_seismic_resolution.insert(std::pair<std::string, std::vector<double> >("Rho", blocked_log));
      else {
        LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type+
                             "\" in BlockedLogs::setLogFromVerticalTrend()\n");
        exit(1);
      }
    }
  }
  else if (type == "WELL_SYNTHETIC_SEISMIC") {
    if (well_synt_seismic_data.size() == 0) {
      well_synt_seismic_data.resize(n_angles_);  //n_angles_ set in CommonData before wavelet is estimated
      for (int i = 0; i < n_angles_; i++) {
        well_synt_seismic_data[i].resize(n_blocks_);
        for (int j = 0; j < static_cast<int>(n_blocks_); j++)
          well_synt_seismic_data[i][j] = RMISSING; //Declare in case the wavelet is not estimated for all angles
      }
    }
    std::vector<fftw_real> trace(vertical_trend.size());
    for (size_t i=0;i<vertical_trend.size();i++)
      trace[i] = static_cast<float>(vertical_trend[i]);
    Utils::ShiftTrace(&(trace[0]), trace.size());
    std::vector<double> v_trend(vertical_trend.size());
    for (size_t i=0;i<vertical_trend.size();i++)
      v_trend[i] = static_cast<double>(trace[i]);
    for (size_t i=0 ; i<n_blocks_ ; i++)
      well_synt_seismic_data[i_angle][i] = v_trend[k_pos_[i]];
  }
}

//------------------------------------------------------------------------------
int BlockedLogsCommon::FindMostProbable(const std::vector<int>  & count,
                                        int                       n_facies,
                                        int                       block_index) const
{
  int  maxIndex     = IMISSING;
  int  maxCount     = 0;
  bool inconclusive = false;

  for (int i=0 ; i < n_facies ; i++ ) {
    if (count[i] > 0 && count[i] > maxCount) {
      maxCount     = count[i];
      maxIndex     = i;
      inconclusive = false;
    }
    else if (count[i] > 0 && count[i] == maxCount) {
      inconclusive = true;
    }
  }

  if (inconclusive) {
    std::vector<int> equal;
    for (int i=0 ; i < n_facies ; i++ ) {
      if (count[i] == maxCount) {
        equal.push_back(i);
      }
    }
    int j = (block_index + 1) % equal.size();
    maxIndex = equal[j];
  }

  return (maxIndex);
}

//--------------------------------------------------------------------------------------
void  BlockedLogsCommon::SetLogFromGrid(FFTGrid    * grid,
                                        int          i_angle,
                                        int          n_angles,
                                        std::string  type)
{
  std::vector<double> blocked_log(n_blocks_);
  //int n_facies = facies_blocked_.size();
  grid->setAccessMode(FFTGrid::RANDOMACCESS);
  for (size_t m = 0 ; m < n_blocks_ ; m++) {
    blocked_log[m] = grid->getRealValue(i_pos_[m], j_pos_[m], k_pos_[m]);
  }
  grid->endAccess();

  if (n_angles_ == 0)
    n_angles_ = n_angles;

  if (type == "SEISMIC_DATA") {
    //Real seismic data is shifted up by half a cell to match inversion. Shift back
    Utils::ShiftTrace(blocked_log,false);
    real_seismic_data_.insert(std::pair<int, std::vector<double> >(i_angle, blocked_log));
  }
  else if (type == "FACIES_PROB") {
    facies_prob_.insert(std::pair<int, std::vector<double> >(i_angle, blocked_log));
  }
  else if (type == "VP_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Vp", blocked_log));
  }
  else if (type == "VS_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Vs", blocked_log));
  }
  else if (type == "RHO_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Rho", blocked_log));
  }
  else {
    LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type
                         +"\" in BlockedLogs::setLogFromGrid()\n");
    exit(1);
  }
}

//--------------------------------------------------------------------------------------
void  BlockedLogsCommon::SetLogFromGrid(SegY         * segy,
                                        const Simbox & simbox,
                                        int            i_angle,
                                        int            n_angles,
                                        std::string    type)
{
  std::vector<double> blocked_log(n_blocks_);

  double x     = 0.0;
  double y     = 0.0;
  double z     = 0.0;
  double value = 0.0;

  double shift = 0.0;
  if (type == "SEISMIC_DATA")
    shift = 0.5*segy->GetDz();

  for (size_t m = 0 ; m < n_blocks_ ; m++) {

    simbox.getCoord(i_pos_[m], j_pos_[m], k_pos_[m], x, y, z);
    value = segy->GetValue(x, y, z+shift);
    blocked_log[m] = value;

    //blocked_log[m] = grid->getRealValue(i_pos_[m], j_pos_[m], k_pos_[m]);
  }

  if (n_angles_ == 0)
    n_angles_ = n_angles;

  if (type == "SEISMIC_DATA") {
    real_seismic_data_.insert(std::pair<int, std::vector<double> >(i_angle, blocked_log));
  }

}

//--------------------------------------------------------------------------------------
void  BlockedLogsCommon::SetLogFromGrid(StormContGrid * storm,
                                        int             i_angle,
                                        int             n_angles,
                                        std::string     type)
{
  std::vector<double> blocked_log(n_blocks_);

  double grid_x = 0.0;
  double grid_y = 0.0;
  double grid_z = 0.0;
  double value  = 0.0;

  for (size_t m = 0 ; m < n_blocks_ ; m++) {
    storm->FindCenterOfCell(i_pos_[m], j_pos_[m], k_pos_[m], grid_x, grid_y, grid_z);
    value = storm->GetValueZInterpolated(grid_x, grid_y, grid_z);
    blocked_log[m] = value;

    //blocked_log[m] = grid->getRealValue(i_pos_[m], j_pos_[m], k_pos_[m]);
  }

  if (n_angles_ == 0)
    n_angles_ = n_angles;

  if (type == "SEISMIC_DATA") {
    real_seismic_data_.insert(std::pair<int, std::vector<double> >(i_angle, blocked_log));
  }
  else if (type == "FACIES_PROB") {
    facies_prob_.insert(std::pair<int, std::vector<double> >(i_angle, blocked_log));
  }
  else if (type == "VP_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Vp", blocked_log));
  }
  else if (type == "VS_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Vs", blocked_log));
  }
  else if (type == "RHO_PREDICTED") {
    continuous_logs_predicted_.insert(std::pair<std::string, std::vector<double> >("Rho", blocked_log));
  }

}


//------------------------------------------------------------------------------
void BlockedLogsCommon::WriteWell(const int                        formats,
                                  const float                      max_hz_background,
                                  const float                      max_hz_seismic,
                                  const std::vector<std::string> & facies_name,
                                  const std::vector<int>         & facies_label) const
{
  if ((formats & IO::RMSWELL) > 0) {
    WriteRMSWell(max_hz_background,
                 max_hz_seismic,
                 facies_name,
                 facies_label);
  }
  if ((formats & IO::NORSARWELL) > 0)
    WriteNorsarWell(max_hz_background,
                    max_hz_seismic);
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::WriteRMSWell(const float                      max_hz_background,
                                     const float                      max_hz_seismic,
                                     const std::vector<std::string> & facies_name,
                                     const std::vector<int>         & facies_label) const
{
  std::string well_name(well_name_);
  NRLib::Substitute(well_name,"/","_");
  NRLib::Substitute(well_name," ","_");
  std::string base_name = IO::PrefixBlockedWells() + well_name.c_str() + IO::SuffixRmsWells();
  std::string file_name = IO::makeFullFileName(IO::PathToWells(), base_name);

  std::ofstream file;
  NRLib::OpenWrite(file, file_name);

  if (!file) {
    LogKit::LogMessage(LogKit::Error,"Error opening "+file_name+" for writing.");
    std::exit(1);
  }

  int n_facies                 = static_cast<int>(facies_map_.size());
  bool got_facies              = (n_facies > 0);
  bool got_actual_synt_seismic = (actual_synt_seismic_data_.size() != 0);
  bool got_well_synt_seismic   = (well_synt_seismic_data_.size() != 0);
  bool got_vp_rho_fac_log      = (vp_facies_filtered_.size() > 0);

  bool got_facies_prob         = (facies_prob_.size() > 0);
  bool got_real_seismic        = (real_seismic_data_.size() > 0);
  bool got_filtered_log        = (cont_logs_seismic_resolution_.size() > 0);
  bool got_predicted           = (continuous_logs_predicted_.size() > 0);

  bool got_interval_log        = (interval_log_.size() > 0);

  int n_logs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz}
  if (got_filtered_log)
    n_logs += 3;
  if (got_vp_rho_fac_log)
    n_logs += 2;
  if (got_facies)
    n_logs += 1;
  if (got_facies_prob)
    n_logs += n_facies;
  if (got_real_seismic)
    n_logs += n_angles_;
  if (got_actual_synt_seismic)
    n_logs += n_angles_;
  if (got_well_synt_seismic)
    n_logs += n_angles_;
  if (got_predicted)
    n_logs += 3;
  if (got_interval_log && got_filtered_log)
    n_logs += 1;
  if (got_interval_log && got_vp_rho_fac_log)
    n_logs += 1;

  std::vector<std::string> params(3);
  params[0] = "Vp";
  params[1] = "Vs";
  params[2] = "Rho";
  file << std::fixed
       << std::setprecision(2);
  //
  // Write HEADER
  //
  file << "1.0\n"
       << "CRAVA\n"
       << IO::PrefixBlockedWells() + well_name_ << " " << x_pos_blocked_[first_B_] << " " << y_pos_blocked_[first_B_] << "\n"
       << n_logs << "\n";

  for (int i = 0; i < 3; i++) {
    file << params[i] << "  UNK lin\n";
    file << params[i] << static_cast<int>(max_hz_background) << "  UNK lin\n";
    file << params[i] << static_cast<int>(max_hz_seismic)    << "  UNK lin\n";
  }
  if (got_filtered_log) {
    for (int i=0;i<3;i++)
      file << params[i] << "_SeismicResolution UNK lin\n";
  }
  if (got_interval_log) {
    file << "Interval Log for SeismicResolution UNK lin\n";
  }
  if (got_predicted) {
    for (int i = 0; i < 3; i++)
      file << params[i] << "_Predicted UNK lin\n";
  }
  if (got_vp_rho_fac_log) {
    file << params[0] << "_ForFacies UNK lin\n";
    file << params[2] << "_ForFacies UNK lin\n";
  }
  if (got_interval_log) {
    file << "Interval Log for ForFacies UNK lin\n";
  }
  if (got_facies) {
    file << "FaciesLog  DISC ";
    for (int i = 0; i < static_cast<int>(facies_name.size()); i++)
      file << " " << facies_label[i] << " " << facies_name[i];
    file << "\n";
  }
  if (got_facies_prob) {
    for (int i=0 ; i<n_facies ; i++)
      file << "FaciesProbabilities" << i << " UNK lin\n";
  }
  if (got_real_seismic) {
    for (int i = 0; i < n_angles_; i++)
      file << "RealSeis" << i << " UNK lin\n";
  }
  if (got_actual_synt_seismic) {
    for (int i = 0; i < n_angles_; i++)
      file << "ActualSyntSeis" << i << " UNK lin\n";
  }
  if (got_well_synt_seismic) {
    for (int i = 0; i < n_angles_; i++)
      file << "WellOptimizedSyntSeis" << i << " UNK lin\n";
  }

  std::vector<int> facies_numbers;
  if (got_facies) {
    for (std::map<int,std::string>::const_iterator it = facies_map_.begin(); it != facies_map_.end(); it++) {
      facies_numbers.push_back(it->first);
    }
  }

  //
  // Write LOGS
  //
  const std::vector<double> & vp  = GetVpBlocked();
  const std::vector<double> & vs  = GetVsBlocked();
  const std::vector<double> & rho = GetRhoBlocked();

  const std::vector<double> & vp_highcut_background  = GetVpHighCutBackground();
  const std::vector<double> & vs_highcut_background  = GetVsHighCutBackground();
  const std::vector<double> & rho_highcut_background = GetRhoHighCutBackground();

  const std::vector<double> & vp_highcut_seismic  = GetVpHighCutSeismic();
  const std::vector<double> & vs_highcut_seismic  = GetVsHighCutSeismic();
  const std::vector<double> & rho_highcut_seismic = GetRhoHighCutSeismic();

  std::vector<double> vp_predicted;
  std::vector<double> vs_predicted;
  std::vector<double> rho_predicted;
  if (got_predicted == true) {
    vp_predicted  = GetVpPredicted();
    vs_predicted  = GetVsPredicted();
    rho_predicted = GetRhoPredicted();
  }

  std::vector<double> vp_seismic_resolution;
  std::vector<double> vs_seismic_resolution;
  std::vector<double> rho_seismic_resolution;
  if (got_filtered_log == true) {
    vp_seismic_resolution  = GetVpSeismicResolution();
    vs_seismic_resolution  = GetVsSeismicResolution();
    rho_seismic_resolution = GetRhoSeismicResolution();
  }

  const std::vector<int> & interval_log = GetIntervalLog();

  for (int i = first_B_; i < last_B_ + 1; i++) {
    file << std::right
         << std::fixed
         << std::setprecision(2)
         << std::setw(9)  << x_pos_blocked_[i] << " "
         << std::setw(10) << y_pos_blocked_[i] << " "
         << std::setw(7)  << z_pos_blocked_[i] << "  "
         << std::setw(7)  << (vp[i]==RMISSING                    ? WELLMISSING : exp(vp[i]))                    << " "
         << std::setw(7)  << (vp_highcut_background[i]==RMISSING ? WELLMISSING : exp(vp_highcut_background[i])) << " "
         << std::setw(7)  << (vp_highcut_seismic[i]==RMISSING    ? WELLMISSING : exp(vp_highcut_seismic[i]))    << " "
         << std::setw(7)  << (vs[i]==RMISSING                    ? WELLMISSING : exp(vs[i]))                    << " "
         << std::setw(7)  << (vs_highcut_background[i]==RMISSING ? WELLMISSING : exp(vs_highcut_background[i])) << " "
         << std::setw(7)  << (vs_highcut_seismic[i]==RMISSING    ? WELLMISSING : exp(vs_highcut_seismic[i]))    << " "
         << std::setprecision(5)
         << std::setw(7) << (rho[i]==RMISSING                    ? WELLMISSING : exp(rho[i]))                    << " "
         << std::setw(7) << (rho_highcut_background[i]==RMISSING ? WELLMISSING : exp(rho_highcut_background[i])) << " "
         << std::setw(7) << (rho_highcut_seismic[i]==RMISSING    ? WELLMISSING : exp(rho_highcut_seismic[i]))    << " ";
    if (got_filtered_log == true && vp_seismic_resolution.size() > 0) {
      file << std::setw(7) << (vp_seismic_resolution[i]==RMISSING  ? WELLMISSING : exp(vp_seismic_resolution[i]))  << "  "
           << std::setw(7) << (vs_seismic_resolution[i]==RMISSING  ? WELLMISSING : exp(vs_seismic_resolution[i]))  << "  "
           << std::setw(7) << (rho_seismic_resolution[i]==RMISSING ? WELLMISSING : exp(rho_seismic_resolution[i])) << "  ";
    }
    if (got_interval_log == true)
      file << std::setw(7) << (interval_log[i]==IMISSING  ? WELLMISSING : interval_log[i])  << "  ";
    if (got_predicted == true && vp_predicted.size() > 0) {
      file << std::setw(7) << (vp_predicted[i]==RMISSING  ? WELLMISSING : exp(vp_predicted[i]))  << "  "
           << std::setw(7) << (vs_predicted[i]==RMISSING  ? WELLMISSING : exp(vs_predicted[i]))  << "  "
           << std::setw(7) << (rho_predicted[i]==RMISSING ? WELLMISSING : exp(rho_predicted[i])) << "  ";
    }
    if (got_vp_rho_fac_log == true) {
      file << std::setw(7) << (vp_facies_filtered_[i]==RMISSING  ? WELLMISSING : exp(vp_facies_filtered_[i]))  << "  "
           << std::setw(7) << (rho_facies_filtered_[i]==RMISSING ? WELLMISSING : exp(rho_facies_filtered_[i])) << "  ";
    }
    if (got_interval_log == true)
      file << std::setw(7) << (interval_log[i]==IMISSING  ? WELLMISSING : interval_log[i])  << "  ";
    if (got_facies)
      file << (facies_blocked_[i]==IMISSING ? static_cast<int>(WELLMISSING) : facies_numbers[facies_blocked_[i]])      << "  ";
    file << std::scientific;
    if (got_facies_prob) {
      for (int a = 0; a < n_facies; a++)
        file << std::setw(12) << (GetFaciesProb(a)[i]==RMISSING ? WELLMISSING : GetFaciesProb(a)[i])          << " ";
      file << " ";
    }
    if (got_real_seismic) {
      for (int a = 0; a < n_angles_; a++)
        file << std::setw(12) << (GetRealSeismicData(a)[i]==RMISSING ? WELLMISSING : GetRealSeismicData(a)[i])          << " ";
      file << " ";
    }
    if (got_actual_synt_seismic) {
      for (int a = 0; a < n_angles_; a++)
        file << std::setw(12) << (actual_synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : actual_synt_seismic_data_[a][i])          << " ";
      file << " ";
    }
    if (got_well_synt_seismic) {
      for (int a = 0; a < n_angles_; a++)
        file << std::setw(12) << (well_synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : well_synt_seismic_data_[a][i])          << " ";
      file << " ";
    }
    file << "\n";
  }

  file.close();
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::WriteNorsarWell(const float max_hz_background,
                                        const float max_hz_seismic) const
{

  double vert_scale = 0.001;
  double hor_scale  = 0.001;

  //Note: At current, only write Vp, Vs and Rho, as others are not supported.

  std::string well_name(well_name_);
  NRLib::Substitute(well_name,"/","_");
  NRLib::Substitute(well_name," ","_");

  //Handle main file.
  std::string base_name = IO::PrefixBlockedWells() + well_name + IO::SuffixNorsarWells();
  std::string file_name = IO::makeFullFileName(IO::PathToWells(), base_name);
  std::ofstream main_file;
  NRLib::OpenWrite(main_file, file_name);
  main_file << std::fixed
           << std::setprecision(2);

  int n_data = last_B_ - first_B_ + 1;

  std::vector<double> md(n_data,0);

  if (continuous_logs_blocked_.count("MD") == 0)
    md[0] = z_pos_blocked_[first_B_]*vert_scale;
  else
    md[0] = GetMDBlocked()[first_B_];
  double dmax = 0;
  double dmin = 1e+30;
  for (int i = first_B_+1; i <= last_B_; i++) {
    double dx = x_pos_blocked_[i]-x_pos_blocked_[i-1];
    double dy = y_pos_blocked_[i]-y_pos_blocked_[i-1];
    double dz = z_pos_blocked_[i]-z_pos_blocked_[i-1];
    double d  = sqrt(dx*dx+dy*dy+dz*dz);
    if (d > dmax)
      dmax = d;
    else if (d < dmin)
      dmin = d;
    if (continuous_logs_blocked_.count("MD") == 0)
      md[i-first_B_] = md[i-first_B_-1] + d*vert_scale;
    else
      md[i-first_B_] = GetMDBlocked()[i];
  }

  main_file << "[Version information]\nVERSION 1000\nFORMAT ASCII\n\n";
  main_file << "[Well information]\n";
  main_file << std::setprecision(5);
  main_file << "MDMIN      km       " << md[0]*vert_scale            << "\n";
  main_file << "MDMAX      km       " << md[n_data-1]*vert_scale     << "\n";
  main_file << "MDMINSTEP  km       " << dmin*vert_scale             << "\n";
  main_file << "MDMAXSTEP  km       " << dmax*vert_scale             << "\n";
  main_file << "UTMX       km       " << x_pos_blocked_[0]*hor_scale << "\n";
  main_file << "UTMY       km       " << y_pos_blocked_[0]*hor_scale << "\n";
  main_file << "EKB        km       " << 0.0f                        << "\n";
  main_file << "UNDEFVAL   no_unit  " << WELLMISSING                 << "\n\n";


  main_file << "[Well track data information]\n";
  main_file << "NUMMD  " << n_data << "\n";
  main_file << "NUMPAR 5\n";
  main_file << "MD      km\n";
  main_file << "TVD     km\n";
  main_file << "TWT     s\n";
  main_file << "UTMX    km\n";
  main_file << "UTMY    km\n";

  std::string log_base_name = IO::PrefixBlockedWells() + well_name + IO::SuffixNorsarLog();
  std::string log_file_name = IO::makeFullFileName(IO::PathToWells(), log_base_name);
  std::string only_name     = NRLib::RemovePath(log_file_name);

  int n_facies = static_cast<int>(facies_blocked_.size());

  bool got_facies              = (n_facies > 0);
  bool got_facies_prob         = (facies_prob_.size() > 0);
  bool got_real_seismic        = (real_seismic_data_.size() > 0);
  bool got_actual_synt_seismic = (actual_synt_seismic_data_.size() > 0);
  bool got_well_synt_seismic   = (well_synt_seismic_data_.size() > 0);
  bool got_filtered_log        = (cont_logs_seismic_resolution_.size() > 0);
  bool got_vp_rho_fac_log      = (vp_facies_filtered_.size() > 0);

  int n_logs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz, seisRes}
  if (got_filtered_log)
    n_logs += 3;
  if (got_facies)
    n_logs += 1;
  if (got_facies_prob)
    n_logs += n_facies;
  if (got_real_seismic)
    n_logs += n_angles_;
  if (got_actual_synt_seismic)
    n_logs += n_angles_;
  if (got_well_synt_seismic)
    n_logs += n_angles_;

  std::vector<std::string> params(3);
  params[0] = "VP";
  params[1] = "VS";
  params[2] = "RHO";

  std::vector<std::string> unit(3);
  unit[0] = "km/s";
  unit[1] = "km/s";
  unit[2] = "tons/m3";

  int n_files = 3;
  if (got_filtered_log)
    n_files += 1;
  if (got_vp_rho_fac_log)
    n_files += 1;

  std::vector<std::string> postfix(5);
  postfix[0] = "orig";
  postfix[1] = NRLib::ToString(max_hz_background)+"Hz";
  postfix[2] = NRLib::ToString(max_hz_seismic)+"Hz";
  postfix[3] = "filtered";
  postfix[4] = "forFacies";

  for (int f=0; f < n_files; f++) {
    main_file << "\n[Well log data information]\n";
    main_file << "LOGNAME log_" << postfix[f] << "\n";
    main_file << "IN_FILE " << only_name << f <<"\n";
    main_file << "NUMMD " << n_data << "\n";
    main_file << "NUMPAR " << 4 << "\n"; //Also count md.
    main_file << "MD      km\n";
    for (int i =0 ; i<3 ; i++)
      main_file << params[i] << " " << unit[i] << "\n";
  }
  main_file.close();

  //Write the track file.
  std::string track_base_name = IO::PrefixBlockedWells() + well_name + IO::SuffixNorsarTrack();
  std::string track_file_name = IO::makeFullFileName(IO::PathToWells(), track_base_name);
  std::ofstream track_file;
  NRLib::OpenWrite(track_file, track_file_name.c_str());
  track_file << std::right
             << std::fixed
             << std::setprecision(2)
             << "[NORSAR Well Track]\n";

  //Note: logFileName created above, needed in mainFile.
  std::vector<std::ofstream *> log_files;
  log_files.resize(n_files);
  for (int f=0; f < n_files; f++) {
    log_files[f] = new std::ofstream();
    file_name = log_file_name+NRLib::ToString(f);
    NRLib::OpenWrite(*(log_files[f]), file_name.c_str());
    *(log_files[f]) << "[NORSAR Well Log]\n";
    *(log_files[f]) << "[See header (.nwh) file for log information]\n";
  }

  const std::vector<double> & vp =  GetVpBlocked();
  const std::vector<double> & vs =  GetVsBlocked();
  const std::vector<double> & rho = GetRhoBlocked();

  const std::vector<double> & vp_highcut_background  = GetVpHighCutBackground();
  const std::vector<double> & vs_highcut_background  = GetVsHighCutBackground();
  const std::vector<double> & rho_highcut_background = GetRhoHighCutBackground();

  const std::vector<double> & vp_highcut_seismic  = GetVpHighCutSeismic();
  const std::vector<double> & vs_highcut_seismic  = GetVsHighCutSeismic();
  const std::vector<double> & rho_highcut_seismic = GetRhoHighCutSeismic();

  const std::vector<double> & vp_seismic_resolution  = GetVpSeismicResolution();
  const std::vector<double> & vs_seismic_resolution  = GetVsSeismicResolution();
  const std::vector<double> & rho_seismic_resolution = GetRhoSeismicResolution();

  for (int i = first_B_; i <= last_B_; i++) {
    track_file << std::setprecision(5) << std::setw(7)
               << md[i-first_B_] << " " << std::setw(7) << z_pos_blocked_[i]*vert_scale << " " << z_pos_blocked_[i]*vert_scale
               << " " << std::setw(10)<< x_pos_blocked_[i]*hor_scale << " " << std::setw(10)<< y_pos_blocked_[i]*hor_scale << "\n";

    *(log_files[0]) << std::right   << std::fixed << std::setprecision(5)
                    << std::setw(7) << md[i-first_B_] << " "
                    << std::setw(7) << (vp[i]==RMISSING                     ? WELLMISSING : exp(vp[i]))                     << " "
                    << std::setw(7) << (vs[i]==RMISSING                     ? WELLMISSING : exp(vs[i]))                     << " "
                    << std::setw(7) << (rho[i]==RMISSING                    ? WELLMISSING : exp(rho[i]))                    << "\n";
    *(log_files[1]) << std::right   << std::fixed << std::setprecision(5)
                    << std::setw(7) << md[i-first_B_] << " "
                    << std::setw(7) << (vp_highcut_background[i]==RMISSING  ? WELLMISSING : exp(vp_highcut_background[i]))  << " "
                    << std::setw(7) << (vs_highcut_background[i]==RMISSING  ? WELLMISSING : exp(vs_highcut_background[i]))  << " "
                    << std::setw(7) << (rho_highcut_background[i]==RMISSING ? WELLMISSING : exp(rho_highcut_background[i])) << "\n";
    *(log_files[2]) << std::right   << std::fixed << std::setprecision(5)
                    << std::setw(7) << md[i-first_B_] << " "
                    << std::setw(7) << (vp_highcut_seismic[i]==RMISSING     ? WELLMISSING : exp(vp_highcut_seismic[i]))     << " "
                    << std::setw(7) << (vs_highcut_seismic[i]==RMISSING     ? WELLMISSING : exp(vs_highcut_seismic[i]))     << " "
                    << std::setw(7) << (rho_highcut_seismic[i]==RMISSING    ? WELLMISSING : exp(rho_highcut_seismic[i]))    << "\n";
    if (got_filtered_log) {
      *(log_files[3]) << std::right   << std::fixed << std::setprecision(5)
                      << std::setw(7) << md[i-first_B_] << " "
                      << std::setw(7) << (vp_seismic_resolution[i]==RMISSING  ? WELLMISSING : exp(vp_seismic_resolution[i]))  << " "
                      << std::setw(7) << (vs_seismic_resolution[i]==RMISSING  ? WELLMISSING : exp(vs_seismic_resolution[i]))  << " "
                      << std::setw(7) << (rho_seismic_resolution[i]==RMISSING ? WELLMISSING : exp(rho_seismic_resolution[i])) << "\n";
    }
    if (got_vp_rho_fac_log) {
      *(log_files[3]) << std::right   << std::fixed << std::setprecision(5)
                      << std::setw(7) << md[i-first_B_] << " "
                      << std::setw(7) << (vp_facies_filtered_[i]==RMISSING  ? WELLMISSING : exp(vp_facies_filtered_[i]))  << " "
                      << std::setw(7) << (rho_facies_filtered_[i]==RMISSING ? WELLMISSING : exp(rho_facies_filtered_[i])) << "\n";
    }
  }
  track_file.close();
  for (int f=0; f < n_files; f++) {
    log_files[f]->close();
    delete log_files[f];
  }
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::SetSpatialFilteredLogs(std::vector<double>       & filtered_log,
                                               int                         n_data,
                                               std::string                 type,
                                               const std::vector<double> & bg)
{
  std::vector<double> blocked_log(n_blocks_);

  assert(n_blocks_ == static_cast<unsigned int>(n_data));
  for (int i = 0; i < n_data; i++)
    blocked_log[i] = filtered_log[i]+bg[i];

  if (type == "VP_SEISMIC_RESOLUTION")
    cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Vp", blocked_log));
  else if (type == "VS_SEISMIC_RESOLUTION")
    cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Vs", blocked_log));
  else if (type == "RHO_SEISMIC_RESOLUTION")
    cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Rho", blocked_log));
  else if (type == "VP_FOR_FACIES")
    vp_facies_filtered_ = blocked_log;
  else if (type == "RHO_FOR_FACIES")
    rho_facies_filtered_ = blocked_log;
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::GenerateSyntheticSeismic(const NRLib::Matrix        & reflection_matrix,
                                                 std::vector<Wavelet *>     & wavelet,
                                                 int                          nz,
                                                 int                          nzp,
                                                 const Simbox               * simbox,
                                                 bool                         well_opt)
{
  int          i, j;
  int          start, length;

  fftw_complex   cAmp;
  Wavelet1D    * local_wavelet;

  int    cnzp = nzp/2+1;
  int    rnzp = 2*cnzp;

  std::vector<double> synt_seis(nz);
  std::vector<double> vp_vert(n_layers_);
  std::vector<double> vs_vert(n_layers_);
  std::vector<double> rho_vert(n_layers_);

  fftw_real    * cpp_r       = new fftw_real[rnzp];
  fftw_complex * cpp_c       = reinterpret_cast<fftw_complex*>(cpp_r);
  fftw_real    * synt_seis_r = new fftw_real[rnzp];
  fftw_complex * synt_seis_c = reinterpret_cast<fftw_complex*>(synt_seis_r);

  GetVerticalTrend(GetVpBlocked(), vp_vert);
  GetVerticalTrend(GetVsBlocked(), vs_vert);
  GetVerticalTrend(GetRhoBlocked(), rho_vert);

  std::vector<bool> has_data(n_layers_);
  for (i=0; i < n_layers_; i++)
    has_data[i] = vp_vert[i] != RMISSING && vs_vert[i] != RMISSING && rho_vert[i] != RMISSING;

  FindContinuousPartOfData(has_data, n_layers_, start, length);

  float scale = static_cast<float>(simbox->getRelThick(i_pos_[0], j_pos_[0]));

  int n_angles = static_cast<int>(wavelet.size());
  for (i=0; i < n_angles; i++) {
    if (wavelet[i] != NULL) {
      for (j=0; j < rnzp; j++) {
        cpp_r[j] = 0;
        synt_seis_r[j] = 0;
      }
      float * refl_coef = new float[3];
      refl_coef[0] = static_cast<float>(reflection_matrix(i,0));
      refl_coef[1] = static_cast<float>(reflection_matrix(i,1));
      refl_coef[2] = static_cast<float>(reflection_matrix(i,2));
      FillInCpp(refl_coef, start, length, cpp_r, nzp);
      Utils::fft(cpp_r,cpp_c,nzp);
      delete [] refl_coef;
      local_wavelet = wavelet[i]->createLocalWavelet1D(i_pos_[0], j_pos_[0]);
      local_wavelet->fft1DInPlace();
     // float sf = wavelet[i]->getLocalStretch(ipos_[0],jpos_[0]);
     // what about relative thickness ????
     // float sf = wavelet[i]->getLocalStretch(ipos_[0],jpos_[0])*Relativethikness... Need simbox;

      for (j=0; j < cnzp; j++) {
        cAmp =  local_wavelet->getCAmp(j, scale);
        synt_seis_c[j].re = cpp_c[j].re*cAmp.re + cpp_c[j].im*cAmp.im;
        synt_seis_c[j].im = cpp_c[j].im*cAmp.re - cpp_c[j].re*cAmp.im;
      }

      Utils::fftInv(synt_seis_c,synt_seis_r,nzp);

      for (j=0; j<nz; j++)
        synt_seis[j] = 0.0; // Do not use RMISSING (fails in setLogFromVerticalTrend())

      for (j=start; j < start+length; j++)
        synt_seis[j] = synt_seis_r[j];

      SetNAngles(n_angles);
      std::string log_type;
      if (well_opt == true)
        log_type = "WELL_SYNTHETIC_SEISMIC";
      else
        log_type = "ACTUAL_SYNTHETIC_SEISMIC";
      SetLogFromVerticalTrend(synt_seis,
                              cont_logs_seismic_resolution_,
                              actual_synt_seismic_data_,
                              well_synt_seismic_data_,
                              log_type,
                              i);

      //localWavelet->fft1DInPlace();
      delete local_wavelet;
    }
  }

  delete [] cpp_r;
  delete [] synt_seis_r;

}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::FindMeanVsVp(const NRLib::Surface<double> & top,
                                     const NRLib::Surface<double> & bot,
                                     double                       & mean_vs_vp,
                                     int                          & n_vs_vp) const
{
  std::vector<bool> active_cell(n_blocks_, true);

  for (size_t i = 0; i < n_blocks_; i++) {
    double z_top  = top.GetZ(x_pos_blocked_[i], y_pos_blocked_[i]);
    double z_base = bot.GetZ(x_pos_blocked_[i], y_pos_blocked_[i]);

    if ((z_pos_blocked_[i] < z_top) || (z_pos_blocked_[i] > z_base))
      active_cell[i] = false;
  }

  mean_vs_vp = 0.0;
  n_vs_vp    = 0;

  for (size_t i = 0; i < n_blocks_; i++) {
    if (GetVpHighCutBackground()[i] != RMISSING && GetVsHighCutBackground()[i] != RMISSING && active_cell[i] == true) {
      mean_vs_vp += GetVsHighCutBackground()[i] / GetVpHighCutBackground()[i];
      n_vs_vp    += 1;
    }
  }

  mean_vs_vp /= n_vs_vp;
}

//--------------------------------------------------------------------------------------
bool
BlockedLogsCommon::VolumeFocus(const NRLib::Volume                            & volume,
                               std::vector<double>                            & x_pos_blocked,
                               std::vector<double>                            & y_pos_blocked,
                               std::vector<double>                            & z_pos_blocked,
                               std::vector<int>                               & facies_blocked,
                               std::vector<int>                               & s_pos,
                               std::vector<int>                               & i_pos,
                               std::vector<int>                               & j_pos,
                               std::vector<int>                               & k_pos,
                               std::map<std::string, std::vector<double> >    & continuous_logs_blocked,
                               std::map<std::string, std::vector<int> >       & discrete_logs_blocked,
                               std::map<std::string, std::vector<double> >    & cont_logs_seismic_resolution,
                               std::map<std::string, std::vector<double> >    & cont_logs_highcut_background,
                               std::map<std::string, std::vector<double> >    & cont_logs_highcut_seismic,
                               std::vector<std::vector<double> >              & actual_synt_seismic_data,
                               std::vector<std::vector<double> >              & well_synt_seismic_data) const
{
  size_t first_index = 0;
  bool   first_found = 0;

  std::vector<std::pair<size_t, size_t> > intervals;
  for (size_t i=0;(i<x_pos_blocked_.size());i++) {
    bool inside = volume.IsInside(x_pos_blocked_[i], y_pos_blocked_[i], z_pos_blocked_[i]);
    if (first_found == false) {
      if (inside == true) {
        first_index = i;
        first_found = true;
      }
    }
    else if (inside == false) {
      std::pair<size_t, size_t> interval(first_index,i);
      intervals.push_back(interval);
      first_found = false;
    }
  }
  if (first_found == true || intervals.size() > 0) {
    if (first_found == true) {
      std::pair<size_t, size_t> interval(first_index,x_pos_blocked_.size());
      intervals.push_back(interval);
    }
    UpdateLog(x_pos_blocked, intervals);
    UpdateLog(y_pos_blocked, intervals);
    UpdateLog(z_pos_blocked, intervals);
    UpdateLog(facies_blocked, intervals);

    UpdateLog(s_pos, intervals);
    UpdateLog(i_pos, intervals);
    UpdateLog(j_pos, intervals);
    UpdateLog(k_pos, intervals);

    std::map<std::string, std::vector<double> >::iterator map_ci;
    for (map_ci = continuous_logs_blocked.begin(); map_ci != continuous_logs_blocked.end(); ++map_ci)
      UpdateLog(map_ci->second, intervals);

    std::map<std::string, std::vector<int> >::iterator map_di;
    for (map_di = discrete_logs_blocked.begin(); map_di != discrete_logs_blocked.end(); ++map_di)
      UpdateLog(map_di->second, intervals);

    for (map_ci = cont_logs_seismic_resolution.begin(); map_ci != cont_logs_seismic_resolution.end(); ++map_ci)
      UpdateLog(map_ci->second, intervals);

    for (map_ci = cont_logs_highcut_background.begin(); map_ci != cont_logs_highcut_background.end(); ++map_ci)
      UpdateLog(map_ci->second, intervals);

    for (map_ci = cont_logs_highcut_seismic.begin(); map_ci != cont_logs_highcut_seismic.end(); ++map_ci)
      UpdateLog(map_ci->second, intervals);

    for (size_t i=0;i<actual_synt_seismic_data_.size();i++)
      UpdateLog(actual_synt_seismic_data[i], intervals);

    for (size_t i=0;i<well_synt_seismic_data_.size();i++)
      UpdateLog(well_synt_seismic_data[i], intervals);
  }
  if (intervals.size() > 0)
    return(true);
  else
    return(false);
}


//--------------------------------------------------------------------------------------
void
BlockedLogsCommon::UpdateLog(std::vector<double>                           & data,
                             const std::vector<std::pair<size_t, size_t> > & intervals) const
{
  if (data.size() > 0) {
    size_t start = 0;
    size_t end;
    for (size_t index=0;index <intervals.size();index++) {
      end = intervals[index].first;
      fill(data.begin()+start, data.begin()+end, RMISSING);
      start = intervals[index].second;
    }
    fill(data.begin()+start, data.end(), RMISSING);
  }
}

//--------------------------------------------------------------------------------------
void
BlockedLogsCommon::UpdateLog(std::vector<int>                              & data,
                             const std::vector<std::pair<size_t, size_t> > & intervals) const
{
  if (data.size() > 0) {
    size_t start = 0;
    size_t end;
    for (size_t index=0;index <intervals.size();index++) {
      end = intervals[index].first;
      fill(data.begin()+start, data.begin()+end, IMISSING);
      start = intervals[index].second;
    }
    fill(data.begin()+start, data.end(), IMISSING);
  }
}


//--------------------------------------------------------------------------------------
void BlockedLogsCommon::FindXYZForVirtualPart(const Simbox              * simbox,
                                              const std::vector<int>    & i_pos,
                                              const std::vector<int>    & j_pos,
                                              const std::vector<int>    & k_pos,
                                              const int                 & n_blocks,
                                              const int                 & first_B,
                                              const int                 & last_B,
                                              std::vector<double>       & x_pos_blocked,
                                              std::vector<double>       & y_pos_blocked,
                                              std::vector<double>       & z_pos_blocked) const
{
  //
  // If the ends have undefined coordinates we use the nearest defined
  // coordinate for x and y and the block cell centre for z
  //
  for (int b = 0 ; b < first_B ; b++) {
    double x,y,z;
    simbox->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    x_pos_blocked[b] = x_pos_blocked[first_B];
    y_pos_blocked[b] = y_pos_blocked[first_B];
    z_pos_blocked[b] = z;
  }

  for (int b = last_B + 1; b < n_blocks; b++) {
    double x,y,z;
    simbox->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    x_pos_blocked[b] = x_pos_blocked[last_B];
    y_pos_blocked[b] = y_pos_blocked[last_B];
    z_pos_blocked[b] = z;
  }
}

//--------------------------------------------------------------------------------------
void BlockedLogsCommon::FindXYZForVirtualPart(const MultiIntervalGrid   * multiple_interval_grid,
                                              const std::vector<int>    & i_pos,
                                              const std::vector<int>    & j_pos,
                                              const std::vector<int>    & k_pos,
                                              int                         n_blocks,
                                              int                         first_B,
                                              int                         last_B,
                                              std::vector<double>       & x_pos_blocked,
                                              std::vector<double>       & y_pos_blocked,
                                              std::vector<double>       & z_pos_blocked) const
{
  //
  // If the ends have undefined coordinates we use the nearest defined
  // coordinate for x and y and the block cell centre for z
  //
  const std::vector<Simbox *> interval_simboxes = multiple_interval_grid->GetIntervalSimboxes();
  int i_interval = 0;
  int n_intervals = multiple_interval_grid->GetNIntervals();

  for (int b = 0 ; b < first_B ; b++) {
    double x,y,z;
    interval_simboxes[i_interval]->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    while (interval_simboxes[i_interval]->IsPointBetweenVisibleSurfaces(x, y, z) == false) {
      i_interval++;
      interval_simboxes[i_interval]->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    }
    x_pos_blocked[b] = x_pos_blocked[first_B];
    y_pos_blocked[b] = y_pos_blocked[first_B];
    z_pos_blocked[b] = z;
  }

  for (int b = last_B + 1; b < n_blocks; b++) {
    double x,y,z;
    interval_simboxes[i_interval]->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    while (i_interval < (n_intervals -1) && interval_simboxes[i_interval]->IsPointBetweenVisibleSurfaces(x, y, z) == false) {
      i_interval++;
      interval_simboxes[i_interval]->getCoord(i_pos[b], j_pos[b], k_pos[b], x, y, z);
    }
    x_pos_blocked[b] = x_pos_blocked[last_B];
    y_pos_blocked[b] = y_pos_blocked[last_B];
    z_pos_blocked[b] = z;
  }
}

//------------------------------------------------------------------------------
void BlockedLogsCommon::FindTrendPositions(const std::vector<int> & i_pos,
                                           const std::vector<int> & j_pos,
                                           const std::vector<int> & k_pos,
                                           const int              & n_blocks,
                                           const CravaTrend       & trend_cubes,
                                           std::vector<double>    & s1,
                                           std::vector<double>    & s2)
{
  s1.resize(n_blocks, RMISSING);
  s2.resize(n_blocks, RMISSING);

  std::vector<double> positions(2, RMISSING);
  for (int i = 0; i < n_blocks; i++) {
    positions = trend_cubes.GetTrendPosition(i_pos[i], j_pos[i], k_pos[i]);
    s1[i] = positions[0];
    s2[i] = positions[1];
  }
}

//------------------------------------------------------------------------------
void BlockedLogsCommon::AssignToFacies(const std::vector<double>         & well_log,
                                       const std::vector<int>            & facies_log,
                                       const std::vector<int>            & facies_numbers,
                                       std::vector<std::vector<double> > & blocked_log) const
{
  if (well_log.size() > 0) {
    int n_facies = static_cast<int>(blocked_log.size());
    int n_blocks = static_cast<int>(blocked_log[0].size());

    for (int m = 0; m < n_blocks; m++) {
      if (well_log[m] != RMISSING && facies_log[m] != IMISSING) {
        for (int j = 0; j < n_facies; j++) {
          if (facies_numbers[j] == facies_log[m])
            blocked_log[j][m] = well_log[m];
        }
      }
    }
  }
}

void BlockedLogsCommon::CalculateBulkShear(const int                         & n_blocks,
                                           const int                         & n_facies,
                                           std::vector<std::vector<double> > & bulk_modulus,
                                           std::vector<std::vector<double> > & shear_modulus)
{
  bulk_modulus.resize(n_facies,  std::vector<double>(n_blocks, RMISSING));
  shear_modulus.resize(n_facies, std::vector<double>(n_blocks, RMISSING));

  for (int i = 0; i < n_facies; i++) {
    for (int j = 0; j < n_blocks; j++) {
      if (vp_for_facies_[i][j] != RMISSING && vs_for_facies_[i][j] != RMISSING && rho_for_facies_[i][j] != RMISSING) {

        double bulk;
        double shear;

        DEMTools::CalcElasticParamsFromSeismicParams(exp(vp_for_facies_[i][j]), exp(vs_for_facies_[i][j]), exp(rho_for_facies_[i][j]), bulk, shear);

        bulk_modulus[i][j]  = static_cast<float>(bulk);
        shear_modulus[i][j] = static_cast<float>(shear);
      }
    }
  }
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetVpForFacies(const std::string & facies_name)
{
  std::vector<double> vp_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      vp_given_facies = vp_for_facies_[i];
  }

  for(size_t i=0;i<vp_given_facies.size();i++) {
    if(vp_given_facies[i] != RMISSING)
      vp_given_facies[i] = exp(vp_given_facies[i]);
  }

  return vp_given_facies;
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetVsForFacies(const std::string & facies_name)
{
  std::vector<double> vs_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      vs_given_facies = vs_for_facies_[i];
  }

  for(size_t i=0;i<vs_given_facies.size();i++) {
    if(vs_given_facies[i] != RMISSING)
      vs_given_facies[i] = exp(vs_given_facies[i]);
  }

  return vs_given_facies;
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetRhoForFacies(const std::string & facies_name)
{
  std::vector<double> rho_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      rho_given_facies = rho_for_facies_[i];
  }

  for(size_t i=0;i<rho_given_facies.size();i++) {
    if(rho_given_facies[i] != RMISSING)
      rho_given_facies[i] = exp(rho_given_facies[i]);
  }

  return rho_given_facies;
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetBulkForFacies(const std::string & facies_name)
{
  std::vector<double> bulk_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      bulk_given_facies = bulk_modulus_[i];
  }

  return bulk_given_facies;
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetShearForFacies(const std::string & facies_name)
{
  std::vector<double> shear_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      shear_given_facies = shear_modulus_[i];
  }

  return shear_given_facies;
}

//------------------------------------------------------------------------------
std::vector<double>
BlockedLogsCommon::GetPorosityForFacies(const std::string & facies_name)
{
  std::vector<double> porosity_given_facies;
  for (size_t i = 0; i < facies_names_.size(); i++) {
    if (facies_name == facies_names_[i])
      porosity_given_facies = porosity_for_facies_[i];
  }

  return porosity_given_facies;
}
