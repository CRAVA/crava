/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "src/commondata.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"


BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well    * const well_data,
                                     const Simbox         * const estimation_simbox,
                                     bool                   interpolate):
n_blocks_(0),
well_name_(""),
n_layers_(estimation_simbox->getnz()){

  BlockWell(well_data, estimation_simbox, continuous_log_names_, discrete_log_names_, continuous_logs_, discrete_logs_, interpolate);
  n_continuous_logs_ = static_cast<int>(continuous_logs_.size());
  n_discrete_logs_ = static_cast<int>(discrete_logs_.size());
}

BlockedLogsCommon::~BlockedLogsCommon(){

}

void BlockedLogsCommon::BlockWell(const NRLib::Well                  * const well_data,
                                  const Simbox                       * const estimation_simbox,
                                  std::map<std::string, int>         & continuous_log_names,
                                  std::map<std::string, int>         & discrete_log_names,
                                  std::vector<std::vector<double> >  & continuous_logs,
                                  std::vector<std::vector<int> >     & discrete_logs,
                                  bool                                 interpolate){

  well_name_ = well_data->GetWellName();
  // Get size of vector
  std::map<std::string,std::vector<double> > all_logs_temp = well_data->GetContLog();
  std::map<std::string,std::vector<double> >::iterator it = all_logs_temp.begin();
  unsigned int nd = it->second.size();

  (void) discrete_log_names;
  (void) discrete_logs;

  std::map<std::string,std::vector<double> > continuous_logs_well = well_data->GetContLog();
  std::map<std::string,std::vector<int> > discrete_logs_well = well_data->GetDiscLog();

  int index_continuous = 0;
  for (std::map<std::string,std::vector<double> >::iterator it=continuous_logs_well.begin(); it!=continuous_logs_well.end(); ++it){
    continuous_log_names.insert(std::pair<std::string, int>(it->first, index_continuous));
    continuous_logs.push_back(continuous_logs_well.find(it->first)->second);
    index_continuous++;
  }

  int index_discrete = 0;
  for (std::map<std::string,std::vector<int> >::iterator it=discrete_logs_well.begin(); it!=discrete_logs_well.end(); ++it){
    discrete_log_names.insert(std::pair<std::string, int>(it->first, index_discrete));
    discrete_logs.push_back(discrete_logs_well.find(it->first)->second);
    index_discrete++;
  }

  // are all continuous logs of the same length?
  std::vector<int> b_ind(nd); // Gives which block each well log entry contributes to

  FindSizeAndBlockPointers(estimation_simbox, b_ind);
  FindBlockIJK(estimation_simbox, b_ind);

  BlockCoordinateLog(b_ind, this->GetXpos(), x_pos_);
  BlockCoordinateLog(b_ind, this->GetYpos(), y_pos_);
  BlockCoordinateLog(b_ind, this->GetZpos(), z_pos_);

  // Could also loop over all continuous logs and fetch the names here
  if(this->HasContLog("MD")){
    std::vector<double> md_blocked;
    BlockCoordinateLog(b_ind, continuous_logs_[continuous_log_names_.find("MD")->second], md_blocked);
    continuous_logs_blocked_.push_back(md_blocked);
    continuous_log_names_blocked_.insert(std::pair<std::string, int>("MD", continuous_logs_blocked_.size()-1));
  }
  if(this->HasContLog("Vp")){
    std::vector<double> vp_blocked;
    BlockCoordinateLog(b_ind, continuous_logs_[continuous_log_names_.find("Vp")->second], vp_blocked);
    continuous_logs_blocked_.push_back(vp_blocked);
    continuous_log_names_blocked_.insert(std::pair<std::string, int>("Vp", continuous_logs_blocked_.size()-1));
  }
  if(this->HasContLog("Vs")){
    std::vector<double> vs_blocked;
    BlockCoordinateLog(b_ind, continuous_logs_[continuous_log_names_.find("Vs")->second], vs_blocked);
    continuous_logs_blocked_.push_back(vs_blocked);
    continuous_log_names_blocked_.insert(std::pair<std::string, int>("Vs", continuous_logs_blocked_.size()-1));
  }
  if(this->HasContLog("Rho")){
    std::vector<double> rho_blocked;
    BlockCoordinateLog(b_ind, continuous_logs_[continuous_log_names_.find("Rho")->second], rho_blocked);
    continuous_logs_blocked_.push_back(rho_blocked);
    continuous_log_names_blocked_.insert(std::pair<std::string, int>("Rho", continuous_logs_blocked_.size()-1));
  }

  // all discrete logs (20130625 EN: No discrete logs for now)

  if(interpolate){
    for(unsigned int i=1;i<nd;i++) {
      if(abs(b_ind[i]-b_ind[i-1]) > 1) {
        int start, end;
        if(b_ind[i] > b_ind[i-1]) {
          start = b_ind[i-1];
          end   = b_ind[i];
        }
        else {
          start = b_ind[i];
          end   = b_ind[i-1];
        }
        for(int j = start+1;j<end;j++) {
          float t = static_cast<float>(j-start)/static_cast<float>(end-start);
          InterpolateContinuousLog(x_pos_, start, end, j, t);
          InterpolateContinuousLog(y_pos_, start, end, j, t);
          InterpolateContinuousLog(z_pos_, start, end, j, t);

          // all blocked continuous logs
          for(std::map<std::string,int >::iterator it=continuous_log_names_blocked_.begin(); it!=continuous_log_names_blocked_.end(); ++it){
            InterpolateContinuousLog(continuous_logs_blocked_[it->second], start, end, j, t);
          }

        }
      }
    }
  }

}


void  BlockedLogsCommon::FindSizeAndBlockPointers(const Simbox         * const estimation_simbox,
                                                  std::vector<int>     & b_ind){
  int   nd = static_cast<int>(b_ind.size());
  const std::vector<double> x_pos = GetXpos();
  const std::vector<double> y_pos = GetYpos();
  const std::vector<double> tvd   = GetTVD();

  //
  // Find first cell in Simbox that the well hits
  //
  int first_I(IMISSING);
  int first_J(IMISSING);
  int first_K(IMISSING);
  for (int m = 0 ; m < nd ; m++) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], tvd[m], first_I, first_J, first_K);
    if (first_I != IMISSING && first_J != IMISSING && first_K != IMISSING) {
      first_M_ = m;
      break;
    }
  }
  //
  // Find last cell in Simbox that the well hits
  //
  int last_I(IMISSING);
  int last_J(IMISSING);
  int last_K(IMISSING);
  for (int m = nd - 1 ; m > 0 ; m--) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], tvd[m], last_I, last_J, last_K);
    if (last_I != IMISSING && last_J != IMISSING && last_K != IMISSING) {
      last_M_ = m;
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
  b_ind[first_M_] = first_K; // The first defined well log entry contributes to this block.

  //
  // The well positions used to be given in float rather than double. Unfortunately, this
  // allowed a well to oscillate between two or more cells, leading to a breakdown of the
  // algorithm below. To remedy for this we introduced array simboxInd which records the
  // indices of the simbox cells that are already accounted for, so that these are not
  // enlisted more than one time.
  //
  int * simbox_ind = new int[nd];                                     // help hack
  const int nx    = estimation_simbox->getnx();                                 // help hack
  const int ny    = estimation_simbox->getny();                                 // help hack
  simbox_ind[0] = nx*ny*old_K + nx*old_J + old_I;                        // help hack

  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], tvd[m], new_I ,new_J, new_K);

    if (new_I != old_I || new_J != old_J || new_K != old_K) {

      int  this_ind = nx*ny*new_K + nx*new_J + new_I;                    // help hack
      bool block_not_listed = true;                                    // help hack
      for (int l = 0 ; l < n_defined_blocks ; l++) {                   // help hack
        if (this_ind == simbox_ind[l]) {                               // help hack
          block_not_listed = false;                                    // help hack
          break;                                                     // help hack
        }                                                            // help hack
      }                                                              // help hack
      if (block_not_listed) {                                          // help hack
        simbox_ind[n_defined_blocks+1] = this_ind;                       // help hack
        old_I = new_I;
        old_J = new_J;
        old_K = new_K;
        n_defined_blocks++;
      }
      else {
      }
    }
    b_ind[m] = first_K + n_defined_blocks;
  }
  n_defined_blocks++;
  //
  // Why we cannot use nBlocks_ = nDefined:
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
  n_blocks_ = first_K + n_defined_blocks + (n_layers_ - last_K - 1);

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"first_M_, last_M_          = %d, %d    \n",first_M_,last_M_);
    LogKit::LogFormatted(LogKit::Low,"n_layers_                  = %d        \n",n_layers_);
    LogKit::LogFormatted(LogKit::Low,"first_I,first_J,first_K     = %d, %d, %d\n",first_I,first_J,first_K);
    LogKit::LogFormatted(LogKit::Low,"last_I,last_J,last_K        = %d, %d, %d\n",last_I,last_J,last_K);
    LogKit::LogFormatted(LogKit::Low,"n_defined_blocks, n_blocks_ = %d, %d    \n",n_defined_blocks,n_blocks_);
  }
  delete [] simbox_ind;
}

void    BlockedLogsCommon::FindBlockIJK(const Simbox             * const estimation_simbox,
                                        const std::vector<int>   & bInd){
  i_pos_.resize(n_blocks_);
  j_pos_.resize(n_blocks_);
  k_pos_.resize(n_blocks_);

  const std::vector<double> x_pos = this->GetXpos();
  const std::vector<double> y_pos = this->GetYpos();
  const std::vector<double> z_pos = this->GetZpos();

  //
  // Set IJK for virtual part of well in upper part of simbox
  //
  int b = -1; // block counter;
  int first_I, first_J, first_K;
  estimation_simbox->getIndexes(x_pos[first_M_], y_pos[first_M_], z_pos[first_M_], first_I, first_J, first_K);
  for (int k = 0 ; k < first_K ; k++) {
    b++;
    i_pos_[b] = first_I;
    j_pos_[b] = first_J;
    k_pos_[b] = k;
  }

  //
  // Set IJK for the defined part of the well
  //
  b = first_K;
  i_pos_[b] = first_I;
  j_pos_[b] = first_J;
  k_pos_[b] = first_K;
  int i, j, k;
  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    if (bInd[m] != bInd[m - 1]) {
      b++;
      estimation_simbox->getIndexes(x_pos[m], y_pos[m], z_pos[m], i, j, k);
      i_pos_[b] = i;
      j_pos_[b] = j;
      k_pos_[b] = k;
    }
  }
  first_B_ = first_K;
  last_B_  = b;

  //
  // Set IJK for the virtual part of well in lower part of simbox
  //
  int last_I,  last_J,  last_K;
  estimation_simbox->getIndexes(x_pos[last_M_], y_pos[last_M_], z_pos[last_M_], last_I, last_J, last_K);
  for (int k = last_K + 1 ; k < n_layers_ ; k++) {
    b++;
    i_pos_[b] = last_I;
    j_pos_[b] = last_J;
    k_pos_[b] = k;
  }

  dz_ = static_cast<float>(estimation_simbox->getRelThick(i_pos_[0],j_pos_[0])*estimation_simbox->getdz());

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"firstB_, lastB_        = %d, %d    \n",first_B_,last_B_);
    LogKit::LogFormatted(LogKit::Low,"firstI, firstJ, firstK = %d, %d, %d\n",first_I, first_J, first_K);
    LogKit::LogFormatted(LogKit::Low,"lastI,  lastJ,  lastK  = %d, %d, %d\n",last_I, last_J, last_K);
    for (int b = 0 ; b < n_blocks_ ; b++)
      LogKit::LogFormatted(LogKit::Low,"b=%d   i,j,k=%d,%d,%d\n",b,i_pos_[b],j_pos_[b],k_pos_[b]);
  }
}

void BlockedLogsCommon::BlockCoordinateLog(const std::vector<int>    &  b_ind,
                                           const std::vector<double> &  coord,
                                           std::vector<double>       &  blocked_coord)
{
  //
  // Initialise arrays
  //
  blocked_coord.resize(n_blocks_, 0.0f);
  std::vector<int> count(n_blocks_, 0);

  //
  // Block log
  //
  for (int m = first_M_ ; m < last_M_ + 1 ; m++) {
    blocked_coord[b_ind[m]] += coord[m];
    count[b_ind[m]]++;
  }
  for (int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0)
      blocked_coord[l] /= count[l];
    else
      blocked_coord[l]  = RMISSING;
  }

}

//------------------------------------------------------------------------------
void BlockedLogsCommon::BlockContinuousLog(const int                 *  b_ind,
                                           const std::vector<double> &  well_log,
                                           std::vector<double>       &  blocked_log){
  //
  // Initialise arrays
  //

  blocked_log.resize(n_blocks_, 0.0f);
  std::vector<int> count(n_blocks_, 0);

  //
  // Block log
  //
  for (int m = first_M_ ; m < last_M_ + 1 ; m++) {
    if (well_log[m] != RMISSING) {
      blocked_log[b_ind[m]] += log(well_log[m]); //NBNB-PAL: Flytt denne logaritmen nedover...
      count[b_ind[m]]++;
      //LogKit::LogFormatted(LogKit::Low,"m=%d bInd[m]  log(wellLog[m])  %d  %.5f \n",m,bInd[m],log(wellLog[m]));
    }
  }
  for (int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0) {
      blocked_log[l] /= count[l];
      //LogKit::LogFormatted(LogKit::Low,"l=%d   count[l]=%d  sum=%.3f  blockedLog[l]=%.4f \n",l,count[l],sum, blockedLog[l]);
    }
    else
      blocked_log[l] = RMISSING;
  }

}

//------------------------------------------------------------------------------
void  BlockedLogsCommon::InterpolateContinuousLog(std::vector<double>   & blocked_log,
                                                  int                     start,
                                                  int                     end,
                                                  int                     index,
                                                  float                   rel)
{
  if(blocked_log[start] != RMISSING && blocked_log[end] != RMISSING && blocked_log[index] == RMISSING)
    blocked_log[index] = rel*blocked_log[end]+(1-rel)*blocked_log[start];
}
