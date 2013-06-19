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


BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well   * const well_data,
                                     const Simbox        * const estimation_simbox):
n_blocks_(0),
well_name_(""),
n_layers_(estimation_simbox->getnz()){
  BlockWell(well_data, estimation_simbox, continuous_log_names_, discrete_log_names_, continuous_logs_, discrete_logs_);
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
                                  std::vector<std::vector<int> >     & discrete_logs){
                                    
  well_name_ = well_data->GetWellName();

  std::map<std::string,std::vector<double> > continuous_logs_well = well_data->GetContLog();
  std::map<std::string,std::vector<int> > discrete_logs_well = well_data->GetDiscLog();

  int index_continuous = 0;
  for (std::map<std::string,std::vector<double> >::iterator it=continuous_logs_well.begin(); it!=continuous_logs_well.end(); ++it){
    continuous_log_names_.insert(std::pair<std::string, int>(it->first, index_continuous));
    continuous_logs.push_back(continuous_logs_well.find(it->first)->second);
    index_continuous++;
  }

  int index_discrete = 0;
  for (std::map<std::string,std::vector<int> >::iterator it=discrete_logs_well.begin(); it!=discrete_logs_well.end(); ++it){
    discrete_log_names_.insert(std::pair<std::string, int>(it->first, index_discrete));
    discrete_logs.push_back(discrete_logs_well.find(it->first)->second);
    index_discrete++;
  }

  // are all continuous logs of the same length?
  std::vector<int> bInd(continuous_logs[0].size()); // Gives which block each well log entry contributes to

  FindSizeAndBlockPointers(well_data, estimation_simbox, bInd);
  //findBlockIJK(well, stormgrid, bInd);

}


void  BlockedLogsCommon::FindSizeAndBlockPointers(const NRLib::Well    * const well_data,
                                                  const Simbox         * const estimation_simbox,
                                                  std::vector<int>     & bInd){
  int   nd = static_cast<int>(bInd.size());
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
    bInd[m] = IMISSING;
  }
  int new_I, new_J, new_K;
  int old_I = first_I;
  int old_J = first_J;
  int old_K = first_K;

  int n_defined_blocks = 0;
  bInd[first_M_] = first_K; // The first defined well log entry contributes to this block.

  //
  // The well positions used to be given in float rather than double. Unfortunately, this
  // allowed a well to oscillate between two or more cells, leading to a breakdown of the
  // algorithm below. To remedy for this we introduced array simboxInd which records the
  // indices of the simbox cells that are already accounted for, so that these are not
  // enlisted more than one time.
  //
  int * simboxInd = new int[nd];                                     // help hack
  const int nx    = estimation_simbox->getnx();                                 // help hack
  const int ny    = estimation_simbox->getny();                                 // help hack
  simboxInd[0] = nx*ny*old_K + nx*old_J + old_I;                        // help hack

  for (int m = first_M_ + 1 ; m < last_M_ + 1 ; m++) {
    estimation_simbox->getIndexes(x_pos[m], y_pos[m], tvd[m], new_I ,new_J, new_K);

    if (new_I != old_I || new_J != old_J || new_K != old_K) {

      int  thisInd = nx*ny*new_K + nx*new_J + new_I;                    // help hack
      bool blockNotListed = true;                                    // help hack
      for (int l = 0 ; l < n_defined_blocks ; l++) {                   // help hack
        if (thisInd == simboxInd[l]) {                               // help hack
          blockNotListed = false;                                    // help hack
          break;                                                     // help hack
        }                                                            // help hack
      }                                                              // help hack
      if (blockNotListed) {                                          // help hack
        simboxInd[n_defined_blocks+1] = thisInd;                       // help hack
        old_I = new_I;
        old_J = new_J;
        old_K = new_K;
        n_defined_blocks++;
      }
      else {
      }
    }
    bInd[m] = first_K + n_defined_blocks;
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
    LogKit::LogFormatted(LogKit::Low,"n_layers_                 = %d        \n",n_layers_);
    LogKit::LogFormatted(LogKit::Low,"first_I,first_J,first_K     = %d, %d, %d\n",first_I,first_J,first_K);
    LogKit::LogFormatted(LogKit::Low,"last_I,last_J,last_K        = %d, %d, %d\n",last_I,last_J,last_K);
    LogKit::LogFormatted(LogKit::Low,"n_defined_blocks, n_blocks_ = %d, %d    \n",n_defined_blocks,n_blocks_);
  }
  delete [] simboxInd;
}