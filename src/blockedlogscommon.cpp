/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#include "src/blockedlogscommon.h"
#include "src/simbox.h"
#include "fftw.h"
#include "lib/utils.h"

BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well                * well_data,
                                     const Simbox                     * const estimation_simbox,
                                     bool                               interpolate,
                                     bool                             & failed,
                                     std::string                      & err_text){
  n_angles_ = 0;
  well_name_ = well_data->GetWellName();
  n_layers_ = estimation_simbox->getnz();
  n_blocks_ = 0;
  interpolate_ = interpolate;

  // Get all continuous and discrete logs
  std::vector<std::string> cont_logs_to_be_blocked;
  std::vector<std::string> disc_logs_to_be_blocked;
  const std::map<std::string,std::vector<double> > cont_logs = well_data->GetContLog();
  const std::map<std::string,std::vector<int> > disc_logs = well_data->GetDiscLog();
  for(std::map<std::string,std::vector<double> >::const_iterator it = cont_logs.begin(); it!=cont_logs.end(); it++){
    cont_logs_to_be_blocked.push_back(it->first);
  }
  for(std::map<std::string,std::vector<int> >::const_iterator it = disc_logs.begin(); it!=disc_logs.end(); it++){
    disc_logs_to_be_blocked.push_back(it->first);
  }

  // FACIES
  if(well_data->HasFaciesLog()){
    facies_log_defined_ = true;
    facies_map_ = well_data->GetFaciesMap();
  }

  // Remove missing values

  RemoveMissingLogValues(well_data, x_pos_unblocked_, y_pos_unblocked_, z_pos_unblocked_, twt_unblocked_,
                         facies_unblocked_, continuous_logs_unblocked_, discrete_logs_unblocked_, cont_logs_to_be_blocked, 
                         disc_logs_to_be_blocked, n_data_, failed, err_text);

  if(failed)
    err_text += "Logs were not successfully read from well " + well_name_ +".\n";

  if (!failed)
    BlockWell(estimation_simbox, continuous_logs_unblocked_, discrete_logs_unblocked_, continuous_logs_blocked_, 
              discrete_logs_blocked_, n_data_, facies_log_defined_, facies_map_, interpolate, failed, err_text);

  n_continuous_logs_ = static_cast<int>(continuous_logs_blocked_.size());
  n_discrete_logs_ = static_cast<int>(discrete_logs_blocked_.size());

}

BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well                * well_data,
                                     const std::vector<std::string>   & cont_logs_to_be_blocked,
                                     const std::vector<std::string>   & disc_logs_to_be_blocked,
                                     const Simbox                     * const estimation_simbox,
                                     bool                               interpolate,
                                     bool                             & failed,
                                     std::string                      & err_text){
  n_angles_ = 0;
  well_name_ = well_data->GetWellName();
  n_layers_ = estimation_simbox->getnz();
  n_blocks_ = 0;
  interpolate_ = interpolate;

  // Get well name
  well_name_ = well_data->GetWellName();

  // 20130627 EN: Missing data are removed upon construction of a well_data object, whereas
  // NRLib::Well objects, which are used here, keep the logs as they are in the input files.
  RemoveMissingLogValues(well_data, x_pos_unblocked_, y_pos_unblocked_, z_pos_unblocked_, twt_unblocked_,
                         facies_unblocked_, continuous_logs_unblocked_, discrete_logs_unblocked_, cont_logs_to_be_blocked, 
                         disc_logs_to_be_blocked, n_data_, failed, err_text);
  if(failed)
    err_text += "Logs were not successfully read from well " + well_name_ +".\n";

  if (!failed)
    BlockWell(estimation_simbox, continuous_logs_unblocked_, discrete_logs_unblocked_, continuous_logs_blocked_, 
              discrete_logs_blocked_, n_data_, interpolate,facies_map_, facies_log_defined_, failed, err_text);

  n_continuous_logs_ = static_cast<int>(continuous_logs_blocked_.size());
  n_discrete_logs_ = static_cast<int>(discrete_logs_blocked_.size());
}

BlockedLogsCommon::~BlockedLogsCommon(){
  
}

void BlockedLogsCommon::BlockWell(const Simbox                                        * const estimation_simbox,
                                  const std::map<std::string, std::vector<double> >   & continuous_logs_unblocked,
                                  const std::map<std::string, std::vector<int> >      & discrete_logs_unblocked,
                                  std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                                  std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                                  unsigned int                                          n_data,
                                  bool                                                  facies_log_defined,
                                  const std::map<int, std::string>                    & facies_map,
                                  bool                                                  interpolate,
                                  bool                                                & failed,
                                  std::string                                         & err_text){

  std::vector<int> b_ind(n_data); // Gives which block each well log entry contributes to

  try{

    FindSizeAndBlockPointers(estimation_simbox, b_ind);
    FindBlockIJK(estimation_simbox, b_ind);

    // Coordinate logs and necessary logs

    BlockCoordinateLog(b_ind, x_pos_unblocked_, x_pos_blocked_);
    BlockCoordinateLog(b_ind, y_pos_unblocked_, y_pos_blocked_);
    BlockCoordinateLog(b_ind, z_pos_unblocked_, z_pos_blocked_);
    BlockContinuousLog(b_ind, twt_unblocked_, twt_blocked_    );

    // Continuous logs

    for(std::map<std::string, std::vector<double> >::const_iterator it = continuous_logs_unblocked.begin(); it!=continuous_logs_unblocked.end(); it++){
      std::vector<double> temp_vector_blocked;
      BlockContinuousLog(b_ind, it->second, temp_vector_blocked);
      continuous_logs_blocked.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector_blocked));
    }

    // Discrete logs
    if (facies_log_defined)
      BlockFaciesLog(b_ind, facies_unblocked_, facies_map, facies_map.size(), facies_blocked_);

    (void) discrete_logs_unblocked;
    (void) discrete_logs_blocked;

    if(interpolate){
      for(unsigned int i=1;i<n_data;i++) {
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
            // Coordinate logs
            InterpolateContinuousLog(x_pos_blocked_, start, end, j, t);
            InterpolateContinuousLog(y_pos_blocked_, start, end, j, t);
            InterpolateContinuousLog(z_pos_blocked_, start, end, j, t);
            InterpolateContinuousLog(twt_blocked_, start, end, j, t);

            // all blocked continuous logs
            for(std::map<std::string,std::vector<double> >::iterator it=continuous_logs_blocked_.begin(); it!=continuous_logs_blocked_.end(); ++it){
              InterpolateContinuousLog(it->second, start, end, j, t);
            }

          }
        }
      }
    }
  }catch(NRLib::Exception & e) {
    err_text += "Error: "+std::string(e.what());
    failed = true;
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
    for (unsigned int b = 0 ; b < n_blocks_ ; b++)
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
  for (unsigned int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0)
      blocked_coord[l] /= count[l];
    else
      blocked_coord[l]  = RMISSING;
  }

}

//------------------------------------------------------------------------------
void BlockedLogsCommon::BlockFaciesLog(const std::vector<int>              & b_ind,
                                       const std::vector<int>              & well_log,
                                       const std::map<int,std::string>     & facies_map,
                                       int                                   n_facies,
                                       std::vector<int>                    &  blocked_log)
{
  if (well_log.size() > 0) {
    //
    // Set undefined
    //
    std::vector<int> facies_numbers;
    for(std::map<int,std::string>::const_iterator it = facies_map.begin(); it != facies_map.end(); it++){
      facies_numbers.push_back(it->first);
    }
    //facies_numbers_ = facies_numbers;
    //n_facies_ = n_facies;
    blocked_log.resize(n_blocks_);
    for (unsigned int m = 0 ; m < n_blocks_ ; m++)
      blocked_log[m] = IMISSING;

    int   max_allowed_value = 100;  // Largest allowed value (facies number).
    std::vector<int> count(n_facies);
    std::vector<int> table(max_allowed_value);

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
    if(value!=IMISSING)
      count[table[value]]++;

    for (int m = first_M_+1 ; m < last_M_ + 1 ; m++) {
      if (b_ind[m] != b_ind[m - 1]) { // bInd[m] is the block number which sample 'm' lies in
        blocked_log[b_ind[m-1]] = FindMostProbable(count, n_facies, b_ind[m-1]);
        for (int i = 0 ; i < n_facies ; i++)
          count[i] = 0;
      }
    value = well_log[m];
    if(value!=IMISSING)
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
                                           const std::vector<double>  &  well_log,
                                           std::vector<double>        &  blocked_log){
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
  for (unsigned int l = 0 ; l < n_blocks_ ; l++) {
    if (count[l] > 0) {
      blocked_log[l] /= count[l];
      //LogKit::LogFormatted(LogKit::Low,"l=%d   count[l]=%d  sum=%.3f  blocked_log[l]=%.4f \n",l,count[l],sum, blocked_log[l]);
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

//------------------------------------------------------------------------------
void    BlockedLogsCommon::RemoveMissingLogValues(const NRLib::Well                            * well_data,
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
                                                  std::string                                  & err_text){

  // Get size of data vector including WELLMISSING data points
  unsigned int n_data_with_wellmissing = well_data->GetNData();
  // The number of legal data excluding WELLMISSING is initially set to 0
  n_data = 0;
  // Value from Well_data class
  double OPENWORKS_MISSING = -999.25;

  // Get the continuous and discrete logs from the Well object

  try{
    std::map<std::string,std::vector<double> > continuous_logs_well = well_data->GetContLog();
    std::map<std::string,std::vector<int> > discrete_logs_well = well_data->GetDiscLog();

    std::map<std::string, std::vector<double> >   continuous_logs_unblocked_temp;
    std::map<std::string, std::vector<int> >      discrete_logs_unblocked_temp;

    // Pick only the variables that are requested in the constructor --------------------------------------------------------

    // Find the continuous vectors in the wells that are to be blocked
    for(unsigned int i=0; i<cont_logs_to_be_blocked.size(); i++){
      std::map<std::string,std::vector<double> >::iterator it = continuous_logs_well.find(cont_logs_to_be_blocked[i]);
      // if the well log contains this continuous log
      if(it != continuous_logs_well.end()){
        continuous_logs_unblocked_temp.insert(std::pair<std::string, std::vector<double> >(it->first, it->second));
      }
    }

    // Find the discrete vectors in the wells that are to be blocked
    for(unsigned int i=0; i<disc_logs_to_be_blocked.size(); i++){
      std::map<std::string,std::vector<int> >::iterator it = discrete_logs_well.find(disc_logs_to_be_blocked[i]);
      // if the well log contains this discrete log
      if(it != discrete_logs_well.end()){
        discrete_logs_unblocked_temp.insert(std::pair<std::string, std::vector<int> >(it->first, it->second));
      }
    }

    for(std::map<std::string,std::vector<double> >::iterator it = continuous_logs_unblocked_temp.begin(); it!=continuous_logs_unblocked_temp.end(); it++){
      std::vector<double> temp_vector;
      continuous_logs_unblocked.insert(std::pair<std::string, std::vector<double> >(it->first, temp_vector));
    }

    for(std::map<std::string,std::vector<int> >::iterator it = discrete_logs_unblocked_temp.begin(); it!=discrete_logs_unblocked_temp.end(); it++){
      std::vector<int> temp_vector;
      discrete_logs_unblocked.insert(std::pair<std::string, std::vector<int> >(it->first, temp_vector));
    }


    // Remove WELLMISSING data from wells
    for(unsigned int i=0; i<n_data_with_wellmissing; i++){
      // TWT is the variable we are testing for WELLMISSING values
      double dummy = continuous_logs_unblocked_temp.find("TWT")->second[i];
      if(dummy != WELLMISSING && dummy != OPENWORKS_MISSING){
        x_pos_unblocked.push_back(continuous_logs_well.find("X_pos")->second[i]);
        y_pos_unblocked.push_back(continuous_logs_well.find("Y_pos")->second[i]);
        z_pos_unblocked.push_back(continuous_logs_well.find("TVD")->second[i]);
        twt_unblocked.push_back(continuous_logs_well.find("TWT")->second[i]);
        if(facies_log_defined_)
          facies_unblocked.push_back(discrete_logs_well.find("Facies")->second[i]);
      
        // Loop over continuous variables and push back this element
        for(std::map<std::string,std::vector<double> >::iterator it = continuous_logs_unblocked_temp.begin(); it!=continuous_logs_unblocked_temp.end(); it++){
          continuous_logs_unblocked.find(it->first)->second.push_back(it->second[i]);
        }

        // Loop over discrete variables and push back this element
        for(std::map<std::string,std::vector<int> >::iterator it = discrete_logs_unblocked_temp.begin(); it!=discrete_logs_unblocked_temp.end(); it++){
          discrete_logs_unblocked.find(it->first)->second.push_back(it->second[i]);
        }

        n_data++;
      }
    }
  }catch(NRLib::Exception & e) {
    err_text += "Error: "+std::string(e.what());
    failed = true;
  }

}

void BlockedLogsCommon::FindOptimalWellLocation(std::vector<SeismicStorage>   & seismic_data,
                                                const Simbox                                        * estimation_simbox,
                                                float                                        ** refl_coef,
                                                int                                             n_angles,
                                                const std::vector<float>                      & angle_weight,
                                                float                                           max_shift,
                                                int                                             i_max_offset,
                                                int                                             j_max_offset,
                                                const std::vector<Surface *>                    limits,
                                                int                                           & i_move,
                                                int                                           & j_move,
                                                float                                         & k_move){
  int   polarity;
  int   i,j,k,l,m;
  int   start, length;
  float sum;
  float shift_F;
  float max_tot;
  float f1,f2,f3;

  int nx            = seismic_data[0].GetNx();
  int ny            = seismic_data[0].GetNy();
  int nzp           = seismic_data[0].GetNz();
  int cnzp          = nzp/2+1;
  int rnzp          = 2*cnzp;
  int i_tot_offset    = 2*i_max_offset+1;
  int j_tot_offset    = 2*j_max_offset+1;
  int polarityMax   = 0;
  float shift       = 0.0f;
  float max_value_tot = 0;
  float total_weight = 0;
  float dz          = static_cast<float>(estimation_simbox->getdz());

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

  for( i=0; i<n_angles; i++ ){
    max_value_max[i] = 0.0f;
    shift_I_max[i]   = 0;
  }

  // make offset vectors
  for(i=-i_max_offset; i<i_max_offset+1; i++){
    i_offset[i+i_max_offset]=i;
  }
  for(j=-j_max_offset; j<j_max_offset+1; j++){
    j_offset[j+j_max_offset]=j;
  }

  GetVerticalTrendLimited(GetVpUnblocked(), vp_vert, limits);
  GetVerticalTrendLimited(GetVsUnblocked(), vs_vert, limits);
  GetVerticalTrendLimited(GetRhoUnblocked(), rho_vert, limits);

  std::vector<bool> has_data(n_layers_);
  for(i = 0 ; i < n_layers_ ; i++) {
    has_data[i] = vp_vert[i] != RMISSING && vs_vert[i] != RMISSING && rho_vert[i] != RMISSING;
  }
  FindContinuousPartOfData(has_data,n_layers_,start,length);

  for( j=0; j<n_angles; j++ ){
    seis_r[j]              = new fftw_real[rnzp];
    cpp_r[j]               = new fftw_real[rnzp];
    cor_cpp_r[j]           = new fftw_real[rnzp];
    ccor_seis_cpp_r[j]     = new fftw_real[rnzp];
    ccor_seis_cpp_Max_r[j] = new fftw_real[rnzp];
  }

  // Calculate reflection coefficients
  for( j=0; j<n_angles; j++ ){
    for(i=0; i<rnzp; i++){
      cpp_r[j][i] = 0;
    }
    FillInCpp(refl_coef[j],start,length,cpp_r[j],nzp);
    Utils::fft(cpp_r[j],cpp_c[j],nzp);
    EstimateCor(cpp_c[j],cpp_c[j],cor_cpp_c[j],cnzp);
    Utils::fftInv(cor_cpp_c[j],cor_cpp_r[j],nzp);
  }

  std::vector<NRLib::Grid<float> > seis_cube_small(n_angles,NRLib::Grid<float> (i_tot_offset,j_tot_offset,n_blocks_));

  for (j = 0 ; j < n_angles ; j++)
  {
    //seismic_data[j]->setAccessMode(FFTGrid::RANDOMACCESS);
    for (k = 0; k < i_tot_offset; k++)
    {
      for (l = 0; l < j_tot_offset; l++)
      {
        GetBlockedGrid(&seismic_data[j], estimation_simbox, seis_log,i_offset[k],j_offset[l]);
        for (m = 0; m < static_cast<int>(n_blocks_); m++)
        {
          seis_cube_small[j](k, l, m) = static_cast<float>(seis_log[m]);
        }
      }
    }
    //seis_cube[j]->endAccess();
  }

  // Loop through possible well locations
  for(k=0; k<i_tot_offset; k++){
    if(i_pos_[0]+i_offset[k]<0 || i_pos_[0]+i_offset[k]>nx-1) //Check if position is within seismic range
      continue;

    for(l=0; l<j_tot_offset; l++){
      if(j_pos_[0]+j_offset[l]<0 || j_pos_[0]+j_offset[l]>ny-1) //Check if position is within seismic range
        continue;

      for( j=0; j<n_angles; j++ ){

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
      dz = static_cast<float>(estimation_simbox->getRelThick(i_pos_[0]+i_offset[k],j_pos_[0]+j_offset[l])*estimation_simbox->getdz());
      sum = 0;
      for( j=0; j<n_angles; j++ ){
        if(angle_weight[j] > 0){
          for(i=0;i<ceil(max_shift/dz);i++)//zero included
            sum+=ccor_seis_cpp_r[j][i];
          for(i=0;i<floor(max_shift/dz);i++)
            sum+=ccor_seis_cpp_r[j][nzp-i-1];
        }
      }
      polarity=-1;
      if(sum > 0)
        polarity=1;

      // Find maximum correlation and corresponding shift for each angle
      max_tot = 0.0;
      for( j=0; j<n_angles; j++ ){
        if(angle_weight[j]>0){
          max_value[j] = 0.0f;
          shift_I[j]=0;
          for(i=0;i<ceil(max_shift/dz);i++){
            if(ccor_seis_cpp_r[j][i]*polarity > max_value[j]){
              max_value[j] = ccor_seis_cpp_r[j][i]*polarity;
              shift_I[j] = i;
            }
          }
          for(i=0;i<floor(max_shift/dz);i++){
            if(ccor_seis_cpp_r[j][nzp-1-i]*polarity > max_value[j]){
              max_value[j] = ccor_seis_cpp_r[j][nzp-1-i]*polarity;
              shift_I[j] = -1-i;
            }
          }
          max_tot += angle_weight[j]*max_value[j]; //Find weighted total maximum correlation
        }
      }

      if(max_tot > max_value_tot){
        max_value_tot = max_tot;
        polarityMax = polarity;
        i_move       = i_offset[k];
        j_move       = j_offset[l];
        for(m=0; m<n_angles; m++){
          shift_I_max[m]   = shift_I[m];
          max_value_max[m] = max_value[m];
          for(i=0;i<rnzp;i++)
            ccor_seis_cpp_Max_r[m][i] = ccor_seis_cpp_r[m][i];
        }
      }
    }
  }

  for(i=0; i<n_angles; i++){
    shift_I[i] = shift_I_max[i];
    max_value[i] = max_value_max[i];
    for(j=0;j<rnzp;j++)
      ccor_seis_cpp_r[i][j] = ccor_seis_cpp_Max_r[i][j];
  }
  polarity = polarityMax;

  // Find kMove in optimal location
  for(j=0; j<n_angles; j++){
    if(angle_weight[j]>0){
      if(shift_I[j] < 0){
        if(ccor_seis_cpp_r[j][nzp+shift_I[j]-1]*polarity < max_value[j]) //then local max
        {
          f1 = ccor_seis_cpp_r[j][nzp+shift_I[j]-1];
          f2 = ccor_seis_cpp_r[j][nzp+shift_I[j]];
          int ind3;
          if(shift_I[j]==-1)
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
        if(ccor_seis_cpp_r[j][shift_I[j]+1]*polarity < max_value[j]) //then local max
        {
          f3 = ccor_seis_cpp_r[j][shift_I[j]+1];
          f2 = ccor_seis_cpp_r[j][shift_I[j]];
          int ind1;
          if(shift_I[j]==0)
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

  for( j=0; j<n_angles; j++ ){
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

void BlockedLogsCommon::GetVerticalTrendLimited(const std::vector<double>          & blocked_log, 
                                                std::vector<double>                & trend, 
                                                const std::vector<Surface *>       & limits){
  if (blocked_log.size() > 0 && trend.size() > 0) {
    std::vector<int> count(n_layers_);
    for (int k = 0 ; k < n_layers_ ; k++) {
      trend[k] = 0.0;
      count[k] = 0;
    }
    for (int m = 0 ; m < static_cast<int>(n_blocks_) ; m++) {
      if (blocked_log[m] != RMISSING) {
        if(limits.size() == 0 ||
           (limits[0]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) <= z_pos_blocked_[m] &&
            limits[1]->GetZ(x_pos_blocked_[m],y_pos_blocked_[m]) >= z_pos_blocked_[m])) {
          trend[k_pos_[m]] += blocked_log[m];
          count[k_pos_[m]]++;
        }
      }
    }
    for (int k = 0 ; k < n_layers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
    if(interpolate_)
      InterpolateTrend(blocked_log,trend,limits);

  }
  else {
    if (blocked_log.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrendLimited(): Trying to use an undefined blocked log\n");
    if (trend.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrendLimited(): Trying to use an undefined trend\n");
    exit(1);
  }
}

void
BlockedLogsCommon::InterpolateTrend(const double  * blocked_log, 
                                    double        * trend)
{
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if(abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if(k_pos_[m] < k_pos_[m-1])
        delta = -1;
      double step_mult = static_cast<double>(delta)/static_cast<double>(k_pos_[m]-k_pos_[m-1]);
      double t = step_mult;
      for(int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if(trend[j] == RMISSING)
          trend[j] = t*blocked_log[m]+(1-t)*blocked_log[m-1];
        t += step_mult;
      }
    }
  }
}

void         BlockedLogsCommon::InterpolateTrend(const std::vector<double>    & blocked_log, 
                                                 std::vector<double>          & trend){
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if(abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if(k_pos_[m] < k_pos_[m-1])
        delta = -1;
      for(int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if(trend[j] == RMISSING)
          trend[j] = blocked_log[m-1];
      }
    }
  }
}

void  BlockedLogsCommon::InterpolateTrend(const std::vector<double>      & blocked_log, 
                                          std::vector<double>            & trend, 
                                          const std::vector<Surface *>   & limits){
  for (int m = 1 ; m < static_cast<int>(n_blocks_) ; m++) {
    if(abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if(k_pos_[m] < k_pos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(k_pos_[m]-k_pos_[m-1]);
      float t = step_mult;
      for(int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if(trend[j] == RMISSING) {
          if(limits.size() == 0 ||
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
                                                 int                     & length) const{
  int  i;
  int  l_pice            =  0;
  int  length_max_pice    = -1;
  int  start_longest_pice =  0;
  bool previous_had_data  = false;

  for(i = 0; i < nz ;i++){
    if(hasData[i]){
      if(! previous_had_data)
        l_pice=1;
      else
        l_pice++;
      previous_had_data = true;
    }
    else{
      if(previous_had_data){
        if(length_max_pice < l_pice){
          length_max_pice  = l_pice;
          start_longest_pice = i-l_pice;
        }
      }
      previous_had_data=false;
    }
  }

  if(previous_had_data){
    if(length_max_pice < l_pice){
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
                                  int           nzp){
  int i;

  for(i=0;i<nzp;i++)
    cpp_r[i]=0;

  std::vector<double> vp_vert(n_layers_);
  std::vector<double> vs_vert(n_layers_);
  std::vector<double> rho_vert(n_layers_);

  GetVerticalTrend(GetVpBlocked(), vp_vert);
  GetVerticalTrend(GetVsBlocked(), vs_vert);
  GetVerticalTrend(GetRhoBlocked(), rho_vert);

  for(i=start;i < start+length-1;i++)
  {
    double ei1 = ComputeElasticImpedance(vp_vert[i], static_cast<float>(vs_vert[i]),static_cast<float>(rho_vert[i]),coeff);
    double ei2 = ComputeElasticImpedance(vp_vert[i+1], static_cast<float>(vs_vert[i+1]),static_cast<float>(rho_vert[i+1]),coeff);
    cpp_r[i] =  static_cast<fftw_real>(ei2-ei1);
  }

}

double BlockedLogsCommon::ComputeElasticImpedance(double         vp,
                                                 float         vs,
                                                 float         rho,
                                                 const float * coeff) const
{
  // vp, vs, rho are logtransformed
  double ang_imp;

  ang_imp = coeff[0]*vp+coeff[1]*vs+coeff[2]*rho;

  return ang_imp;
}



void     BlockedLogsCommon::GetVerticalTrend(const std::vector<double>  & blocked_log, 
                                             std::vector<double>        & trend){
  if (blocked_log.size() > 0 && trend.size() > 0) {
    std::vector<double> count(n_layers_);
    for (int k = 0 ; k < n_layers_ ; k++) {
      trend[k] = 0.0;
      count[k] = 0;
    }
    for (int m = 0 ; m < static_cast<int>(n_blocks_) ; m++) {
      if (blocked_log[m] != RMISSING) {
        trend[k_pos_[m]] += blocked_log[m];
        count[k_pos_[m]]++;
      }
    }
    for (int k = 0 ; k < n_layers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
    if(interpolate_ == true)
      InterpolateTrend(blocked_log,trend);

  }
  else {
    if (blocked_log.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrend(): Trying to use an undefined blocked log\n");
    if (trend.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsCommon::GetVerticalTrend(): Trying to use an undefined trend\n");
    exit(1);
  }
}

void BlockedLogsCommon::EstimateCor(fftw_complex * var1_c,
                                    fftw_complex * var2_c,
                                    fftw_complex * ccor_1_2_c,
                                    int            cnzp) const{
  for(int i=0;i<cnzp;i++){
    ccor_1_2_c[i].re =  var1_c[i].re*var2_c[i].re + var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }
}

void BlockedLogsCommon::GetBlockedGrid(const SeismicStorage   * grid,
                                       const Simbox           * estimation_simbox,
                                       std::vector<double>    & blocked_log,
                                       int                      i_offset,
                                       int                      j_offset){
  for (int m = 0 ; m < static_cast<float>(n_blocks_) ; m++) {
    //LogKit::LogFormatted(LogKit::Low,"m=%d  ipos_[m], jpos_[m], kpos_[m] = %d %d %d\n",m,ipos_[m], jpos_[m], kpos_[m]);
    blocked_log[m] = grid->GetRealTraceValue(estimation_simbox, i_pos_[m]+i_offset, j_pos_[m]+j_offset, k_pos_[m]);

  }
}


void BlockedLogsCommon::FillInSeismic(std::vector<double>   & seismic_data,
                                      int                     start,
                                      int                     length,
                                      fftw_real             * seis_r,
                                      int                     nzp) const{
  int i;
  for(i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for(i=start; i<start+length; i++)
    seis_r[i] = static_cast<fftw_real>(seismic_data[i]);
/*
  int lTregion = 3;
  int* modify  = getIndexPrior(start,lTregion,nzp);
  int* conditionto = getIndexPost(start-1,lTregion,nzp);
  //NBNB Odd: interpolate endpoints?
*/

}

void BlockedLogsCommon::SetLogFromVerticalTrend(std::vector<double>   & blocked_log,
                                                std::vector<double>   & z_pos,
                                                int                     n_blocks,
                                                std::vector<double>   & vertical_trend,
                                                double                  z0,
                                                double                  dzVt,
                                                int                     nz)
{
  //
  // Initialise as undefined
  //
  for (int i=0 ; i<n_blocks ; i++)
    blocked_log[i] = RMISSING;

  //
  // Aritmethic mean of values in overlapping cells
  //
  for (int i=0 ; i<n_blocks ; i++) {
    double dz;
    if (i==n_blocks-1)
      dz = z_pos[i]-z_pos[i-1];
    else
      dz = z_pos[i+1]-z_pos[i];
    double zi = z_pos[i];
    double a  = zi - 0.5*dz;     // Top of blocked log cell
    double b  = z0 + 0.5*dzVt;   // Base of first vertical trend cell

    int j=0;
    while (b<a && j<nz) {
      b += dzVt;
      j++;
    }
    // Now 'j' is the first vertical trend cell overlapping blocked log cell 'i'

    float value;
    if (j==nz) {
      // We have come to the end of the end-of-vertical-trend
      value = static_cast<float>(vertical_trend[j-1]);
    }
    else if (b >= a+dz) {
      // One single vertical-trend-cell covers a blocked log cell completely.
      value = static_cast<float>(vertical_trend[j]);
    }
    else {
      double zj = b + 0.5*dzVt; // Center of vertical trend cell
      value = static_cast<float>( vertical_trend[j]* (zj+dzVt-zi)/dzVt + vertical_trend[j]*(zi-zj)/dzVt);
    }
    blocked_log[i] = value;

    //LogKit::LogFormatted(LogKit::Error,"i j  blockedLog[i]   %d %d  %7.3f\n",i,j,blockedLog[i]);
  }
}

void  BlockedLogsCommon::SetLogFromVerticalTrend(std::vector<double>    & vertical_trend,
                                                 double                   z0,              // z-value of center in top layer
                                                 double                   dz,              // dz in vertical trend
                                                 int                      nz,              // layers in vertical trend
                                                 std::string              type,
                                                 int                      i_angle)
{
  if (type != "WELL_SYNTHETIC_SEISMIC")
  {
    std::vector<double> blocked_log(n_blocks_);

    SetLogFromVerticalTrend(blocked_log, z_pos_blocked_, n_blocks_,
                            vertical_trend, z0, dz, nz);

    if (type == "ALPHA_SEISMIC_RESOLUTION")
      cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Vp", blocked_log));
    else if (type == "BETA_SEISMIC_RESOLUTION")
      cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Vs", blocked_log));
    else if (type == "RHO_SEISMIC_RESOLUTION")
      cont_logs_seismic_resolution_.insert(std::pair<std::string, std::vector<double> >("Rho", blocked_log));
    else if (type == "ACTUAL_SYNTHETIC_SEISMIC") {
      if (actual_synt_seismic_data_.size() == 0)
        actual_synt_seismic_data_.resize(n_angles_); // nAngles is set along with real_seismic_data_
      actual_synt_seismic_data_[i_angle] = blocked_log;
    }
    else {
      LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type+
                           "\" in BlockedLogs::setLogFromVerticalTrend()\n");
      exit(1);
    }
  }
  else if (type == "WELL_SYNTHETIC_SEISMIC") {
    if (well_synt_seismic_data_.size() == 0)
    {
      well_synt_seismic_data_.resize(n_angles_);
      for (int i=0; i<n_angles_; i++)
      {
        well_synt_seismic_data_[i].resize(n_blocks_);
        for (int j=0; j<static_cast<int>(n_blocks_); j++)
          well_synt_seismic_data_[i][j] = RMISSING; //Declare in case the wavelet is not estimated for all angles
      }
    }
    SetLogFromVerticalTrend(well_synt_seismic_data_[i_angle], z_pos_blocked_, n_blocks_,
                            vertical_trend, z0, dz, nz);
  }
}

//------------------------------------------------------------------------------
int BlockedLogsCommon::FindMostProbable(const std::vector<int>  & count,
                                        int                       n_facies,
                                        int                       block_index)
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

/*
void  BlockedLogs::SetLogFromGrid(SeismicStorage    * grid,
                                  int                 i_angle,
                                  int                 n_angles,
                                  std::string         type)
{
  std::vector<double> blocked_log(n_blocks_);

  for (int m = 0 ; m < n_blocks_ ; m++) {
    grid->GetRealTraceValue(
    blocked_log[m] = grid->getRealValue(ipos_[m], jpos_[m], kpos_[m]);
  }

  if (nAngles_ == 0)
    nAngles_ = nAngles;

  if (type == "REFLECTION_COEFFICIENT") {
    if (cpp_ == NULL)
      cpp_ = new float * [nAngles_];
    cpp_[iAngle] = blockedLog;
  }
  else if (type == "SEISMIC_DATA") {
    if (real_seismic_data_ == NULL)
      real_seismic_data_ = new float * [nAngles_];
    real_seismic_data_[iAngle] = blockedLog;
  }
  else if (type == "FACIES_PROB") {
    if (facies_prob_ == NULL)
      facies_prob_ = new float * [nFacies_];
    facies_prob_[iAngle] = blockedLog;
  }
  else if (type == "ALPHA_PREDICTED")
    alpha_predicted_ = blockedLog;
  else if (type == "BETA_PREDICTED")
    beta_predicted_ = blockedLog;
  else if (type == "RHO_PREDICTED")
    rho_predicted_ = blockedLog;

  else {
    LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type
                         +"\" in BlockedLogs::setLogFromGrid()\n");
    exit(1);
  }
}
*/
