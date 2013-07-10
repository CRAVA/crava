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


BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well                * const well_data,
                                     const std::vector<std::string>   & cont_logs_to_be_blocked,
                                     const std::vector<std::string>   & disc_logs_to_be_blocked,
                                     const Simbox                     * const estimation_simbox,
                                     bool                               interpolate,
                                     bool                             & failed,
                                     std::string                      & err_text):
n_blocks_(0),
well_name_(""),
n_layers_(estimation_simbox->getnz()){

  // Get well name
  well_name_ = well_data->GetWellName();

  // 20130627 EN: Missing data are removed upon construction of a well_data object, whereas
  // NRLib::Well objects, which are used here, keep the logs as they are in the input files.
  RemoveMissingLogValues(well_data, x_pos_unblocked_, y_pos_unblocked_, z_pos_unblocked_, twt_unblocked_,
                         continuous_logs_unblocked_, discrete_logs_unblocked_, cont_logs_to_be_blocked, 
                         disc_logs_to_be_blocked, n_data_, failed, err_text);
  if(failed)
    err_text += "Logs were not successfully read from well " + well_name_ +".\n";

  if (!failed)
    BlockWell(estimation_simbox, continuous_logs_unblocked_, discrete_logs_unblocked_,
             continuous_logs_blocked_, discrete_logs_blocked_, n_data_, interpolate, failed, err_text);

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

    // Discrete logs (20130625 EN: No blocking of discrete logs for now)

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

//------------------------------------------------------------------------------
void    BlockedLogsCommon::RemoveMissingLogValues(const NRLib::Well                            * const well_data,
                                                  std::vector<double>                          & x_pos_unblocked,
                                                  std::vector<double>                          & y_pos_unblocked,
                                                  std::vector<double>                          & z_pos_unblocked,
                                                  std::vector<double>                          & twt_unblocked,
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

void BlockedLogsCommon::FindOptimalWellLocation(std::map<int, std::vector<SeismicStorage> >   seismic_data,
                                                Simbox                     * timeSimbox,
                                                float                     ** reflCoef,
                                                int                          nAngles,
                                                const std::vector<float>   & angleWeight,
                                                float                        maxShift,
                                                int                          iMaxOffset,
                                                int                          jMaxOffset,
                                                const std::vector<Surface *> limits,
                                                int                        & iMove,
                                                int                        & jMove,
                                                float                      & kMove){
  int   polarity;
  int   i,j,k,l,m;
  int   start,length;
  float sum;
  float shiftF;
  float maxTot;
  float f1,f2,f3;
  seismic_data.find(0)->second[0];
  int nx            = seisCube[0]->getNx();
  int ny            = seisCube[0]->getNy();
  int nzp           = seisCube[0]->getNzp();
  int cnzp          = nzp/2+1;
  int rnzp          = 2*cnzp;
  int iTotOffset    = 2*iMaxOffset+1;
  int jTotOffset    = 2*jMaxOffset+1;
  int polarityMax   = 0;
  float shift       = 0.0f;
  float maxValueTot = 0;
  float totalWeight = 0;
  float dz          = static_cast<float>(timeSimbox->getdz());

  float  * seisLog   = new float[n_blocks_];
  float  * seisData  = new float[n_layers_];
  float  * alphaVert = new float[n_layers_];
  float  * betaVert  = new float[n_layers_];
  float  * rhoVert   = new float[n_layers_];


  std::vector<int>   iOffset(iTotOffset);
  std::vector<int>   jOffset(jTotOffset);
  std::vector<int>   shiftI(nAngles);
  std::vector<int>   shiftIMax(nAngles);
  std::vector<float> maxValue(nAngles);
  std::vector<float> maxValueMax(nAngles);

  fftw_real    ** cpp_r               = new fftw_real*[nAngles];
  fftw_complex ** cpp_c               = reinterpret_cast<fftw_complex**>(cpp_r);

  fftw_real    ** cor_cpp_r           = new fftw_real*[nAngles];
  fftw_complex ** cor_cpp_c           = reinterpret_cast<fftw_complex**>(cor_cpp_r);

  fftw_real    ** seis_r              = new fftw_real*[nAngles];
  fftw_complex ** seis_c              = reinterpret_cast<fftw_complex**>(seis_r);

  fftw_real    ** ccor_seis_cpp_r     = new fftw_real*[nAngles];
  fftw_complex ** ccor_seis_cpp_c     = reinterpret_cast<fftw_complex**>(ccor_seis_cpp_r);

  fftw_real    ** ccor_seis_cpp_Max_r = new fftw_real*[nAngles];

  for( i=0; i<nAngles; i++ ){
    maxValueMax[i] = 0.0f;
    shiftIMax[i]   = 0;
  }

  // make offset vectors
  for(i=-iMaxOffset; i<iMaxOffset+1; i++){
    iOffset[i+iMaxOffset]=i;
  }
  for(j=-jMaxOffset; j<jMaxOffset+1; j++){
    jOffset[j+jMaxOffset]=j;
  }

  getVerticalTrendLimited(alpha_, alphaVert, limits);
  getVerticalTrendLimited(beta_, betaVert, limits);
  getVerticalTrendLimited(rho_, rhoVert, limits);

  std::vector<bool> hasData(nLayers_);
  for(i = 0 ; i < nLayers_ ; i++) {
    hasData[i] = alphaVert[i] != RMISSING && betaVert[i] != RMISSING && rhoVert[i] != RMISSING;
  }
  findContiniousPartOfData(hasData,nLayers_,start,length);

  for( j=0; j<nAngles; j++ ){
    seis_r[j]              = new fftw_real[rnzp];
    cpp_r[j]               = new fftw_real[rnzp];
    cor_cpp_r[j]           = new fftw_real[rnzp];
    ccor_seis_cpp_r[j]     = new fftw_real[rnzp];
    ccor_seis_cpp_Max_r[j] = new fftw_real[rnzp];
  }

  // Calculate reflection coefficients
  for( j=0; j<nAngles; j++ ){
    for(i=0; i<rnzp; i++){
      cpp_r[j][i] = 0;
    }
    fillInCpp(reflCoef[j],start,length,cpp_r[j],nzp);
    Utils::fft(cpp_r[j],cpp_c[j],nzp);
    estimateCor(cpp_c[j],cpp_c[j],cor_cpp_c[j],cnzp);
    Utils::fftInv(cor_cpp_c[j],cor_cpp_r[j],nzp);
  }

  std::vector<NRLib::Grid<float> > seisCubeSmall(nAngles,NRLib::Grid<float> (iTotOffset,jTotOffset,nBlocks_));

  for (j = 0 ; j < nAngles ; j++)
  {
    seisCube[j]->setAccessMode(FFTGrid::RANDOMACCESS);
    for (k = 0; k < iTotOffset; k++)
    {
      for (l = 0; l < jTotOffset; l++)
      {
        getBlockedGrid(seisCube[j],seisLog,iOffset[k],jOffset[l]);
        for (m = 0; m < nBlocks_; m++)
        {
          seisCubeSmall[j](k, l, m) = seisLog[m];
        }
      }
    }
    seisCube[j]->endAccess();
  }

  // Loop through possible well locations
  for(k=0; k<iTotOffset; k++){
    if(ipos_[0]+iOffset[k]<0 || ipos_[0]+iOffset[k]>nx-1) //Check if position is within seismic range
      continue;

    for(l=0; l<jTotOffset; l++){
      if(jpos_[0]+jOffset[l]<0 || jpos_[0]+jOffset[l]>ny-1) //Check if position is within seismic range
        continue;

      for( j=0; j<nAngles; j++ ){

        for (m=0; m<nBlocks_; m++)
          seisLog[m] = seisCubeSmall[j](k,l,m);

        getVerticalTrend(seisLog, seisData);
        fillInSeismic(seisData,start,length,seis_r[j],nzp);

        Utils::fft(seis_r[j],seis_c[j],nzp);
        estimateCor(seis_c[j],cpp_c[j],ccor_seis_cpp_c[j],cnzp);
        Utils::fftInv(ccor_seis_cpp_c[j],ccor_seis_cpp_r[j],nzp);
      }

      // if the sum from -maxShift to maxShift ms is
      // positive then polarity is positive
      dz = static_cast<float>(timeSimbox->getRelThick(ipos_[0]+iOffset[k],jpos_[0]+jOffset[l])*timeSimbox->getdz());
      sum = 0;
      for( j=0; j<nAngles; j++ ){
        if(angleWeight[j] > 0){
          for(i=0;i<ceil(maxShift/dz);i++)//zero included
            sum+=ccor_seis_cpp_r[j][i];
          for(i=0;i<floor(maxShift/dz);i++)
            sum+=ccor_seis_cpp_r[j][nzp-i-1];
        }
      }
      polarity=-1;
      if(sum > 0)
        polarity=1;

      // Find maximum correlation and corresponding shift for each angle
      maxTot = 0.0;
      for( j=0; j<nAngles; j++ ){
        if(angleWeight[j]>0){
          maxValue[j] = 0.0f;
          shiftI[j]=0;
          for(i=0;i<ceil(maxShift/dz);i++){
            if(ccor_seis_cpp_r[j][i]*polarity > maxValue[j]){
              maxValue[j] = ccor_seis_cpp_r[j][i]*polarity;
              shiftI[j] = i;
            }
          }
          for(i=0;i<floor(maxShift/dz);i++){
            if(ccor_seis_cpp_r[j][nzp-1-i]*polarity > maxValue[j]){
              maxValue[j] = ccor_seis_cpp_r[j][nzp-1-i]*polarity;
              shiftI[j] = -1-i;
            }
          }
          maxTot += angleWeight[j]*maxValue[j]; //Find weighted total maximum correlation
        }
      }

      if(maxTot > maxValueTot){
        maxValueTot = maxTot;
        polarityMax = polarity;
        iMove       = iOffset[k];
        jMove       = jOffset[l];
        for(m=0; m<nAngles; m++){
          shiftIMax[m]   = shiftI[m];
          maxValueMax[m] = maxValue[m];
          for(i=0;i<rnzp;i++)
            ccor_seis_cpp_Max_r[m][i] = ccor_seis_cpp_r[m][i];
        }
      }
    }
  }

  for(i=0; i<nAngles; i++){
    shiftI[i] = shiftIMax[i];
    maxValue[i] = maxValueMax[i];
    for(j=0;j<rnzp;j++)
      ccor_seis_cpp_r[i][j] = ccor_seis_cpp_Max_r[i][j];
  }
  polarity = polarityMax;

  // Find kMove in optimal location
  for(j=0; j<nAngles; j++){
    if(angleWeight[j]>0){
      if(shiftI[j] < 0){
        if(ccor_seis_cpp_r[j][nzp+shiftI[j]-1]*polarity < maxValue[j]) //then local max
        {
          f1 = ccor_seis_cpp_r[j][nzp+shiftI[j]-1];
          f2 = ccor_seis_cpp_r[j][nzp+shiftI[j]];
          int ind3;
          if(shiftI[j]==-1)
            ind3 = 0;
          else
            ind3=nzp+shiftI[j]+1;
          f3 = ccor_seis_cpp_r[j][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI[j]+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI[j]);
      }
      else //positive or zero shift
      {
        if(ccor_seis_cpp_r[j][shiftI[j]+1]*polarity < maxValue[j]) //then local max
        {
          f3 = ccor_seis_cpp_r[j][shiftI[j]+1];
          f2 = ccor_seis_cpp_r[j][shiftI[j]];
          int ind1;
          if(shiftI[j]==0)
            ind1 = nzp-1;
          else
            ind1=shiftI[j]-1;
          f1 = ccor_seis_cpp_r[j][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI[j]+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI[j]);
      }
      shift += angleWeight[j]*shiftF*dz;//weigthing shift according to wellWeight
      totalWeight += angleWeight[j];
    }
  }

  shift/=totalWeight;
  kMove = shift;

  for( j=0; j<nAngles; j++ ){
    delete [] ccor_seis_cpp_Max_r[j];
    delete [] ccor_seis_cpp_r[j];
    delete [] cor_cpp_r[j];
    delete [] seis_r[j];
    delete [] cpp_r[j];
  }
  delete [] alphaVert;
  delete [] betaVert;
  delete [] rhoVert;
  delete [] seisLog;
  delete [] seisData;

  delete [] ccor_seis_cpp_Max_r;
  delete [] ccor_seis_cpp_r;
  delete [] cor_cpp_r;
  delete [] seis_r;
  delete [] cpp_r;
}
