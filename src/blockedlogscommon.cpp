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
//#include "src/fftgrid.h"
#include "nrlib/flens/nrlib_flens.hpp"
//#include "src/seismicstorage.h"


BlockedLogsCommon::BlockedLogsCommon(const NRLib::Well    * const well_data,
                                     const Simbox         * const estimation_simbox,
                                     bool                   interpolate):
n_blocks_(0),
well_name_(""),
interpolate_(interpolate),
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


void BlockedLogsCommon::SetSeismicGradient(double                            v0,
                                           const NRLib::Grid2D<float>   &    structure_depth_grad_x,
                                           const NRLib::Grid2D<float>   &    structure_depth_grad_y,
                                           const NRLib::Grid2D<float>   &    ref_time_grad_x,
                                           const NRLib::Grid2D<float>   &    ref_time_grad_y,
                                           std::vector<double>          &    x_gradient,
                                           std::vector<double>          &    y_gradient)
{
  x_gradient.resize(n_blocks_);
  y_gradient.resize(n_blocks_);

  double mp= 2.0/(v0*0.001);
  for(int k = 0; k < n_blocks_; k++){
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

void BlockedLogsCommon::FindSeismicGradient(const std::vector<SeismicStorage> & seismic_data,
                                            const Simbox                * const estimation_simbox,
                                            int                                 n_angles,
                                            std::vector<double>               & x_gradient,
                                            std::vector<double>               & y_gradient,
                                            std::vector<std::vector<double> > & sigma_gradient)
{
  int i, j, k, l;
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
  //NBNB marita Må testes om det fungerer for forskjellige brønner
  int nx = estimation_simbox->getnx();
  int ny = estimation_simbox->getny();
  int di_neg = 0; int di_pos = 0; int dj_neg = 0; int dj_pos = 0; //The max replacement in well in x and y direction.
  for(k = 1; k < n_blocks_; k++){
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
  if(di != 0 || dj != 0){
    //adjust the well location
    for(k = 0; k < n_blocks_; k++){
      i_pos_[k] += di;
      j_pos_[k] += dj;
    }
  }

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/traces");
  std::ofstream out(buffer);
  delete [] buffer;

  //Temp
  //SegY * segy = seismic_data[0].GetSegY();
  //segy->GetAllValues

  double dz, ztop, dzW, ztopW;
  for(l = 0; l < n_angles; l++){
    for(j = -yEx; j <= yEx; j++){
      for(i = -xEx; i <= xEx; i++){


        seis_trace = seismic_data[l].GetRealTrace(estimation_simbox, i0, j0); ///H Correct values returned?
        //seis_trace = seisCube[l]->getRealTrace2(i0, j0);

        SmoothTrace(seis_trace);
        //if(j == 0 ){
        //  for(int s = 0; s < seisTrace.size(); s++)
        //    out << seisTrace[s] << std::endl;
        //}

        dzW =  estimation_simbox->getdz(i0,j0);
        ztopW =  estimation_simbox->getTop(i0,j0);
        FindPeakTrace(seis_trace, z_peak_well, peak_well, b_well, dzW, ztopW);

        //seis_trace = seisCube[l]->getRealTrace2(i0+i, j0+j);
        seis_trace = seismic_data[l].GetRealTrace(estimation_simbox, i0, j0);
        SmoothTrace(seis_trace);
        if(i==0){
          for(size_t s = 0; s < seis_trace.size(); s++)
            out << seis_trace[s] << std::endl;
        }

        dz =  estimation_simbox->getdz(i0+i, j0+j);
        ztop =  estimation_simbox->getTop(i0+i, j0+j);
        FindPeakTrace(seis_trace, z_peak, peak, b, dz, ztop);

        PeakMatch(z_peak,peak,b,z_peak_well,peak_well,b_well);//Finds the matching peaks in two traces
        z_shift[(i+2) + (j+2)*nZx] = ComputeShift(z_peak,z_peak_well,z_pos_[0]);

        for(k = 1; k < n_blocks_; k++){
          //Check if well changes lateral position
          if((i_pos_[k]- i_pos_[k-1] == 0) && (j_pos_[k] - j_pos_[k-1] == 0))
            z_shift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = ComputeShift(z_peak,z_peak_well,z_pos_[k]);
          else{
            //well has changed lateral position and we adapt to the new well position
            //seis_trace = seisCube[l]->getRealTrace2(i_pos_[k],j_pos_[k]);
            seis_trace = seismic_data[l].GetRealTrace(estimation_simbox, i_pos_[k], j_pos_[k]);
            SmoothTrace(seis_trace);
            dzW = estimation_simbox->getdz(i_pos_[k],j_pos_[k]);
            ztopW = estimation_simbox->getTop(i_pos_[k],j_pos_[k]);
            FindPeakTrace(seis_trace, z_peak_well, peak_well, b_well, dzW, ztopW);

            //seis_trace = seisCube[l]->getRealTrace2(i_pos_[k]+i, j_pos_[k]+j);
            seis_trace = seismic_data[l].GetRealTrace(estimation_simbox, i_pos_[k]+i, j_pos_[k]+j);
            SmoothTrace(seis_trace);
            dz = estimation_simbox->getdz(i_pos_[k]+i, j_pos_[k]+j);
            ztop = estimation_simbox->getTop(i_pos_[k]+i, j_pos_[k]+j);
            FindPeakTrace(seis_trace, z_peak, peak, b, dz, ztop);

            PeakMatch(z_peak,peak,b,z_peak_well,peak_well,b_well);

            z_shift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = ComputeShift(z_peak, z_peak_well, z_pos_[k]);
          }
        }
      }
    }

    double dx = estimation_simbox->getdx();
    double dy = estimation_simbox->getdy();

    ComputeGradient(q_epsilon, q_epsilon_data, z_shift, nZx, nZx, dx, dy);
  }

  SmoothGradient(x_gradient, y_gradient, q_epsilon, q_epsilon_data, sigma_gradient);
  // NBNB Odd slår av estimeringen for å teste om det gir bedre resultat
 /* for(k = 0; k < nBlocks_; k++){
    xGradient[k]=0.0;
    yGradient[k]=0.0;
  }*/

}

void BlockedLogsCommon::SmoothTrace(std::vector<float> &trace)
{
  float smoothing_distance = 40; //ms in each direction
  int  L = static_cast<int>(ceil(smoothing_distance/dz_)); //number of lags in the gauss kernel
  float sigma = 10 /dz_; // ms / (ms/cell)

  unsigned int n_trace = static_cast<unsigned int>(trace.size());
  std::vector<float> gk(2*L+1);
  std::vector<float> s_trace(n_trace);
  unsigned int i;
  int j;
  float tmp;
  for(j = -L; j <= L; j++){
    tmp = (j*j)/(2*sigma*sigma);
    gk[j+L] = exp(-tmp);
  }

  float N;
  for(i = 0; i < n_trace; i++){
    N = 0;
    for(j = -L; j <= L; j++){
      if(i+j >= 0 && i+j < n_trace){
        s_trace[i] += gk[j+L]*trace[i+j];
        N += gk[j+L];
      }
    }
    s_trace[i] /= N;
  }

  for(i = 0; i < n_trace; i++)
    trace[i] = s_trace[i];
}

void BlockedLogsCommon::FindPeakTrace(std::vector<float> &trace, std::vector<double> &z_peak, std::vector<double> &peak,
                                      std::vector<double> &b, double dz, double z_top)
{
  int k;
  double x1, x2, x3, y1, y2, y3, y11, y12, y21;
  int N = static_cast<int>(trace.size());
  z_peak.resize(N); peak.resize(N); b.resize(N);

  double a;
  double c;
  int counter = 0;
  for(k = 1; k < N-1; k++){
    if((trace[k] >= trace[k-1] && trace[k] > trace[k+1]) || (trace[k] <= trace[k-1] && trace[k] < trace[k+1])){
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

void BlockedLogsCommon::PeakMatch(std::vector<double> &z_peak, std::vector<double> &peak, std::vector<double> &b,
                                  std::vector<double> &z_peak_w, std::vector<double> &peak_w, std::vector<double> &bW)
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
  for(i = 0; i < z_peak_w.size(); i++){
    for(j = lim; j < z_peak.size(); j++){
      diffz = fabs(z_peak_w[i] - z_peak[j]);
      if(diffz < maxdiffz){
        //Check if the peaks point in the same direction
        if((bW[i] < 0 && b[j] < 0)||(bW[i] >= 0 && b[j] >= 0)){
          // Check for difference in peak size
          if((fabs(peak_w[i] - peak[j]))/(fabs(peak_w[i]) + fabs(peak[j])) < diffp){
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
  for(i = 0; i < static_cast<unsigned int>(counter); i++){
    z_peak_w[i] = pW[i];
    z_peak[i] =  p[i];
  }

}

double BlockedLogsCommon::ComputeShift(std::vector<double> &z_peak, std::vector<double> &z_peak_w, double z0)
{
  //This routine computes the position of z0 between two peaks in the well and finds the corresponding distance in
  //the other trace. Then zShift is the difference in z between the two.
  unsigned int N = static_cast<unsigned int>(z_peak.size());
  if(N == 0)
    return RMISSING; //The case of no match in the traces
  else{
    unsigned int i;
    int pos = 0;
    double zShift;
    if(z0 < z_peak_w[0])
      zShift = z_peak_w[0] - z_peak[0];
    else if(z0 >= z_peak_w[N-1])
      zShift = z_peak_w[N-1] - z_peak[N-1];
    else{
      for(i = 0; i < N-1; i++){
        if(z0 >= z_peak_w[i] && z0 < z_peak_w[i+1]){
          pos = i;
          i = N;
        }
      }
      zShift = z0 - (z_peak[pos] + (z0 - z_peak_w[pos])/(z_peak_w[pos+1]-z_peak_w[pos])*(z_peak[pos+1]-z_peak[pos]));
    }

    return zShift;
  }
}

void BlockedLogsCommon::ComputeGradient(std::vector<double> &q_epsilon, std::vector<double> &q_epsilon_data,
                                        std::vector<double> &z_shift, int nx, int ny, double dx, double dy)
{
  //This fit the model zshift(x,y)= beta0 + beta1*x + beta2*y  ==> beta1 is x-gradient and beta2 is y-gradient
  int i, j, k, l;
  std::vector<double> Z(3*nx*ny);
  std::vector<double> Y(nx*ny);
  std::vector<double> cov(9);
  std::vector<double> invcov(9);
  std::vector<double> regM(3*nx*ny);

  static bool append = false;

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/gradNoSmooth");
  std::ofstream out;
  if(append){
    out.open(buffer, std::ios::app|std::ios::out);
  }
  else{
    out.open(buffer);
    append = true;
  }
  delete [] buffer;

  int ndata;
  double data;

  int cy, cz;
  int counter1 = 0;
  int counter2 = 0;
  for(l = 0; l < n_blocks_; l++){
    cy = 0; cz = 0;
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        data = z_shift[i + j*nx + l*nx*ny];
        if(data != RMISSING){
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
    for(j = 0; j < 3; j++){
      for(i = 0; i < 3; i++){
        tmp = 0;
        for(k = 0; k < ndata; k++)
          tmp += Z[i + 3*k] * Z[j + 3*k];
        cov[i + 3*j] = tmp;
      }
    }
    double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  +   cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

    if(det != 0){
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
      for(j = 0; j < 3; j++){
        for(i = 0; i < ndata; i++){
          tmp = 0;
          for(k = 0; k < 3; k++)
            tmp += invcov[k + 3*j]*Z[k + 3*i];
          regM[i + j*ndata] = tmp;
        }
      }


      //Compute beta_1(gradientx) og beta_2(gradienty), beta_0 not necessary
      double beta0 = 0;
      double beta1 = 0;
      double beta2 = 0;
      for(j = 0; j < ndata; j++){
        beta0 += regM[j]*Y[j];
        beta1 += regM[j + ndata]*Y[j];
        beta2 += regM[j + 2*ndata]*Y[j];}
      double sigma2 = 0;
      double sigmatmp;
      for(j = 0; j < ndata; j++){
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
                                       std::vector<std::vector<double> > & sigma_gradient)
{
  int i, j;
  int n_beta = n_blocks_* 2;

  NRLib::SymmetricMatrix q_beta_data = NRLib::SymmetricZeroMatrix(n_beta);

  //Set the prior precicion values
  double a, b, c;
  ComputePrecisionMatrix(a,b,c);

  //Initialize the Qm|d matrix
  for(i = 0; i < n_beta - 2; i++){
    q_beta_data(i,i)   = c + q_epsilon[2*i];
    q_beta_data(i,i+2) = b;
    if(i % 2 == 0)
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

  for(i = 0; i < n_beta; i++){
    double tmp = 0;
    for(j = 0; j < n_beta; j++)
      tmp += I(i,j)*q_epsilon_data[j];
    res(i) = tmp;
  }

  // Return Sigma_gradient
  sigma_gradient.resize(n_beta);
  std::vector<double> tmp_vec(n_beta);
  for(i = 0; i < n_beta; i++){
    for(j = 0; j < n_beta; j++)
      tmp_vec[j] = I(i,j);
    sigma_gradient[i] = tmp_vec;
  }


  int counter = 0;
  for(i = 0; i < n_blocks_; i++){
    x_gradient[i] = res(counter);
    y_gradient[i] = res(counter+1);
    counter += 2;
  }

  /*
   char* buffer2 = new char[1000];
   sprintf(buffer2,"%s.txt", "C:/Outputfiles/gradients");
   std::ofstream out2(buffer2);
   for(i = 0; i < nBlocks_; i++){
       out2 << xGradient[i] << " " << yGradient[i];
     out2 << std::endl;
   }
  delete [] buffer2;
  */
}

void BlockedLogsCommon::ComputePrecisionMatrix(double &a, double &b, double &c)
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

void BlockedLogsCommon::GetBlockedGrid(const Simbox         * estimation_simbox,
                                       const SeismicStorage * seismic_data,
                                       double               * blockedLog,
                                       int                    i_offset,
                                       int                    j_offset)
{
  for (int m = 0 ; m < n_blocks_ ; m++) {
    blockedLog[m] = seismic_data->GetRealTraceValue(estimation_simbox, i_pos_[m]+i_offset, j_pos_[m]+j_offset, k_pos_[m]);
  }
}

void
BlockedLogsCommon::GetVerticalTrend(const std::vector<double> & blocked_log,
                                    double * trend)
{
  if (blocked_log.size() > 0 && trend != NULL) {
    int * count = new int[n_layers_];
    for (int k = 0 ; k < n_layers_ ; k++) {
      trend[k] = 0.0f;
      count[k] = 0;
    }
    for (int m = 0 ; m < n_blocks_ ; m++) {
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
      InterpolateTrend(blocked_log, trend);
    delete [] count;
  }
  else {
    if (blocked_log.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }
}

void
BlockedLogsCommon::InterpolateTrend(const std::vector<double> & blocked_log,
                                    double * trend)
{
  for (int m = 1 ; m < n_blocks_ ; m++) {
    if(abs(k_pos_[m]-k_pos_[m-1]) > 1) {
      int delta = 1;
      if(k_pos_[m] < k_pos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(k_pos_[m]-k_pos_[m-1]);
      float t = step_mult;
      for(int j = k_pos_[m-1]+delta; j != k_pos_[m];j++) {
        if(trend[j] == RMISSING)
          trend[j] = t*blocked_log[m]+(1-t)*blocked_log[m-1];
        t += step_mult;
      }
    }
  }
}

void BlockedLogsCommon::FindContiniousPartOfData(const std::vector<bool> & has_data,
                                                 int                       nz,
                                                 int                     & start,
                                                 int                     & length)
{
  int  i;
  int  l_pice             =  0;
  int  length_max_pice    = -1;
  int  start_longest_pice =  0;
  bool previous_had_data  = false;

  for(i = 0; i < nz ;i++){
    if(has_data[i]){
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
                                  int           nzp)
{
  int i;

  for(i=0;i<nzp;i++)
    cpp_r[i]=0;

  double * alpha_vert = new double[n_layers_];
  double * beta_vert  = new double[n_layers_];
  double * rho_vert   = new double[n_layers_];

  //std::vector<double> alpha_vert(n_layers_);
  //std::vector<double> beta_vert(n_layers_);
  //std::vector<double> rho_vert(n_layers_);

  GetVerticalTrend(GetVp(), alpha_vert);
  GetVerticalTrend(GetVs(), beta_vert);
  GetVerticalTrend(GetRho(), rho_vert);

  for(i=start;i < start+length-1;i++)
  {
    float ei1 = ComputeElasticImpedance(alpha_vert[i],beta_vert[i],rho_vert[i],coeff);
    float ei2 = ComputeElasticImpedance(alpha_vert[i+1],beta_vert[i+1],rho_vert[i+1],coeff);
    cpp_r[i] =  ei2-ei1;
  }
  delete [] alpha_vert;
  delete [] beta_vert;
  delete [] rho_vert;

}

float BlockedLogsCommon::ComputeElasticImpedance(double         alpha,
                                                 double         beta,
                                                 double         rho,
                                                 const float  * coeff) const
{
  // vp, vs, rho are logtransformed
  float ang_imp;

  ang_imp = float(coeff[0]*alpha+coeff[1]*beta+coeff[2]*rho );

  return(ang_imp);
}

void BlockedLogsCommon::FillInSeismic(double     * seismicData,
                                      int         start,
                                      int         length,
                                      fftw_real * seis_r,
                                      int         nzp) const
{
  int i;
  for(i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for(i=start; i<start+length; i++)
    seis_r[i] = seismicData[i];
/*
  int lTregion = 3;
  int* modify  = getIndexPrior(start,lTregion,nzp);
  int* conditionto = getIndexPost(start-1,lTregion,nzp);
  //NBNB Odd: interpolate endpoints?
*/

}

void BlockedLogsCommon::EstimateCor(fftw_complex * var1_c,
                                    fftw_complex * var2_c,
                                    fftw_complex * ccor_1_2_c,
                                    int            cnzp) const
{
  for(int i=0;i<cnzp;i++){
    ccor_1_2_c[i].re =  var1_c[i].re*var2_c[i].re + var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }
}

void BlockedLogsCommon::SetLogFromVerticalTrend(float      * vertical_trend,
                                                double       z0,              // z-value of center in top layer
                                                double       dz,              // dz in vertical trend
                                                int          nz,              // layers in vertical trend
                                                std::string  type,
                                                int          i_angle)
{
  if (type != "WELL_SYNTHETIC_SEISMIC")
  {
    float * blocked_log = new float[n_blocks_];

    SetLogFromVerticalTrend(blocked_log, GetZpos(), n_blocks_,
                            vertical_trend, z0, dz, nz);

    if (type == "ALPHA_SEISMIC_RESOLUTION")
      alpha_seismic_resolution_ = blocked_log;
    else if (type == "BETA_SEISMIC_RESOLUTION")
      beta_seismic_resolution_ = blocked_log;
    else if (type == "RHO_SEISMIC_RESOLUTION")
      rho_seismic_resolution_ = blocked_log;
    else if (type == "ACTUAL_SYNTHETIC_SEISMIC") {
      if (actual_synt_seismic_data_ == NULL)
        actual_synt_seismic_data_ = new float * [n_angles_]; // "nAngles is set along with real_seismic_data_" Now: Currently not set.
      actual_synt_seismic_data_[i_angle] = blocked_log;
    }
    else {
      LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type+
                           "\" in BlockedLogs::setLogFromVerticalTrend()\n");
      exit(1);
    }
  }
  else if (type == "WELL_SYNTHETIC_SEISMIC") {
    if (well_synt_seismic_data_ == NULL)
    {
      well_synt_seismic_data_ = new float * [n_angles_];
      for (int i=0; i<n_angles_; i++)
      {
        well_synt_seismic_data_[i] = new float[n_blocks_];
        for (int j=0; j<n_blocks_; j++)
          well_synt_seismic_data_[i][j] = RMISSING; //Declare in case the wavelet is not estimated for all angles
      }
    }
    SetLogFromVerticalTrend(well_synt_seismic_data_[i_angle], GetZpos(), n_blocks_,
                            vertical_trend, z0, dz, nz);
  }
}

void BlockedLogsCommon::SetLogFromVerticalTrend(float     *& blocked_log,
                                                const std::vector<double> & zpos,
                                                int          n_blocks,
                                                float      * vertical_trend,
                                                double       z0,
                                                double       dzVt,
                                                int          nz)
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
      dz = zpos[i]-zpos[i-1];
    else
      dz = zpos[i+1]-zpos[i];
    double zi = zpos[i];
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
      value = vertical_trend[j-1];
    }
    else if (b >= a+dz) {
      // One single vertical-trend-cell covers a blocked log cell completely.
      value = vertical_trend[j];
    }
    else {
      double zj = b + 0.5*dzVt; // Center of vertical trend cell
      value = vertical_trend[j]* static_cast<float>((zj+dzVt-zi)/dzVt) + vertical_trend[j]*static_cast<float>((zi-zj)/dzVt);
    }
    blocked_log[i] = value;

    //LogKit::LogFormatted(LogKit::Error,"i j  blockedLog[i]   %d %d  %7.3f\n",i,j,blockedLog[i]);
  }
}