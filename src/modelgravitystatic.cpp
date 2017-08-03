/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelgravitystatic.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/analyzelog.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/waveletfilter.h"
#include "src/tasklist.h"
#include "src/parameteroutput.h"
#include "nrlib/surface/surface.hpp"

ModelGravityStatic::ModelGravityStatic(ModelSettings        *& modelSettings,
                                       ModelGeneral         *& modelGeneral,
                                       const InputFiles      * inputFiles)
{
  modelGeneral_           = modelGeneral;

  failed_                 = false;
  before_injection_start_ = false;

  bool failedLoadingModel = false;
  bool failedReadingFile  = false;
  std::string errText("");

  bool doGravityInversion = true;
  int numberGravityFiles  = 0;
  for(int i = 0; i<modelSettings->getNumberOfVintages(); i++){
    if(modelSettings->getGravityTimeLapse(i))
      numberGravityFiles++;
  }

  if(numberGravityFiles == 0){
    // Everything is ok - we do not need gravity inversion
    failedLoadingModel = false;
    doGravityInversion = false;
  }

  if(numberGravityFiles == 1){
    failedLoadingModel = true;
    doGravityInversion = false;
    errText+="Need at least two gravity surveys for inversion.";
  }

  // Set up gravimetric baseline
  if(doGravityInversion){

    LogKit::WriteHeader("Setting up gravimetric baseline");

    // Find first gravity data file
    std::string fileName = inputFiles->getGravimetricData(0);

    int nObs = 30;     //user input
    int nColumns = 5;  // We require data files to have five columns

    observation_location_utmx_.resize(nObs);
    observation_location_utmy_.resize(nObs);
    observation_location_depth_.resize(nObs);
    gravity_response_.resize(nObs);
    gravity_std_dev_.resize(nObs);

    ReadGravityDataFile(fileName, "gravimetric base survey",
                        nObs, nColumns,
                        observation_location_utmx_,
                        observation_location_utmy_,
                        observation_location_depth_,
                        gravity_response_,
                        gravity_std_dev_,
                        failedReadingFile,
                        errText);
    failedLoadingModel = failedReadingFile;

    const Simbox * fullTimeSimbox = modelGeneral->GetSimbox();

    x_upscaling_factor_ = fullTimeSimbox->GetNXpad()/5 + 1;   // should be user input...
    y_upscaling_factor_ = fullTimeSimbox->GetNYpad()/5 + 1;
    z_upscaling_factor_ = fullTimeSimbox->GetNZpad()/5 + 1;

    SetUpscaledPaddingSize(fullTimeSimbox);  // NB: Changes upscaling factors!

    dx_upscaled_ = fullTimeSimbox->GetLX()/nx_upscaled_;
    dy_upscaled_ = fullTimeSimbox->GetLY()/ny_upscaled_;
    dz_upscaled_ = fullTimeSimbox->GetLZ()/nz_upscaled_;

    LogKit::LogFormatted(LogKit::Low, "Generating smoothing kernel ...");
    MakeUpscalingKernel(modelSettings, fullTimeSimbox);
    LogKit::LogFormatted(LogKit::Low, "ok.\n");

    LogKit::LogFormatted(LogKit::Low, "Generating lag index table (size " + NRLib::ToString(nxp_upscaled_) + " x "
                                                                          + NRLib::ToString(nyp_upscaled_) + " x "
                                                                          + NRLib::ToString(nzp_upscaled_) + ") ...");
    MakeLagIndex(nxp_upscaled_, nyp_upscaled_, nzp_upscaled_); // Including padded region!
    LogKit::LogFormatted(LogKit::Low, "ok.\n");
  }

  if (failedLoadingModel) {
    LogKit::WriteHeader("Error(s) with gravimetric surveys");
    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel || failedReadingFile;
  failed_details_.push_back(failedReadingFile);
}


ModelGravityStatic::~ModelGravityStatic(void)
{
}

void
ModelGravityStatic::ReadGravityDataFile(const std::string   & fileName,
                                        const std::string   & readReason,
                                        int                   nObs,
                                        int                   nColumns,
                                        std::vector <float> & obs_loc_utmx,
                                        std::vector <float> & obs_loc_utmy,
                                        std::vector <float> & obs_loc_depth,
                                        std::vector <float> & gravity_response,
                                        std::vector <float> & gravity_std_dev,
                                        bool                  failed,
                                        std::string         & errText)
{
  float * tmpRes = new float[nObs*nColumns+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,fileName);
  std::string text = "Reading "+readReason+" from file "+fileName+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index = 0;
  failed = false;

  while(failed == false && inFile >> storage) {
    if(index < nObs*nColumns) {
      try {
        tmpRes[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        errText += "Error in "+fileName+"\n";
        errText += e.what();
        failed = true;
      }
    }
    index++;
  }
  if(failed == false) {
    if(index != nObs*nColumns) {
      failed = true;
      errText += "Found "+NRLib::ToString(index)+" in file "+fileName+", expected "+NRLib::ToString(nObs*nColumns)+".\n";
    }
  }

  if(failed == false) {
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
    index = 0;
    for(int i=0;i<nObs;i++) {
      obs_loc_utmx[i] = tmpRes[index];
      index++;
      obs_loc_utmy[i] = tmpRes[index];
      index++;
      obs_loc_depth[i] = tmpRes[index];
      index++;
      gravity_response[i] = tmpRes[index];
      index++;
      gravity_std_dev[i] = tmpRes[index];
      index++;
    }
  }
  else
    LogKit::LogFormatted(LogKit::Low,"failed.\n");

  delete [] tmpRes;
}

void
ModelGravityStatic::MakeUpscalingKernel(ModelSettings * modelSettings, const Simbox * fullTimeSimbox)
{
  int nx = fullTimeSimbox->getnx();
  int ny = fullTimeSimbox->getny();
  int nz = fullTimeSimbox->getnz();

  int nxp = fullTimeSimbox->GetNXpad();
  int nyp = fullTimeSimbox->GetNYpad();
  int nzp = fullTimeSimbox->GetNZpad();

  upscaling_kernel_ = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);
  upscaling_kernel_->setType(FFTGrid::PARAMETER);
  upscaling_kernel_->fillInConstant(0.0);

  upscaling_kernel_->setAccessMode(FFTGrid::RANDOMACCESS);

  for(int k = 0; k < nz_upscaled_; k++)
    for(int j = 0; j < ny_upscaled_; j++)
      for(int i = 0; i < nx_upscaled_; i++)
        upscaling_kernel_->setRealValue(i, j, k, 1.0);

  upscaling_kernel_->endAccess();

  upscaling_kernel_->multiplyByScalar(static_cast<float>(nxp_upscaled_*nyp_upscaled_*nzp_upscaled_)/static_cast<float>(nxp*nyp*nzp));
}

void ModelGravityStatic::MakeLagIndex(int nx_upscaled, int ny_upscaled, int nz_upscaled)
{
  int N_up = nx_upscaled*ny_upscaled*nz_upscaled;
  lag_index_.resize(N_up);
  for (int i = 0; i < N_up; ++i) {
    lag_index_[i].resize(N_up);

    for (int j = 0; j < N_up; ++j)
      lag_index_[i][j].resize(3);
  }

  int I, J;
  for(int k1 = 1; k1 <= nz_upscaled; k1++)
    for(int j1 = 1; j1 <= ny_upscaled; j1++)
      for(int i1 = 1; i1 <= nx_upscaled; i1++){
        I =  i1 + (j1-1)*nx_upscaled + (k1-1)*nx_upscaled*ny_upscaled;

        for(int k2 = 1; k2 <= nz_upscaled; k2++)
          for(int j2 = 1; j2 <= ny_upscaled; j2++)
            for(int i2 = 1; i2 <= nx_upscaled; i2++){
              J = i2 + (j2-1)*nx_upscaled + (k2-1)*nx_upscaled*ny_upscaled;

              int lag_i = i2 - i1;
              int lag_j = j2 - j1;
              int lag_k = k2 - k1;

              int ind1, ind2, ind3;
              if(abs(lag_i) <= nx_upscaled/2 && abs(lag_j) <= ny_upscaled/2 && abs(lag_k) <= nz_upscaled/2) {
                if(lag_i >= 0)
                  ind1 = lag_i + 1;
                else
                  ind1 = nx_upscaled + lag_i + 1;

                if(lag_j >= 0)
                  ind2 = lag_j + 1;
                else
                  ind2 = ny_upscaled + lag_j + 1;

                if(lag_k >= 0)
                  ind3 = lag_k + 1;
                else
                  ind3 = nz_upscaled + lag_k + 1;

                lag_index_[I-1][J-1][0] = ind1 - 1;   // NB: -1
                lag_index_[I-1][J-1][1] = ind2 - 1;
                lag_index_[I-1][J-1][2] = ind3 - 1;
              }
              else
              {
                lag_index_[I-1][J-1][0] = -1;
                lag_index_[I-1][J-1][1] = -1;
                lag_index_[I-1][J-1][2] = -1;
              }
            }
      }
}
void
ModelGravityStatic::SetUpscaledPaddingSize(const Simbox * fullTimeSimbox)
{
  // Find original nxp, nyp, nzp
  int nxpad = fullTimeSimbox->GetNXpad();
  int nypad = fullTimeSimbox->GetNXpad();
  int nzpad = fullTimeSimbox->GetNXpad();

  int nxpad_up = SetPaddingSize(nxpad, x_upscaling_factor_);
  int nypad_up = SetPaddingSize(nypad, y_upscaling_factor_);
  int nzpad_up = SetPaddingSize(nzpad, z_upscaling_factor_);

  // Initilizing!
  nxp_upscaled_ = nxpad_up;
  nyp_upscaled_ = nypad_up;
  nzp_upscaled_ = nzpad_up;

  nx_upscaled_ = nxpad_up;
  ny_upscaled_ = nypad_up;
  nz_upscaled_ = nzpad_up;

  // Set true upscaling factors
  x_upscaling_factor_ = nxpad/nxp_upscaled_;
  y_upscaling_factor_ = nypad/nyp_upscaled_;
  z_upscaling_factor_ = nzpad/nzp_upscaled_;
}


int
ModelGravityStatic::SetPaddingSize(int original_nxp, int upscaling_factor)
{
  int leastint = static_cast<int>(ceil(static_cast<double>(original_nxp)/static_cast<double>(upscaling_factor)));

  std::vector<int> exp_list = findClosestFactorableNumber(original_nxp);

  int closestprod = original_nxp;

  int factor = 1;

  for(int i=0;i<exp_list[0]+1;i++)
    for(int j=0;j<exp_list[1]+1;j++)
      for(int k=0;k<exp_list[2]+1;k++)
        for(int l=0;l<exp_list[3]+1;l++)
          for(int m=0;m<exp_list[4]+1;m++)
            for(int n=exp_list[4];n<exp_list[5]+1;n++)
            {
              factor = static_cast<int>(pow(2.0f,i)*pow(3.0f,j)*pow(5.0f,k)*pow(7.0f,l)*pow(11.0f,m)*pow(13.0f,n));
              if ((factor >=  leastint) &&  (factor <  closestprod))
              {
                closestprod=factor;
              }
            }
  return closestprod;
}

// Same as in FFTGrid-class, however, this one returns list of exponents
std::vector<int> ModelGravityStatic::findClosestFactorableNumber(int leastint)
{
  int i,j,k,l,m,n;
  int factor   =       1;

  std::vector<int> exp_list(6);

  int maxant2    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(2.0f) ));
  int maxant3    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(3.0f) ));
  int maxant5    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(5.0f) ));
  int maxant7    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(7.0f) ));
  int maxant11   = 0;
  int maxant13   = 0;

  int closestprod= static_cast<int>(pow(2.0f,maxant2));
  exp_list[0] = maxant2;
  exp_list[1] = 0;
  exp_list[2] = 0;
  exp_list[3] = 0;
  exp_list[4] = 0;
  exp_list[5] = 0;

  for(i=0;i<maxant2+1;i++)
    for(j=0;j<maxant3+1;j++)
      for(k=0;k<maxant5+1;k++)
        for(l=0;l<maxant7+1;l++)
          for(m=0;m<maxant11+1;m++)
            for(n=maxant11;n<maxant13+1;n++)
            {
              factor = static_cast<int>(pow(2.0f,i)*pow(3.0f,j)*pow(5.0f,k)*
                pow(7.0f,l)*pow(11.0f,m)*pow(13.0f,n));
              if ((factor >=  leastint) &&  (factor <  closestprod))
              {
                exp_list[0] = i;
                exp_list[1] = j;
                exp_list[2] = k;
                exp_list[3] = l;
                exp_list[4] = m;
                exp_list[5] = n;
                closestprod=factor;
              }
            }
  return exp_list;
}
