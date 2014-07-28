/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/definitions.h"
#include "src/modelgravitystatic.h"
#include "src/modelgravitydynamic.h"
#include "src/gravimetricinversion.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/tasklist.h"
#include "src/seismicparametersholder.h"
//#include "src/commondata.h"

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/volume/volume.hpp"

ModelGravityDynamic::ModelGravityDynamic(const ModelSettings          * modelSettings,
                                         ModelGravityStatic           * modelGravityStatic,
                                         const Simbox                 * simbox,
                                         CommonData                   * commonData,
                                         int                            t,
                                         SeismicParametersHolder      & seismicParameters)

{
  // BuildGMatrix: For each obs (=30). Create x,y,z based on simbox, and create localDistanceSquared = distance x and x0(1 of 30 obs from file).
  // Change observation_location_depth_ to be per interval? Only use the x0, y0, z0 which are inside simbox?

  //H Changes to multiple intervals: Removed modelGeneral, only use the new interval simbox
  //modelGeneral_           = modelGeneral;    // For å bruke full size time simbox. Løse dette bedre...

  debug_                  = true;
  failed_                 = false;
  thisTimeLapse_          = t;

  bool failedLoadingModel = false;
  //bool failedReadingFile  = false;
  std::string errText("");

 // int nObs = 30;     // should this be given in input file
  //int nColumns = 5;  // We require data files to have five columns

  // Check that timeLapse is ok
  if(thisTimeLapse_ < 1 && thisTimeLapse_ > modelSettings->getNumberOfVintages()){
    errText += "Not valid time lapse";
    failedLoadingModel = true;
  }

  if(failedLoadingModel == false){
    LogKit::WriteHeader("Setting up gravimetric time lapse");

    observation_location_utmx_  = commonData->GetGravityObservationUtmxTimeLapse(thisTimeLapse_);
    observation_location_utmy_  = commonData->GetGravityObservationUtmyTimeLapse(thisTimeLapse_);
    observation_location_depth_ = commonData->GetGravityObservationDepthTimeLapse(thisTimeLapse_);
    gravity_response_           = commonData->GetGravityResponseTimeLapse(thisTimeLapse_);
    gravity_std_dev_            = commonData->GetGravityStdDevTimeLapse(thisTimeLapse_);


    LogKit::LogFormatted(LogKit::Low, "Setting up forward model matrix ...");
    BuildGMatrix(modelGravityStatic, seismicParameters, simbox);
    LogKit::LogFormatted(LogKit::Low, "ok.\n");
  }

}

ModelGravityDynamic::~ModelGravityDynamic(void)
{
}

void ModelGravityDynamic::BuildGMatrix(ModelGravityStatic      * modelGravityStatic,
                                       SeismicParametersHolder & seismicParameters,
                                       const Simbox            * simbox)
{
  // Building gravity matrix for each time vintage, using updated mean Vp in generating the grid.
  double gamma = 6.67384e-11; // units: m^3/(kg*s^2)

  //const Simbox * fullSizeTimeSimbox = NULL;
  //if(simbox == NULL)
  //  fullSizeTimeSimbox = modelGeneral_->GetTimeSimbox();
  //else
  //  fullSizeTimeSimbox = simbox;

    // Use vp_current, found in Seismic parameters holder here.
  FFTGrid * expMeanVp      = new FFTGrid(seismicParameters.GetMeanVp());  // for upscaling
  FFTGrid * meanVpFullSize = new FFTGrid(expMeanVp);                    // for full size matrix

  int nx = simbox->getnx();
  int ny = simbox->getny();
  int nz = simbox->getnz();

  double dx = simbox->getdx();
  double dy = simbox->getdy();
  //double dz = fullSizeTimeSimbox->getdz();  //Tvilsomt?

  int nxp = expMeanVp->getNxp();
  int nyp = expMeanVp->getNyp();
  int nzp = expMeanVp->getNzp();

  int nxp_upscaled = modelGravityStatic->GetNxp_upscaled();
  int nyp_upscaled = modelGravityStatic->GetNyp_upscaled();
  int nzp_upscaled = modelGravityStatic->GetNzp_upscaled();

  int upscaling_factor_x = nxp/nxp_upscaled;
  int upscaling_factor_y = nyp/nyp_upscaled;
  int upscaling_factor_z = nzp/nzp_upscaled;

  // dimensions of one grid cell
  double dx_upscaled = dx*upscaling_factor_x;
  double dy_upscaled = dy*upscaling_factor_y;

  int nx_upscaled = modelGravityStatic->GetNx_upscaled();
  int ny_upscaled = modelGravityStatic->GetNy_upscaled();
  int nz_upscaled = modelGravityStatic->GetNz_upscaled();

  int N_upscaled = nx_upscaled*ny_upscaled*nz_upscaled;
  int N_fullsize = nx*ny*nz;

  int nObs = 30;

  G_         .resize(nObs, N_upscaled);
  G_fullsize_.resize(nObs, N_fullsize);

  // Need to be in real domain for transforming from log domain
  if(expMeanVp->getIsTransformed())
    expMeanVp->invFFTInPlace();

  if(meanVpFullSize->getIsTransformed())
    meanVpFullSize->invFFTInPlace();

  float sigma_squared = GravimetricInversion::GetSigmaForTransformation(seismicParameters.GetCovVp());
  GravimetricInversion::MeanExpTransform(expMeanVp,      sigma_squared);
  GravimetricInversion::MeanExpTransform(meanVpFullSize, sigma_squared);


  //Smooth (convolve) and subsample
  FFTGrid * upscalingKernel_conj = modelGravityStatic->GetUpscalingKernel();
   if(upscalingKernel_conj->getIsTransformed() == false)
    upscalingKernel_conj->fftInPlace();
  upscalingKernel_conj->conjugate();  // Conjugate only in FFT domain.

  // Need to be in FFT domain for convolution and subsampling
  if(expMeanVp->getIsTransformed() == false)
    expMeanVp->fftInPlace();

  expMeanVp->multiply(upscalingKernel_conj);  // Now is expMeanAlpha smoothed

  FFTGrid * upscaledMeanVp;
  GravimetricInversion::Subsample(upscaledMeanVp, expMeanVp, nx_upscaled, ny_upscaled, nz_upscaled, nxp_upscaled, nyp_upscaled, nzp_upscaled);

  upscaledMeanVp->invFFTInPlace();

  float x0, y0, z0; // Coordinates for the observations points
  int J = 0;        // Index in matrix counting cell number
  int I = 0;

  float  vp;
  double dt;
  double localMass;
  double localDistanceSquared;

  NRLib::Matrix UpscaledCoord(nx_upscaled*ny_upscaled*nz_upscaled, 3);
  NRLib::Matrix FullSizeCoord(nx*ny*nz, 3);

  for(int i = 0; i < nObs; i++){
    x0 = observation_location_utmx_[i];
    y0 = observation_location_utmy_[i];
    z0 = observation_location_depth_[i];

     J = 0; I = 0;

    // Loop through upscaled simbox to get x, y, z for each grid cell
    for(int ii = 0; ii < nx_upscaled; ii++){   // eller k først
      for(int jj = 0; jj < ny_upscaled; jj++){
        for(int kk = 0; kk < nz_upscaled; kk++){   // eller i sist??
          double x, y, z;
          vp = upscaledMeanVp->getRealValue(ii,jj,kk);

          int istart = ii*upscaling_factor_x;
          int istop  = (ii+1)*upscaling_factor_x;
          int jstart = jj*upscaling_factor_y;
          int jstop  = (jj+1)*upscaling_factor_y;
          int kstart = kk*upscaling_factor_z;
          int kstop = (kk+1)*upscaling_factor_z;

          double x_local  = 0;
          double y_local  = 0;
          double z_local  = 0;
          double dt_local = 0;
          //Find centerposition of coarse grid cell using indices of fine grid and averaging their cell centers.
          for(int iii=istart; iii<istop; iii++){
            for(int jjj=jstart; jjj<jstop; jjj++){
              for(int kkk=kstart; kkk<kstop; kkk++){
                simbox->getCoord(iii, jjj, kkk, x, y, z);
                x_local += x;
                y_local += y;
                z_local += z;   // NB NB Need time depth mapping!!

                dt_local += simbox->getdz(iii, jjj); // also average dt?
              }
            }
          }
          x_local  /= (upscaling_factor_x*upscaling_factor_y*upscaling_factor_z);
          y_local  /= (upscaling_factor_x*upscaling_factor_y*upscaling_factor_z);
          z_local  /= (upscaling_factor_x*upscaling_factor_y*upscaling_factor_z);
          dt_local /= (upscaling_factor_x*upscaling_factor_y*upscaling_factor_z);

          // Find fraction of dx_upscaled and dy_upscaled according to indicies
          double xfactor = 1;
          if(istop <= nx){  // inside
            xfactor = 1;
          }
          else if(istart > nx){ //outside
            xfactor = 0;
          }
          else{
            xfactor = (nx - istart)/upscaling_factor_x;
          }

          double yfactor = 1;
          if(jstop <= ny){
            yfactor = 1;
          }
          else if(jstart > ny){
            yfactor = 0;
          }
          else{
            yfactor = (ny - jstart)/upscaling_factor_y;
          }
          if(debug_){ // OBS
            dt_local = 1;
            vp = 160;
          }

          if(debug_){
            UpscaledCoord(J, 0) = x_local;
            UpscaledCoord(J, 1) = y_local;
            UpscaledCoord(J, 2) = z_local;
          }

          localMass = (xfactor*dx_upscaled)*(yfactor*dy_upscaled)*dt_local*vp*0.5*1000; // units kg
          localDistanceSquared = pow((x_local-x0),2) + pow((y_local-y0),2) + pow((z_local-z0),2); //units m^2

          G_(i,J) = localMass/localDistanceSquared;
          J++;
        }
      }
    }
    // Loop through full size simbox to get x, y, z for each grid cell
    for(int ii = 0; ii < nx; ii++){
      for(int jj = 0; jj < ny; jj++){
        for(int kk = 0; kk < nz; kk++){
          double x, y, z;
          simbox->getCoord(ii, jj, kk, x, y, z); // assuming these are center positions...
          vp = meanVpFullSize->getRealValue(ii, jj, kk);
          dt = simbox->getdz(ii, jj);
          if(debug_){
            dt = 1;  // such that vp*dt = 80
            vp= 160; // = 2*80, dz i matlab er 80   //
          }
         // int simbox_i = expMeanAlpha->getXSimboxIndex(ii);
         // int simbox_j = expMeanAlpha->getYSimboxIndex(jj);

          if(debug_){
            FullSizeCoord(I, 0) = x;
            FullSizeCoord(I, 1) = y;
            FullSizeCoord(I, 2) = z;
          }

          localMass = dx*dy*dt*vp*0.5*1000; // units kg
          localDistanceSquared = pow((x-x0),2) + pow((y-y0),2) + pow((z-z0),2); //units m^2

          G_fullsize_(i,I) = localMass/localDistanceSquared;
          I++;
        }
      }
    }
  }
  G_ = G_*gamma;
  G_fullsize_ = G_fullsize_*gamma;

  //delete expMeanAlpha;
  //delete meanAlphaFullSize;
  //delete upscaledMeanAlpha;

  if(debug_){
    NRLib::WriteMatrixToFile("UpscaledCoords.txt", UpscaledCoord);
    NRLib::WriteMatrixToFile("FullSizeCoords.txt", FullSizeCoord);
    // Dump G_ matrix
    NRLib::WriteMatrixToFile("UpscaledGMatrix_dynamicClass.txt", G_);
    // Dump G_fullsize_ matrix
    NRLib::WriteMatrixToFile("FullSizeGMatrix_dynamicClass.txt", G_fullsize_);
  }
}


