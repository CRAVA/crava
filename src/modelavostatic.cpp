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
#include "src/modelgeneral.h"
#include "src/modelavostatic.h"
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

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


//ModelAVOStatic::ModelAVOStatic(ModelSettings        *& modelSettings,
//                               ModelGeneral         *& modelGeneral,
//                               const InputFiles      * inputFiles,
//                               GridMapping           * timeCutMapping,
//                               Simbox                * timeSimbox,
//                               Simbox               *& timeBGSimbox,
//                               Simbox                * timeSimboxConstThick,
//                               std::vector<WellData *> wells)
//{
//  forwardModeling_        = modelSettings->getForwardModeling();
//
//  bool failedModelGeneral = modelGeneral->getFailed();
//
//  failed_                 = false;
//  bool failedExtraSurf    = false;
//  bool failedPriorFacies  = false;
//
//  bool failedLoadingModel = false;
//
//  std::string errText("");
//
//  if(!failedModelGeneral)
//  {
//    if (modelSettings->getForwardModeling() == false)
//    {
//      //
//      // INVERSION/ESTIMATION
//      //
//
//      loadExtraSurfaces(waveletEstimInterval_, faciesEstimInterval_, wellMoveInterval_,
//                        timeSimbox, inputFiles, errText, failedExtraSurf);
//
//      blockLogs(wells, timeSimbox, timeBGSimbox, timeSimboxConstThick, modelSettings);
//
//      checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
//      bool estimationMode = modelSettings->getEstimationMode();
//      if (estimationMode == false && !failedExtraSurf)
//      {
//        Simbox * timeCutSimbox = NULL;
//        if (timeCutMapping != NULL)
//          timeCutSimbox = timeCutMapping->getSimbox(); // For the got-enough-data test
//        else
//          timeCutSimbox = timeSimbox;
//
//        modelGeneral->processPriorFaciesProb(faciesEstimInterval_,
//                                             wells,
//                                             timeSimbox,
//                                             timeCutSimbox,
//                                             modelSettings,
//                                             failedPriorFacies,
//                                             errText,
//                                             inputFiles);
//      }
//    }
//    else // forward modeling
//      checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
//  }
//  failedLoadingModel = failedExtraSurf || failedPriorFacies;
//
//  if (failedLoadingModel) {
//    LogKit::WriteHeader("Error(s) while loading data");
//    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
//    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
//  }
//
//  failed_ = failedLoadingModel;
//  failed_details_.push_back(failedExtraSurf);
//  failed_details_.push_back(failedPriorFacies);
//}

ModelAVOStatic::ModelAVOStatic(ModelSettings        *& model_settings,
                               const InputFiles      * input_files,
                               CommonData            * common_data,
                               const Simbox          * simbox,
                               int                     i_interval)
{
  forward_modeling_        = model_settings->getForwardModeling();

  if (model_settings->getForwardModeling() == false)
  {
    //
    // INVERSION/ESTIMATION
    //
    CheckAvailableMemory(simbox, model_settings, input_files);
    bool estimationMode = model_settings->getEstimationMode();
    if (estimationMode == false)
      facies_estim_interval_ = common_data->GetFaciesEstimInterval(); //Read in in CommonData under SetupPriorFaciesProb based on estimation_simbox. Should this have been per interval?

    //Maybe load in from modelSettings here, things that are needed in doAVOInversion


    //Set up errCorr here
    int nx  = common_data->GetMultipleIntervalGrid()->GetBackgroundParametersInterval(i_interval)[0]->GetNI();
    int ny  = common_data->GetMultipleIntervalGrid()->GetBackgroundParametersInterval(i_interval)[0]->GetNJ();
    int nz  = common_data->GetMultipleIntervalGrid()->GetBackgroundParametersInterval(i_interval)[0]->GetNK();
    int nxp = simbox->GetNXpad();
    int nyp = simbox->GetNYpad();
    int nzp = simbox->GetNZpad();

    err_corr_ = CreateFFTGrid(nx, ny, nz,
                              nxp, nyp, nzp,
                              model_settings->getFileGrid());
    err_corr_ ->setType(FFTGrid::COVARIANCE);
    err_corr_ ->createRealGrid();

    float corr_grad_I = 0.0f;
    float corr_grad_J = 0.0f;
    common_data->GetCorrGradIJ(corr_grad_I, corr_grad_J, simbox);

    err_corr_->fillInErrCorr(common_data->GetPriorCorrXY(i_interval), corr_grad_I, corr_grad_J);

    //CheckAvailableMemory(simbox, model_settings, input_files);

  }
  else // forward modeling
    CheckAvailableMemory(simbox, model_settings, input_files);
}


ModelAVOStatic::~ModelAVOStatic(void)
{

  if (wavelet_estim_interval_.size() == 2) {
    if (wavelet_estim_interval_[0] != NULL)
      delete wavelet_estim_interval_[0];
    if (wavelet_estim_interval_[1] != NULL)
      delete wavelet_estim_interval_[1];
  }

  if (facies_estim_interval_.size() == 2) {
    if (facies_estim_interval_[0] != NULL)
      delete facies_estim_interval_[0];
    if (facies_estim_interval_[1] != NULL)
      delete facies_estim_interval_[1];
  }

  if (well_move_interval_.size() == 2) {
    if (well_move_interval_[0] != NULL)
      delete well_move_interval_[0];
    if (well_move_interval_[1] != NULL)
      delete well_move_interval_[1];
  }
}

//void
//ModelAVOStatic::blockLogs(std::vector<WellData *> & wells,
//                          Simbox                  * timeSimbox,
//                          Simbox                  * timeBGSimbox,
//                          Simbox                  * timeSimboxConstThick,
//                          ModelSettings          *& modelSettings)
//{
//  int     nWells         = modelSettings->getNumberOfWells();
//
//  if(nWells > 0) {
//    for (int i=0 ; i<nWells ; i++)
//    {
//      wells[i]->findMeanVsVp(waveletEstimInterval_);
//
//      wells[i]->setBlockedLogsOrigThick( new BlockedLogs(wells[i], timeSimbox, modelSettings->getRunFromPanel()) );
//      wells[i]->setBlockedLogsConstThick( new BlockedLogs(wells[i], timeSimboxConstThick) );
//      if (timeBGSimbox==NULL)
//        wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeSimbox) ); // Need a copy constructor?
//      else
//        wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeBGSimbox) );
//    }
//  }
//}

void
ModelAVOStatic::CheckAvailableMemory(const Simbox     * time_simbox,
                                     ModelSettings    * model_settings,
                                     const InputFiles * input_files)
{
  LogKit::WriteHeader("Estimating amount of memory needed");
  //
  // Find the size of first seismic volume
  //
  float mem_one_seis = 0.0f;
  if (input_files->getNumberOfSeismicFiles(0) > 0 && input_files->getSeismicFile(0,0) != "") {
    mem_one_seis = static_cast<float> (NRLib::FindFileSize(input_files->getSeismicFile(0,0)));
  }

  //
  // Find the size of one grid
  //
  FFTGrid * dummy_grid = new FFTGrid(time_simbox->getnx(),
                                     time_simbox->getny(),
                                     time_simbox->getnz(),
                                     time_simbox->GetNXpad(),
                                     time_simbox->GetNYpad(),
                                     time_simbox->GetNZpad());
  long long int grid_size_pad = static_cast<long long int>(4)*dummy_grid->getrsize();

  delete dummy_grid;
  dummy_grid = new FFTGrid(time_simbox->getnx(),
                           time_simbox->getny(),
                           time_simbox->getnz(),
                           time_simbox->getnx(),
                           time_simbox->getny(),
                           time_simbox->getnz());
  long long int grid_size_base = 4*dummy_grid->getrsize();
  delete dummy_grid;
  int n_grid_parameters   = 3;                                      // Vp + Vs + Rho, padded
  int n_grid_background   = 3;                                      // Vp + Vs + Rho, padded
  int n_grid_covariances  = 6;                                      // Covariances, padded
  int n_grid_seismic_data = model_settings->getNumberOfAngles(0);     // One for each angle stack, padded

  std::map<std::string, float> facies_prob = model_settings->getPriorFaciesProb(); //Used to find number of facies grids needed

  int n_grid_facies       = static_cast<int>(facies_prob.size())+1; // One for each facies, one for undef, unpadded.
  int n_grid_histograms   = static_cast<int>(facies_prob.size());   // One for each facies, 2MB.
  int n_grid_kriging      = 1;                                      // One grid for kriging, unpadded.
  int n_grid_compute      = 1;                                      // Computation grid, padded (for convenience)
  int n_grid_file_mode    = 1;                                      // One grid for intermediate file storage

  int n_grids;
  long long int grid_mem;
  if (model_settings->getForwardModeling() == true) {
    if (model_settings->getFileGrid())  // Use disk buffering
      n_grids = n_grid_file_mode;
    else
      n_grids = n_grid_parameters + 1;

    grid_mem = n_grids*grid_size_pad;
  }
  else {
    if (model_settings->getFileGrid()) { // Use disk buffering
      n_grids = n_grid_file_mode;
      if (model_settings->getKrigingParameter() > 0) {
        n_grids += n_grid_kriging;
      }
      if (model_settings->getNumberOfSimulations() > 0)
        n_grids = n_grid_parameters;
      if (model_settings->getUseLocalNoise(0)) {
        n_grids = 2*n_grid_parameters;
      }

      grid_mem = n_grids*grid_size_pad;
    }
    else {
      //baseP and baseU are the padded and unpadde grids allocated at each peak.
      int base_P = n_grid_parameters + n_grid_covariances;
      if (model_settings->getUseLocalNoise(0) == true || (model_settings->getEstimateFaciesProb() && model_settings->getFaciesProbRelative()))
        base_P += n_grid_background;
      int base_U = 0;
      if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        base_U += static_cast<int>(facies_prob.size());

      //First peak: At inversion
      int peak_1P = base_P + n_grid_seismic_data; //Need seismic data as well here.
      int peak_1U = base_U;

      long long int peak_grid_mem = peak_1P*grid_size_pad + peak_1U*grid_size_base; //First peak must be currently largest.
      int peak_n_grid = peak_1P;                                             //Also in number of padded grids

      if (model_settings->getNumberOfSimulations() > 0) { //Second possible peak when simulating.
        int peak_2P = base_P + 3; //Three extra parameter grids for simulated parameters.
        if (model_settings->getUseLocalNoise(0) == true &&
           (model_settings->getEstimateFaciesProb() == false || model_settings->getFaciesProbRelative() == false))
          peak_2P -= n_grid_background; //Background grids are released before simulation in this case.
        int peak_2U = base_U;     //Base level is the same, but may increase.
        bool compute_grid_used = ((model_settings->getOutputGridsElastic() & (IO::AI + IO::LAMBDARHO + IO::LAMELAMBDA + IO::LAMEMU + IO::MURHO + IO::POISSONRATIO + IO::SI + IO::VPVSRATIO)) > 0);
        if (compute_grid_used == true)
          peak_2P += n_grid_compute;
        else if (model_settings->getKrigingParameter() > 0) //Note the else, since this grid will use same memory as computation grid if both are active.
          peak_2U += n_grid_kriging;

        if (peak_2P > peak_n_grid)
          peak_n_grid = peak_2P;

        long long int peak_2_mem = peak_2P*grid_size_pad + peak_2U*grid_size_base;
        if (peak_2_mem > peak_grid_mem)
          peak_grid_mem = peak_2_mem;
      }

      if (model_settings->getEstimateFaciesProb() == true) {//Third possible peak when computing facies prob.
        int peak_3P = base_P;                //No extra padded grids, so this one can not peak here.
        int peak_3U = base_U + n_grid_facies;  //But this one will, and may trigger new memory max.
        if ((model_settings->getOtherOutputFlag() & IO::FACIES_LIKELIHOOD) > 0)
          peak_3U += 1; //Also needs to store seismic likelihood.

        long long int peak_3_mem = peak_3P*grid_size_pad + peak_3U*grid_size_base + 2000000*n_grid_histograms; //These are 2MB when Vs is used.
        if (peak_3_mem > peak_grid_mem)
          peak_grid_mem = peak_3_mem;
      }
      n_grids  = peak_n_grid;
      grid_mem = peak_grid_mem;
    }
  }
  FFTGrid::setMaxAllowedGrids(n_grids);
  if (model_settings->getDebugFlag()>0)
    FFTGrid::setTerminateOnMaxGrid(true);

  int   work_size   = 2500 + static_cast<int>( 0.65*grid_size_pad); //Size of memory used beyond grids.

  float mem0        = 4.0f * work_size;
  float mem1        = static_cast<float>(grid_mem);
  float mem2        = static_cast<float>(model_settings->getNumberOfAngles(0))*grid_size_pad + mem_one_seis; //Peak memory when reading seismic, overestimated.

  float needed_mem   = mem0 + std::max(mem1, mem2);

  float mega_bytes   = needed_mem/(1024.f*1024.f);
  float giga_bytes   = mega_bytes/1024.f;

  LogKit::LogFormatted(LogKit::High,"\nMemory needed for reading seismic data       : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding internal grids (%2d): %10.2f MB\n",n_grids, mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding other entities     : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (mega_bytes > 1000.0f)
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f gigaBytes\n",giga_bytes);
  else
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f megaBytes\n",mega_bytes);

  if (mem2>mem1)
    LogKit::LogFormatted(LogKit::Low,"\n This estimate is too high because seismic data are cut to fit the internal grid\n");
  if (!model_settings->getFileGrid()) {
    //
    // Check if we can hold everything in memory.
    //
    model_settings->setFileGrid(false);
    char ** memchunk  = new char*[n_grids];

    int i = 0;
    try {
      for(i = 0 ; i < n_grids ; i++)
        memchunk[i] = new char[static_cast<size_t>(grid_size_pad)];
    }
    catch (std::bad_alloc& ) //Could not allocate memory
    {
      model_settings->setFileGrid(true);
      LogKit::LogFormatted(LogKit::Low,"Not enough memory to hold all grids. Using file storage.\n");
    }

    for(int j=0 ; j<i ; j++)
      delete [] memchunk[j];
    delete [] memchunk;
  }
}

//void
//ModelAVOStatic::loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
//                                  std::vector<Surface *> & faciesEstimInterval,
//                                  std::vector<Surface *> & wellMoveInterval,
//                                  Simbox                 * timeSimbox,
//                                  const InputFiles       * inputFiles,
//                                  std::string            & errText,
//                                  bool                   & failed)
//{
//  const double x0 = timeSimbox->getx0();
//  const double y0 = timeSimbox->gety0();
//  const double lx = timeSimbox->getlx();
//  const double ly = timeSimbox->getly();
//  const int    nx = timeSimbox->getnx();
//  const int    ny = timeSimbox->getny();
//  //
//  // Get wavelet estimation interval
//  //
//  const std::string & topWEI  = inputFiles->getWaveletEstIntFileTop(0); //Same for all time lapses
//  const std::string & baseWEI = inputFiles->getWaveletEstIntFileBase(0);//Same for all time lapses
//
//  if (topWEI != "" && baseWEI != "") {
//    waveletEstimInterval.resize(2);
//    try {
//      if (NRLib::IsNumber(topWEI))
//        waveletEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWEI.c_str()));
//      else {
//        Surface tmpSurf(topWEI);
//        waveletEstimInterval[0] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//
//    try {
//      if (NRLib::IsNumber(baseWEI))
//        waveletEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWEI.c_str()));
//      else {
//        Surface tmpSurf(baseWEI);
//        waveletEstimInterval[1] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//  }
//  //
//  // Get facies estimation interval
//  //
//  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
//  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);
//
//  if (topFEI != "" && baseFEI != "") {
//    faciesEstimInterval.resize(2);
//    try {
//      if (NRLib::IsNumber(topFEI))
//        faciesEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topFEI.c_str()));
//      else {
//        Surface tmpSurf(topFEI);
//        faciesEstimInterval[0] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//
//    try {
//      if (NRLib::IsNumber(baseFEI))
//        faciesEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseFEI.c_str()));
//      else {
//        Surface tmpSurf(baseFEI);
//        faciesEstimInterval[1] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//  }
//  //
//  // Get well move interval
//  //
//  const std::string & topWMI  = inputFiles->getWellMoveIntFile(0);
//  const std::string & baseWMI = inputFiles->getWellMoveIntFile(1);
//
//  if (topWMI != "" && baseWMI != "") {
//    wellMoveInterval.resize(2);
//    try {
//      if (NRLib::IsNumber(topWMI))
//        wellMoveInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWMI.c_str()));
//      else {
//        Surface tmpSurf(topWMI);
//        wellMoveInterval[0] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//
//    try {
//      if (NRLib::IsNumber(baseWMI))
//        wellMoveInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWMI.c_str()));
//      else {
//        Surface tmpSurf(baseWMI);
//        wellMoveInterval[1] = new Surface(tmpSurf);
//      }
//    }
//    catch (NRLib::Exception & e) {
//      errText += e.what();
//      failed = true;
//    }
//  }
//}

//void ModelAVOStatic::writeWells(std::vector<WellData *> wells, ModelSettings * modelSettings) const
//{
//  int nWells  = modelSettings->getNumberOfWells();
//  for(int i=0;i<nWells;i++)
//    wells[i]->writeWell(modelSettings->getWellFormatFlag());
//}

//void ModelAVOStatic::writeBlockedWells(std::vector<WellData *> wells, ModelSettings * modelSettings, std::vector<std::string> facies_name, std::vector<int> facies_label) const
//{
//  int nWells  = modelSettings->getNumberOfWells();
//  for(int i=0;i<nWells;i++)
//    wells[i]->getBlockedLogsOrigThick()->writeWell(modelSettings, facies_name, facies_label);
//}

//void ModelAVOStatic::WriteBlockedWells(std::map<std::string, BlockedLogsCommon *> blocked_wells,
//                                       ModelSettings                            * model_settings,
//                                       std::vector<std::string>                   facies_name,
//                                       std::vector<int>                           facies_label) const
//{
//  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
//    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
//    BlockedLogsCommon * blocked_log = iter->second;
//
//    blocked_log->WriteWell(model_settings->getWellFormatFlag(),
//                           model_settings->getMaxHzBackground(),
//                           model_settings->getMaxHzSeismic(),
//                           facies_name,
//                           facies_label);
//
//  }
//}

//void ModelAVOStatic::addSeismicLogs(std::vector<WellData *> wells,
//                                    FFTGrid              ** seisCube,
//                                    const ModelSettings   * modelSettings,
//                                    int                     nAngles)
//{
//  int nWells  = modelSettings->getNumberOfWells();
//
//    for (int iAngle = 0 ; iAngle < nAngles ; iAngle++)
//    {
//      seisCube[iAngle]->setAccessMode(FFTGrid::RANDOMACCESS);
//      for(int i=0;i<nWells;i++)
//        wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
//      seisCube[iAngle]->endAccess();
//  }
//}

void ModelAVOStatic::AddSeismicLogs(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                    std::vector<FFTGrid *>                     seis_cube,
                                    int                                        n_angles)
{
  for (int i_angle = 0; i_angle < n_angles; i_angle++) {
    seis_cube[i_angle]->setAccessMode(FFTGrid::RANDOMACCESS);

    for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
      std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
      BlockedLogsCommon * blocked_log = iter->second;

      blocked_log->SetLogFromGrid(seis_cube[i_angle], i_angle, n_angles, "SEISMIC_DATA");
    }
    //for(int i=0;i<nWells;i++)
    //  wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
    seis_cube[i_angle]->endAccess();
  }
}

//void ModelAVOStatic::generateSyntheticSeismic(Wavelet              ** wavelet,
//                                              std::vector<WellData *> wells,
//                                              const float *   const * reflectionMatrix,
//                                              const Simbox          * timeSimbox,
//                                              const ModelSettings   * modelSettings,
//                                              int                     nAngles)
//{
//  int nWells  = modelSettings->getNumberOfWells();
//  int nzp     = modelSettings->getNZpad();
//  int nz      = timeSimbox->getnz();
//
//  int i;
//
//  for( i=0; i<nWells; i++ )
//  {
//    if( wells[i]->isDeviated() == false )
//      wells[i]->getBlockedLogsOrigThick()->generateSyntheticSeismic(reflectionMatrix,nAngles,wavelet,nz,nzp,timeSimbox);
//  }
//}

void ModelAVOStatic::GenerateSyntheticSeismic(std::vector<Wavelet *>                   & wavelet,
                                              std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                              const float *                      const * reflection_matrix,
                                              const Simbox                             * time_simbox,
                                              int                                        n_angles)
{
  int nzp = time_simbox->GetNZpad();
  int nz  = time_simbox->getnz();

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log->GetIsDeviated() == true)
      blocked_log->GenerateSyntheticSeismic(reflection_matrix, n_angles, wavelet, nz, nzp, time_simbox);
  }
}

//-------------------------------------------------------------------
FFTGrid *
ModelAVOStatic::CreateFFTGrid(int nx,  int ny,  int nz,
                              int nxp, int nyp, int nzp,
                              bool file_grid)
{
  FFTGrid * fft_grid;
  if (file_grid)
    fft_grid = new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fft_grid = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);
  return(fft_grid);
}


//void ModelAVOStatic::deleteDynamicWells(std::vector<WellData *> wells,
//                                        int                     nWells)
//{
//  for(int i=0; i<nWells; i++){
//    wells[i]->getBlockedLogsConstThick()->deleteDynamicBlockedLogs();
//    wells[i]->getBlockedLogsExtendedBG()->deleteDynamicBlockedLogs();
//    wells[i]->getBlockedLogsOrigThick()->deleteDynamicBlockedLogs();
//  }
//}
