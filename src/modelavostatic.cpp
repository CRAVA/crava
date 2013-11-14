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
//#include "src/welldata.h"
#include "src/blockedlogs.h"
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

ModelAVOStatic::ModelAVOStatic(ModelSettings        *& modelSettings,
                               ModelGeneral         *& modelGeneral,
                               const InputFiles      * inputFiles,
                               CommonData            * commonData,
                               Simbox                * timeSimbox,
                               int                     i_interval)
{
  forwardModeling_        = modelSettings->getForwardModeling();

  if (modelSettings->getForwardModeling() == false)
  {
    //
    // INVERSION/ESTIMATION
    //

    faciesEstimInterval_ = commonData->GetFaciesEstimInterval();

    checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
  }
  else // forward modeling
    checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
}


ModelAVOStatic::~ModelAVOStatic(void)
{

  if(waveletEstimInterval_.size() == 2) {
    if (waveletEstimInterval_[0] != NULL)
      delete waveletEstimInterval_[0];
    if (waveletEstimInterval_[1] != NULL)
      delete waveletEstimInterval_[1];
  }

  if(faciesEstimInterval_.size() == 2) {
    if (faciesEstimInterval_[0] != NULL)
      delete faciesEstimInterval_[0];
    if (faciesEstimInterval_[1] != NULL)
      delete faciesEstimInterval_[1];
  }

  if(wellMoveInterval_.size() == 2) {
    if (wellMoveInterval_[0] != NULL)
      delete wellMoveInterval_[0];
    if (wellMoveInterval_[1] != NULL)
      delete wellMoveInterval_[1];
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
ModelAVOStatic::checkAvailableMemory(Simbox           * timeSimbox,
                                     ModelSettings    * modelSettings,
                                     const InputFiles * inputFiles)
{
  LogKit::WriteHeader("Estimating amount of memory needed");
  //
  // Find the size of first seismic volume
  //
  float memOneSeis = 0.0f;
  if (inputFiles->getNumberOfSeismicFiles(0) > 0 && inputFiles->getSeismicFile(0,0) != "") {
    memOneSeis = static_cast<float> (NRLib::FindFileSize(inputFiles->getSeismicFile(0,0)));
  }

  //
  // Find the size of one grid
  //
  FFTGrid * dummyGrid = new FFTGrid(timeSimbox->getnx(),
                                    timeSimbox->getny(),
                                    timeSimbox->getnz(),
                                    modelSettings->getNXpad(),
                                    modelSettings->getNYpad(),
                                    modelSettings->getNZpad());
  long long int gridSizePad = static_cast<long long int>(4)*dummyGrid->getrsize();

  delete dummyGrid;
  dummyGrid = new FFTGrid(timeSimbox->getnx(),
                          timeSimbox->getny(),
                          timeSimbox->getnz(),
                          timeSimbox->getnx(),
                          timeSimbox->getny(),
                          timeSimbox->getnz());
  long long int gridSizeBase = 4*dummyGrid->getrsize();
  delete dummyGrid;
  int nGridParameters  = 3;                                      // Vp + Vs + Rho, padded
  int nGridBackground  = 3;                                      // Vp + Vs + Rho, padded
  int nGridCovariances = 6;                                      // Covariances, padded
  int nGridSeismicData = modelSettings->getNumberOfAngles(0);     // One for each angle stack, padded

  std::map<std::string, float> facies_prob = modelSettings->getPriorFaciesProb(); //Used to find number of facies grids needed

  int nGridFacies       = static_cast<int>(facies_prob.size())+1; // One for each facies, one for undef, unpadded.
  int nGridHistograms   = static_cast<int>(facies_prob.size());   // One for each facies, 2MB.
  int nGridKriging      = 1;                                      // One grid for kriging, unpadded.
  int nGridCompute      = 1;                                      // Computation grid, padded (for convenience)
  int nGridFileMode     = 1;                                      // One grid for intermediate file storage

  int nGrids;
  long long int gridMem;
  if(modelSettings->getForwardModeling() == true) {
    if (modelSettings->getFileGrid())  // Use disk buffering
      nGrids = nGridFileMode;
    else
      nGrids = nGridParameters + 1;

    gridMem = nGrids*gridSizePad;
  }
  else {
    if (modelSettings->getFileGrid()) { // Use disk buffering
      nGrids = nGridFileMode;
      if(modelSettings->getKrigingParameter() > 0) {
        nGrids += nGridKriging;
      }
      if(modelSettings->getNumberOfSimulations() > 0)
        nGrids = nGridParameters;
      if(modelSettings->getUseLocalNoise(0)) {
        nGrids = 2*nGridParameters;
      }

      gridMem = nGrids*gridSizePad;
    }
    else {
      //baseP and baseU are the padded and unpadde grids allocated at each peak.
      int baseP = nGridParameters + nGridCovariances;
      if(modelSettings->getUseLocalNoise(0) == true || (modelSettings->getEstimateFaciesProb() && modelSettings->getFaciesProbRelative()))
        baseP += nGridBackground;
      int baseU = 0;
      if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        baseU += static_cast<int>(facies_prob.size());

      //First peak: At inversion
      int peak1P = baseP + nGridSeismicData; //Need seismic data as well here.
      int peak1U = baseU;

      long long int peakGridMem = peak1P*gridSizePad + peak1U*gridSizeBase; //First peak must be currently largest.
      int peakNGrid   = peak1P;                                             //Also in number of padded grids

      if(modelSettings->getNumberOfSimulations() > 0) { //Second possible peak when simulating.
        int peak2P = baseP + 3; //Three extra parameter grids for simulated parameters.
        if(modelSettings->getUseLocalNoise(0) == true &&
           (modelSettings->getEstimateFaciesProb() == false || modelSettings->getFaciesProbRelative() == false))
          peak2P -= nGridBackground; //Background grids are released before simulation in this case.
        int peak2U = baseU;     //Base level is the same, but may increase.
        bool computeGridUsed = ((modelSettings->getOutputGridsElastic() & (IO::AI + IO::LAMBDARHO + IO::LAMELAMBDA + IO::LAMEMU + IO::MURHO + IO::POISSONRATIO + IO::SI + IO::VPVSRATIO)) > 0);
        if(computeGridUsed == true)
          peak2P += nGridCompute;
        else if(modelSettings->getKrigingParameter() > 0) //Note the else, since this grid will use same memory as computation grid if both are active.
          peak2U += nGridKriging;

        if(peak2P > peakNGrid)
          peakNGrid = peak2P;

        long long int peak2Mem = peak2P*gridSizePad + peak2U*gridSizeBase;
        if(peak2Mem > peakGridMem)
          peakGridMem = peak2Mem;
      }

      if(modelSettings->getEstimateFaciesProb() == true) {//Third possible peak when computing facies prob.
        int peak3P = baseP;                //No extra padded grids, so this one can not peak here.
        int peak3U = baseU + nGridFacies;  //But this one will, and may trigger new memory max.
        if((modelSettings->getOtherOutputFlag() & IO::FACIES_LIKELIHOOD) > 0)
          peak3U += 1; //Also needs to store seismic likelihood.

        long long int peak3Mem = peak3P*gridSizePad + peak3U*gridSizeBase + 2000000*nGridHistograms; //These are 2MB when Vs is used.
        if(peak3Mem > peakGridMem)
          peakGridMem = peak3Mem;
      }
      nGrids  = peakNGrid;
      gridMem = peakGridMem;
    }
  }
  FFTGrid::setMaxAllowedGrids(nGrids);
  if(modelSettings->getDebugFlag()>0)
    FFTGrid::setTerminateOnMaxGrid(true);

  int   workSize    = 2500 + static_cast<int>( 0.65*gridSizePad); //Size of memory used beyond grids.

  float mem0        = 4.0f * workSize;
  float mem1        = static_cast<float>(gridMem);
  float mem2        = static_cast<float>(modelSettings->getNumberOfAngles(0))*gridSizePad + memOneSeis; //Peak memory when reading seismic, overestimated.

  float neededMem   = mem0 + std::max(mem1, mem2);

  float megaBytes   = neededMem/(1024.f*1024.f);
  float gigaBytes   = megaBytes/1024.f;

  LogKit::LogFormatted(LogKit::High,"\nMemory needed for reading seismic data       : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding internal grids (%2d): %10.2f MB\n",nGrids, mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding other entities     : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (megaBytes > 1000.0f)
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f gigaBytes\n",gigaBytes);
  else
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f megaBytes\n",megaBytes);

  if(mem2>mem1)
    LogKit::LogFormatted(LogKit::Low,"\n This estimate is too high because seismic data are cut to fit the internal grid\n");
  if (!modelSettings->getFileGrid()) {
    //
    // Check if we can hold everything in memory.
    //
    modelSettings->setFileGrid(false);
    char ** memchunk  = new char*[nGrids];

    int i = 0;
    try {
      for(i = 0 ; i < nGrids ; i++)
        memchunk[i] = new char[static_cast<size_t>(gridSizePad)];
    }
    catch (std::bad_alloc& ) //Could not allocate memory
    {
      modelSettings->setFileGrid(true);
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

void ModelAVOStatic::writeBlockedWells(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                       ModelSettings                            * modelSettings,
                                       std::vector<std::string>                   facies_name,
                                       std::vector<int>                           facies_label) const
{
  //int nWells  = modelSettings->getNumberOfWells();
  //for(int i=0;i<nWells;i++)
  //  wells[i]->getBlockedLogsOrigThick()->writeWell(modelSettings, facies_name, facies_label);

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    blocked_log->WriteWell(modelSettings->getWellFormatFlag(),
                           modelSettings->getMaxHzBackground(),
                           modelSettings->getMaxHzSeismic(),
                           facies_name,
                           facies_label);

  }
}

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

void ModelAVOStatic::addSeismicLogs(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                    FFTGrid                                 ** seisCube,
                                    const ModelSettings                      * modelSettings,
                                    int                                        nAngles)
{
  int nWells  = modelSettings->getNumberOfWells();

    for (int iAngle = 0 ; iAngle < nAngles ; iAngle++)
    {
      seisCube[iAngle]->setAccessMode(FFTGrid::RANDOMACCESS);

      for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
        std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
        BlockedLogsCommon * blocked_log = iter->second;

        blocked_log->SetLogFromGrid(seisCube[iAngle], iAngle, nAngles, "SEISMIC_DATA");
      }
      //for(int i=0;i<nWells;i++)
      //  wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
      seisCube[iAngle]->endAccess();
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

void ModelAVOStatic::generateSyntheticSeismic(Wavelet                                 ** wavelet,
                                              std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                              const float *                      const * reflectionMatrix,
                                              const Simbox                             * timeSimbox,
                                              const ModelSettings                      * modelSettings,
                                              int                                        nAngles)
{
  int nWells  = modelSettings->getNumberOfWells();
  int nzp     = modelSettings->getNZpad();
  int nz      = timeSimbox->getnz();

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if(blocked_log->GetIsDeviated() == true)
      blocked_log->GenerateSyntheticSeismic(reflectionMatrix, nAngles, wavelet, nz, nzp, timeSimbox);
  }
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
