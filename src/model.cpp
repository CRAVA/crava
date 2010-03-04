#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/definitions.h"
#include "src/model.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/corr.h"
#include "src/analyzelog.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/welldata.h"
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
#include "lib/lib_matr.h"
#include "lib/global_def.h"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


Model::Model(const std::string & fileName)
{
  modelSettings_          = NULL;
  timeSimbox_             = new Simbox();
  timeSimboxConstThick_   = NULL;
  wells_                  = NULL;
  background_             = NULL;
  correlations_           = NULL;
  priorFacies_            = NULL;
  seisCube_               = NULL;
  wavelet_                = NULL;
 
  correlationDirection_   = NULL;     
  reflectionMatrix_       = NULL;
  randomGen_              = NULL;
  failed_                 = false;
  gradX_                  = 0.0;
  gradY_                  = 0.0;
  
  priorFaciesProbCubes_   = NULL;
 
  timeDepthMapping_       = NULL;
  timeCutMapping_         = NULL;
  velocityFromInversion_  = false;

  bool failedWavelet      = false;
  bool failedSeismic      = false;
  bool failedSimbox       = false;
  bool failedWells        = false;
  bool failedReflMat      = false;
  bool failedExtraSurf    = false;
  bool failedBackground   = false;
  bool failedPriorCorr    = false;
  bool failedDepthConv    = false;
  bool failedPriorFacies  = false;

  bool failedModelFile    = false;
  bool failedInputFiles   = false;
  bool failedLoadingModel = false;

  Simbox * timeCutSimbox  = NULL;
  Simbox * timeBGSimbox   = NULL;

  InputFiles * inputFiles;

  XmlModelFile modelFile(fileName);
  inputFiles     = modelFile.getInputFiles();
  modelSettings_ = modelFile.getModelSettings();

  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
 
  std::string errTxt = inputFiles->addInputPathAndCheckFiles();
  if(errTxt != "") {
    Utils::writeHeader("Error opening files");
    LogKit::LogMessage(LogKit::ERROR, "\n"+errTxt);
    LogKit::LogFormatted(LogKit::ERROR,"\nAborting\n");
    failedInputFiles = true;
  }

  if(!failedModelFile && !failedInputFiles)
  {
    LogKit::SetScreenLog(modelSettings_->getLogLevel());

    std::string logFileName = IO::makeFullFileName("",IO::FileLog()+IO::SuffixTextFiles());
    LogKit::SetFileLog(logFileName,modelSettings_->getLogLevel());

    if(modelSettings_->getDebugFlag() > 0)
    {
      std::string fName = IO::makeFullFileName("",IO::FileDebug()+IO::SuffixTextFiles());
      LogKit::SetFileLog(fName, LogKit::DEBUGHIGH+LogKit::DEBUGLOW);
    }
    LogKit::EndBuffering();
    

    if(inputFiles->getSeedFile() == "")
      randomGen_ = new RandomGen(modelSettings_->getSeed());
    else
      randomGen_ = new RandomGen(inputFiles->getSeedFile().c_str());

    if(modelSettings_->getNumberOfSimulations() == 0)
      modelSettings_->setWritePrediction(true); //write predicted grids. 
    
    printSettings(modelSettings_, inputFiles);
    
    Utils::writeHeader("Defining modelling grid");

    std::string errText("");

    makeTimeSimboxes(timeSimbox_, timeCutSimbox, timeBGSimbox, timeSimboxConstThick_,  //Handles correlation direction too.
                     correlationDirection_, modelSettings_, inputFiles, 
                     errText, failedSimbox);

    if(!failedSimbox)
    { 
      //
      // FORWARD MODELLING
      //
      if (modelSettings_->getForwardModeling() == true)
      {
        checkAvailableMemory(timeSimbox_, modelSettings_, inputFiles); 
        processBackground(background_, wells_, timeSimbox_, timeBGSimbox,
                          modelSettings_, inputFiles,
                          errText, failedBackground);
        if (!failedBackground)
        {
          processReflectionMatrixFromBackground(reflectionMatrix_, background_, 
                                                modelSettings_, inputFiles,
                                                errText, failedReflMat);  
          if (!failedReflMat)
          {
            processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                            timeSimbox_, waveletEstimInterval_, 
                            modelSettings_, inputFiles, errText, failedWavelet);
          }              
        }
      }
      else
      {
        //
        // INVERSION/ESTIMATION
        //
        if(timeCutSimbox!=NULL)  {
          timeCutMapping_ = new GridMapping();
          timeCutMapping_->makeTimeTimeMapping(timeCutSimbox);
        }

        processWells(wells_, timeSimbox_, timeBGSimbox, timeSimboxConstThick_, 
                     randomGen_, modelSettings_, inputFiles, errText, failedWells);

        checkAvailableMemory(timeSimbox_, modelSettings_, inputFiles); 

        loadExtraSurfaces(waveletEstimInterval_, faciesEstimInterval_, wellMoveInterval_, 
                          timeSimbox_, inputFiles, errText, failedExtraSurf);


        bool estimationMode = modelSettings_->getEstimationMode();

        if (!failedWells && !failedDepthConv)
        {
          bool backgroundDone = false;

          if(!(modelSettings_->getOptimizeWellLocation() == true &&
               modelSettings_->getGenerateBackground() == true)) 
          { 
            if(estimationMode == false || 
               modelSettings_->getEstimateBackground() == true ||
               modelSettings_->getEstimateCorrelations() == true || 
               modelSettings_->getEstimateWaveletNoise() == true ||
               modelSettings_->getOptimizeWellLocation() == true) 
            {
              processBackground(background_, wells_, timeSimbox_, timeBGSimbox,
                                modelSettings_, inputFiles,
                                errText, failedBackground);
              backgroundDone = true;
            }

            if(failedBackground == false && backgroundDone == true &&
              (estimationMode == false || modelSettings_->getEstimateWaveletNoise() || 
               modelSettings_->getOptimizeWellLocation() == true))
            {
              processReflectionMatrixFromBackground(reflectionMatrix_, background_, modelSettings_, 
                                                    inputFiles, errText, failedReflMat);
            }

            if(failedBackground == false && backgroundDone == true && 
               failedReflMat == false && failedExtraSurf == false &&
               (estimationMode == false || modelSettings_->getEstimateWaveletNoise() || 
                modelSettings_->getOptimizeWellLocation() == true))
            {
              processSeismic(seisCube_, timeSimbox_, modelSettings_, 
                             inputFiles, errText, failedSeismic);
              if(failedSeismic == false && modelSettings_->getOptimizeWellLocation() == true)
              {
                for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
                  seisCube_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
                processWellLocation(seisCube_, wells_, reflectionMatrix_,
                                    timeSimbox_, modelSettings_, wellMoveInterval_, randomGen_);
                for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
                  seisCube_[i]->endAccess();
              }
            }
          }
          else 
            //
            // When well locations are to be optimized, the reflection matrix is estimated from 
            // the wells instead of the background as the background is dependent of the well 
            // locations. Need to alternate the order of the process-functions, as the well 
            // locations need to be estimated before the background model is processed.
            //
          {
            processReflectionMatrixFromWells(reflectionMatrix_, wells_, modelSettings_, 
                                             inputFiles, errText, failedReflMat);
            if (failedReflMat == false && failedExtraSurf == false)
            {
              processSeismic(seisCube_, timeSimbox_, modelSettings_, 
                             inputFiles, errText, failedSeismic);
              if(failedSeismic == false) 
              {
                processWellLocation(seisCube_, wells_, reflectionMatrix_,
                                    timeSimbox_, modelSettings_, wellMoveInterval_, randomGen_);       
              }
            }
            if(estimationMode == false || 
               modelSettings_->getEstimateBackground() == true ||
               modelSettings_->getEstimateCorrelations() == true || 
               modelSettings_->getEstimateWaveletNoise() == true) 
            {
              processBackground(background_, wells_, timeSimbox_, timeBGSimbox, 
                                modelSettings_, inputFiles, errText, failedBackground);
              backgroundDone = true;
            }
          }

          if(failedBackground == false && backgroundDone == true && 
             failedSeismic == false && failedReflMat == false &&
             (estimationMode == false || modelSettings_->getEstimateCorrelations() == true))
          {
            processPriorCorrelations(correlations_, background_, wells_, timeSimbox_, modelSettings_, 
                                     inputFiles, errText, failedPriorCorr);
          }

          if(failedSeismic == false && failedBackground == false && 
            (estimationMode == false || modelSettings_->getEstimateWaveletNoise() || modelSettings_->getOptimizeWellLocation()))
          {
            addSeismicLogs(wells_, seisCube_, modelSettings_); 
            if(failedReflMat == false && failedExtraSurf == false) 
            {
              processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                              timeSimbox_, waveletEstimInterval_,
                              modelSettings_, inputFiles, errText, failedWavelet);
            }
          }
          if(failedSeismic == false && failedWavelet == false && failedReflMat == false &&
            (modelSettings_->getOptimizeWellLocation() || modelSettings_->getEstimateWaveletNoise() ))
          {
            generateSyntheticSeismic(wavelet_, wells_, reflectionMatrix_, timeSimbox_, modelSettings_);
          }
        }

        if(estimationMode == false && modelSettings_->getDoDepthConversion() == true)
        {
          processDepthConversion(timeCutSimbox, timeSimbox_, modelSettings_, 
                                 inputFiles, errText, failedDepthConv);
        }

        if (estimationMode == false && !failedWells && !failedExtraSurf)
        {
          processPriorFaciesProb(faciesEstimInterval_,
                                 priorFacies_,
                                 wells_,
                                 randomGen_,
                                 timeSimbox_->getnz(),
                                 static_cast<float> (timeSimbox_->getdz()),
                                 modelSettings_,
                                 failedPriorFacies,
                                 errText,
                                 inputFiles);
        }

        if (!failedWells)
        {
          if(((modelSettings_->getWellOutputFlag() & IO::WELLS) > 0) ||
             (estimationMode == true && modelSettings_->getEstimateBackground() == true))
          {
            writeWells(wells_, modelSettings_);
          }
          if(estimationMode  == true && (modelSettings_->getEstimateWaveletNoise() == true ||
                                   (modelSettings_->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0))
            writeBlockedWells(wells_, modelSettings_);
        }
      }
    }
    failedLoadingModel = failedSimbox  || failedSeismic   || failedPriorCorr  ||
                         failedWells   || failedReflMat   || failedBackground ||
                         failedWavelet || failedDepthConv || failedExtraSurf  || failedPriorFacies;

    if (failedLoadingModel) {
      Utils::writeHeader("Error(s) while loading data");
      LogKit::LogFormatted(LogKit::ERROR,"\n"+errText);
      LogKit::LogFormatted(LogKit::ERROR,"\nAborting\n");
    }

    delete inputFiles;
  }
  
  failed_ = failedModelFile || failedInputFiles || failedLoadingModel;

  if(timeCutSimbox != NULL)
    delete timeCutSimbox;
  if(timeBGSimbox != NULL)
    delete timeBGSimbox;
}


Model::~Model(void)
{
  if(!modelSettings_->getForwardModeling()) 
  {
    for(int i=0 ; i<modelSettings_->getNumberOfWells() ; i++)
      if(wells_[i] != NULL)
        delete wells_[i];
    delete [] wells_;
  }

  if (wavelet_ != NULL) 
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(wavelet_[i] != NULL)
        delete wavelet_[i];
    delete [] wavelet_;
  }

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

  if (priorFacies_ != NULL) 
    delete [] priorFacies_;

  if(wellMoveInterval_.size() == 2) {
    if (wellMoveInterval_[0] != NULL)
      delete wellMoveInterval_[0];
    if (wellMoveInterval_[1] != NULL)
      delete wellMoveInterval_[1];
  }

  if(reflectionMatrix_ != NULL) {
    for(int i = 0;i<modelSettings_->getNumberOfAngles();i++)
      delete [] reflectionMatrix_[i] ;
    delete [] reflectionMatrix_ ;
  }

  for (int i = 0 ; i < static_cast<int>(localNoiseScale_.size()) ; i++) {
    if (localNoiseScale_[i] != NULL)
      delete localNoiseScale_[i];
  }
  
  if (correlations_ != NULL)
    delete correlations_;

  if(timeDepthMapping_!=NULL)
    delete timeDepthMapping_;

  if(timeCutMapping_!=NULL)
    delete timeCutMapping_;

  if(correlationDirection_ !=NULL)
    delete correlationDirection_;

  delete randomGen_;
  delete modelSettings_;
  delete timeSimbox_;
  delete timeSimboxConstThick_;

}

void
Model::releaseGrids(void)
{
  delete background_;
  seisCube_ = NULL;
}

float **
Model::readMatrix(const std::string & fileName, int n1, int n2, 
                  const std::string & readReason, 
                  std::string       & errText)
{
  float * tmpRes = new float[n1*n2+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,fileName);
  std::string text = "Reading "+readReason+" from file "+fileName+" ... ";
  LogKit::LogFormatted(LogKit::LOW,text);
  std::string storage;
  int index = 0;
  int error = 0;
  
  while(error == 0 && inFile >> storage) {
    if(index < n1*n2) {
      try {
        tmpRes[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        errText += "Error in "+fileName+"\n"; 
        errText += e.what();
        error = 1;
      }
    }
    index++;
  }
  if(error == 0) {
    if(index != n1*n2) {
      error = 1;
      errText += "Found "+NRLib::ToString(index)+" in file "+fileName+", expected "+NRLib::ToString(n1*n2)+".\n";
    }
  }

  float ** result = NULL;
  if(error == 0) {
    LogKit::LogFormatted(LogKit::LOW,"ok.\n");
    result = new float * [n1];
    int i, j;
    index = 0;
    for(i=0;i<n1;i++) {
      result[i] = new float[n2];
      for(j=0;j<n2;j++) {
        result[i][j] = tmpRes[index];
        index++;
      }
    }
  }
  else
    LogKit::LogFormatted(LogKit::LOW,"failed.\n");
  delete [] tmpRes;
  return(result);
}

void
Model::checkAvailableMemory(Simbox        * timeSimbox,
                            ModelSettings * modelSettings,
                            InputFiles    * inputFiles)
{
  Utils::writeHeader("Estimating amount of memory needed");
  //
  // Find the size of first seismic volume
  //
  float memOneSeis = 0.0f;
  if (inputFiles->getNumberOfSeismicFiles() > 0 && inputFiles->getSeismicFile(0) != "") {
    memOneSeis = static_cast<float> (NRLib::FindFileSize(inputFiles->getSeismicFile(0)));
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
  long long int gridSizePad = 4*dummyGrid->getrsize();
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
  int nGridSeismicData = modelSettings->getNumberOfAngles();     // One for each angle stack, padded

  int nGridFacies       = modelSettings->getNumberOfFacies()+1;   // One for each facies, one for undef, unpadded.
  int nGridHistograms   = modelSettings->getNumberOfFacies();     // One for each facies, 2MB.
  int nGridKriging      = 1;                                      // One grid for kriging, unpadded.
  int nGridCompute      = 1;                                      // Computation grid, padded (for convenience)
  int nGridFileMode     = 1;                                      // One grid for intermediate file storage

  int nGrids;
  long long int gridMem;
  if(modelSettings->getForwardModeling() == true) {
    if (modelSettings->getFileGrid())  // Use disk buffering
      nGrids = nGridFileMode;
    else
      nGrids = nGridParameters + nGridSeismicData;

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
      if(modelSettings->getUseLocalNoise()) {
        nGrids = 2*nGridParameters;
      }

      gridMem = nGrids*gridSizePad;
    }
    else {
      //baseP and baseU are the padded and unpadde grids allocated at each peak.
      int baseP = nGridParameters + nGridCovariances;
      if(modelSettings->getUseLocalNoise() == true || modelSettings->getFaciesProbRelative() == true)
        baseP += nGridBackground;
      int baseU = 0;
      if(modelSettings->getIsPriorFaciesProbGiven()==2)
        baseU += modelSettings->getNumberOfFacies();

      //First peak: At inversion
      int peak1P = baseP + nGridSeismicData; //Need seismic data as well here.
      int peak1U = baseU;
      
      long long int peakGridMem = peak1P*gridSizePad + peak1U*gridSizeBase; //First peak must be currently largest.
      int peakNGrid   = peak1P;                                             //Also in number of padded grids

      if(modelSettings->getNumberOfSimulations() > 0) { //Second possible peak when simulating.
        int peak2P = baseP + 3; //Need three extra parameter grids for simulated parameters.
        int peak2U = baseU;     //Base level is the same, but may increase.
        bool computeGridUsed = ((modelSettings->getGridOutputFlag() & (IO::AI + IO::LAMBDARHO + IO::LAMELAMBDA + IO::LAMEMU + IO::MURHO + IO::POISSONRATIO + IO::SI + IO::VPVSRATIO)) > 0);
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
  float mem2        = static_cast<float>(modelSettings->getNumberOfAngles())*gridSizePad + memOneSeis; //Peak memory when reading seismic, overestimated.
 
  float neededMem   = mem0 + std::max(mem1, mem2);

  float megaBytes   = neededMem/(1024.f*1024.f);
  float gigaBytes   = megaBytes/1024.f;

  LogKit::LogFormatted(LogKit::HIGH,"\nMemory needed for reading seismic data       : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::HIGH,  "Memory needed for holding internal grids (%2d): %10.2f MB\n",nGrids, mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::HIGH,  "Memory needed for holding other entities     : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (gigaBytes < 0.01f)
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f megaBytes\n",megaBytes);
  else
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f gigaBytes\n",gigaBytes);

  if(mem2>mem1)
    LogKit::LogFormatted(LogKit::LOW,"\n This estimate is too high because seismic data are cut to fit the internal grid\n");
  if (!modelSettings->getFileGrid()) {
    //
    // Check if we can hold everything in memory.
    //
    float ** memchunk  = new float*[nGrids];

    for(int i = 0 ; i < nGrids ; i++)
      memchunk[i] = new float[gridSizePad];

    if(memchunk[nGrids-1] == NULL)  //Could not allocate memory
    {
      modelSettings->setFileGrid(true);
      LogKit::LogFormatted(LogKit::LOW,"Not enough memory to hold all grids. Using file storage.\n");
    }
    else
    {
      modelSettings->setFileGrid(false);
    }
    
    for(int i=0 ; i<nGrids ; i++)
      if(memchunk[i] != NULL) delete [] memchunk[i];
    if(memchunk != NULL) delete [] memchunk;
  }
}

int
Model::readSegyFile(const std::string       & fileName, 
                    FFTGrid                *& target, 
                    Simbox                 *& timeSimbox, 
                    ModelSettings          *& modelSettings, 
                    const SegyGeometry     *& geometry,
                    int                       gridType,
                    float                     offset,
                    const TraceHeaderFormat * format,
                    std::string             & errText,
                    bool                      nopadding)
{
  SegY * segy = NULL;
  bool failed = false;
  target = NULL;

  try
  {
    //
    // Currently we have only one optional TraceHeaderFormat, but this can 
    // be augmented to a list with several formats ...
    //
    if(format == NULL) { //Unknown format
      std::vector<TraceHeaderFormat*> traceHeaderFormats(0);
      if (modelSettings->getTraceHeaderFormat() != NULL)
      {
        traceHeaderFormats.push_back(modelSettings->getTraceHeaderFormat());
      }

      segy = new SegY(fileName, 
                      offset, 
                      traceHeaderFormats, 
                      true); // Add standard formats to format search
    }
    else //Known format, read directly.
      segy = new SegY(fileName, offset, *format);

    bool onlyVolume = modelSettings->getAreaParameters() != NULL; // This is now always true
    segy->ReadAllTraces(timeSimbox, 
                        modelSettings->getZPadFac(),
                        onlyVolume);
    segy->CreateRegularGrid();
  }
  catch (NRLib::Exception & e)
  {
    errText += e.what();
    failed = true;
  }
  
  int outsideTraces = 0;
  if (failed == false)
  {
    const SegyGeometry * geo;
    geo = segy->GetGeometry();
    geo->WriteGeometry();
    if (gridType == FFTGrid::DATA) 
      geometry = new SegyGeometry(geo);
    
    int xpad, ypad, zpad;
    if(nopadding==false)
    {
      xpad = modelSettings->getNXpad();
      ypad = modelSettings->getNYpad();
      zpad = modelSettings->getNZpad();
    }
    else
    {
      xpad = timeSimbox->getnx();
      ypad = timeSimbox->getny();
      zpad = timeSimbox->getnz();
    }
    target = createFFTGrid(timeSimbox->getnx(), 
                           timeSimbox->getny(), 
                           timeSimbox->getnz(), 
                           xpad, 
                           ypad, 
                           zpad, 
                           modelSettings->getFileGrid());
    target->setType(gridType);
    outsideTraces = target->fillInFromSegY(segy, timeSimbox, nopadding);
  }
  if (segy != NULL)
    delete segy;
  return(outsideTraces);
}


int
Model::readStormFile(const std::string  & fName, 
                     FFTGrid           *& target, 
                     const int            gridType,
                     const std::string  & parName, 
                     Simbox             * timeSimbox, 
                     ModelSettings     *& modelSettings, 
                     std::string        & errText,
                     bool                 scale,
                     bool                 nopadding)
{
  StormContGrid * stormgrid = NULL;
  bool failed = false;
  
  try
  {   
    stormgrid = new StormContGrid(0,0,0);
    stormgrid->ReadFromFile(fName);
    stormgrid->SetMissingCode(RMISSING);
  }
  catch (NRLib::Exception & e) 
  {
    errText += e.what();
    failed = true;
  }
  int xpad, ypad, zpad;
  if(nopadding==false)
  {
    xpad = modelSettings->getNXpad();
    ypad = modelSettings->getNYpad();
    zpad = modelSettings->getNZpad();
  }
  else
  {
    xpad = timeSimbox->getnx();
    ypad = timeSimbox->getny();
    zpad = timeSimbox->getnz();
  }
  int outsideTraces = 0;
  if(failed == false)
  {
    target = createFFTGrid(timeSimbox->getnx(), 
                           timeSimbox->getny(), 
                           timeSimbox->getnz(), 
                           xpad,
                           ypad,
                           zpad,
                           modelSettings->getFileGrid());
    target->setType(gridType);
    outsideTraces = target->fillInFromStorm(timeSimbox,stormgrid, parName, scale);
  }  

  if (stormgrid != NULL)
    delete stormgrid;

  return(outsideTraces);
}

int 
Model::setPaddingSize(int nx, double px)
{
  int leastint    =  static_cast<int>(ceil(nx*(1.0f+px)));
  int closestprod =  FFTGrid::findClosestFactorableNumber(leastint);
  return(closestprod);
}

void 
Model::makeTimeSimboxes(Simbox        *& timeSimbox,
                        Simbox        *& timeCutSimbox,
                        Simbox        *& timeBGSimbox,
                        Simbox        *& timeSimboxConstThick,
                        Surface       *& correlationDirection,
                        ModelSettings *& modelSettings, 
                        InputFiles     * inputFiles,
                        std::string    & errText,
                        bool           & failed)
{
  std::string gridFile("");

  int  areaSpecification      = modelSettings->getAreaSpecification();

  bool estimationModeNeedILXL = modelSettings->getEstimationMode() && 
                                (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA ||
                                 (modelSettings->getGridOutputFlag() & IO::ORIGINAL_SEISMIC_DATA) > 0);
                                 
  if(modelSettings->getForwardModeling())
    gridFile = inputFiles->getBackFile(0);    // Get geometry from earth model (Vp)
  else {
    if (modelSettings->getEstimationMode() == false || estimationModeNeedILXL)
      gridFile = inputFiles->getSeismicFile(0); // Get area from first seismic data volume
  }
  SegyGeometry * ILXLGeometry = NULL; //Geometry with correct XL and IL settings.
  //
  // Set area geometry information
  // -----------------------------
  //
  std::string areaType;
  if (areaSpecification == ModelSettings::AREA_FROM_UTM) 
  {
    // The geometry is already present in modelSettings (read from model file).
    LogKit::LogFormatted(LogKit::HIGH,"\nArea information has been taken from model file\n");
    areaType = "Model file";
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_SURFACE) 
  {
    LogKit::LogFormatted(LogKit::HIGH,"\nFinding area information from surface \'"+inputFiles->getAreaSurfaceFile()+"\'\n");
    areaType = "Surface";
    Surface surf(inputFiles->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    modelSettings->setAreaParameters(&geometry);
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_GRID_DATA) 
  {
    LogKit::LogFormatted(LogKit::HIGH,"\nFinding inversion area from grid data in file \'"+gridFile+"\'\n");
    areaType = "Grid data";
    std::string tmpErrText;
    SegyGeometry * geometry;
    getGeometryFromGridOnFile(gridFile,
                              modelSettings->getTraceHeaderFormat(0),
                              geometry,
                              tmpErrText);

    if(geometry != NULL) {
      if(modelSettings->getAreaILXL().size() > 0) {
        SegyGeometry * fullGeometry = geometry;
        bool interpolated, aligned;
        try {
          geometry = fullGeometry->GetILXLSubGeometry(modelSettings->getAreaILXL(), interpolated, aligned);
          std::string text;
          if(interpolated == true) {
            if(aligned == true) {
              text  = "Check IL/XL specification: Specified IL- or XL-step is not an integer multiple\n";
              text += "   of those found in the seismic data. Furthermore, the distance between first\n";
              text += "   and last XL and/or IL does not match the step size.\n";
              TaskList::addTask(text);
            }
            else {
              text  = "Check IL/XL specifiaction: Specified IL- or XL-step is not an integer multiple\n";
              text  = "   of those found in the seismic data.\n";
              TaskList::addTask(text);
            }
          }
          else if(aligned == true) {
            text  = "Check IL/XL specification: Either start or end of IL and/or XL interval does not\n";
            text += "   align with IL/XL in seismic data, or end IL and/or XL is not an integer multiple\n";
            text += "   of steps away from the start.\n";
            TaskList::addTask(text);
          }
        }
        catch (NRLib::Exception & e) {
          errText += "Error: "+std::string(e.what());
          geometry->WriteGeometry(true);
          geometry = NULL;
          failed = true;
        }
        delete fullGeometry;
      }
      if(!failed) {
        modelSettings->setAreaParameters(geometry);
        ILXLGeometry = geometry;
      }
    }
    else {
      errText += tmpErrText;
      failed = true;
    }
  }

  if(!failed)
  {
    const SegyGeometry * areaParams = modelSettings->getAreaParameters(); 
    failed = timeSimbox->setArea(areaParams, errText);

    if(failed)
    {
      writeAreas(areaParams,timeSimbox,areaType);
      errText += "The specified AREA extends outside the surface(s).\n";
    }
    else
    {
      LogKit::LogFormatted(LogKit::LOW,"\nResolution                x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::LOW,"-------------------------------------------------------------------------------------------------\n");
      double azimuth = (-1)*timeSimbox->getAngle()*(180.0/M_PI);
      if (azimuth < 0)
        azimuth += 360.0;
      LogKit::LogFormatted(LogKit::LOW,"%-12s     %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n", 
                           areaType.c_str(),
                           timeSimbox->getx0(), timeSimbox->gety0(), 
                           timeSimbox->getlx(), timeSimbox->getly(), azimuth, 
                           timeSimbox->getdx(), timeSimbox->getdy());
    }

    if(!failed)
    {
      //
      // Set IL/XL information in geometry
      // ---------------------------------
      //
      // Skip for estimation mode if possible:
      //   a) For speed 
      //   b) Grid data may not be available.
      if (modelSettings->getEstimationMode() == false || estimationModeNeedILXL == true) {
        if(ILXLGeometry == NULL) {
          int gridType = IO::findGridType(gridFile);
          bool ilxl_info_available = ((gridType == IO::SEGY) || (gridType == IO::CRAVA));
          if (ilxl_info_available) {
            LogKit::LogFormatted(LogKit::HIGH,"\nFinding IL/XL information from grid data file \'"+gridFile+"\'\n");
            std::string tmpErrText;
            getGeometryFromGridOnFile(gridFile,
                                      modelSettings->getTraceHeaderFormat(0),
                                      ILXLGeometry,
                                      tmpErrText);
            if(ILXLGeometry == NULL) {
              errText += tmpErrText;
              failed = true;
            }
          }
          else {
            LogKit::LogFormatted(LogKit::HIGH,"\nCannot extract IL/XL information from non-SEGY grid data file \'"+gridFile+"\'\n");
          }
        }
        if(ILXLGeometry != NULL) {
          if(timeSimbox->isAligned(ILXLGeometry))
            timeSimbox->setILXL(ILXLGeometry);
          delete ILXLGeometry;
        }
      }
      
      // Rotate variograms relative to simbox
      modelSettings->rotateVariograms(static_cast<float> (timeSimbox_->getAngle()));
      
      //
      // Set SURFACES
      //
      int outputFormat = modelSettings->getGridOutputFormat();
      int outputDomain = modelSettings->getGridOutputDomain();
      int outputFlag   = modelSettings->getGridOutputFlag();
      int otherOutput  = modelSettings->getOtherOutputFlag();
      
      setSimboxSurfaces(timeSimbox, 
                        inputFiles->getTimeSurfFiles(), 
                        modelSettings->getForwardModeling(),
                        modelSettings->getEstimationMode(),
                        modelSettings->getGenerateBackground(),
                        modelSettings->getParallelTimeSurfaces(), 
                        modelSettings->getTimeDTop(), 
                        modelSettings->getTimeLz(), 
                        modelSettings->getTimeDz(), 
                        modelSettings->getTimeNz(),
                        outputFormat,
                        outputDomain,
                        outputFlag,
                        errText,
                        failed);
 
      if(!failed)
      {
        if(modelSettings->getUseLocalWavelet() && timeSimbox->getIsConstantThick())
        {
          LogKit::LogFormatted(LogKit::WARNING,"\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");
          TaskList::addTask("If local wavelet is to be used, constant thickness in depth should be removed.");
        }


        int status = timeSimbox->calculateDz(modelSettings->getLzLimit(),errText);
        estimateZPaddingSize(timeSimbox, modelSettings);   
        
        if(status == Simbox::BOXOK)
        {
          logIntervalInformation(timeSimbox, "Time output interval:","Two-way-time");
          //
          // Make extended time simbox
          //
          if(inputFiles->getCorrDirFile() != "") {
            //
            // Get correlation direction
            //
            try {
              Surface tmpSurf(inputFiles->getCorrDirFile());
              correlationDirection = new Surface(tmpSurf);
            }
            catch (NRLib::Exception & e) {
              errText += e.what();
              failed = true;
            }

            if(failed == false && modelSettings->getForwardModeling() == false) {
              //Extends timeSimbox for correlation coverage. Original stored in timeCutSimbox
              setupExtendedTimeSimbox(timeSimbox, correlationDirection, timeCutSimbox, 
                                      outputFormat, outputDomain, modelSettings->getOtherOutputFlag()); 
            }      

            estimateZPaddingSize(timeSimbox, modelSettings);   

            status = timeSimbox->calculateDz(modelSettings->getLzLimit(),errText);

            if(status == Simbox::BOXOK) 
              logIntervalInformation(timeSimbox, "Time inversion interval (extended relative to output interval due to correlation):","Two-way-time");
            else
            {
              errText += "Could not make the time simulation grid.\n";
              failed = true;
            }

            if(modelSettings->getForwardModeling() == false) {
              setupExtendedBackgroundSimbox(timeSimbox, correlationDirection, timeBGSimbox, 
                                            outputFormat, outputDomain, modelSettings->getOtherOutputFlag());
              status = timeBGSimbox->calculateDz(modelSettings->getLzLimit(),errText);
              if(status == Simbox::BOXOK)
                logIntervalInformation(timeBGSimbox, "Time interval used for background modelling:","Two-way-time");
              else
              {
                errText += "Could not make the grid for background model.\n";
                failed = true;
              }
            }
          }
          
          if(failed == false) {
            estimateXYPaddingSizes(timeSimbox, modelSettings);
            
            LogKit::LogFormatted(LogKit::LOW,"\nTime simulation grids:\n");
            LogKit::LogFormatted(LogKit::LOW,"  Output grid         %4i * %4i * %4i   : %10i\n",
                                 timeSimbox->getnx(),timeSimbox->getny(),timeSimbox->getnz(),
                                 timeSimbox->getnx()*timeSimbox->getny()*timeSimbox->getnz()); 
            LogKit::LogFormatted(LogKit::LOW,"  FFT grid            %4i * %4i * %4i   : %10i\n",
                                 modelSettings->getNXpad(),modelSettings->getNYpad(),modelSettings->getNZpad(),
                                 modelSettings->getNXpad()*modelSettings->getNYpad()*modelSettings->getNZpad());
          }
          
          //
          // Make time simbox with constant thicknesses (needed for log filtering and facies probabilities)
          //
          timeSimboxConstThick = new Simbox(timeSimbox);
          Surface * tsurf = new Surface(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface()));
          timeSimboxConstThick->setDepth(tsurf, 0, timeSimbox->getlz(), timeSimbox->getdz());
          
          if((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
            std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_ConstThick";
            std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_ConstThick";
            timeSimboxConstThick->writeTopBotGrids(topSurf, 
                                                   baseSurf,
                                                   IO::PathToInversionResults(),
                                                   outputFormat);
          }
        }
        else
        {
          errText += "Could not make time simulation grid.\n";
          failed = true;
        }
      }
      else
      {
        timeSimbox->externalFailure();
        failed = true;
      }
    }
  }
}
  
void 
Model::logIntervalInformation(const Simbox      * simbox, 
                              const std::string & header_text1, 
                              const std::string & header_text2)
{
  LogKit::LogFormatted(LogKit::LOW,"\n"+header_text1+"\n");
  double zmin, zmax;
  simbox->getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::LOW," %13s          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       header_text2.c_str(),
                       zmin+simbox->getlz()*simbox->getAvgRelThick()*0.5,
                       zmin,zmax); 
  LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
                       simbox->getlz()*simbox->getAvgRelThick(),
                       simbox->getlz()*simbox->getMinRelThick(),
                       simbox->getlz());
  LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
                       simbox->getdz()*simbox->getAvgRelThick(),
                       simbox->getdz(),
                       simbox->getdz()*simbox->getMinRelThick());
}
  
void 
Model::setSimboxSurfaces(Simbox                        *& simbox, 
                         const std::vector<std::string> & surfFile, 
                         bool                             generateSeismic,
                         bool                             estimationMode,
                         bool                             generateBackground,
                         bool                             parallelSurfaces, 
                         double                           dTop,
                         double                           lz, 
                         double                           dz, 
                         int                              nz,
                         int                              outputFormat,
                         int                              outputDomain,
                         int                              outputFlag,
                         std::string                    & errText,
                         bool                           & failed)
{
  const std::string & topName = surfFile[0]; 

  Surface * z0Grid = NULL;
  Surface * z1Grid = NULL;
  try {
    if (NRLib::IsNumber(topName)) {
      // Find the smallest surface that covers the simbox. For simplicity 
      // we use only four nodes (nx=ny=2).
      double xMin, xMax;
      double yMin, yMax;
      findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                  simbox->getlx(), simbox->getly(),
                                  simbox->getAngle(), 
                                  xMin,yMin,xMax,yMax);
      z0Grid = new Surface(xMin-100, yMin-100, xMax-xMin+200, yMax-yMin+200, 2, 2, atof(topName.c_str()));
    } 
    else {
      Surface tmpSurf(topName);
      z0Grid = new Surface(tmpSurf);
    }
  }
  catch (NRLib::Exception & e) {
    errText += e.what();
    failed = true;
  }

  if(!failed) {
    if(parallelSurfaces) { //Only one reference surface
      simbox->setDepth(z0Grid, dTop, lz, dz);
    }
    else {
      const std::string & baseName = surfFile[1]; 
      try {
        if (NRLib::IsNumber(baseName)) {
          // Find the smallest surface that covers the simbox. For simplicity 
          // we use only four nodes (nx=ny=2).
          double xMin, xMax;
          double yMin, yMax;
          findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                      simbox->getlx(), simbox->getly(),
                                      simbox->getAngle(), 
                                      xMin,yMin,xMax,yMax);
          z1Grid = new Surface(xMin-100, yMin-100, xMax-xMin+200, yMax-yMin+200, 2, 2, atof(baseName.c_str()));
        }
        else {
          Surface tmpSurf(baseName);
          z1Grid = new Surface(tmpSurf);
        }
      }
      catch (NRLib::Exception & e) {
        errText += e.what();
        failed = true;
      }
      if(!failed) {
        try {
          simbox->setDepth(z0Grid, z1Grid, nz);
        }
        catch (NRLib::Exception & e) {
          errText += e.what();
          std::string text(std::string("Seismic data"));
          writeAreas(modelSettings_->getAreaParameters(),simbox,text);
          failed = true;
        }
      }
    }
    if (!failed) {
      if((outputDomain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        if (generateSeismic) {
          simbox->writeTopBotGrids(topSurf, 
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   outputFormat);
        }
        else if (!estimationMode){
          simbox->writeTopBotGrids(topSurf, 
                                   baseSurf,
                                   IO::PathToInversionResults(),
                                   outputFormat);
        }
        if((outputFormat & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((outputFlag & IO::BACKGROUND) > 0 || 
              (outputFlag & IO::BACKGROUND_TREND) > 0 || 
              estimationMode && generateBackground) {
            simbox->writeTopBotGrids(topSurf, 
                                     baseSurf,
                                     IO::PathToBackground(),
                                     outputFormat);
          }
          if ((outputFlag & IO::CORRELATION) > 0) {
            simbox->writeTopBotGrids(topSurf, 
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     outputFormat);
          }
          if ((outputFlag & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            simbox->writeTopBotGrids(topSurf, 
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     outputFormat);
          }
          if ((outputFlag & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            simbox->writeTopBotGrids(topSurf, 
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     outputFormat);          
          }
        }
      }
    }
  }
}

void
Model::setupExtendedTimeSimbox(Simbox   * timeSimbox, 
                               Surface  * corrSurf, 
                               Simbox  *& timeCutSimbox,
                               int        outputFormat,
                               int        outputDomain,
                               int        otherOutput)
{
  timeCutSimbox = new Simbox(timeSimbox);
  double * corrPlanePars = findPlane(corrSurf);
  Surface * meanSurf = new Surface(*corrSurf);
  int i;
  for(i=0;i<static_cast<int>(meanSurf->GetN());i++)
    (*meanSurf)(i) = 0;

  meanSurf->AddNonConform(&(timeSimbox->GetTopSurface()));
  meanSurf->AddNonConform(&(timeSimbox->GetBotSurface()));
  meanSurf->Multiply(0.5);
  double * refPlanePars = findPlane(meanSurf);
  delete meanSurf;

  for(i=0;i<3;i++)
    refPlanePars[i] -= corrPlanePars[i];
  gradX_ = refPlanePars[1];
  gradY_ = refPlanePars[2];

  Surface * refPlane = createPlaneSurface(refPlanePars, corrSurf);

  std::string fileName = "Correlation_Rotation_Plane";
  IO::writeSurfaceToFile(*refPlane, fileName, IO::PathToCorrelations(), outputFormat);

  refPlane->AddNonConform(corrSurf);
  delete [] corrPlanePars;
  delete [] refPlanePars;

  Surface * topSurf = new Surface(*refPlane);
  topSurf->SubtractNonConform(&(timeSimbox->GetTopSurface()));
  double shiftTop = topSurf->Max();
  shiftTop *= -1.0;
  topSurf->Add(shiftTop);
  topSurf->AddNonConform(&(timeSimbox->GetTopSurface()));

  Surface * botSurf = new Surface(*refPlane);
  botSurf->SubtractNonConform(&(timeSimbox->GetBotSurface()));
  double shiftBot = botSurf->Min();
  shiftBot *= -1.0;
  double thick    = shiftBot-shiftTop;
  double dz       = timeCutSimbox->getdz();
  int    nz       = int(thick/dz);
  double residual = thick - nz*dz;
  if (residual > 0.0) {
    shiftBot += dz-residual;
    nz++;
  }
  if (nz != timeCutSimbox->getnz()) {
    LogKit::LogFormatted(LogKit::HIGH,"\nNumber of layers in inversion increased from %d",timeCutSimbox->getnz());
    LogKit::LogFormatted(LogKit::HIGH," to %d in grid created using correlation direction.\n",nz);
  }
  botSurf->Add(shiftBot);
  botSurf->AddNonConform(&(timeSimbox->GetBotSurface()));

  timeSimbox->setDepth(topSurf, botSurf, nz);
  
  if((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Extended";
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Extended";
    timeSimbox->writeTopBotGrids(topSurf, 
                                 baseSurf,
                                 IO::PathToInversionResults(),
                                 outputFormat);
  }

  delete refPlane;
}

void
Model::setupExtendedBackgroundSimbox(Simbox   * timeSimbox, 
                                     Surface  * corrSurf, 
                                     Simbox  *& timeBGSimbox,
                                     int        outputFormat,
                                     int        outputDomain,
                                     int        otherOutput)
{
  //
  // Move correlation surface for easier handling.
  //
  Surface * tmpSurf = new Surface(*corrSurf);
  double avg = tmpSurf->Avg();
  if (avg > 0)
    tmpSurf->Subtract(avg);
  else
    tmpSurf->Add(avg); // This situation is not very likely, but ...

  //
  // Find top surface of background simbox. 
  //
  // The funny/strange dTop->Multiply(-1.0) is due to NRLIB's current
  // inability to set dTop equal to Simbox top surface.
  //
  Surface * dTop = new Surface(*tmpSurf);
  dTop->SubtractNonConform(&(timeSimbox->GetTopSurface()));
  dTop->Multiply(-1.0);
  double shiftTop = dTop->Min();
  delete dTop;
  Surface * topSurf = new Surface(*tmpSurf);
  topSurf->Add(shiftTop);

  //
  // Find base surface of background simbox
  //
  Surface * dBot = new Surface(*tmpSurf);
  dBot->SubtractNonConform(&(timeSimbox->GetBotSurface()));
  dBot->Multiply(-1.0);
  double shiftBot = dBot->Max();
  delete dBot;
  Surface * botSurf = new Surface(*tmpSurf);
  botSurf->Add(shiftBot);

  //
  // Calculate number of layers of background simbox 
  //
  tmpSurf->Assign(0.0);
  tmpSurf->AddNonConform(botSurf);
  tmpSurf->SubtractNonConform(topSurf);
  double dMax = tmpSurf->Max();
  double dt = timeSimbox->getdz();
  int nz;
  //
  // NBNB-PAL: I think it is a good idea to use a maximum dt of 10ms.
  //
  //if (dt < 10.0) {
  //  LogKit::LogFormatted(LogKit::HIGH,"\nReducing sampling density for background",dt);
  //  LogKit::LogFormatted(LogKit::HIGH," modelling from %.2fms to 10.0ms\n");
  //  dt = 10.0;  // A sampling density of 10.0ms is good enough for BG model
  // }
  nz = static_cast<int>(ceil(dMax/dt));
  delete tmpSurf;

  //
  // Make new simbox
  //
  timeBGSimbox = new Simbox(timeSimbox);
  timeBGSimbox->setDepth(topSurf, botSurf, nz);

  if((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_BG";
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_BG";
    timeBGSimbox->writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToBackground(),
                                   outputFormat);
  }
}

double *
Model::findPlane(Surface * surf)
{
  double ** A = new double * [3];
  double * b = new double[3];
  int i, j;
  for(i=0;i<3;i++) {
    A[i] = new double[3];
    for(j=0;j<3;j++)
      A[i][j] = 0;
    b[i] = 0;
  }

  double x, y, z;
  int nData = 0;
  for(i=0;i<static_cast<int>(surf->GetN());i++) {
    surf->GetXY(i, x, y);
    z = (*surf)(i);
    if(!surf->IsMissing(z)) {
      nData++;
      A[0][1] += x;
      A[0][2] += y;
      A[1][1] += x*x;
      A[1][2] += x*y;
      A[2][2] += y*y;
      b[0] += z;
      b[1] += x*z;
      b[2] += y*z;
    }
  }

  A[0][0] = nData;
  A[1][0] = A[0][1];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];

  lib_matrCholR(3, A);
  lib_matrAxeqbR(3, A, b);

  for(i=0;i<3;i++)
    delete [] A[i];
  delete [] A;

  return(b);
}


Surface *
Model::createPlaneSurface(double * planeParams, Surface * templateSurf)
{
  Surface * result = new Surface(*templateSurf);
  double x,y;
  int i;
  for(i=0;i<static_cast<int>(result->GetN());i++) {
    result->GetXY(i,x,y);
    (*result)(i) = planeParams[0]+planeParams[1]*x+planeParams[2]*y;
  }
  return(result);
}


void
Model::estimateXYPaddingSizes(Simbox         * timeSimbox, 
                              ModelSettings *& modelSettings)
{
  double dx      = timeSimbox->getdx();
  double dy      = timeSimbox->getdy();
  double lx      = timeSimbox->getlx();
  double ly      = timeSimbox->getly();
  int    nx      = timeSimbox->getnx();
  int    ny      = timeSimbox->getny();

  double xPadFac = modelSettings->getXPadFac();
  double yPadFac = modelSettings->getYPadFac();
  double xPad    = xPadFac*lx;
  double yPad    = yPadFac*ly;

  if (modelSettings->getEstimateXYPadding())
  {
    float  range1 = modelSettings->getLateralCorr()->getRange();
    float  range2 = modelSettings->getLateralCorr()->getSubRange();
    float  angle  = modelSettings->getLateralCorr()->getAngle();
    double factor = 0.5;  // Lateral correlation is not very important. Half a range is probably more than enough

    xPad          = factor * std::max(fabs(range1*cos(angle)), fabs(range2*sin(angle)));
    yPad          = factor * std::max(fabs(range1*sin(angle)), fabs(range2*cos(angle)));
    xPad          = std::max(xPad, dx);     // Always require at least on grid cell
    yPad          = std::max(yPad, dy);     // Always require at least one grid cell
    xPadFac       = std::min(1.0, xPad/lx); // A padding of more than 100% is insensible
    yPadFac       = std::min(1.0, yPad/ly);
  }

  int nxPad = setPaddingSize(nx, xPadFac);
  int nyPad = setPaddingSize(ny, yPadFac);
  xPadFac   = static_cast<double>(nxPad - nx)/static_cast<double>(nx); // Update xPadFac
  yPadFac   = static_cast<double>(nyPad - ny)/static_cast<double>(ny); // Update yPadFac
  xPad      = xPadFac*lx;                                              // Update xPad
  yPad      = yPadFac*ly;                                              // Update yPad

  modelSettings->setNXpad(nxPad);
  modelSettings->setNYpad(nyPad);
  modelSettings->setXPadFac(xPadFac);
  modelSettings->setYPadFac(yPadFac);

  std::string text1;
  std::string text2;
  int logLevel = LogKit::MEDIUM;
  if (modelSettings->getEstimateXYPadding()) {
    text1 = " estimated from lateral correlation ranges in internal grid";
    logLevel = LogKit::LOW;
  }  
  if (modelSettings->getEstimateZPadding()) {
    text2 = " estimated from an assumed wavelet length"; 
    logLevel = LogKit::LOW;
  }  

  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text1+":\n");
  LogKit::LogFormatted(logLevel,"  xPad, xPadFac, nx, nxPad                 : %6.fm, %4.2f, %5d, %4d\n", 
                       xPad, xPadFac, timeSimbox->getnx(), nxPad);
  LogKit::LogFormatted(logLevel,"  yPad, yPadFac, ny, nyPad                 : %6.fm, %4.2f, %5d, %4d", 
                       yPad, yPadFac, timeSimbox->getny(), nyPad);
  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text2+":\n");
  LogKit::LogFormatted(logLevel,"  zPad, zPadFac, nz, nzPad                 : %5.fms, %4.2f, %5d, %4d\n", 
                       modelSettings->getZPadFac()*(timeSimbox->getlz()*timeSimbox->getMinRelThick()), 
                       modelSettings->getZPadFac(), timeSimbox->getnz(), modelSettings->getNZpad());
}

void
Model::estimateZPaddingSize(Simbox         * timeSimbox,
                            ModelSettings *& modelSettings)
{
  int    nz          = timeSimbox->getnz();
  double minLz       = timeSimbox->getlz()*timeSimbox->getMinRelThick();
  double zPadFac     = modelSettings->getZPadFac();
  double zPad        = zPadFac*minLz;

  if (modelSettings->getEstimateZPadding())
  {
    double wLength = 200.0;                     // Assume a wavelet is approx 200ms.
    zPad           = wLength/2.0;               // Use half a wavelet as padding
    zPadFac        = std::min(1.0, zPad/minLz); // More than 100% padding is not sensible
  }

  int nzPad = setPaddingSize(nz, zPadFac);
  zPadFac   = static_cast<double>(nzPad - nz)/static_cast<double>(nz); // Update yPadFac
 
  modelSettings->setNZpad(nzPad);
  modelSettings->setZPadFac(zPadFac);
}

void
Model::processSeismic(FFTGrid      **& seisCube,
                      Simbox        *& timeSimbox,
                      ModelSettings *& modelSettings, 
                      InputFiles     * inputFiles,
                      std::string    & errText,
                      bool           & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  if(inputFiles->getNumberOfSeismicFiles() > 0)
  {
    Utils::writeHeader("Reading seismic data");

    int nAngles = modelSettings->getNumberOfAngles();
    const SegyGeometry ** geometry = new const SegyGeometry * [nAngles];
    seisCube = new FFTGrid * [nAngles];

    bool outsideWarning = false;
    for (int i = 0 ; i < nAngles ; i++) {
      geometry[i] = NULL;
      std::string tmpErrText("");
      std::string fileName = inputFiles->getSeismicFile(i);
      std::string dataName = "Seismic data angle stack"+NRLib::ToString(i);
      float       offset = modelSettings->getLocalSegyOffset(i);
      if(offset < 0)
        offset = modelSettings->getSegyOffset();

      int outsideTraces = readGridFromFile(fileName,
                                           dataName,
                                           offset,
                                           seisCube[i],
                                           geometry[i],
                                           modelSettings->getTraceHeaderFormat(i),
                                           FFTGrid::DATA,
                                           timeSimbox,
                                           modelSettings,
                                           tmpErrText);
      if(tmpErrText != "")
      {
        tmpErrText += "Reading of file \'"+fileName+"\' for "+dataName+" failed.\n";
        errText += tmpErrText;
        failed = true;
      } 
      else { 
        seisCube[i]->setAngle(modelSettings->getAngle(i));
        if(outsideTraces > 0) {
          if(outsideTraces == seisCube[i]->getNxp()*seisCube[i]->getNyp()) {
            errText += "Error: Data in file "+fileName+" was completely outside the inversion area.\n";
            failed = true;
          }
          else {
            LogKit::LogMessage(LogKit::WARNING, "Warning: "+NRLib::ToString(outsideTraces)+" traces in the grid were outside the data area in file "
              +fileName+". Note that this includes traces in the padding.\n");
            outsideWarning = true;
          }
        }
      }
    }
    LogKit::LogFormatted(LogKit::LOW,"\n");
    if(outsideWarning == true)
      TaskList::addTask("Check seismic volumes and inversion area: One or more of the seismic input files did not have data enough for the entire area (including padding).\n");

    if(failed == false)
    {
      bool segyVolumesRead = false;
      for (int i = 0 ; i < nAngles ; i++) 
      {
        if (geometry[i] != NULL)
          segyVolumesRead = true;
      }
      if (segyVolumesRead) 
      {
        LogKit::LogFormatted(LogKit::LOW,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
        LogKit::LogFormatted(LogKit::LOW,"-------------------------------------------------------------------------------------------------\n");
        for (int i = 0 ; i < nAngles ; i++)
        {   
          if (geometry[i] != NULL) {
            double geoAngle = (-1)*timeSimbox->getAngle()*(180/M_PI);
            if (geoAngle < 0)
              geoAngle += 360.0f;
            LogKit::LogFormatted(LogKit::LOW,"Seismic data %d   %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",i,
                                 geometry[i]->GetX0(), geometry[i]->GetY0(), 
                                 geometry[i]->Getlx(), geometry[i]->Getly(), geoAngle,
                                 geometry[i]->GetDx(), geometry[i]->GetDy());
          }
        }
      }

      if((modelSettings->getGridOutputFlag() & IO::ORIGINAL_SEISMIC_DATA) > 0) {
        for(int i=0;i<nAngles;i++) {
          std::string angle    = NRLib::ToString(modelSettings->getAngle(i)*(180/M_PI), 1);
          std::string baseName = IO::PrefixOriginalSeismicData() + angle;
          std::string sgriLabel = std::string("Original seismic data for angle stack ") + angle;
          float       offset = modelSettings->getLocalSegyOffset(i);
          if(offset < 0)
            offset = modelSettings->getSegyOffset();
          seisCube[i]->writeFile(baseName, 
                                 IO::PathToSeismicData(), 
                                 timeSimbox, 
                                 sgriLabel, 
                                 offset,
                                 timeDepthMapping_, 
                                 timeCutMapping_, *modelSettings->getTraceHeaderFormatOutput());

        }
      }

      for (int i =0 ; i < nAngles ; i++)
        if (geometry[i] != NULL)
          delete geometry[i];
      delete [] geometry;
    }
  }
  Timings::setTimeSeismic(wall,cpu);
}

void 
Model::processWells(WellData     **& wells,
                    Simbox         * timeSimbox,
                    Simbox         * timeBGSimbox,
                    Simbox         * timeSimboxConstThick,
                    RandomGen      * randomGen,
                    ModelSettings *& modelSettings, 
                    InputFiles     * inputFiles,
                    std::string    & errText,
                    bool           & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  Utils::writeHeader("Reading and processing wells");

  bool    faciesLogGiven = modelSettings->getFaciesLogGiven();
  int     nWells         = modelSettings->getNumberOfWells();
  int     nFacies        = 0;

  int error = 0;

  std::string tmpErrText("");
  wells = new WellData *[nWells];
  for(int i=0 ; i<nWells ; i++) {
    wells[i] = new WellData(inputFiles->getWellFile(i), 
                            modelSettings->getLogNames(),
                            modelSettings->getInverseVelocity(),
                            modelSettings, 
                            modelSettings->getIndicatorFacies(i),
                            modelSettings->getIndicatorWavelet(i),
                            modelSettings->getIndicatorBGTrend(i),
                            modelSettings->getIndicatorRealVs(i),
                            faciesLogGiven);
    if(wells[i]->checkError(tmpErrText) != 0) {
      errText += tmpErrText;
      error = 1;
    }
  }

  if (error == 0) {
    if(modelSettings->getFaciesLogGiven()) { 
      checkFaciesNames(wells, modelSettings, inputFiles, tmpErrText, error);
      nFacies = modelSettings->getNumberOfFacies(); // nFacies is set in checkFaciesNames()
    }
    if (error>0)
      errText += "Prior facies probabilities failed.\n"+tmpErrText;

    int   * validWells    = new int[nWells];
    bool  * validIndex    = new bool[nWells];
    int   * nMerges       = new int[nWells];
    int   * nInvalidAlpha = new int[nWells];
    int   * nInvalidBeta  = new int[nWells];
    int   * nInvalidRho   = new int[nWells];
    float * rankCorr      = new float[nWells];
    float * devAngle      = new float[nWells];
    int  ** faciesCount   = NULL;
    
    if(nFacies > 0) { 
      faciesCount = new int * [nWells]; 
      for (int i = 0 ; i < nWells ; i++)
        faciesCount[i] = new int[nFacies];
    }
    
    int count = 0;
    int nohit=0;
    int empty=0;
    int facieslognotok = 0;
    LogKit::LogFormatted(LogKit::LOW,"\n");
    for (int i=0 ; i<nWells ; i++)
    {
      bool skip = false;
      LogKit::LogFormatted(LogKit::LOW,wells[i]->getWellname()+" : \n");
      if(wells[i]!=NULL) {
        if(wells[i]->checkSimbox(timeSimbox) == 1) {
          skip = true;
          nohit++;
          TaskList::addTask("Consider increasing the inversion volume such that well "+wells[i]->getWellname()+ " can be included");
        }
        if(wells[i]->getNd() == 0) {
          LogKit::LogFormatted(LogKit::LOW,"  IGNORED (no log entries found)\n");
          skip = true;
          empty++;
          TaskList::addTask("Check the log entries in well "+wells[i]->getWellname()+".");
        }
        if(wells[i]->isFaciesOk()==0) {
          LogKit::LogFormatted(LogKit::LOW,"   IGNORED (facies log has wrong entries)\n");
          skip = true;
          facieslognotok++;
          TaskList::addTask("Check the facies logs in well "+wells[i]->getWellname()+".\n       The facies logs in this well are wrong and the well is ignored");
        }
        if(skip)
          validIndex[i] = false;
        else {
          validIndex[i] = true;
          wells[i]->removeDuplicateLogEntries(nMerges[i]);
          wells[i]->setWrongLogEntriesUndefined(nInvalidAlpha[i], nInvalidBeta[i], nInvalidRho[i]);
          wells[i]->filterLogs();
          wells[i]->lookForSyntheticVsLog(rankCorr[i]);
          wells[i]->calculateDeviation(devAngle[i], timeSimbox);
          wells[i]->setBlockedLogsOrigThick( new BlockedLogs(wells[i], timeSimbox, randomGen) );
          wells[i]->setBlockedLogsConstThick( new BlockedLogs(wells[i], timeSimboxConstThick, randomGen) );
          if (timeBGSimbox==NULL)
            wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeSimbox, randomGen) ); // Need a copy constructor?
          else
            wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeBGSimbox, randomGen) );
          if (nFacies > 0)
            wells[i]->countFacies(timeSimbox,faciesCount[i]);
          validWells[count] = i;
          count++;      
        }
      }
    }
    //
    // Write summary.
    //
    LogKit::LogFormatted(LogKit::LOW,"\n");
    LogKit::LogFormatted(LogKit::LOW,"                                      Invalid                                    \n");
    LogKit::LogFormatted(LogKit::LOW,"Well                    Merges      Vp   Vs  Rho  synthVs/Corr    Deviated/Angle \n");
    LogKit::LogFormatted(LogKit::LOW,"---------------------------------------------------------------------------------\n");
    for(int i=0 ; i<nWells ; i++) {
      if (validIndex[i]) 
        LogKit::LogFormatted(LogKit::LOW,"%-23s %6d    %4d %4d %4d     %3s / %5.3f      %3s / %4.1f\n",
                             wells[i]->getWellname().c_str(),
                             nMerges[i],
                             nInvalidAlpha[i], 
                             nInvalidBeta[i], 
                             nInvalidRho[i],
                             (wells[i]->hasSyntheticVsLog() ? "yes" : " no"),
                             rankCorr[i],
                             (devAngle[i] > modelSettings->getMaxDevAngle() ? "yes" : " no"),
                             devAngle[i]);
      else  
        LogKit::LogFormatted(LogKit::LOW,"%-23s      -       -    -    -       - /     -       -  /    -\n",
                             wells[i]->getWellname().c_str());
    }
    
    //
    // Print facies count for each well
    //
    if(nFacies > 0) { 
      //
      // Probabilities
      //
      LogKit::LogFormatted(LogKit::LOW,"\nFacies distributions for each well: \n");
      LogKit::LogFormatted(LogKit::LOW,"\nWell                    ");
      for (int i = 0 ; i < nFacies ; i++)
        LogKit::LogFormatted(LogKit::LOW,"%12s ",modelSettings->getFaciesName(i).c_str());
      LogKit::LogFormatted(LogKit::LOW,"\n");
      for (int i = 0 ; i < 24+13*nFacies ; i++)
        LogKit::LogFormatted(LogKit::LOW,"-");
      LogKit::LogFormatted(LogKit::LOW,"\n");
      for (int i = 0 ; i < nWells ; i++) {
        if (validIndex[i]) {
          float tot = 0.0;
          for (int f = 0 ; f < nFacies ; f++)
            tot += static_cast<float>(faciesCount[i][f]);
          LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[i]->getWellname().c_str());
          for (int f = 0 ; f < nFacies ; f++) {
            if (tot > 0) {
              float faciesProb = static_cast<float>(faciesCount[i][f])/tot;
              LogKit::LogFormatted(LogKit::LOW,"%12.4f ",faciesProb);
            }
            else
              LogKit::LogFormatted(LogKit::LOW,"         -   ");
          }
          LogKit::LogFormatted(LogKit::LOW,"\n");
        } 
        else {
          LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[i]->getWellname().c_str());
          for (int f = 0 ; f < nFacies ; f++)
            LogKit::LogFormatted(LogKit::LOW,"         -   ");
          LogKit::LogFormatted(LogKit::LOW,"\n");

        }
      }
      LogKit::LogFormatted(LogKit::LOW,"\n");
      //
      // Counts
      //
      LogKit::LogFormatted(LogKit::MEDIUM,"\nFacies counts for each well: \n");
      LogKit::LogFormatted(LogKit::MEDIUM,"\nWell                    ");
      for (int i = 0 ; i < nFacies ; i++)
        LogKit::LogFormatted(LogKit::MEDIUM,"%12s ",modelSettings->getFaciesName(i).c_str());
      LogKit::LogFormatted(LogKit::MEDIUM,"\n");
      for (int i = 0 ; i < 24+13*nFacies ; i++)
        LogKit::LogFormatted(LogKit::MEDIUM,"-");
      LogKit::LogFormatted(LogKit::MEDIUM,"\n");
      for (int i = 0 ; i < nWells ; i++) {
        if (validIndex[i]) {
          float tot = 0.0;
          for (int f = 0 ; f < nFacies ; f++)
            tot += static_cast<float>(faciesCount[i][f]);
          LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[i]->getWellname().c_str());
          for (int f = 0 ; f < nFacies ; f++) {
            LogKit::LogFormatted(LogKit::MEDIUM,"%12d ",faciesCount[i][f]);
          }
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");
        } 
        else {
          LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[i]->getWellname().c_str());
          for (int f = 0 ; f < nFacies ; f++)
            LogKit::LogFormatted(LogKit::MEDIUM,"         -   ");
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");

        }
      }
      LogKit::LogFormatted(LogKit::MEDIUM,"\n");
    }
    
    //
    // Remove invalid wells
    //
    for(int i=0 ; i<nWells ; i++)
      if (!validIndex[i]) 
        delete wells[i];
    for(int i=0 ; i<count ; i++)
      wells[i] = wells[validWells[i]];
    for(int i=count ; i<nWells ; i++)
      wells[i] = NULL;
    nWells = count;
    modelSettings->setNumberOfWells(nWells);
    
    delete [] validWells;
    delete [] validIndex;
    delete [] nMerges;
    delete [] nInvalidAlpha;
    delete [] nInvalidBeta;
    delete [] nInvalidRho;
    delete [] rankCorr;
    delete [] devAngle;
    if(nFacies > 0) { 
      for (int i = 0 ; i<nWells ; i++)
        delete [] faciesCount[i];
      delete [] faciesCount;
    }
    if (nohit>0)
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: %d well(s) do not hit the inversion volume and will be ignored.\n",nohit);
    if (empty>0)
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: %d well(s) contain no log entries and will be ignored.\n",empty);
    if(facieslognotok>0)
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: %d well(s) have wrong facies logs and will be ignored.\n",facieslognotok);
    if (nWells==0 && modelSettings->getNoWellNedded()==false) {
      LogKit::LogFormatted(LogKit::LOW,"\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
      LogKit::LogFormatted(LogKit::LOW,"\n       below is correct. If it is not, you probably have problems with coordinate scaling.");
      LogKit::LogFormatted(LogKit::LOW,"\n                                   X0          Y0        DeltaX      DeltaY      Angle");
      LogKit::LogFormatted(LogKit::LOW,"\n       -------------------------------------------------------------------------------");
      LogKit::LogFormatted(LogKit::LOW,"\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n", 
                           timeSimbox->getx0(), timeSimbox->gety0(), 
                           timeSimbox->getlx(), timeSimbox->getly(), 
                           (timeSimbox->getAngle()*180)/M_PI);
      errText += "No wells available for estimation.";
      error = 1;
    }
  }
  failed = error > 0;
  Timings::setTimeWells(wall,cpu);
}


void Model::addSeismicLogs(WellData     ** wells, 
                           FFTGrid      ** seisCube, 
                           ModelSettings * modelSettings)
{
  int nWells  = modelSettings->getNumberOfWells();
  int nAngles = modelSettings->getNumberOfAngles();
  
    for (int iAngle = 0 ; iAngle < nAngles ; iAngle++)
    {
      seisCube[iAngle]->setAccessMode(FFTGrid::RANDOMACCESS);
      for(int i=0;i<nWells;i++) 
        wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
      seisCube[iAngle]->endAccess();
  }
}
   
void Model::generateSyntheticSeismic(Wavelet      ** wavelet,
                                     WellData     ** wells,
                                     float        ** reflectionMatrix,
                                     Simbox        * timeSimbox,
                                     ModelSettings * modelSettings) const
{
  int nWells  = modelSettings->getNumberOfWells();
  int nAngles = modelSettings->getNumberOfAngles();
  int nzp     = modelSettings->getNZpad();
  int nz      = timeSimbox->getnz();
  int i;

  for( i=0; i<nWells; i++ )
  {
    if( wells[i]->isDeviated() == false )
      wells[i]->getBlockedLogsOrigThick()->generateSyntheticSeismic(reflectionMatrix,nAngles,wavelet,nz,nzp);
  }
}

void Model::writeWells(WellData ** wells, ModelSettings * modelSettings)
{
  int nWells  = modelSettings->getNumberOfWells();
  for(int i=0;i<nWells;i++)
    wells[i]->writeWell(modelSettings->getWellFormatFlag());
}

void Model::writeBlockedWells(WellData ** wells, ModelSettings * modelSettings)
{
  int nWells  = modelSettings->getNumberOfWells();
  for(int i=0;i<nWells;i++)
    wells[i]->getBlockedLogsOrigThick()->writeWell(modelSettings);
}

void Model::checkFaciesNames(WellData      ** wells,
                             ModelSettings *& modelSettings,
                             InputFiles     * inputFiles,
                             std::string    & tmpErrText,
                             int            & error)
{
  int min,max;
  int globalmin = 0;
  int globalmax = 0;
  bool first = true;
  for (int w = 0; w < modelSettings->getNumberOfWells(); w++) {
    if(wells[w]->isFaciesLogDefined())
    {
      wells[w]->getMinMaxFnr(min,max);
      if(first==true)
      {
        globalmin = min;
        globalmax = max;
        first = false;
      }
      else
      {
        if(min<globalmin)
          globalmin = min;
        if(max>globalmax)
          globalmax = max;
      }
    }
  }

  int nnames = globalmax - globalmin + 1;
  std::vector<std::string> names(nnames);
  
  for(int w=0 ; w<modelSettings->getNumberOfWells() ; w++)
  {
    if(wells[w]->isFaciesLogDefined())
    {
      for(int i=0 ; i < wells[w]->getNFacies() ; i++)
      {
        std::string name = wells[w]->getFaciesName(i);
        int         fnr  = wells[w]->getFaciesNr(i) - globalmin;

        if(names[fnr] == "") {
          names[fnr] = name;
        }
        else if(names[fnr] != name)
        {
          tmpErrText += "Problem with facies logs. Facies names and numbers are not uniquely defined.\n";
          error++;
        }
      }
    }
  }

  LogKit::LogFormatted(LogKit::LOW,"\nFaciesLabel      FaciesName           ");
  LogKit::LogFormatted(LogKit::LOW,"\n--------------------------------------\n");
  for(int i=0 ; i<nnames ; i++)
    if(names[i] != "") 
      LogKit::LogFormatted(LogKit::LOW,"    %2d           %-20s\n",i+globalmin,names[i].c_str());
  
  int nFacies = 0;
  for(int i=0 ; i<nnames ; i++)
    if(names[i] != "") 
      nFacies++;
  
  for(int i=0 ; i<nnames ; i++)
  {
    if(names[i] != "")
    {
      modelSettings->addFaciesName(names[i]);
      modelSettings->addFaciesLabel(globalmin + i);
    }
  }
  
  // Compare names in wells with names given in .xml-file
  if(modelSettings->getIsPriorFaciesProbGiven()==1)
  {
    typedef std::map<std::string,float> mapType;
    mapType myMap = modelSettings->getPriorFaciesProb();

    for(int i=0;i<nFacies;i++)
    {
      mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
      if (iter==myMap.end())
      {
        tmpErrText += "Problem with facies logs. Facies "+modelSettings->getFaciesName(i)+" is not one of the facies given in the xml-file.\n";
        error++;
      }
    }
  }
  else if(modelSettings->getIsPriorFaciesProbGiven()==2)
  {
    typedef std::map<std::string,std::string> mapType;
    mapType myMap = inputFiles->getPriorFaciesProbFile();

    for(int i=0;i<nFacies;i++)
    {
      mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
      if (iter==myMap.end())
      {
        tmpErrText += "Problem with facies logs. Facies "+modelSettings->getFaciesName(i)+" is not one of the facies given in the xml-file.\n";
        error++;
      }
    }
  }
}

void 
Model::processBackground(Background   *& background,
                         WellData     ** wells,
                         Simbox        * timeSimbox,
                         Simbox        * timeBGSimbox,
                         ModelSettings * modelSettings, 
                         InputFiles    * inputFiles,
                         std::string   & errText,
                         bool          & failed)
{
  if (modelSettings->getForwardModeling())
    Utils::writeHeader("Earth Model");
  else
    Utils::writeHeader("Prior Expectations / Background Model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
    
  FFTGrid * backModel[3];
  const int nx    = timeSimbox->getnx();
  const int ny    = timeSimbox->getny();
  const int nz    = timeSimbox->getnz();
  const int nxPad = modelSettings->getNXpad();
  const int nyPad = modelSettings->getNYpad();
  const int nzPad = modelSettings->getNZpad();
  if (modelSettings->getGenerateBackground()) 
  {
    FFTGrid * velocity = NULL;
    std::string backVelFile = inputFiles->getBackVelFile();
    if (backVelFile != "")
      loadVelocity(velocity, timeSimbox, modelSettings, backVelFile, errText, failed);
    if (!failed) 
    {
      if(modelSettings->getBackgroundVario() == NULL)
      {
        errText += "There is no variogram available for the background modelling.\n";
        failed = true;
      }
      for (int i=0 ; i<3 ; i++)
      {
        //backModel[i] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);    
        backModel[i] = createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());
        backModel[i]->setType(FFTGrid::PARAMETER);
      }
      background = new Background(backModel, wells, velocity, timeSimbox, timeBGSimbox, modelSettings);
    }
    if(velocity != NULL)
      delete velocity;
  }
  else 
  {
    std::vector<std::string> parName;
    parName.push_back("Vp "+modelSettings->getBackgroundType());
    parName.push_back("Vs "+modelSettings->getBackgroundType());
    parName.push_back("Rho "+modelSettings->getBackgroundType());

    for(int i=0 ; i<3 ; i++)
    {
      float constBackValue = modelSettings->getConstBackValue(i);

      const std::string & backFile = inputFiles->getBackFile(i);

      if(constBackValue < 0)
      {
        if(backFile.size() > 0)
        {
          const SegyGeometry      * dummy1 = NULL;
          const TraceHeaderFormat * dummy2 = NULL;
          const float               offset = modelSettings->getSegyOffset();
          std::string errorText("");
          int outsideTraces = readGridFromFile(backFile,
                                               parName[i],
                                               offset,
                                               backModel[i],
                                               dummy1,
                                               dummy2,
                                               FFTGrid::PARAMETER,
                                               timeSimbox,
                                               modelSettings,
                                               errorText);
          if(errorText != "")
          {
            errorText += "Reading of file '"+backFile+"' for parameter '"+parName[i]+"' failed\n";
            errText += errorText;
            failed = true;
          } 
          else {
            backModel[i]->logTransf();
            if(outsideTraces > 0) {
                errText += "Error: Background model in file "+backFile+" does not cover the inversion area.\n";
                failed = true;
            }
          }
        }
        else
        {
          errText += "Reading of file for parameter "+parName[i]+" failed. No file name is given.\n";
          failed = true;
        }
      }
      else if(constBackValue > 0)
      {
        backModel[i] = createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());            
        backModel[i]->setType(FFTGrid::PARAMETER);
        backModel[i]->fillInConstant(float( log( constBackValue )));
      }
      else
      {
        errText += "Trying to set background model to 0 for parameter "+parName[i]+"\n";
        failed = true;
      }
    }
    if (failed == false) {
      background = new Background(backModel);
    }
  }

  if (failed == false) {
    if(((modelSettings->getGridOutputFlag() & IO::BACKGROUND) > 0) || 
       ( modelSettings->getEstimationMode() && modelSettings->getGenerateBackground())) {
      background->writeBackgrounds(timeSimbox, 
                                   timeDepthMapping_, 
                                   timeCutMapping_, 
                                   *modelSettings->getTraceHeaderFormatOutput()); 
    }
  }
    
  Timings::setTimePriorExpectation(wall,cpu);
}

int
Model::readGridFromFile(const std::string       & fileName,
                        const std::string       & parName,
                        const float               offset,
                        FFTGrid                *& grid,
                        const SegyGeometry     *& geometry,   
                        const TraceHeaderFormat * format,
                        int                       gridType,
                        Simbox                  * timeSimbox,
                        ModelSettings           * modelSettings,
                        std::string             & errText,
                        bool                      nopadding)
{
  int fileType = IO::findGridType(fileName);

  int outsideTraces = 0;
  if(fileType == IO::CRAVA) 
  {
    int xpad, ypad, zpad;
    if(nopadding==false)
    {
      xpad = modelSettings->getNXpad();
      ypad = modelSettings->getNYpad();
      zpad = modelSettings->getNZpad();
    }
    else
    {
      xpad = timeSimbox->getnx();
      ypad = timeSimbox->getny();
      zpad = timeSimbox->getnz();
    }
    LogKit::LogFormatted(LogKit::LOW,"\nReading grid \'"+parName+"\' from file "+fileName);
    grid = createFFTGrid(timeSimbox->getnx(),
                             timeSimbox->getny(), 
                             timeSimbox->getnz(),
                             xpad, 
                             ypad, 
                             zpad,
                             modelSettings->getFileGrid());
    
    grid->setType(gridType);
    grid->readCravaFile(fileName, errText, nopadding);
  }
  else if(fileType == IO::SEGY) 
    outsideTraces  = readSegyFile(fileName, grid, timeSimbox, modelSettings, geometry, 
                                  gridType, offset, format, errText, nopadding);
  else if(fileType == IO::STORM) 
    outsideTraces = readStormFile(fileName, grid, gridType, parName, timeSimbox, modelSettings, errText, false, nopadding);
  else if(fileType == IO::SGRI)
    outsideTraces = readStormFile(fileName, grid, gridType, parName, timeSimbox, modelSettings, errText, true, nopadding);
  else 
  {
    errText += "\nReading of file \'"+fileName+"\' for grid type \'"
               +parName+"\'failed. File type not recognized.\n";
  }
  return(outsideTraces);
}


void 
Model::processPriorCorrelations(Corr         *& correlations,
                                Background    * background,
                                WellData     ** wells,
                                Simbox        * timeSimbox,
                                ModelSettings * modelSettings, 
                                InputFiles    * inputFiles,
                                std::string   & errText,
                                bool          & failed)
{
  bool printResult = ((modelSettings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0 ||
                      modelSettings->getEstimationMode() == true);
  if (modelSettings->getDoInversion() || printResult)
  {
    Utils::writeHeader("Prior Covariance");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    const std::string & paramCorrFile = inputFiles->getParamCorrFile();
    const std::string & corrTFile     = inputFiles->getTempCorrFile();

    bool estimateParamCorr = paramCorrFile == "";
    bool estimateTempCorr  = corrTFile     == "";

    //
    // Read parameter correlation (Var0) from file 
    //
    float ** paramCorr = NULL;
    bool failedParamCorr = false;
    if(!estimateParamCorr) 
    {
      std::string tmpErrText("");
      paramCorr = readMatrix(paramCorrFile, 3, 3, "parameter correlation", tmpErrText);
      if(paramCorr == NULL) 
      {
        errText += "Reading of file "+paramCorrFile+" for parameter correlation matrix failed\n";
        errText += tmpErrText;
        failedParamCorr = true;
      }
    }

    //
    // Estimate lateral correlation from seismic data
    //
    Surface * CorrXY = findCorrXYGrid(modelSettings);

    if(modelSettings->getLateralCorr()==NULL) // NBNB-PAL: this will never be true (default lateral corr)
    {
      estimateCorrXYFromSeismic(CorrXY,
                                seisCube_,
                                modelSettings->getNumberOfAngles());
    }

    int nCorrT = modelSettings->getNZpad();
    if((nCorrT % 2) == 0)
      nCorrT = nCorrT/2+1;
    else
      nCorrT = nCorrT/2+1;
    float * corrT = NULL;

    bool failedTempCorr = false;
    if(!estimateTempCorr) 
    {
      std::string tmpErrText("");
      float ** corrMat = readMatrix(corrTFile, 1, nCorrT+1, "temporal correlation", tmpErrText);
      if(corrMat == NULL) 
      {
        errText += "Reading of file '"+corrTFile+"' for temporal correlation failed\n";
        errText += tmpErrText;
        failedTempCorr = true;
      }
      corrT = new float[nCorrT];
      if (!failedTempCorr) 
      {
        for(int i=0;i<nCorrT;i++)
          corrT[i] = corrMat[0][i+1];
        delete [] corrMat[0];
        delete [] corrMat;
      }
    }

    float ** pointVar0 = NULL;
    if (estimateParamCorr || estimateTempCorr) //Need well estimation
    {
      Analyzelog * analyze = new Analyzelog(wells, 
                                            background,
                                            timeSimbox, 
                                            modelSettings);

      if(estimateParamCorr)
        paramCorr = analyze->getVar0();
      else
        delete [] analyze->getVar0();

      pointVar0 = analyze->getPointVar0();

      float * estCorrT = analyze->getCorrT();
      if(estimateTempCorr) {
        corrT = new float[nCorrT];
        int nEst = analyze->getNumberOfLags();
        int i, max = nEst;
        if(max > nCorrT)
          max = nCorrT;
        for(i=0;i<max;i++)
          corrT[i] = estCorrT[i];
        if(i<nCorrT) {
          LogKit::LogFormatted(LogKit::WARNING, 
            "Warning: Only able to estimate %d of %d lags needed in temporal correlation. The rest are set to 0.\n", nEst, nCorrT);
          for(;i<nCorrT;i++)
            corrT[i] = 0.0f;
        }
      }
      delete [] estCorrT;

      delete analyze;
    }

    if (failedParamCorr || failedTempCorr)
      failed = true;
        
    if (!failed) {
      correlations = new Corr(pointVar0, 
                              paramCorr, 
                              corrT, 
                              nCorrT,
                              static_cast<float>(timeSimbox->getdz()), 
                              CorrXY);
      if(printResult)
        correlations->writeFilePriorVariances(modelSettings);
      correlations->printPriorVariances();
    }

    if(failedTempCorr == false && failedParamCorr == false && correlations == NULL)
    {
      errText += "Could not construct prior covariance. Unknown why...\n";
      failed = true;
    }

    Timings::setTimePriorCorrelation(wall,cpu);
  }
}

Surface * 
Model::findCorrXYGrid(ModelSettings * modelSettings)
{
  float dx  = static_cast<float>(timeSimbox_->getdx());
  float dy  = static_cast<float>(timeSimbox_->getdy());

  int   nx  = modelSettings->getNXpad();
  int   ny  = modelSettings->getNYpad();

  Surface * grid = new Surface(0, 0, dx*nx, dy*ny, nx, ny, RMISSING);

  if(modelSettings->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
  {
    int refi,refj;
    for(int j=0;j<ny;j++)
    {
      for(int i=0;i<nx;i++)
      {
        if(i<(nx/2+1))
        {
          refi = i;
        }
        else
        {
          refi = i-nx;
        }
        if(j< (ny/2+1))
        {
          refj = j;
        }
        else
        {
          refj = j-ny;
        }
        (*grid)(j*nx+i) = modelSettings->getLateralCorr()->corr(refi*dx, refj*dy);
      }
    }
  }
  return(grid);
}

void 
Model::estimateCorrXYFromSeismic(Surface *& corrXY,
                                 FFTGrid ** seisCube,
                                 int        nAngles)
{
  FFTGrid * transf;
  float   * grid;

  int n = corrXY->GetNI()*corrXY->GetNJ();
  grid = new float[n];

  for(int i=0 ; i<n ; i++)
    grid[i] = 0.0;

  for(int i=0 ; i<nAngles ; i++)
  {
    if(seisCube[i]->isFile())
      transf = new FFTFileGrid(static_cast<FFTFileGrid *>(seisCube[i])); //move new out of loop? Copy grid instead
    else
      transf = new FFTGrid(seisCube[i]); //move new out of loop? Copy grid instead

    transf->setAccessMode(FFTGrid::RANDOMACCESS);
    transf->fftInPlace();
    transf->square();
    transf->invFFTInPlace();
    transf->collapseAndAdd( grid ); //the result of the collapse (the result for z=0) is is added to grid
    transf->endAccess();
    delete transf;
  }
  float sill = grid[0];
  for(int i=0;i<n;i++)
    (*corrXY)(i) = grid[i]/sill;
  delete [] grid;
}

void 
Model::processReflectionMatrixFromWells(float       **& reflectionMatrix,
                                        WellData     ** wells,
                                        ModelSettings * modelSettings, 
                                        InputFiles    * inputFiles,
                                        std::string   & errText,
                                        bool          & failed)
{
  Utils::writeHeader("Making reflection matrix from well information");  
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from wells.
  //
  const std::string & reflMatrFile = inputFiles->getReflMatrFile();

  if(reflMatrFile != "") {
    std::string tmpErrText("");
    reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", tmpErrText);
    if(reflectionMatrix == NULL) {
      errText += "Reading of file "+reflMatrFile+ " for reflection matrix failed\n";
      errText += tmpErrText;
      failed = true;
    }
    LogKit::LogFormatted(LogKit::LOW,"Reflection parameters read from file.\n\n");
  }
  else {

    double vsvp = vsvpFromWells(wells, modelSettings);

    setupDefaultReflectionMatrix(reflectionMatrix, vsvp, modelSettings);
  }
}

void 
Model::processReflectionMatrixFromBackground(float       **& reflectionMatrix,
                                             Background    * background,
                                             ModelSettings * modelSettings, 
                                             InputFiles    * inputFiles,
                                             std::string   & errText,
                                             bool          & failed)
{
  if (modelSettings->getForwardModeling())
    Utils::writeHeader("Making reflection matrix from earth model");  
  else
    Utils::writeHeader("Making reflection matrix from background model");  
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model.
  //
  const std::string & reflMatrFile = inputFiles->getReflMatrFile();

  if(reflMatrFile != "") {
    std::string tmpErrText("");
    reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", tmpErrText);
    if(reflectionMatrix == NULL) {
      errText += "Reading of file "+reflMatrFile+" for reflection matrix failed\n";
      errText += tmpErrText;
      failed = true;
    }
    LogKit::LogFormatted(LogKit::LOW,"Reflection parameters read from file.\n\n");
  }
  else {

    if (background != NULL) {
      double vsvp  = background->getMeanVsVp(); 
      setupDefaultReflectionMatrix(reflectionMatrix, vsvp, modelSettings);
    }
    else {
      errText += "\nFailed to set up reflection matrix. Background model is empty.\n";
      failed = true;
    }
  }
}

void
Model::setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                    double          vsvp,
                                    ModelSettings * modelSettings)
{
  int      i;
  float ** A      = new float * [modelSettings->getNumberOfAngles()];

  // For debugging
  //background->setClassicVsVp();

  double vsvp2 = vsvp*vsvp;
  for(i = 0; i < modelSettings->getNumberOfAngles(); i++)
  {
    double angle = static_cast<double>(modelSettings->getAngle(i));
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(modelSettings->getSeismicType(i) == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float( (1.0 +tan2t )/2.0 ) ; 
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0 );
    }
    else if(modelSettings->getSeismicType(i) == ModelSettings::PSSEIS) {
      double cost  = cos(angle);
      double cosp  = sqrt(1-vsvp2*sint2);
      double fac   = 0.5*sint/cosp;

      A[i][0] = 0;
      A[i][1] = float(4.0*fac*(vsvp2*sint2-vsvp*cost*cosp));
      A[i][2] = float(fac*(-1.0+2*vsvp2*sint2+2*vsvp*cost*cosp));
    }
  }
  reflectionMatrix = A;
  LogKit::LogFormatted(LogKit::LOW,"\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",1.0f/vsvp);
}

double Model::vsvpFromWells(WellData     ** wells,
                            ModelSettings * modelSettings)
{
  int    i;
  int    nWells = modelSettings->getNumberOfWells();
  double vsvp;
  float  muA;
  float  muB;
  float  vp = 0;
  float  vs = 0;
  
  for( i=0; i<nWells; i++ )
  {
    wells[i]->getMeanVsVp(muA, muB);
    vp += muA;
    vs += muB;
  }

  vsvp = vs/vp;
  return vsvp;
}



void
Model::processWellLocation(FFTGrid                     ** seisCube, 
                           WellData                    ** wells, 
                           float                       ** reflectionMatrix,
                           Simbox                       * timeSimbox,
                           ModelSettings                * modelSettings,
                           const std::vector<Surface *> & interval, 
                           RandomGen                    * randomGen)
{
  Utils::writeHeader("Estimating optimized well location");
  
  double  deltaX, deltaY;
  float   sum;
  float   kMove;
  float   moveAngle;
  int     iMove;
  int     jMove;
  int     i,j,w;
  int     iMaxOffset;
  int     jMaxOffset;
  int     nMoveAngles = 0;
  int     nWells      = modelSettings->getNumberOfWells();
  int     nAngles     = modelSettings->getNumberOfAngles();
  float   maxShift    = modelSettings->getMaxWellShift();
  float   maxOffset   = modelSettings->getMaxWellOffset();
  double  angle       = timeSimbox->getAngle();
  double  dx          = timeSimbox->getdx();
  double  dy          = timeSimbox->getdx();

  std::vector<float> angleWeight(nAngles); 
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"  Well             Shift[ms]       DeltaI   DeltaX[m]   DeltaJ   DeltaY[m] \n");
  LogKit::LogFormatted(LogKit::LOW,"  ----------------------------------------------------------------------------------\n");

  for (w = 0 ; w < nWells ; w++) {
    if( wells[w]->isDeviated()==true )
      continue;

    BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
    nMoveAngles = modelSettings->getNumberOfWellAngles(w);
    
    if( nMoveAngles==0 )
      continue;

    for( i=0; i<nAngles; i++ )
      angleWeight[i] = 0;

    for( i=0; i<nMoveAngles; i++ ){
      moveAngle   = modelSettings->getWellMoveAngle(w,i);

      for( j=0; j<nAngles; j++ ){
        if( moveAngle==modelSettings->getAngle(j)){
          angleWeight[j] = modelSettings->getWellMoveWeight(w,i);
          break;
        }
      }
    }

    sum = 0;
    for( i=0; i<nAngles; i++ )
      sum += angleWeight[i];
    if( sum == 0 )
      continue;

    iMaxOffset = static_cast<int>(std::ceil(maxOffset/dx));
    jMaxOffset = static_cast<int>(std::ceil(maxOffset/dy));

    bl->findOptimalWellLocation(seisCube,timeSimbox,reflectionMatrix,nAngles,angleWeight,maxShift,iMaxOffset,jMaxOffset,interval,iMove,jMove,kMove);

    deltaX = iMove*dx*cos(angle) - jMove*dy*sin(angle);
    deltaY = iMove*dx*sin(angle) + jMove*dy*cos(angle);
    wells[w]->moveWell(timeSimbox,deltaX,deltaY,kMove);
    wells[w]->deleteBlockedLogsOrigThick();
    wells[w]->setBlockedLogsOrigThick( new BlockedLogs(wells[w], timeSimbox, randomGen) );
    LogKit::LogFormatted(LogKit::LOW,"  %-13s %11.2f %12d %11.2f %8d %11.2f \n", 
    wells[w]->getWellname().c_str(), kMove, iMove, deltaX, jMove, deltaY);
  }

   for (w = 0 ; w < nWells ; w++){
     nMoveAngles = modelSettings->getNumberOfWellAngles(w);

    if( wells[w]->isDeviated()==true && nMoveAngles > 0 )
    {
      LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: Well %7s is treated as deviated and can not be moved.\n",
          wells[w]->getWellname().c_str());
      TaskList::addTask("Well "+NRLib::ToString(wells[w]->getWellname())+" can not be moved. Remove <optimize-location-to> for this well");
    }
   }
}

void
Model::processWavelets(Wavelet                    **& wavelet,
                       FFTGrid                     ** seisCube,
                       WellData                    ** wells,
                       float                       ** reflectionMatrix,
                       Simbox                       * timeSimbox,
                       const std::vector<Surface *> & waveletEstimInterval,    
                       ModelSettings                * modelSettings, 
                       InputFiles                   * inputFiles,
                       std::string                  & errText,
                       bool                         & failed)
{
  int error = 0;
  Utils::writeHeader("Processing/generating wavelets");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  bool estimateStuff = false;
  for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++) {  
    estimateStuff = estimateStuff || modelSettings->getEstimateWavelet(i); 
    estimateStuff = estimateStuff || modelSettings->getEstimateSNRatio(i);
    if(modelSettings->getEstimateWavelet(i) == true)
      modelSettings->setWaveletScale(i,1.0);
  }
  if (estimateStuff) {
    LogKit::LogFormatted(LogKit::HIGH,"\nWells that cannot be used in wavelet generation or noise estimation:");
    LogKit::LogFormatted(LogKit::HIGH,"\n  Deviated wells.");
    LogKit::LogFormatted(LogKit::HIGH,"\n  Wells with too little data.\n");
  }

  unsigned int nAngles = modelSettings->getNumberOfAngles();

  wavelet = new Wavelet * [nAngles];

  std::vector<Grid2D *> shiftGrids(nAngles);
  std::vector<Grid2D *> gainGrids(nAngles);
  localNoiseScale_.resize(nAngles);

  bool has3Dwavelet = false;
  for(unsigned int i=0 ; i < nAngles ; i++)
  {
    localNoiseScale_[i] = NULL;
    shiftGrids[i]       = NULL;
    gainGrids[i]        = NULL;
    if (modelSettings->getWaveletDim(i) == Wavelet::THREE_D)
      has3Dwavelet = true;
  }

  NRLib::Grid2D<float>      refTimeGradX;          ///< Time gradient in x-direction for reference time surface (t0)
  NRLib::Grid2D<float>      refTimeGradY;          ///< Time gradient in x-direction for reference time surface (t0)
  unsigned int nWells = modelSettings->getNumberOfWells();
  std::vector<std::vector<double> > tGradX(nWells);
  std::vector<std::vector<double> > tGradY(nWells);

  if (has3Dwavelet) {
    if (inputFiles->getRefSurfaceFile() != "") {
      if (findTimeGradientSurface(inputFiles->getRefSurfaceFile(), timeSimbox, refTimeGradX, refTimeGradY) == false) {
        errText += "Simbox is not completely inside reference time surface in (x,y).\n";
        error = 1;
      }
    }
    for (unsigned int w=0; w<nWells; w++) {
      BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick(); 
      bl->findSeismicGradient(seisCube[0], timeSimbox, /*nAngles,*/tGradX[w], tGradY[w]);
    }
  }
  
  if (timeSimbox->getdz() > 4.01f && modelSettings->getEstimateNumberOfWavelets() > 0) 
  { // Require this density for wavelet estimation
    LogKit::LogFormatted(LogKit::LOW,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
    LogKit::LogFormatted(LogKit::LOW,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
    LogKit::LogFormatted(LogKit::LOW,"         The number of layers must be increased.                                  \n");
    std::string text("");
    text += "Increase the number of layers to improve the quality of the wavelet estimation.\n";
    text += "   The minimum sampling density is "+NRLib::ToString(timeSimbox->getdz())+", and it should be lower than 4.0.\n";
    text += "   To obtain the desired density, the number of layers should be at least "+NRLib::ToString(timeSimbox->GetLZ()/4);
    TaskList::addTask(text);
  }

  float globalScale = 1.0;
  for(unsigned int i=0 ; i < nAngles ; i++) { 
    if(modelSettings_->getForwardModeling()==false)
      seisCube[i]->setAccessMode(FFTGrid::RANDOMACCESS);
    if(modelSettings->getUseLocalWavelet()==true) {
      if(inputFiles->getScaleFile(i)!="") {
        Surface help(inputFiles->getScaleFile(i));
        gainGrids[i] = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
        resampleSurfaceToGrid2D(timeSimbox, &help, gainGrids[i]);
      }
    }
    
    globalScale = modelSettings->getWaveletScale(i);
    
    float angle = float(modelSettings->getAngle(i)*180.0/M_PI);
    LogKit::LogFormatted(LogKit::LOW,"\nAngle stack : %.1f deg",angle);
    if (modelSettings->getEstimateWavelet(i)) 
    {
      if (modelSettings->getWaveletScale(i) != 1.0) {
        LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: The wavelet scale specified in the model file ("
                             +NRLib::ToString(modelSettings->getWaveletScale(i),2)
                             +") has no effect when the wavelet\n         is estimated and not read from file\n\n");
      }
      if (modelSettings->getWaveletDim(i) == Wavelet::ONE_D)
        wavelet[i] = new Wavelet1D(timeSimbox, 
                                   seisCube[i], 
                                   wells, 
                                   waveletEstimInterval,                                   
                                   modelSettings, 
                                   reflectionMatrix[i],
                                   i);
      else {
        wavelet[i] = new Wavelet3D(inputFiles->getWaveletFilterFile(i),
                                   waveletEstimInterval,
                                   refTimeGradX,
                                   refTimeGradY,
                                   tGradX,
                                   tGradY,
                                   seisCube[i],
                                   modelSettings,
                                   wells,
                                   timeSimbox,
                                   reflectionMatrix[i],
                                   i,
                                   error,
                                   errText);
        modelSettings->setEstimateWaveletNoise(false);
      }
    }
    else { //Not estimation modus
      const std::string & waveletFile = inputFiles->getWaveletFile(i);
      int fileFormat = getWaveletFileFormat(waveletFile,errText);
      if(fileFormat < 0) {
        errText += "Unknown file format of file '"+waveletFile+"'.\n";
        error++;
      }
      else {
        bool flip = false;
        if (fileFormat == Wavelet::OLD)
          flip = true;
        if (modelSettings->getWaveletDim(i) == Wavelet::ONE_D) {
          wavelet[i] = new Wavelet1D(waveletFile, 
                                     fileFormat, 
                                     modelSettings, 
                                     reflectionMatrix[i],
                                     modelSettings->getAngle(i),
                                     error, 
                                     errText);
          // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
          if(seisCube!=NULL) // If forward modeling, we have no seismic, can not prescale wavelet.
          {
            float prescale = wavelet[i]->findGlobalScaleForGivenWavelet(modelSettings, timeSimbox, seisCube[i], wells);
            double limHigh = 3.0;
            double limLow  = 0.33; 
            float  scale   = 1.0;
            if(!modelSettings->getEstimateGlobalWaveletScale(i) && modelSettings->getWaveletScale(i)!=scale && modelSettings->getEstimateLocalScale(i) && (prescale>limHigh || prescale<limLow))
            {
              errText += "The wavelet given for angle no "+NRLib::ToString(i)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
              error++;
            }
            else if(!modelSettings->getEstimateGlobalWaveletScale(i) && modelSettings->getWaveletScale(i)!=scale && (prescale>limHigh || prescale<limLow))
            {
              LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: The wavelet given for angle no "+ NRLib::ToString(i) +"is badly scaled. Ask Crava to estimate global wavelet scale.\n");
              TaskList::addTask("Ask Crava to estimate global wavelet scale");
            }
            else if(modelSettings->getEstimateGlobalWaveletScale(i)) // prescale, then we have correct size order, and later scale estimation will be ok.
              wavelet[i]->multiplyRAmpByConstant(prescale);
          }
          if (error == 0)
            wavelet[i]->resample(static_cast<float>(timeSimbox->getdz()), 
                                 timeSimbox->getnz(), 
                                 static_cast<float>(modelSettings->getZPadFac()),
                                 flip);
        }
        else { //3D-wavelet constructed by 1D wavelet and filter
          wavelet[i] = new Wavelet3D(waveletFile, 
                                     fileFormat, 
                                     modelSettings, 
                                     reflectionMatrix[i],
                                     modelSettings->getAngle(i),
                                     error, 
                                     errText,
                                     inputFiles->getWaveletFilterFile(i));
          if (error == 0) {
            wavelet[i]->resample(static_cast<float>(timeSimbox->getdz()), 
                                 timeSimbox->getnz(), 
                                 static_cast<float>(modelSettings->getZPadFac()),
                                 flip);
          }
        }
      }
    }

    if ((modelSettings->getEstimationMode() == false) && (wavelet[i]->getDim() == 3) && !timeSimbox->getIsConstantThick()) {
      errText += "Simbox must have constant thickness if forward modelling or inversion when 3D wavelet.\n";
      error++;
    }

    if (error == 0) {
      wavelet[i]->scale(modelSettings->getWaveletScale(i));
      bool localEst = (modelSettings->getEstimateLocalScale(i) || modelSettings->getEstimateLocalShift(i) ||
                       modelSettings->getEstimateLocalNoise(i) || modelSettings->getEstimateGlobalWaveletScale(i));

      if ((modelSettings->getEstimateSNRatio(i) || localEst) && modelSettings->getForwardModeling() == false)
      {
        if (wavelet[i]->getDim() == 3) { //Not possible to estimate signal-to-noise ratio for 3D wavelets
          errText += "Estimation of signal-to-noise ratio is not possible for 3D wavelets.\n";
          errText += "The s/n ratio must be specified in the model file\n";
          error++;
        }
        else {
          float SNRatio = wavelet[i]->calculateSNRatioAndLocalWavelet(timeSimbox, seisCube[i], wells, 
                                                                      shiftGrids[i], gainGrids[i],
                                                                      modelSettings, errText, error, 
                                                                      localNoiseScale_[i], i, globalScale);
          if(modelSettings->getEstimateSNRatio(i))
            modelSettings->setSNRatio(i,SNRatio);
        }
      }

      if (error == 0) {
        if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 || 
          (modelSettings->getEstimationMode() && modelSettings->getEstimateWavelet(i)))
        {
          std::string type;
          if (modelSettings->getEstimateWavelet(i))
            type = "Estimated_";
          else if (modelSettings->getWaveletScale(i) == 1.00)
            type = "";
          else
            type = "Scaled_";
          wavelet[i]->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0); // dt_max = 1.0;
        }
        
        float SNRatio = modelSettings->getSNRatio(i);
        float SNLow  = 1.0;
        float SNHigh = 10.0;
        if ((SNRatio <=SNLow  || SNRatio > SNHigh) && modelSettings->getForwardModeling()==false) {
          errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio)+" for cube "+NRLib::ToString(i+1)+".\n";
          errText += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
          error++;
        }
        
        bool useLocalNoise = modelSettings->getEstimateLocalNoise(i) || inputFiles->getLocalNoiseFile(i) != "";
        bool useLocalShift = modelSettings->getEstimateLocalShift(i) || inputFiles->getShiftFile(i) != "";
        bool useLocalGain  = modelSettings->getEstimateLocalScale(i) || inputFiles->getScaleFile(i) != "";

        if (useLocalNoise) {
          writeLocalGridsToFile(inputFiles->getLocalNoiseFile(i),                              
                                IO::PrefixLocalNoise(),
                                1.0,  // Scale map with this factor before writing to disk
                                modelSettings,
                                i,
                                timeSimbox,
                                localNoiseScale_[i]);
        }

        if (useLocalShift) {
          writeLocalGridsToFile(inputFiles->getShiftFile(i),
                                IO::PrefixLocalWaveletShift(),
                                1.0,
                                modelSettings,
                                i,
                                timeSimbox,
                                shiftGrids[i]);
          wavelet[i]->setShiftGrid(shiftGrids[i]);
        }

        if (useLocalGain) {
          writeLocalGridsToFile("", // Gain grids have already been read.
                                IO::PrefixLocalWaveletGain(),
                                1.0,
                                modelSettings,
                                i,
                                timeSimbox,
                                gainGrids[i]);
          wavelet[i]->setGainGrid(gainGrids[i]);
        }
      }
    }
    if(modelSettings_->getForwardModeling()==false) // else, no seismic data
      seisCube[i]->endAccess();
  } // end i (angles)
  
  // check if local noise is set for some angles.
  bool localNoiseSet = false;
  for(unsigned int i=0;i<nAngles;i++)
    if(localNoiseScale_[i]!=NULL)
      localNoiseSet = true;
  if(localNoiseSet==true)
  {
    for(unsigned int i=0;i<nAngles;i++)
      if(localNoiseScale_[i]==NULL)
        localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(), 
                                         timeSimbox->getny(), 
                                         1.0);
  }
 
  Timings::setTimeWavelets(wall,cpu);
  failed = error > 0;
  if(estimateStuff == true && modelSettings->getEstimationMode() == true) {
    WellData ** wells = getWells();
    for (int i=0 ; i<modelSettings->getNumberOfWells() ; i++)
      wells[i]->getBlockedLogsOrigThick()->writeWell(modelSettings);
  }
}

int
Model::getWaveletFileFormat(const std::string & fileName, std::string & errText)
{
  int fileformat = -1;
  int line = 0;
  std::string dummyStr;
  // test for old file format
  std::ifstream file;
  NRLib::OpenRead(file,fileName);
  for(int i = 0; i < 5; i++) {
    NRLib::GetNextToken(file,dummyStr,line);
    if (NRLib::CheckEndOfFile(file)) {
      errText += "End of wavelet file '"+fileName+"' is premature.\n";
      return 0;
    } 
  }  
  file.close();
  file.clear();
  std::string targetString = "CMX";
  int  pos = Utils::findEnd(dummyStr, 0, targetString);
  if(pos>=0)
    fileformat= Wavelet::OLD;

  targetString = "pulse";
  pos = Utils::findEnd(dummyStr, 0, targetString);
  if(pos>=0)
    fileformat= Wavelet::NORSAR;

  if(fileformat<0) { // not old or norsar format
      // test for jason file format
    NRLib::OpenRead(file,fileName);
    line         = 0;
    int thisLine = 0;
    bool lineIsComment = true; 
    while (lineIsComment == true) {
      NRLib::GetNextToken(file,dummyStr,line);
      if (NRLib::CheckEndOfFile(file)) {
        errText += "End of wavelet file "+fileName+" is premature\n";
        return 0;
      } 
      else {
        if (thisLine == line) {
          NRLib::DiscardRestOfLine(file,line,false);
          thisLine = line;
        }
        if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
          lineIsComment = false;
      }
    }  
    file.close();
    if (NRLib::ParseType<double>(dummyStr)!=0) // not convertable number
      fileformat= Wavelet::JASON;
  }
  return fileformat;
}

void Model::processPriorFaciesProb(const std::vector<Surface *> & faciesEstimInterval,
                                   float                       *& priorFacies,
                                   WellData                    ** wells,
                                   RandomGen                    * randomGen,
                                   int                            nz,
                                   float                          dz,
                                   ModelSettings                * modelSettings,
                                   bool                         & failed,
                                   std::string                  & errTxt,
                                   InputFiles                   * inputFiles)
{
  if (modelSettings->getEstimateFaciesProb())
  {
    Utils::writeHeader("Prior Facies Probabilities");
    int nFacies = modelSettings->getNumberOfFacies();

    if(modelSettings->getIsPriorFaciesProbGiven()==0)
    {
      if (nFacies > 0) 
      {
        int nWells  = modelSettings->getNumberOfWells();
        int nFacies = modelSettings->getNumberOfFacies();
        int ndata   = nWells*nz;

        int ** faciesCount = new int * [nWells]; 
        for (int w = 0 ; w < nWells ; w++)
          faciesCount[w] = new int[nFacies];

        for (int w = 0 ; w < nWells ; w++)
          for (int i = 0 ; i < nFacies ; i++)
            faciesCount[w][i] = 0;

        int * faciesLog = new int[ndata];   // NB! *internal* log numbering (0, 1, 2, ...)
        for (int i = 0 ; i < ndata ; i++)
          faciesLog[i] = IMISSING;

        float * vtAlpha   = new float[nz];  // vt = vertical trend
        float * vtBeta    = new float[nz];
        float * vtRho     = new float[nz];
        int   * vtFacies  = new int[nz];

        int nUsedWells = 0;

        for (int w = 0 ; w < nWells ; w++)
        {
          if(wells[w]->getNFacies() > 0) // Well has facies log
          { 
            //
            // Note that we use timeSimbox to calculate prior facies probabilities
            // instead of the simbox with parallel top and base surfaces. This
            // will make the prior probabilities slightly different, but that
            // should not be a problem.
            //
            BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
            int nBlocks = bl->getNumberOfBlocks();
            //
            // Set facies data outside facies estimation interval IMISSING
            //
            int * blFaciesLog = new int[nBlocks];
            Utils::copyVector(bl->getFacies(), blFaciesLog, nBlocks);

            if (faciesEstimInterval.size() > 0) {
              const double * xPos  = bl->getXpos();
              const double * yPos  = bl->getYpos();
              const double * zPos  = bl->getZpos();
              for (int i = 0 ; i < nBlocks ; i++) {
                const double zTop  = faciesEstimInterval[0]->GetZ(xPos[i],yPos[i]);
                const double zBase = faciesEstimInterval[1]->GetZ(xPos[i],yPos[i]);
                if ( (zPos[i] - 0.5*dz) < zTop || (zPos[i] + 0.5*dz) > zBase)
                  blFaciesLog[i] = IMISSING;
              }
            }

            bl->getVerticalTrend(bl->getAlpha(),vtAlpha);
            bl->getVerticalTrend(bl->getBeta(),vtBeta);
            bl->getVerticalTrend(bl->getRho(),vtRho);
            bl->getVerticalTrend(blFaciesLog,vtFacies,randomGen);
            delete [] blFaciesLog;

            for(int i=0 ; i<nz ; i++)
            {
              int facies;
              if(vtAlpha[i] != RMISSING && vtBeta[i] != RMISSING && vtRho[i] != RMISSING)
                facies = vtFacies[i];
              else
                facies = IMISSING;

              faciesLog[w*nz + i] = facies;

              if(facies != IMISSING)
                faciesCount[w][facies]++;
            }
            nUsedWells++;
          }
        }
        delete [] vtAlpha;
        delete [] vtBeta; 
        delete [] vtRho;  
        delete [] vtFacies; 

        if (nUsedWells > 0) {
          //
          // Probabilities
          //
          LogKit::LogFormatted(LogKit::LOW,"\nFacies distributions for each blocked well: \n");
          LogKit::LogFormatted(LogKit::LOW,"\nBlockedWell              ");
          for (int i = 0 ; i < nFacies ; i++)
            LogKit::LogFormatted(LogKit::LOW,"%12s ",modelSettings->getFaciesName(i).c_str());
          LogKit::LogFormatted(LogKit::LOW,"\n");
          for (int i = 0 ; i < 24+13*nFacies ; i++)
            LogKit::LogFormatted(LogKit::LOW,"-");
          LogKit::LogFormatted(LogKit::LOW,"\n");
          for (int w = 0 ; w < nWells ; w++)
          {
            if(wells[w]->getNFacies() > 0) // Well has facies log
            { 
              float tot = 0.0;
              for (int i = 0 ; i < nFacies ; i++)
                tot += static_cast<float>(faciesCount[w][i]);
              LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[w]->getWellname().c_str());
              for (int i = 0 ; i < nFacies ; i++) {
                float faciesProb = static_cast<float>(faciesCount[w][i])/tot;
                LogKit::LogFormatted(LogKit::LOW," %12.4f",faciesProb);
              }
              LogKit::LogFormatted(LogKit::LOW,"\n");
            }
          }
          LogKit::LogFormatted(LogKit::LOW,"\n");
          //
          // Counts
          //
          LogKit::LogFormatted(LogKit::MEDIUM,"\nFacies counts for each well: \n");

          LogKit::LogFormatted(LogKit::MEDIUM,"\nBlockedWell              ");
          for (int i = 0 ; i < nFacies ; i++)
            LogKit::LogFormatted(LogKit::MEDIUM,"%12s ",modelSettings->getFaciesName(i).c_str());
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");
          for (int i = 0 ; i < 24+13*nFacies ; i++)
            LogKit::LogFormatted(LogKit::MEDIUM,"-");
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");
          for (int w = 0 ; w < nWells ; w++)
          {
            if(wells[w]->getUseForFaciesProbabilities())
            { 
              float tot = 0.0;
              for (int i = 0 ; i < nFacies ; i++)
                tot += static_cast<float>(faciesCount[w][i]);
              LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[w]->getWellname().c_str());
              for (int i = 0 ; i < nFacies ; i++) {
                LogKit::LogFormatted(LogKit::MEDIUM,"%12d ",faciesCount[w][i]);
              }
              LogKit::LogFormatted(LogKit::MEDIUM,"\n");
            }
          }
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");

          for (int w = 0 ; w < nWells ; w++)
            delete [] faciesCount[w];
          delete [] faciesCount;

          //
          // Make prior facies probabilities
          //
          float sum = 0.0f;
          int * nData = new int[nFacies];
          for(int i=0 ; i<nFacies ; i++)
            nData[i] = 0;

          for(int i=0 ; i<ndata ; i++) {
            if(faciesLog[i] != IMISSING)
              nData[faciesLog[i]]++;
          }
          delete [] faciesLog;

          for(int i=0 ; i<nFacies ; i++)
            sum += nData[i];

          if (sum > 0) {
            LogKit::LogFormatted(LogKit::LOW,"Facies         Probability\n");
            LogKit::LogFormatted(LogKit::LOW,"--------------------------\n");
            priorFacies = new float[nFacies];
            for(int i=0 ; i<nFacies ; i++) {
              priorFacies[i] = float(nData[i])/sum;
              LogKit::LogFormatted(LogKit::LOW,"%-15s %10.4f\n",modelSettings->getFaciesName(i).c_str(),priorFacies[i]);
            }
          }
          else { 
            LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: No valid facies log entries have been found\n");
            modelSettings->setEstimateFaciesProb(false);
            TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");

          }
          delete [] nData;
        }
        else
        {
          LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: Estimation of facies probabilites have been requested, but there");
          LogKit::LogFormatted(LogKit::WARNING,"\n         are no wells with facies available and CRAVA will therefore not");
          LogKit::LogFormatted(LogKit::WARNING,"\n         be able to estimate these probabilities...\n");
          modelSettings->setEstimateFaciesProb(false);

          TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
        }
      }
      else 
      {
        LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: Estimation of facies probabilites have been requested, but no facies");
        LogKit::LogFormatted(LogKit::WARNING,"\n         have been found and CRAVA will therefore not be able to estimate");
        LogKit::LogFormatted(LogKit::WARNING,"\n         these probabilities...\n");
        modelSettings->setEstimateFaciesProb(false);
        TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
      }
    }
    else if(modelSettings->getIsPriorFaciesProbGiven()==1)
    {
      priorFacies = new float[nFacies];
      typedef std::map<std::string,float> mapType;
      mapType myMap = modelSettings->getPriorFaciesProb();
      
      for(int i=0;i<nFacies;i++)
      {
        mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
        if(iter!=myMap.end())
          priorFacies[i] = iter->second;
        else
        {
          LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: No prior facies probability found for facies %12s\n",modelSettings->getFaciesName(i).c_str());
          modelSettings->setEstimateFaciesProb(false);
          TaskList::addTask("Check that facies " +NRLib::ToString(modelSettings->getFaciesName(i).c_str())+" is given a prior probability in the xml-file");
        }
      }
      LogKit::LogFormatted(LogKit::LOW,"Facies         Probability\n");
      LogKit::LogFormatted(LogKit::LOW,"--------------------------\n");
      for(int i=0 ; i<nFacies ; i++) {
        LogKit::LogFormatted(LogKit::LOW,"%-15s %10.4f\n",modelSettings->getFaciesName(i).c_str(),priorFacies[i]);
      }

    }
    else if(modelSettings->getIsPriorFaciesProbGiven()==2)
    {
      readPriorFaciesProbCubes(inputFiles, 
                                modelSettings, 
                                priorFaciesProbCubes_,
                                timeSimbox_,
                                errTxt,
                                failed);
       
       typedef std::map<std::string,std::string> mapType;
       mapType myMap = inputFiles->getPriorFaciesProbFile();

       LogKit::LogFormatted(LogKit::LOW,"Facies         Probability in file\n");
       LogKit::LogFormatted(LogKit::LOW,"----------------------------------\n");
       for(mapType::iterator it=myMap.begin();it!=myMap.end();it++)
         LogKit::LogFormatted(LogKit::LOW,"%-15s %10s\n",(it->first).c_str(),(it->second).c_str());
       
    }
  }
}

void Model::readPriorFaciesProbCubes(InputFiles      * inputFiles, 
                                     ModelSettings   * modelSettings, 
                                     FFTGrid       **& priorFaciesProbCubes,
                                     Simbox          * timeSimbox,
                                     std::string     & errTxt,
                                     bool            & failed)
{
  int nFacies = modelSettings->getNumberOfFacies();
  priorFaciesProbCubes_ = new FFTGrid*[nFacies];

  typedef std::map<std::string,std::string> mapType;
  mapType myMap = inputFiles->getPriorFaciesProbFile();
  for(int i=0;i<nFacies;i++)
  {
    mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
    
    if(iter!=myMap.end())
    {
      const std::string & faciesProbFile = iter->second;
      const SegyGeometry      * dummy1 = NULL;
      const TraceHeaderFormat * dummy2 = NULL;
      const float               offset = modelSettings->getSegyOffset();
      std::string errorText("");
      readGridFromFile(faciesProbFile,
                       "priorfaciesprob",
                       offset,
                       priorFaciesProbCubes[i],
                       dummy1,
                       dummy2,
                       FFTGrid::PARAMETER,
                       timeSimbox,
                       modelSettings,
                       errorText, true);
      if(errorText != "")
      {
        errorText += "Reading of file \'"+faciesProbFile+"\' for prior facies probability for facies \'"
                     +modelSettings->getFaciesName(i)+"\' failed\n";
        errTxt += errorText;
        failed = true;
      } 
    }
    else
    {
      LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: No prior facies probability found for facies %12s\n",
                           modelSettings->getFaciesName(i).c_str());
      TaskList::addTask("Check that facies "+NRLib::ToString(modelSettings->getFaciesName(i).c_str())+" is given prior probability in the xml-file");
      modelSettings->setEstimateFaciesProb(false);
      break;
    }  
  }
}


void
Model::loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
                         std::vector<Surface *> & faciesEstimInterval,
                         std::vector<Surface *> & wellMoveInterval,
                         Simbox      * timeSimbox,
                         InputFiles  * inputFiles,
                         std::string & errText,
                         bool        & failed)
{
  const double x0 = timeSimbox->getx0();
  const double y0 = timeSimbox->gety0();
  const double lx = timeSimbox->getlx();
  const double ly = timeSimbox->getly();
  const int    nx = timeSimbox->getnx();
  const int    ny = timeSimbox->getny();
  //
  // Get wavelet estimation interval
  //
  const std::string & topWEI  = inputFiles->getWaveletEstIntFile(0);
  const std::string & baseWEI = inputFiles->getWaveletEstIntFile(1);

  if (topWEI != "" && baseWEI != "") {  
    waveletEstimInterval.resize(2);
    try {
      if (NRLib::IsNumber(topWEI)) 
        waveletEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWEI.c_str()));
      else { 
        Surface tmpSurf(topWEI);
        waveletEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
    
    try {
      if (NRLib::IsNumber(baseWEI)) 
        waveletEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWEI.c_str()));
      else { 
        Surface tmpSurf(baseWEI);
        waveletEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
  //
  // Get facies estimation interval
  //
  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" && baseFEI != "") {  
    faciesEstimInterval.resize(2);
    try {
      if (NRLib::IsNumber(topFEI)) 
        faciesEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topFEI.c_str()));
      else { 
        Surface tmpSurf(topFEI);
        faciesEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {      
      errText += e.what();
      failed = true;
    }

    try {
      if (NRLib::IsNumber(baseFEI)) 
        faciesEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseFEI.c_str()));
      else { 
        Surface tmpSurf(baseFEI);
        faciesEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
  //
  // Get well move interval
  //
  const std::string & topWMI  = inputFiles->getWellMoveIntFile(0);
  const std::string & baseWMI = inputFiles->getWellMoveIntFile(1);

  if (topWMI != "" && baseWMI != "") {  
    wellMoveInterval.resize(2);
    try {
      if (NRLib::IsNumber(topWMI)) 
        wellMoveInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWMI.c_str()));
      else { 
        Surface tmpSurf(topWMI);
        wellMoveInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }

    try {
      if (NRLib::IsNumber(baseWMI)) 
        wellMoveInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWMI.c_str()));
      else { 
        Surface tmpSurf(baseWMI);
        wellMoveInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
}

void
Model::printSettings(ModelSettings * modelSettings,
                     InputFiles    * inputFiles)
{
  Utils::writeHeader("Model settings");

  LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings:\n");
  if(modelSettings->getForwardModeling()==true)
    LogKit::LogFormatted(LogKit::LOW,"  Modelling mode                           : forward\n");
  else if (modelSettings->getNumberOfSimulations() == 0)
    LogKit::LogFormatted(LogKit::LOW,"  Modelling mode                           : prediction\n"); 
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"  Modelling mode                           : simulation\n");
    if(inputFiles->getSeedFile()=="") {
      if (modelSettings->getSeed() == 0)
        LogKit::LogFormatted(LogKit::LOW,"  Seed                                     :          0 (default seed)\n");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Seed                                     : %10d\n",modelSettings->getSeed());
    }
    else
      LogKit::LogFormatted(LogKit::LOW,"  Seed read from file                      : %10s\n",inputFiles->getSeedFile().c_str());


    LogKit::LogFormatted(LogKit::LOW,"  Number of realisations                   : %10d\n",modelSettings->getNumberOfSimulations());
  }
  if(modelSettings->getForwardModeling()==false)
  {
    LogKit::LogFormatted(LogKit::LOW,"  Kriging                                  : %10s\n",(modelSettings->getKrigingParameter()>0 ? "yes" : "no"));
    if (modelSettings->getEstimateFaciesProb()) {
      if (modelSettings->getFaciesProbRelative())
        LogKit::LogFormatted(LogKit::LOW,"  Facies probabilities                     : %10s\n","relative");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Facies probabilities                     : %10s\n","absolute");
    }
    LogKit::LogFormatted(LogKit::LOW,"  Synthetic seismic                        : %10s\n",(modelSettings->getGenerateSeismicAfterInv() ? "yes" : "no" ));
  }

  LogKit::LogFormatted(LogKit::LOW,"\nInput/Output settings:\n");
  std::string logText("*NONE*");
  int logLevel = modelSettings->getLogLevel();
  if (logLevel == LogKit::L_ERROR)
    logText = "ERROR";
  else if (logLevel == LogKit::L_WARNING)
    logText = "WARNING";
  else if (logLevel == LogKit::L_LOW)
    logText = "LOW";
  else if (logLevel == LogKit::L_MEDIUM)
    logText = "MEDIUM";
  else if (logLevel == LogKit::L_HIGH)
    logText = "HIGH";
  else if (logLevel == LogKit::L_DEBUGLOW)
     logText = "DEBUGLOW";
  else if (logLevel == LogKit::L_DEBUGHIGH)
    logText = "DEBUGHIGH";
  LogKit::LogFormatted(LogKit::LOW, "  Log level                                : %10s\n",logText.c_str());
  LogKit::LogFormatted(LogKit::HIGH,"  Input directory                          : %10s\n",inputFiles->getInputDirectory().c_str());
  LogKit::LogFormatted(LogKit::HIGH,"  Output directory                         : %10s\n",IO::getOutputPath().c_str());

  int outputFlag = modelSettings->getGridOutputFlag();
  int gridFormat = modelSettings->getGridOutputFormat();

  if (outputFlag > 0) {
    LogKit::LogFormatted(LogKit::MEDIUM,"\nGrid output formats:\n");
    if (gridFormat & IO::SEGY) {
      const std::string & formatName = modelSettings->getTraceHeaderFormatOutput()->GetFormatName();
      LogKit::LogFormatted(LogKit::MEDIUM,"  Segy - %-10s                        :        yes\n",formatName.c_str());
    }
    if (gridFormat & IO::STORM)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Storm                                    :        yes\n");
    if (gridFormat & IO::ASCII)
      LogKit::LogFormatted(LogKit::MEDIUM,"  ASCII                                    :        yes\n");
    if (gridFormat & IO::SGRI)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Norsar                                   :        yes\n");
    if (gridFormat & IO::CRAVA)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Crava                                    :        yes\n");

    LogKit::LogFormatted(LogKit::MEDIUM,"\nGrid output domains:\n");
    if (gridFormat & IO::TIMEDOMAIN)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Time                                     :        yes\n");
    if (gridFormat & IO::DEPTHDOMAIN)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Depth                                    :        yes\n");
  }

  if (modelSettings->getDoInversion()) {
    LogKit::LogFormatted(LogKit::MEDIUM,"\nOutput of elastic parameters:\n");
    if ((outputFlag & IO::VP) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Pressure-wave velocity  (Vp)             :        yes\n");
    if ((outputFlag & IO::VS) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Shear-wave velocity  (Vs)                :        yes\n");
    if ((outputFlag & IO::RHO) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Density  (Rho)                           :        yes\n");
    if ((outputFlag & IO::AI) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Acoustic impedance  (AI)                 :        yes\n");
    if ((outputFlag & IO::VPVSRATIO) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Vp/Vs ratio                              :        yes\n");
    if ((outputFlag & IO::SI) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Shear impedance  (SI)                    :        yes\n");
    if ((outputFlag & IO::MURHO) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  MuRho  (SI*SI)                           :        yes\n");
    if ((outputFlag & IO::LAMBDARHO) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  LambdaRho  (AI*AI - 2*SI*SI)             :        yes\n");
    if ((outputFlag & IO::LAMELAMBDA) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Lame's dirst parameter                   :        yes\n");
    if ((outputFlag & IO::LAMEMU) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Lame's second parameter (shear modulus)  :        yes\n");
    if ((outputFlag & IO::POISSONRATIO) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Poisson ratio                            :        yes\n");
    if ((outputFlag & IO::BACKGROUND) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Background (Vp, Vs, Rho)                 :        yes\n");
    if ((outputFlag & IO::BACKGROUND_TREND) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Background trend (Vp, Vs, Rho)           :        yes\n");
  }    

  if (modelSettings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::MEDIUM,"\nOutput of facies probability volumes:\n");
    if (modelSettings->getFaciesProbRelative())
      LogKit::LogFormatted(LogKit::MEDIUM,"  Use rel. amplitudes for elastic param.   :        yes\n");
    else
      LogKit::LogFormatted(LogKit::MEDIUM,"  Use abs. amplitudes for elastic param.   :        yes\n");
  }

  if (modelSettings->getForwardModeling() ||
      (outputFlag & IO::SYNTHETIC_SEISMIC_DATA) > 0 ||
      (outputFlag & IO::ORIGINAL_SEISMIC_DATA) > 0 ||
      (outputFlag & IO::RESIDUAL) > 0) {
    LogKit::LogFormatted(LogKit::MEDIUM,"\nOutput of seismic data:\n");
    if ((outputFlag & IO::SYNTHETIC_SEISMIC_DATA) > 0 || modelSettings->getForwardModeling())
      LogKit::LogFormatted(LogKit::MEDIUM,"  Synthetic seismic data (forward modelled):        yes\n");
    if ((outputFlag & IO::ORIGINAL_SEISMIC_DATA) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Original seismic data (in output grid)   :        yes\n");
    if ((outputFlag & IO::RESIDUAL) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Seismic data residuals                   :        yes\n");
  }

  bool otherGridOutput 
    = (outputFlag & IO::CORRELATION) > 0
    ||(outputFlag & IO::EXTRA_GRIDS) > 0
    ||(outputFlag & IO::TIME_TO_DEPTH_VELOCITY) > 0;
  
  if (otherGridOutput) {
    if ((outputFlag & IO::CORRELATION) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Posterior correlations                   :        yes\n");
    if ((outputFlag & IO::EXTRA_GRIDS) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Help grids (see use manual)              :        yes\n");
    if ((outputFlag & IO::TIME_TO_DEPTH_VELOCITY) > 0)
      LogKit::LogFormatted(LogKit::MEDIUM,"  Time-to-depth velocity                   :        yes\n");
  }

  // NBNB-PAL: Vi får utvide testen nedenfor etter hvert...
  if (modelSettings->getFileGrid()) {
    LogKit::LogFormatted(LogKit::MEDIUM,"\nAdvanced settings:\n");
    LogKit::LogFormatted(LogKit::MEDIUM, "  Use intermediate disk storage for grids  :        yes\n");
  }

  LogKit::LogFormatted(LogKit::HIGH,"\nUnit settings/assumptions:\n");
  LogKit::LogFormatted(LogKit::HIGH,"  Time                                     : %10s\n","ms TWT");
  LogKit::LogFormatted(LogKit::HIGH,"  Frequency                                : %10s\n","Hz");
  LogKit::LogFormatted(LogKit::HIGH,"  Length                                   : %10s\n","m");
  LogKit::LogFormatted(LogKit::HIGH,"  Velocities                               : %10s\n","m/s");
  LogKit::LogFormatted(LogKit::HIGH,"  Density                                  : %10s\n","g/cm3");
  LogKit::LogFormatted(LogKit::HIGH,"  Angles                                   : %10s\n","   degrees (clockwise relative to north when applicable)");

  //
  // WELL PROCESSING
  //
  if (modelSettings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::HIGH,"\nSettings for well processing:\n");
    LogKit::LogFormatted(LogKit::HIGH,"  Threshold for merging log entries        : %10.2f ms\n",modelSettings->getMaxMergeDist());
    LogKit::LogFormatted(LogKit::HIGH,"  Threshold for Vp-Vs rank correlation     : %10.2f\n",modelSettings->getMaxRankCorr());
    LogKit::LogFormatted(LogKit::HIGH,"  Threshold for deviation angle            : %10.1f (=%.2fm/ms TWT)\n",
                         modelSettings->getMaxDevAngle(),tan(modelSettings->getMaxDevAngle()*M_PI/180.0));
    LogKit::LogFormatted(LogKit::HIGH,"  High cut for background modelling        : %10.1f\n",modelSettings->getMaxHzBackground());
    LogKit::LogFormatted(LogKit::HIGH,"  High cut for seismic resolution          : %10.1f\n",modelSettings->getMaxHzSeismic());
  }
  LogKit::LogFormatted(LogKit::HIGH,"\nRange of allowed parameter values:\n");
  LogKit::LogFormatted(LogKit::HIGH,"  Vp  - min                                : %10.0f\n",modelSettings->getAlphaMin());
  LogKit::LogFormatted(LogKit::HIGH,"  Vp  - max                                : %10.0f\n",modelSettings->getAlphaMax());
  LogKit::LogFormatted(LogKit::HIGH,"  Vs  - min                                : %10.0f\n",modelSettings->getBetaMin());
  LogKit::LogFormatted(LogKit::HIGH,"  Vs  - max                                : %10.0f\n",modelSettings->getBetaMax());
  LogKit::LogFormatted(LogKit::HIGH,"  Rho - min                                : %10.1f\n",modelSettings->getRhoMin());
  LogKit::LogFormatted(LogKit::HIGH,"  Rho - max                                : %10.1f\n",modelSettings->getRhoMax());  

  //
  // WELL DATA
  //
  if (modelSettings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nWell logs:\n");
    const std::vector<std::string> & logNames = modelSettings->getLogNames();

    if (logNames.size() > 0)
    {
      LogKit::LogFormatted(LogKit::LOW,"  Time                                     : %10s\n",  logNames[0].c_str());
      if(NRLib::Uppercase(logNames[1])=="VP" || 
         NRLib::Uppercase(logNames[1])=="LFP_VP")
        LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity                          : %10s\n",logNames[1].c_str());
      else
        LogKit::LogFormatted(LogKit::LOW,"  Sonic                                    : %10s\n",logNames[1].c_str());
      if(NRLib::Uppercase(logNames[3])=="VS" || 
         NRLib::Uppercase(logNames[3])=="LFP_VS")
        LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity                          : %10s\n",logNames[3].c_str());
      else
        LogKit::LogFormatted(LogKit::LOW,"  Shear sonic                              : %10s\n",logNames[3].c_str());
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10s\n",  logNames[2].c_str());
      if (modelSettings->getFaciesLogGiven())
        LogKit::LogFormatted(LogKit::LOW,"  Facies                                   : %10s\n",logNames[4].c_str());
    }
    else
    {
      LogKit::LogFormatted(LogKit::LOW,"  Time                                     : %10s\n","TWT");
      LogKit::LogFormatted(LogKit::LOW,"  Sonic                                    : %10s\n","DT");
      LogKit::LogFormatted(LogKit::LOW,"  Shear sonic                              : %10s\n","DTS");
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10s\n","RHOB");
      LogKit::LogFormatted(LogKit::LOW,"  Facies                                   : %10s\n","FACIES");
    }
    LogKit::LogFormatted(LogKit::LOW,"\nWell files:\n");
    for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++) 
    {
      LogKit::LogFormatted(LogKit::LOW,"  %-2d                                       : %s\n",i+1,inputFiles->getWellFile(i).c_str());
    }
    bool generateBackground = modelSettings->getGenerateBackground();
    bool estimateFaciesProb = modelSettings->getFaciesLogGiven();
    bool estimateWavelet    = false;
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
      estimateWavelet = estimateWavelet || modelSettings->getEstimateWavelet(i);
    if (generateBackground || estimateFaciesProb || estimateWavelet) 
    {
      LogKit::LogFormatted(LogKit::LOW,"\nUse well in estimation of:                   ");
      if (generateBackground) LogKit::LogFormatted(LogKit::LOW,"BackgroundTrend  ");
      if (estimateWavelet)    LogKit::LogFormatted(LogKit::LOW,"WaveletEstimation  ");
      if (estimateFaciesProb) LogKit::LogFormatted(LogKit::LOW,"FaciesProbabilities");
      LogKit::LogFormatted(LogKit::LOW,"\n");
      for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++) 
      {
        LogKit::LogFormatted(LogKit::LOW,"  %-2d                                       : ",i+1);
        if (generateBackground) {
          if (modelSettings->getIndicatorBGTrend(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::LOW,"    %-11s  ","yes");
          else if (modelSettings->getIndicatorBGTrend(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::LOW,"    %-11s  ","no");
          else
            LogKit::LogFormatted(LogKit::LOW,"    %-11s  ","yes");            
        }
        if (estimateWavelet) {
          if (modelSettings->getIndicatorWavelet(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::LOW,"    %-13s  ","yes");
          else if (modelSettings->getIndicatorWavelet(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::LOW,"    %-13s  ","no");
          else
            LogKit::LogFormatted(LogKit::LOW,"    %-13s  ","if possible");
        }
        if (estimateFaciesProb) {
          if (modelSettings->getIndicatorFacies(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::LOW,"    %-12s","yes");
          else if (modelSettings->getIndicatorFacies(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::LOW,"    %-12s","no");
          else
            LogKit::LogFormatted(LogKit::LOW,"    %-12s","if possible");
        }
        LogKit::LogFormatted(LogKit::LOW,"\n");
      }
    }
    if ( modelSettings->getOptimizeWellLocation() ) 
    {
      LogKit::LogFormatted(LogKit::LOW,"\nFor well, optimize location to             : Angle with Weight\n");
      for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++) 
      {
        int nMoveAngles = modelSettings->getNumberOfWellAngles(i);
        if( nMoveAngles > 0 )
        {
          LogKit::LogFormatted(LogKit::LOW," %2d %46.1f %10.1f\n",i+1,(modelSettings->getWellMoveAngle(i,0)*180/M_PI),modelSettings->getWellMoveWeight(i,0));
          for (int j=1; j<nMoveAngles; j++)
            LogKit::LogFormatted(LogKit::LOW," %49.1f %10.1f\n",(modelSettings->getWellMoveAngle(i,j)*180/M_PI),modelSettings->getWellMoveWeight(i,j));
        }
        LogKit::LogFormatted(LogKit::LOW,"\n");
      }
    }
  }

  //
  // AREA
  // 
  std::string gridFile;
  int areaSpecification = modelSettings->getAreaSpecification();
  if(modelSettings->getForwardModeling()) {
    LogKit::LogFormatted(LogKit::LOW,"\nSeismic area:\n");
    gridFile = inputFiles->getBackFile(0);    // Get geometry from earth model (Vp)
  }
  else {
    LogKit::LogFormatted(LogKit::LOW,"\nInversion area");
    if(areaSpecification == ModelSettings::AREA_FROM_GRID_DATA)
      gridFile = inputFiles->getSeismicFile(0); // Get area from first seismic data volume
  }
  if (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA) {
    const std::vector<int> & areaILXL = modelSettings->getAreaILXL();
    LogKit::LogFormatted(LogKit::LOW," taken from grid\n");
    LogKit::LogFormatted(LogKit::LOW,"  Grid                                     : "+gridFile+"\n");
    if(areaILXL.size()>0)
    {
    if (areaILXL[0] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  In-line start                            : %10d\n", areaILXL[0]);
    if (areaILXL[1] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  In-line end                              : %10d\n", areaILXL[1]);
    if (areaILXL[4] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  In-line step                             : %10d\n", areaILXL[4]);
    if (areaILXL[2] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  Cross-line start                         : %10d\n", areaILXL[2]);
    if (areaILXL[3] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  Cross-line end                           : %10d\n", areaILXL[3]);
    if (areaILXL[5] != IMISSING)
      LogKit::LogFormatted(LogKit::LOW,"  Cross-line step                          : %10d\n", areaILXL[5]);
    }
  }
  else if (areaSpecification == ModelSettings::AREA_FROM_UTM) {
    LogKit::LogFormatted(LogKit::LOW," given as UTM coordinates\n");
    const SegyGeometry * geometry = modelSettings->getAreaParameters();
    LogKit::LogFormatted(LogKit::LOW,"  Reference point x                        : %10.1f\n", geometry->GetX0());
    LogKit::LogFormatted(LogKit::LOW,"  Reference point y                        : %10.1f\n", geometry->GetY0());
    LogKit::LogFormatted(LogKit::LOW,"  Length x                                 : %10.1f\n", geometry->Getlx());
    LogKit::LogFormatted(LogKit::LOW,"  Length y                                 : %10.1f\n", geometry->Getly());
    LogKit::LogFormatted(LogKit::LOW,"  Sample density x                         : %10.1f\n", geometry->GetDx());
    LogKit::LogFormatted(LogKit::LOW,"  Sample density y                         : %10.1f\n", geometry->GetDy());
    LogKit::LogFormatted(LogKit::LOW,"  Azimuth                                  : %10.4f\n", geometry->GetAngle());
  }
  else if (areaSpecification == ModelSettings::AREA_FROM_SURFACE) {
    LogKit::LogFormatted(LogKit::LOW," taken from surface\n");
    LogKit::LogFormatted(LogKit::LOW,"  Reference surface                        : "+inputFiles->getAreaSurfaceFile()+"\n");
  }

  //
  // SURFACES
  // 
  LogKit::LogFormatted(LogKit::LOW,"\nTime surfaces:\n");
  if (modelSettings->getParallelTimeSurfaces())
  {
    LogKit::LogFormatted(LogKit::LOW,"  Surface                                  : "+inputFiles->getTimeSurfFile(0)+"\n");
    LogKit::LogFormatted(LogKit::LOW,"  Shift to top surface                     : %10.1f\n", modelSettings->getTimeDTop());
    LogKit::LogFormatted(LogKit::LOW,"  Time slice                               : %10.1f\n", modelSettings->getTimeLz());
    LogKit::LogFormatted(LogKit::LOW,"  Sampling density                         : %10.1f\n", modelSettings->getTimeDz());
    LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n",   int(modelSettings->getTimeLz()/modelSettings->getTimeDz()+0.5));
  }
  else
  {
    const std::string & topName  = inputFiles->getTimeSurfFile(0); 
    const std::string & baseName = inputFiles->getTimeSurfFile(1); 

    if (NRLib::IsNumber(topName))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topName.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Top surface                              : "+topName+"\n");

    if (NRLib::IsNumber(baseName))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n", atof(baseName.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Base surface                             : "+baseName+"\n");
      LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n", modelSettings->getTimeNz());

    LogKit::LogFormatted(LogKit::LOW,"  Minimum allowed value for lmin/lmax      : %10.2f\n", modelSettings->getLzLimit());
  }
  if (inputFiles->getCorrDirFile() != "")
    LogKit::LogFormatted(LogKit::LOW,"\n  Correlation direction                    : "+inputFiles->getCorrDirFile()+"\n");

  if (modelSettings->getDoDepthConversion())
  {
    LogKit::LogFormatted(LogKit::LOW,"\nDepth conversion:\n");
    if (inputFiles->getDepthSurfFile(0) != "")
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : "+inputFiles->getDepthSurfFile(0)+"\n");
    else
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : %s\n", "Made from base depth surface and velocity field");
    if (inputFiles->getDepthSurfFile(1) != "")
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : "+inputFiles->getDepthSurfFile(1)+"\n");
    else
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : %s\n", "Made from top depth surface and velocity field");
    std::string velocityField = inputFiles->getVelocityField();
    if (modelSettings->getVelocityFromInversion()) {
      velocityField = "Use Vp from inversion";
    }
     LogKit::LogFormatted(LogKit::LOW,"  Velocity field                           : "+velocityField+"\n");
  }

  const std::string & topWEI  = inputFiles->getWaveletEstIntFile(0);
  const std::string & baseWEI = inputFiles->getWaveletEstIntFile(1);

  if (topWEI != "" || baseWEI != "") {
    LogKit::LogFormatted(LogKit::LOW,"\nWavelet estimation interval:\n");
    if (NRLib::IsNumber(topWEI))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : "+topWEI+"\n");
    
    if (NRLib::IsNumber(baseWEI))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : "+baseWEI+"\n");
  }

  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" || baseFEI != "") {
    LogKit::LogFormatted(LogKit::LOW,"\nFacies estimation interval:\n");
    if (NRLib::IsNumber(topFEI))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : "+topFEI+"\n");
    
    if (NRLib::IsNumber(baseFEI))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : "+baseFEI+"\n");
  }

  //
  // BACKGROUND
  //
  if (modelSettings->getGenerateBackground()) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBackground model (estimated):\n");
    if (inputFiles->getBackVelFile() != "")
      LogKit::LogFormatted(LogKit::LOW,"  Trend for p-wave velocity                : "+inputFiles->getBackVelFile()+"\n");
    Vario       * vario  = modelSettings->getBackgroundVario();
    GenExpVario * pVario = dynamic_cast<GenExpVario*>(vario);
    LogKit::LogFormatted(LogKit::LOW,"  Variogram\n");
    LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",(vario->getType()).c_str());
    if (pVario != NULL)
    LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pVario->getPower());
    LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",vario->getRange());
    if (vario->getAnisotropic()) 
    {
      LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",vario->getSubRange());
      LogKit::LogFormatted(LogKit::LOW,"    Azimuth                                : %10.1f\n",90.0 - vario->getAngle()*(180/M_PI));
    }
    LogKit::LogFormatted(LogKit::LOW,"  High cut frequency for well logs         : %10.1f\n",modelSettings->getMaxHzBackground());
  }
  else
  {
    if(modelSettings->getForwardModeling()==true)
      LogKit::LogFormatted(LogKit::LOW,"\nEarth model:\n");
    else
      LogKit::LogFormatted(LogKit::LOW,"\nBackground model:\n");
    if (modelSettings->getConstBackValue(0) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  P-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(0));
    else
      LogKit::LogFormatted(LogKit::LOW,"  P-wave velocity read from file           : "+inputFiles->getBackFile(0)+"\n");
    
    if (modelSettings->getConstBackValue(1) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  S-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(1));
    else
      LogKit::LogFormatted(LogKit::LOW,"  S-wave velocity read from file           : "+inputFiles->getBackFile(1)+"\n");
      
    if (modelSettings->getConstBackValue(2) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10.1f\n",modelSettings->getConstBackValue(2));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Density read from file                   : "+inputFiles->getBackFile(2)+"\n");
  }

  TraceHeaderFormat * thf_old = modelSettings->getTraceHeaderFormat();
  if (thf_old != NULL) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nAdditional SegY trace header format:\n");
    if (thf_old != NULL) {
      LogKit::LogFormatted(LogKit::LOW,"  Format name                              : "+thf_old->GetFormatName()+"\n");
      if (thf_old->GetBypassCoordScaling())
        LogKit::LogFormatted(LogKit::LOW,"  Bypass coordinate scaling                :        yes\n");
      if (!thf_old->GetStandardType()) 
      {
        LogKit::LogFormatted(LogKit::LOW,"  Start pos coordinate scaling             : %10d\n",thf_old->GetScalCoLoc());
        LogKit::LogFormatted(LogKit::LOW,"  Start pos trace x coordinate             : %10d\n",thf_old->GetUtmxLoc());
        LogKit::LogFormatted(LogKit::LOW,"  Start pos trace y coordinate             : %10d\n",thf_old->GetUtmyLoc());
        LogKit::LogFormatted(LogKit::LOW,"  Start pos inline index                   : %10d\n",thf_old->GetInlineLoc());
        LogKit::LogFormatted(LogKit::LOW,"  Start pos crossline index                : %10d\n",thf_old->GetCrosslineLoc());
        LogKit::LogFormatted(LogKit::LOW,"  Coordinate system                        : %10s\n",thf_old->GetCoordSys()==0 ? "UTM" : "ILXL" );
      }
    }
  }

  if (modelSettings->getForwardModeling())
  {
    //
    // SEISMIC
    //
    LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for seismic:\n");
    LogKit::LogFormatted(LogKit::LOW,"  Generating seismic                       : %10s\n","yes");
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {
      LogKit::LogFormatted(LogKit::LOW,"\nSettings for AVO stack %d:\n",i+1);
      LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10.1f\n",(modelSettings->getAngle(i)*180/M_PI));
      LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : "+inputFiles->getWaveletFile(i)+"\n");
    }
  }
  else
  {
    //
    // PRIOR CORRELATION
    //
    Vario * corr = modelSettings->getLateralCorr();
    if (corr != NULL) {
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::LogFormatted(LogKit::LOW,"\nPrior correlation (of residuals):\n");
      LogKit::LogFormatted(LogKit::LOW,"  Range of allowed parameter values:\n");
      LogKit::LogFormatted(LogKit::LOW,"    Var{Vp}  - min                         : %10.1e\n",modelSettings->getVarAlphaMin());
      LogKit::LogFormatted(LogKit::LOW,"    Var{Vp}  - max                         : %10.1e\n",modelSettings->getVarAlphaMax());
      LogKit::LogFormatted(LogKit::LOW,"    Var{Vs}  - min                         : %10.1e\n",modelSettings->getVarBetaMin());
      LogKit::LogFormatted(LogKit::LOW,"    Var{Vs}  - max                         : %10.1e\n",modelSettings->getVarBetaMax());
      LogKit::LogFormatted(LogKit::LOW,"    Var{Rho} - min                         : %10.1e\n",modelSettings->getVarRhoMin());
      LogKit::LogFormatted(LogKit::LOW,"    Var{Rho} - max                         : %10.1e\n",modelSettings->getVarRhoMax());  
      LogKit::LogFormatted(LogKit::LOW,"  Lateral correlation:\n");
      LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",corr->getRange());
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",corr->getSubRange());
        LogKit::LogFormatted(LogKit::LOW,"    Azimuth                                : %10.1f\n",90.0 - corr->getAngle()*(180/M_PI));
      }
    }
    //
    // PRIOR FACIES
    //
    if (modelSettings->getIsPriorFaciesProbGiven() > 0)
    {
      LogKit::LogFormatted(LogKit::LOW,"\nPrior facies probabilities:\n");      
      if(modelSettings->getIsPriorFaciesProbGiven()==1)
      {
        typedef std::map<std::string,float> mapType;
        mapType myMap = modelSettings->getPriorFaciesProb();
        
        for(mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::LOW,"   %-12s                            : %10.2f\n",(i->first).c_str(),i->second);
      }
      else if (modelSettings->getIsPriorFaciesProbGiven()==2)
      {
        typedef std::map<std::string,std::string> mapType;
        mapType myMap = inputFiles->getPriorFaciesProbFile();

        for(mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::LOW,"   %-12s                            : %10s\n",(i->first).c_str(),(i->second).c_str());
      }
    }
    //
    // SEISMIC
    //
    if (modelSettings->getNoSeismicNeeded()==false)
    {
      LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for seismic:\n");
      LogKit::LogFormatted(LogKit::LOW,"  White noise component                    : %10.2f\n",modelSettings->getWNC());
      LogKit::LogFormatted(LogKit::LOW,"  Low cut for inversion                    : %10.1f\n",modelSettings->getLowCut());
      LogKit::LogFormatted(LogKit::LOW,"  High cut for inversion                   : %10.1f\n",modelSettings->getHighCut());
      corr  = modelSettings->getAngularCorr();
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::LogFormatted(LogKit::LOW,"  Angular correlation:\n");
      LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",corr->getRange()*180.0/M_PI);
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",corr->getSubRange()*180.0/M_PI);
        LogKit::LogFormatted(LogKit::LOW,"    Angle                                  : %10.1f\n",corr->getAngle());
      }
      bool estimateNoise = false;
      for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++) {
        estimateNoise = estimateNoise || modelSettings->getEstimateSNRatio(i); 
      }
      LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for wavelet:\n");
      if (estimateNoise)
        LogKit::LogFormatted(LogKit::LOW,"  Maximum shift in noise estimation        : %10.1f\n",modelSettings->getMaxWaveletShift());
      LogKit::LogFormatted(LogKit::LOW,"  Minimum relative amplitude               : %10.3f\n",modelSettings->getMinRelWaveletAmp());
      LogKit::LogFormatted(LogKit::LOW,"  Wavelet tapering length                  : %10.1f\n",modelSettings->getWaveletTaperingL());
      if (modelSettings->getOptimizeWellLocation()) {
        LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for well locations:\n");
        LogKit::LogFormatted(LogKit::LOW,"  Maximum offset                           : %10.1f\n",modelSettings->getMaxWellOffset());
        LogKit::LogFormatted(LogKit::LOW,"  Maximum vertical shift                   : %10.1f\n",modelSettings->getMaxWellShift());
      }
      for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nSettings for AVO stack %d:\n",i+1);
        LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10.1f\n",(modelSettings->getAngle(i)*180/M_PI));
        LogKit::LogFormatted(LogKit::LOW,"  SegY start time                          : %10.1f\n",modelSettings->getSegyOffset());
        TraceHeaderFormat * thf = modelSettings->getTraceHeaderFormat(i);
        if (thf != NULL) 
        {
          LogKit::LogFormatted(LogKit::LOW,"  SegY trace header format:\n");
          LogKit::LogFormatted(LogKit::LOW,"    Format name                            : "+thf->GetFormatName()+"\n");
          if (thf->GetBypassCoordScaling())
            LogKit::LogFormatted(LogKit::LOW,"    Bypass coordinate scaling              :        yes\n");
          if (!thf->GetStandardType()) 
          {
            LogKit::LogFormatted(LogKit::LOW,"    Start pos coordinate scaling           : %10d\n",thf->GetScalCoLoc());
            LogKit::LogFormatted(LogKit::LOW,"    Start pos trace x coordinate           : %10d\n",thf->GetUtmxLoc());
            LogKit::LogFormatted(LogKit::LOW,"    Start pos trace y coordinate           : %10d\n",thf->GetUtmyLoc());
            LogKit::LogFormatted(LogKit::LOW,"    Start pos inline index                 : %10d\n",thf->GetInlineLoc());
            LogKit::LogFormatted(LogKit::LOW,"    Start pos crossline index              : %10d\n",thf->GetCrosslineLoc());
            LogKit::LogFormatted(LogKit::LOW,"    Coordinate system                      : %10s\n",thf->GetCoordSys()==0 ? "UTM" : "ILXL" );
          }
        }
        LogKit::LogFormatted(LogKit::LOW,"  Data                                     : "+inputFiles->getSeismicFile(i)+"\n");
        if (modelSettings->getEstimateWavelet(i))
          LogKit::LogFormatted(LogKit::LOW,"  Estimate wavelet                         : %10s\n", "yes");
        else
          LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : "+inputFiles->getWaveletFile(i)+"\n");
        if (modelSettings->getEstimateLocalShift(i))
          LogKit::LogFormatted(LogKit::LOW,"  Estimate local shift map                 : %10s\n", "yes");
        else if (inputFiles->getShiftFile(i) != "")
          LogKit::LogFormatted(LogKit::LOW,"  Local shift map                          : "+inputFiles->getShiftFile(i)+"\n");  
        if (modelSettings->getEstimateLocalScale(i))
          LogKit::LogFormatted(LogKit::LOW,"  Estimate local scale map                 : %10s\n", "yes");
        else if (inputFiles->getScaleFile(i) != "")
          LogKit::LogFormatted(LogKit::LOW,"  Local scale map                          : "+inputFiles->getScaleFile(i)+"\n");      
        if (modelSettings->getMatchEnergies(i)) 
          LogKit::LogFormatted(LogKit::LOW,"  Match empirical and theoretical energies : %10s\n", "yes");
        if (!modelSettings->getEstimateWavelet(i) && !modelSettings->getMatchEnergies(i)) {
          if (modelSettings->getEstimateGlobalWaveletScale(i))
            LogKit::LogFormatted(LogKit::LOW,"  Estimate global wavelet scale            : %10s\n","yes");
          else
            LogKit::LogFormatted(LogKit::LOW,"  Global wavelet scale                     : %10.2e\n",modelSettings->getWaveletScale(i));
        }
        if (modelSettings->getEstimateSNRatio(i)) 
          LogKit::LogFormatted(LogKit::LOW,"  Estimate signal-to-noise ratio           : %10s\n", "yes");
        else
          LogKit::LogFormatted(LogKit::LOW,"  Signal-to-noise ratio                    : %10.1f\n",modelSettings->getSNRatio(i));
        if (modelSettings->getEstimateLocalNoise(i)) { 
          if (inputFiles->getLocalNoiseFile(i) == "") 
            LogKit::LogFormatted(LogKit::LOW,"  Estimate local signal-to-noise ratio map : %10s\n", "yes");
          else
            LogKit::LogFormatted(LogKit::LOW,"  Estimate local signal-to-noise ratio map : "+inputFiles->getLocalNoiseFile(i)+"\n");
        }
        if (modelSettings->getEstimateLocalNoise(i))
          LogKit::LogFormatted(LogKit::LOW,"  Estimate local noise                     : %10s\n", "yes");
        if (inputFiles->getLocalNoiseFile(i) != "")
          LogKit::LogFormatted(LogKit::LOW,"  Local noise                              : "+inputFiles->getLocalNoiseFile(i)+"\n");
      }
    }
  }
}

void
Model::getCorrGradIJ(float & corrGradI, float &corrGradJ) const
{
  double angle  = timeSimbox_->getAngle();
  double cosrot = cos(angle);
  double sinrot = sin(angle);
  double dx     = timeSimbox_->getdx();  
  double dy     = timeSimbox_->getdy();

  double cI =  dx*cosrot*gradX_ + dy*sinrot*gradY_;
  double cJ = -dx*sinrot*gradX_ + dy*cosrot*gradY_;

  corrGradI = float(cI/timeSimbox_->getdz());
  corrGradJ = float(cJ/timeSimbox_->getdz());
}

void
Model::processDepthConversion(Simbox        * timeCutSimbox, 
                              Simbox        * timeSimbox,
                              ModelSettings * modelSettings, 
                              InputFiles    * inputFiles,
                              std::string   & errText, 
                              bool          & failed)
{
  FFTGrid * velocity = NULL;
  if(timeCutSimbox != NULL)
    loadVelocity(velocity, timeCutSimbox, modelSettings, 
                 inputFiles->getVelocityField(), 
                 errText, failed);
  else
    loadVelocity(velocity, timeSimbox, modelSettings, 
                 inputFiles->getVelocityField(), 
                 errText, failed);
  
  if(!failed) 
  {
    timeDepthMapping_ = new GridMapping();
    timeDepthMapping_->setDepthSurfaces(inputFiles->getDepthSurfFiles(), failed, errText);
    if(velocity != NULL) 
    {
      velocity->setAccessMode(FFTGrid::RANDOMACCESS);
      timeDepthMapping_->calculateSurfaceFromVelocity(velocity, timeSimbox);
      timeDepthMapping_->setDepthSimbox(timeSimbox, timeSimbox->getnz(), 
                                        modelSettings->getGridOutputFormat(),
                                        failed, errText);            // NBNB-PAL: Er dettet riktig nz (timeCut vs time)? 
      timeDepthMapping_->makeTimeDepthMapping(velocity, timeSimbox);
      velocity->endAccess();

      if((modelSettings->getGridOutputFlag() & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
        std::string baseName  = IO::FileTimeToDepthVelocity();
        std::string sgriLabel = std::string("Time-to-depth velocity");
        float       offset    = modelSettings->getSegyOffset();
        velocity->writeFile(baseName, 
                            IO::PathToVelocity(), 
                            timeSimbox, 
                            sgriLabel, 
                            offset,
                            timeDepthMapping_, 
                            timeCutMapping_);
      }
    }
    else if (velocity==NULL && velocityFromInversion_==false)
    {
      timeDepthMapping_->setDepthSimbox(timeSimbox, 
                                        timeSimbox->getnz(), 
                                        modelSettings->getGridOutputFormat(),
                                        failed, 
                                        errText);

    }
  }
  if(velocity != NULL)
    delete velocity;
}

void 
Model::loadVelocity(FFTGrid          *& velocity,
                    Simbox            * timeSimbox,
                    ModelSettings     * modelSettings, 
                    const std::string & velocityField, 
                    std::string       & errText,
                    bool              & failed)
{
  Utils::writeHeader("Setup time-to-depth relationship");

  if(modelSettings->getVelocityFromInversion() == true)
  {
    velocityFromInversion_ = true;
    velocity = NULL;
  }
  else if(velocityField == "")
    velocity = NULL;
  else
  {
    const SegyGeometry      * dummy1 = NULL;
    const TraceHeaderFormat * dummy2 = NULL;
    const float               offset = modelSettings->getSegyOffset();
    std::string errorText("");
    readGridFromFile(velocityField,
                     "velocity field",
                     offset,
                     velocity,
                     dummy1,
                     dummy2,
                     FFTGrid::PARAMETER,
                     timeSimbox,
                     modelSettings,
                     errorText);
    
    if (errorText == "") { // No errors
      //
      // Check that the velocity grid is veldefined.
      //
      float logMin = modelSettings->getAlphaMin();
      float logMax = modelSettings->getAlphaMax();
      const int nzp = velocity->getNzp();
      const int nyp = velocity->getNyp();
      const int nxp = velocity->getNxp();
      const int nz = velocity->getNz();
      const int ny = velocity->getNy();
      const int nx = velocity->getNx();
      int tooLow  = 0; 
      int tooHigh = 0; 
      velocity->setAccessMode(FFTGrid::READ);
      int rnxp = 2*(nxp/2 + 1);
      for (int k = 0; k < nzp; k++) 
        for (int j = 0; j < nyp; j++)
          for (int i = 0; i < rnxp; i++) {
            if(i < nx && j < ny && k < nz) {
              float value = velocity->getNextReal();
              if (value < logMin && value != RMISSING) {
                tooLow++;
              }
              if (value > logMax && value != RMISSING)
                tooHigh++;
            }
          }
      velocity->endAccess();
      
      if (tooLow+tooHigh > 0) {
        std::string text;
        text += "\nThe velocity grid used as trend in the background model of Vp";
        text += "\ncontains too small and/or too high velocities:";
        text += "\n  Minimum Vp = "+NRLib::ToString(logMin,2)+"    Number of too low values  : "+NRLib::ToString(tooLow);
        text += "\n  Maximum Vp = "+NRLib::ToString(logMax,2)+"    Number of too high values : "+NRLib::ToString(tooHigh);
        text += "\nThe range of allowed values can changed using the ALLOWED_PARAMETER_VALUES keyword\n";
        text += "\naborting...\n";
        errText += "Reading of file '"+velocityField+"' for background velocity field failed.\n";
        errText += text;
        failed = true;
      } 
    }
    else {
      errorText += "Reading of file \'"+velocityField+"\' for background velocity field failed.\n";
      errText += errorText;
      failed = true;
    }
  }
}

void 
Model::writeAreas(const SegyGeometry * areaParams,
                  Simbox             * timeSimbox,
                  std::string        & text)
{
  double areaX0   = areaParams->GetX0();
  double areaY0   = areaParams->GetY0();
  double areaLx   = areaParams->Getlx();
  double areaLy   = areaParams->Getly();
  double areaDx   = areaParams->GetDx();
  double areaDy   = areaParams->GetDy();
  double areaRot  = areaParams->GetAngle();
  double areaXmin = RMISSING;
  double areaXmax = RMISSING;
  double areaYmin = RMISSING;
  double areaYmax = RMISSING;

  findSmallestSurfaceGeometry(areaX0, areaY0, areaLx, areaLy, areaRot,
                              areaXmin, areaYmin, areaXmax, areaYmax);

  LogKit::LogFormatted(LogKit::LOW,"\nThe top and/or base time surfaces do not cover the area specified by the "+text);
  LogKit::LogFormatted(LogKit::LOW,"\nPlease extrapolate surfaces or specify a smaller AREA in the model file.\n");
  LogKit::LogFormatted(LogKit::LOW,"\nArea/resolution           x0           y0            lx        ly     azimuth          dx      dy\n");
  LogKit::LogFormatted(LogKit::LOW,"-------------------------------------------------------------------------------------------------\n");
  double azimuth = (-1)*areaRot*(180.0/M_PI);
  if (azimuth < 0)
    azimuth += 360.0;
  LogKit::LogFormatted(LogKit::LOW,"Model area       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n\n", 
                       areaX0, areaY0, areaLx, areaLy, areaDx, areaDy, azimuth);

  LogKit::LogFormatted(LogKit::LOW,"Area                    xmin         xmax           ymin        ymax\n");
  LogKit::LogFormatted(LogKit::LOW,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"%-12s     %11.2f  %11.2f    %11.2f %11.2f\n", 
                       text.c_str(),areaXmin, areaXmax, areaYmin, areaYmax);
  const NRLib::Surface<double> & top  = timeSimbox->GetTopSurface();
  const NRLib::Surface<double> & base = timeSimbox->GetBotSurface();
  LogKit::LogFormatted(LogKit::LOW,"Top surface      %11.2f  %11.2f    %11.2f %11.2f\n", 
                       top.GetXMin(), top.GetXMax(), top.GetYMin(), top.GetYMax()); 
  LogKit::LogFormatted(LogKit::LOW,"Base surface     %11.2f  %11.2f    %11.2f %11.2f\n", 
                       base.GetXMin(), base.GetXMax(), base.GetYMin(), base.GetYMax()); 
}

void   
Model::findSmallestSurfaceGeometry(const double   x0, 
                                   const double   y0, 
                                   const double   lx, 
                                   const double   ly, 
                                   const double   rot,
                                   double       & xMin,
                                   double       & yMin,
                                   double       & xMax,
                                   double       & yMax)
{
  xMin = x0 - ly*sin(rot);
  xMax = x0 + lx*cos(rot);
  yMin = y0;
  yMax = y0 + lx*sin(rot) + ly*cos(rot);
  if (rot < 0) {
    xMin = x0;
    xMax = x0 + lx*cos(rot) - ly*sin(rot);
    yMin = y0 + lx*sin(rot);
    yMax = y0 + ly*cos(rot);
  }
}

void
Model::writeLocalGridsToFile(const std::string   & fileName,
                             const std::string   & type,
                             const float           scaleFactor,
                             const ModelSettings * modelSettings,
                             const unsigned int    i,
                             Simbox              * timeSimbox,
                             Grid2D             *& grid)  
{
  bool   estimationMode   = modelSettings->getEstimationMode();
  int    outputFormat     = modelSettings->getGridOutputFormat();
  double angle            = modelSettings->getAngle(i)*180.0/M_PI;

  Surface * help = NULL; 

  if(fileName != "") {
    help = new Surface(fileName);
    grid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
    resampleSurfaceToGrid2D(timeSimbox, help, grid);
  }
  else {
    if (grid != NULL) {
      resampleGrid2DToSurface(timeSimbox, grid, help);
    }
  }
  if ((estimationMode ||
    ((type==IO::PrefixLocalWaveletGain() || type==IO::PrefixLocalWaveletShift()) && (modelSettings->getWaveletOutputFlag() & IO::LOCAL_WAVELETS)>0) ||
    (type==IO::PrefixLocalNoise() && (modelSettings->getOtherOutputFlag() & IO::LOCAL_NOISE)>0)) && 
     help != NULL)
  {    
    std::string baseName = type + NRLib::ToString(angle, 1);
    help->Multiply(scaleFactor);
    if (type==IO::PrefixLocalNoise())
      IO::writeSurfaceToFile(*help, baseName, IO::PathToNoise(), outputFormat);
    else
      IO::writeSurfaceToFile(*help, baseName, IO::PathToWavelets(), outputFormat);
  }
  if (help != NULL)
    delete help;
}


void 
Model::resampleSurfaceToGrid2D(const Simbox  * simbox, 
                               const Surface * surface, 
                               Grid2D        * outgrid)
{
  for(int i=0;i<simbox->getnx();i++) {
    for(int j=0;j<simbox->getny();j++) {
      double x, y, z;
      simbox->getCoord(i, j, 0, x, y, z);
      (*outgrid)(i,j) = static_cast<float>(surface->GetZ(x,y));
    }
  }
}

void 
Model::resampleGrid2DToSurface(const Simbox   * simbox, 
                               const Grid2D   * grid,
                               Surface       *& surface)
{
  double xmin,xmax,ymin,ymax;
  simbox->getMinAndMaxXY(xmin,xmax,ymin,ymax);
  int nx,ny;
  double angle = simbox->getAngle()*180.0/M_PI;
  if(angle > -45 || angle < 45)
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::cos(simbox->getAngle())+0.5));
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::cos(simbox->getAngle())+0.5));
  }
  else
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::sin(simbox->getAngle())+0.5));
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::sin(simbox->getAngle())+0.5));
  }
  surface = new Surface(xmin,ymin,xmax-xmin,ymax-ymin,nx,ny,0.0);
  double x,y;
  int i1,j1;
  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      surface->GetXY(i,j,x,y);
      simbox->getIndexes(x,y,i1,j1);
      if(i1==IMISSING || j1== IMISSING)
        surface->SetMissing(i,j);
      else
        (*surface)(i,j) = (*grid)(i1,j1);
    }
  }
}

bool 
Model::findTimeGradientSurface(const std::string    & refTimeFile,
                               Simbox               * simbox,
                               NRLib::Grid2D<float> & refTimeGradX,
                               NRLib::Grid2D<float> & refTimeGradY)
{
  double x, y;
  bool inside = true;
  unsigned int nx = static_cast<unsigned int> (simbox->getnx());
  unsigned int ny = static_cast<unsigned int> (simbox->getny());
  float dx = static_cast<float> (simbox->getdx());
  float dy = static_cast<float> (simbox->getdy());

  if (!NRLib::IsNumber(refTimeFile)) {
    NRLib::RegularSurfaceRotated<float> t0surface(refTimeFile, NRLib::SURF_SGRI);

    simbox->getXYCoord(0,0,x,y);
    if (t0surface.IsInsideSurface(x,y)) {
      simbox->getXYCoord(0,ny-1,x,y);
      if (t0surface.IsInsideSurface(x,y)) {
        simbox->getXYCoord(nx-1,0,x,y);
        if(t0surface.IsInsideSurface(x,y)) {
          simbox->getXYCoord(nx-1,ny-1,x,y);
          if (!t0surface.IsInsideSurface(x,y))
            inside = false;
        }
        else
          inside = false;
      }
      else
        inside = false;
    }
    else
      inside = false;

    if (inside) {
      refTimeGradX.Resize(nx, ny, RMISSING);
      refTimeGradY.Resize(nx, ny, RMISSING);
      for (unsigned int i = nx-1; i >= 0; i++) {
        for (unsigned int j = ny-1; j >= 0; j++) {
          simbox->getXYCoord(i,j,x,y);
          float z_high = t0surface.GetZInside(x,y);
          simbox->getXYCoord(i,j-1,x,y); //XYCoord is ok even if j = -1, but point is outside simbox
          float z_low = t0surface.GetZInside(x,y);
          if (!t0surface.IsMissing(z_low))
            refTimeGradY(i,j) = (z_high - z_low) / dy;
          else
            refTimeGradY(i,j) = refTimeGradY(i,j+1);

          simbox->getXYCoord(i-1,j,x,y); //XYCoord is ok even if j = -1, but point is outside simbox
          z_low = t0surface.GetZInside(x,y);
          if (!t0surface.IsMissing(z_low))
            refTimeGradX(i,j) = (z_high - z_low) / dx;
          else
            refTimeGradX(i,j) = refTimeGradX(i+1,j);
        }
      }
    }
  }
  else {
    refTimeGradX.Resize(nx, ny, 0.0);
    refTimeGradY.Resize(nx, ny, 0.0);
  }

  return(inside);
}

void
Model::getGeometryFromGridOnFile(const std::string         gridFile,
                                 const TraceHeaderFormat * thf,
                                 SegyGeometry           *& geometry,
                                 std::string             & errText)
{
  geometry = NULL;
  if(gridFile != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int fileType = IO::findGridType(gridFile);
    if(fileType == IO::CRAVA) {
      geometry = geometryFromCravaFile(gridFile);
    }
    else if(fileType == IO::SEGY) {
      try
      {
        geometry = SegY::FindGridGeometry(gridFile, thf);
      }
      catch (NRLib::Exception & e)
      {
        errText = e.what();
      }
    }
    else if(fileType == IO::STORM)
      geometry = geometryFromStormFile(gridFile, errText);
    else if(fileType==IO::SGRI) {
      bool scale = true;
      geometry = geometryFromStormFile(gridFile, errText, scale);
    }
    else {
      errText = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    errText = "Cannot get geometry from file. The file name is empty.\n";
  }
}

SegyGeometry *
Model::geometryFromCravaFile(const std::string & fileName) 
{
  std::ifstream binFile;
  NRLib::OpenRead(binFile, fileName, std::ios::in | std::ios::binary);
  
  std::string fileType;
  getline(binFile,fileType);         

  double x0      = NRLib::ReadBinaryDouble(binFile);
  double y0      = NRLib::ReadBinaryDouble(binFile);
  double dx      = NRLib::ReadBinaryDouble(binFile);
  double dy      = NRLib::ReadBinaryDouble(binFile);
  int    nx      = NRLib::ReadBinaryInt(binFile);
  int    ny      = NRLib::ReadBinaryInt(binFile);
  double IL0     = NRLib::ReadBinaryDouble(binFile);
  double XL0     = NRLib::ReadBinaryDouble(binFile);
  double ilStepX = NRLib::ReadBinaryDouble(binFile);
  double ilStepY = NRLib::ReadBinaryDouble(binFile);
  double xlStepX = NRLib::ReadBinaryDouble(binFile);
  double xlStepY = NRLib::ReadBinaryDouble(binFile);
  double rot     = NRLib::ReadBinaryDouble(binFile);

  binFile.close();
  
  SegyGeometry * geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                             IL0, XL0, ilStepX, ilStepY, 
                                             xlStepX, xlStepY, rot);
  return(geometry);
}

SegyGeometry *
Model::geometryFromStormFile(const std::string & fileName,
                             std::string       & errText,
                             bool scale) 
{
  SegyGeometry  * geometry  = NULL;
  StormContGrid * stormgrid = NULL;
  std::string     tmpErrText;
  float scalehor;
  if(scale==false)
  {
    scalehor = 1.0;
  }
  else //from sgri file
  {
    LogKit::LogFormatted(LogKit::LOW,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scalehor  = 1000.0;
  }
  try
  {   
    stormgrid = new StormContGrid(0,0,0);
    stormgrid->ReadFromFile(fileName);
    stormgrid->SetMissingCode(RMISSING);
  }
  catch (NRLib::Exception & e) 
  {
    tmpErrText = e.what();
  }

  if (tmpErrText == "") {
    double x0      = stormgrid->GetXMin()*scalehor;
    double y0      = stormgrid->GetYMin()*scalehor;
    double dx      = stormgrid->GetDX()*scalehor;
    double dy      = stormgrid->GetDY()*scalehor;
    int    nx      = stormgrid->GetNI();
    int    ny      = stormgrid->GetNJ();
    double rot     = stormgrid->GetAngle();
    double IL0     = 0.0;  ///< Dummy value since information is not contained in format
    double XL0     = 0.0;  ///< Dummy value since information is not contained in format
    double ilStepX =   1;  ///< Dummy value since information is not contained in format
    double ilStepY =   1;  ///< Dummy value since information is not contained in format
    double xlStepX =   1;  ///< Dummy value since information is not contained in format
    double xlStepY =   1;  ///< Dummy value since information is not contained in format
    geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                IL0, XL0, ilStepX, ilStepY, 
                                xlStepX, xlStepY, rot);
  }
  else {
    errText += tmpErrText;
  }

  if (stormgrid != NULL)
    delete stormgrid;
  
  return(geometry);
}

FFTGrid*            
Model::createFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid)
{
  FFTGrid* fftGrid;

  if(fileGrid)
    fftGrid =  new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fftGrid =  new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

  return(fftGrid);
}
