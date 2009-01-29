
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
#include "src/modelfile.h"
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

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/global_def.h"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"

Model::Model(char * fileName)
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
  shiftGrids_             = NULL;
  gainGrids_              = NULL;
  waveletEstimInterval_   = NULL;
  faciesEstimInterval_    = NULL; 
  correlationDirection_   = NULL;     
  reflectionMatrix_       = NULL;
  randomGen_              = NULL;
  failed_                 = false;
  gradX_                  = 0.0;
  gradY_                  = 0.0;
 
  timeDepthMapping_       = NULL;
  timeCutMapping_         = NULL;
  velocityFromInversion_  = false;

  bool failedModelFile    = false;
  bool failedWavelet      = false;
  bool failedSeismic      = false;
  bool failedSimbox       = false;
  bool failedWells        = false;
  bool failedReflMat      = false;
  bool failedExtraSurf    = false;
  bool failedBackground   = false;
  bool failedPriorCorr    = false;
  bool failedDepthConv    = false;
  bool failedLoadingModel = false;

  Simbox * timeCutSimbox  = NULL;
  Simbox * timeBGSimbox   = NULL;

  ModelFile modelFile(fileName);

  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
  else
  {
    InputFiles * inputFiles;
    inputFiles     = modelFile.getInputFiles();
    modelSettings_ = modelFile.getModelSettings();
    
    LogKit::SetScreenLog(modelSettings_->getLogLevel());

    std::string logFileName = ModelSettings::makeFullFileName("logFile.txt");
    LogKit::SetFileLog(logFileName,modelSettings_->getLogLevel());

    if(modelSettings_->getDebugFlag() > 0)
    {
      std::string fName = ModelSettings::makeFullFileName("debug.txt");
      LogKit::SetFileLog(fName, LogKit::DEBUGHIGH+LogKit::DEBUGLOW);
    }
    LogKit::EndBuffering();
    
    if(modelSettings_->getLocalWaveletVario()==NULL) 
      modelSettings_->copyBackgroundVarioToLocalWaveletVario();

    if(inputFiles->getSeedFile() == "")
      randomGen_ = new RandomGen(modelSettings_->getSeed());
    else
      randomGen_ = new RandomGen(inputFiles->getSeedFile().c_str());

    if(modelSettings_->getNumberOfSimulations() == 0)
      modelSettings_->setOutputFlag(ModelSettings::PREDICTION); //write predicted grids. 
    
    printSettings(modelSettings_, inputFiles);
    
    Utils::writeHeader("Reading input data");

    char errText[MAX_STRING];
    sprintf(errText,"%c",'\0');

    makeTimeSimboxes(timeSimbox_, timeCutSimbox, timeBGSimbox, timeSimboxConstThick_,  //Handles correlation direction too.
                     correlationDirection_, modelSettings_, inputFiles, 
                     errText, failedSimbox);

    if(!failedSimbox)
    { 
      //
      // FORWARD MODELLING
      //
      if (modelSettings_->getGenerateSeismic() == true)
      {
        processBackground(background_, wells_, timeSimbox_, timeBGSimbox,
                          modelSettings_, inputFiles,
                          errText, failedBackground);
        if (!failedBackground)
        {
          processReflectionMatrix(reflectionMatrix_, background_, 
                                  modelSettings_, inputFiles,
                                  errText, failedReflMat);  
          if (!failedReflMat)
          {
            processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                            timeSimbox_, waveletEstimInterval_, shiftGrids_, gainGrids_,
                            modelSettings_, inputFiles, errText, failedWavelet);
          }              
        }
        background_->getAlpha()->setOutputFormat(modelSettings_->getFormatFlag()); //static, controls all grids.
      }
      else
      {
        //
        // INVERSION
        //
        processSeismic(seisCube_, timeSimbox_, 
                       modelSettings_, inputFiles,
                       errText, failedSeismic);
        
        if(failedSeismic==false)
        {
          if(timeCutSimbox!=NULL)  {
            timeCutMapping_ = new GridMapping();
            timeCutMapping_->makeTimeTimeMapping(timeCutSimbox);
          }
          processDepthConversion(timeCutSimbox, timeSimbox_,
                                 modelSettings_, inputFiles,
                                 errText, failedDepthConv);
          processWells(wells_, timeSimbox_, timeBGSimbox, timeSimboxConstThick_, 
                       seisCube_, randomGen_, modelSettings_, inputFiles,
                       errText, failedWells);
          loadExtraSurfaces(waveletEstimInterval_, 
                            faciesEstimInterval_,
                            timeSimbox_, inputFiles,
                            errText, failedExtraSurf);
          if (!failedWells && !failedDepthConv)
          {
            processBackground(background_, wells_, timeSimbox_, timeBGSimbox,
                              modelSettings_, inputFiles,
                              errText, failedBackground);

            if (!failedBackground)
            {
              processPriorCorrelations(correlations_, background_, wells_, timeSimbox_,
                                       modelSettings_, inputFiles,
                                       errText, failedPriorCorr);
              processReflectionMatrix(reflectionMatrix_, background_, 
                                      modelSettings_, inputFiles, 
                                      errText, failedReflMat);
              if (!failedReflMat && !failedExtraSurf)
              {
                processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                                timeSimbox_, waveletEstimInterval_, shiftGrids_, gainGrids_,
                                modelSettings_, inputFiles, errText, failedWavelet);
              }
            }
          }
        }
        if (!failedWells && !failedExtraSurf)
        {
          processPriorFaciesProb(priorFacies_,
                                 wells_,
                                 randomGen_,
                                 timeSimbox_->getnz(),
                                 modelSettings_);
        }
      }
    }

    failedLoadingModel = failedSimbox  || failedSeismic   || failedPriorCorr  ||
                         failedWells   || failedReflMat   || failedBackground ||
                         failedWavelet || failedDepthConv || failedExtraSurf;

    if (failedLoadingModel) 
      LogKit::LogFormatted(LogKit::ERROR,"\nERROR(s) while loading model. \n %s", errText);
  }
  failed_ = failedModelFile || failedLoadingModel;
  
  if(timeCutSimbox != NULL)
    delete timeCutSimbox;
  if(timeBGSimbox != NULL)
    delete timeBGSimbox;
}


Model::~Model(void)
{
  if(!modelSettings_->getGenerateSeismic()) 
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

  if(shiftGrids_ != NULL)
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(shiftGrids_[i] != NULL)
        delete shiftGrids_[i];
  }

  if(gainGrids_ != NULL)
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(gainGrids_[i] != NULL)
        delete gainGrids_[i];
  }

  if (waveletEstimInterval_ != NULL) 
  {
    if (waveletEstimInterval_[0] != NULL)
      delete waveletEstimInterval_[0];
    if (waveletEstimInterval_[1] != NULL)
      delete waveletEstimInterval_[1];
    delete [] waveletEstimInterval_;
  }

  if (faciesEstimInterval_ != NULL) 
  {
    if (faciesEstimInterval_[0] != NULL)
      delete faciesEstimInterval_[0];
    if (faciesEstimInterval_[1] != NULL)
      delete faciesEstimInterval_[1];
    delete [] faciesEstimInterval_;
  }

  for(int i = 0;i<modelSettings_->getNumberOfAngles();i++)
    delete [] reflectionMatrix_[i] ;
  delete [] reflectionMatrix_ ;

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
                  char * errText)
{
  float * tmpRes = new float[n1*n2+1];
  FILE * inFile = fopen(fileName.c_str(),"r");
  std::string text = "Reading "+readReason+" from file "+fileName+" .... \n";
  LogKit::LogFormatted(LogKit::LOW,text);
  char storage[MAX_STRING];
  int index = 0;
  int error = 0;
  while(error == 0 && fscanf(inFile,"%s",storage) != EOF) {
    if(index < n1*n2) {
      if(isNumber(storage) == 1) {
        tmpRes[index] = float(atof(storage));
      }
      else {
        sprintf(errText,"Found '%s' in file %s, expected a number.\n", storage, fileName.c_str());
        error = 1;
      }
    }
    index++;
  }
  if(error == 0) {
    if(index != n1*n2) {
      error = 1;
      sprintf(errText,"Found %d elements in file %s, expected %d.\n",
              index, fileName.c_str(),n1*n2);
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
Model::checkAvailableMemory(Simbox            * timeSimbox,
                            ModelSettings     * modelSettings,
                            const std::string & seismicFile)                           
{
  if (modelSettings->getFileGrid() > 0) // Disk buffering is turn on
    return;

  //
  // Find the size of first seismic volume
  //
  float memOneSeis = 0.0f;
  if (seismicFile != "") {
    memOneSeis = static_cast<float> (NRLib2::FindFileSize(seismicFile));
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
  int gridSize = dummyGrid->getrsize();
  delete dummyGrid;

  int nGrids;
  if((modelSettings->getOutputFlag() & ModelSettings::PREDICTION) == 1)
  {
    nGrids = 10 + modelSettings->getNumberOfAngles();
    if(modelSettings->getNumberOfSimulations() > 0 && nGrids < 13)
      nGrids = 13;
  }
  else
    nGrids = 12;

  int   workSize    = 2500 + int( 0.65*gridSize ); //Size of memory used beyond grids.

  float memOneGrid  = 4.0f * static_cast<float>(gridSize);  
  float mem0        = 4.0f * workSize;
  float mem1        = static_cast<float>(nGrids)*memOneGrid;
  float mem2        = static_cast<float>(modelSettings->getNumberOfAngles())*memOneGrid + memOneSeis;
 
  float neededMem   = mem0 + std::max(mem1, mem2);

  float megaBytes   = neededMem/(1024.f*1024.f);
  float gigaBytes   = megaBytes/1024.f;

  LogKit::LogFormatted(LogKit::HIGH,"\nMemory needed for reading seismic data   : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::HIGH,"Memory needed for holding internal grids : %10.2f MB\n",mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::HIGH,"Memory needed for holding other entities : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (gigaBytes < 0.01f)
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f megaBytes\n",megaBytes);
  else
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f gigaBytes\n",gigaBytes);


  //
  // Check if we can hold everything in memory.
  //
  char   * memchunk0 = new char[workSize];
  float ** memchunk  = new float*[nGrids];

  for(int i = 0 ; i < nGrids ; i++)
    memchunk[i] = new float[gridSize];

  if(memchunk[nGrids-1] == NULL)  //Could not allocate memory
  {
    modelSettings->setFileGrid(1);
    LogKit::LogFormatted(LogKit::LOW,"Not enough memory to hold all grids. Using file storage.\n");
  }
  else
  {
    modelSettings->setFileGrid(0);
  }

  for(int i=0 ; i<nGrids ; i++)
    if(memchunk[i] != NULL) delete [] memchunk[i];
  if(memchunk != NULL) delete [] memchunk;

  if(memchunk0 != NULL) delete[] memchunk0;
}

int
Model::readSegyFile(const std::string   & fileName, 
                    FFTGrid            *& target, 
                    Simbox             *& timeSimbox, 
                    ModelSettings      *& modelSettings, 
                    char                * errText, 
                    const SegyGeometry *& geometry,
                    int                   gridType,
                    int                   i)
{
  strcpy(errText, "");
  SegY * segy = NULL;
  int error = 0;
  target = NULL;

  try
  {
    //
    // Currently we have only one optional TraceHeaderFormat, but this can 
    // be augmented to a list with several formats ...
    //
    std::vector<TraceHeaderFormat*> traceHeaderFormats(0);
    if (modelSettings->getTraceHeaderFormat() != NULL)
    {
      traceHeaderFormats.push_back(modelSettings->getTraceHeaderFormat());
    }


    segy = new SegY(fileName, 
                    modelSettings->getSegyOffset(), 
                    traceHeaderFormats, 
                    true); // Add standard formats to format search


    bool onlyVolume = modelSettings->getAreaParameters() != NULL; // This is now always true
    segy->ReadAllTraces(timeSimbox, 
                        modelSettings->getZpad(),
                        onlyVolume);
    segy->CreateRegularGrid();
  }
  catch (NRLib2::Exception & e)
  {
    sprintf(errText,"%s%s",errText,e.what());
    error++;
  }
  
  if (error == 0)
  {
    const SegyGeometry * geo;
    geo = segy->GetGeometry();
    geo->WriteGeometry();
    if (gridType == FFTGrid::DATA) 
      geometry = new SegyGeometry(geo);
    
    error = timeSimbox->insideRectangle(geo);
    if(error == 0)
    {
      if(modelSettings->getFileGrid() == 1)
        target = new FFTFileGrid(timeSimbox->getnx(),
                                 timeSimbox->getny(), 
                                 timeSimbox->getnz(),
                                 modelSettings->getNXpad(), 
                                 modelSettings->getNYpad(), 
                                 modelSettings->getNZpad());
      else
        target = new FFTGrid(timeSimbox->getnx(), 
                             timeSimbox->getny(), 
                             timeSimbox->getnz(),
                             modelSettings->getNXpad(), 
                             modelSettings->getNYpad(), 
                             modelSettings->getNZpad());
      target->setType(gridType);
      target->fillInFromSegY(segy, timeSimbox);
      if (gridType == FFTGrid::DATA)
        target->setAngle(modelSettings->getAngle(i));
    }
    else
      sprintf(errText,"%sSpecified area in command AREA is larger than the data from SegY file %s\n",
              errText,fileName.c_str());
  }
  if (segy != NULL)
    delete segy;

  return(error);
}

//NBNB Following routine only to be used for parameters!
int
Model::readStormFile(const std::string  & fName, 
                     FFTGrid           *& target, 
                     const char         * parName, 
                     Simbox             * timeSimbox, 
                     ModelSettings     *& modelSettings, 
                     char               * errText)
{
  int error = 0;
  StormContGrid *stormgrid = NULL;
  
  try
  {   
    stormgrid = new StormContGrid(0,0,0);
    stormgrid->ReadFromFile(fName);
    stormgrid->SetMissingCode(RMISSING);
  }
  catch (NRLib2::Exception & e) 
  {
    sprintf(errText,"%s%s",errText,e.what());
    error = 1;
  }

  if(error==0)
  {
    if(modelSettings->getFileGrid() == 1)
      target = new FFTFileGrid(timeSimbox->getnx(), 
                               timeSimbox->getny(), 
                               timeSimbox->getnz(), 
                               modelSettings->getNXpad(), 
                               modelSettings->getNYpad(), 
                               modelSettings->getNZpad());
    else
      target = new FFTGrid(timeSimbox->getnx(),
                           timeSimbox->getny(), 
                           timeSimbox->getnz(), 
                           modelSettings->getNXpad(), 
                           modelSettings->getNYpad(), 
                           modelSettings->getNZpad());
    target->setType(FFTGrid::PARAMETER);
    target->fillInFromStorm(timeSimbox,stormgrid, parName);
  }  

  if (stormgrid != NULL)
    delete stormgrid;
  
  return(error);
}

int 
Model::setPaddingSize(int nx, float px)
{
  int leastint    =  static_cast<int>(ceil(nx*(1.0f+px)));
  int closestprod =  findClosestFactorableNumber(leastint);
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
                        char           * errText,
                        bool           & failed)
{
  //
  // Set AREA (Use SegY geometry from first seismic data volume if needed).
  //
  std::string areaType = "Model file";
  std::string seismicFile("");
  if (modelSettings_->getGenerateSeismic()==false && inputFiles->getNumberOfSeismicFiles() > 0)
    seismicFile = inputFiles->getSeismicFile(0); // Also needed for checkAvailableMemory()

  if (modelSettings->getAreaParameters() == NULL)
  {
    LogKit::LogFormatted(LogKit::HIGH,"\nFinding inversion area from seismic data in file %s\n", 
                         seismicFile.c_str());
    areaType = "Seismic data";
    SegyGeometry * geometry = SegY::FindGridGeometry(seismicFile, modelSettings->getTraceHeaderFormat());
    modelSettings->setAreaParameters(geometry);
    delete geometry;
  }
  const SegyGeometry * areaParams = modelSettings->getAreaParameters(); 

  int error = timeSimbox->setArea(areaParams, errText);
  if(error==1)
  {
  writeAreas(areaParams,timeSimbox,areaType);
  sprintf(errText,"%s The specified AREA extends outside the surface(s).\n",errText);
  failed = true;
  }
  // Rotate variograms relative to simbox
  modelSettings_->rotateVariograms(static_cast<float> (timeSimbox_->getAngle()));

  //
  // Set SURFACES
  //
  int outputFormat = modelSettings->getFormatFlag();
  int outputFlag   = modelSettings->getOutputFlag();

  error = 0;
  setSimboxSurfaces(timeSimbox, 
                    inputFiles->getTimeSurfFiles(), 
                    modelSettings->getParallelTimeSurfaces(), 
                    modelSettings->getTimeDTop(), 
                    modelSettings->getTimeLz(), 
                    modelSettings->getTimeDz(), 
                    modelSettings->getTimeNz(),
                    outputFormat,
                    outputFlag,
                    errText,
                    error);


  if(error == 0)
  {
    if(modelSettings->getUseLocalWavelet() && timeSimbox->getIsConstantThick())
      LogKit::LogFormatted(LogKit::WARNING,"\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");

    estimateZPaddingSize(timeSimbox, modelSettings);   

    sprintf(errText,"%c",'\0');
    error = timeSimbox->checkError(modelSettings->getLzLimit(),errText);

    if(error == 0)
    {
      double zmin, zmax;
      timeSimbox->getMinMaxZ(zmin,zmax);
      LogKit::LogFormatted(LogKit::LOW,"\nTime output interval:\n");
      LogKit::LogFormatted(LogKit::LOW,"  Two-way-time          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
        zmin+timeSimbox->getlz()*timeSimbox->getAvgRelThick()*0.5,
        zmin,zmax); 
      LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
        timeSimbox->getlz()*timeSimbox->getAvgRelThick(),
        timeSimbox->getlz()*timeSimbox->getMinRelThick(),
        timeSimbox->getlz());
      LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
        timeSimbox->getdz()*timeSimbox->getAvgRelThick(),
        timeSimbox->getdz(),
        timeSimbox->getdz()*timeSimbox->getMinRelThick());
    }
    else
    {
      sprintf(errText,"%s. Could not make time simulation grid.\n",errText);
      failed = true;
    }

    //
    // Make extended time simbox
    //
    if(inputFiles->getCorrDirFile() != "") {
      //
      // Get correlation direction
      //
      try {
        Surface tmpSurf = NRLib2::ReadStormSurf(inputFiles->getCorrDirFile());
        correlationDirection = new Surface(tmpSurf);
      }
      catch (NRLib2::Exception & e) {
        sprintf(errText,"%s%s",errText,e.what());
        failed = true;
      }
      if(failed == false && modelSettings->getGenerateSeismic() == false) {
        //Extends timeSimbox for correlation coverage. Original stored in timeCutSimbox
        setupExtendedTimeSimbox(timeSimbox, correlationDirection, timeCutSimbox, outputFormat, outputFlag); 
      }      estimateZPaddingSize(timeSimbox, modelSettings);   
      error = timeSimbox->checkError(modelSettings->getLzLimit(),errText);
      if(error == 0)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nTime inversion interval (extended relative to output interval due to correlation):\n");
        double zmin, zmax;
        timeSimbox->getMinMaxZ(zmin,zmax);
        LogKit::LogFormatted(LogKit::LOW,"  Two-way-time          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
          zmin+timeSimbox->getlz()*timeSimbox->getAvgRelThick()*0.5,
          zmin,zmax); 
        LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
          timeSimbox->getlz()*timeSimbox->getAvgRelThick(),
          timeSimbox->getlz()*timeSimbox->getMinRelThick(),
          timeSimbox->getlz());
        LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
          timeSimbox->getdz()*timeSimbox->getAvgRelThick(),
          timeSimbox->getdz(),
          timeSimbox->getdz()*timeSimbox->getMinRelThick());
      }
      else
      {
        sprintf(errText,"%s Could not make the time simulation grid.\n",errText);
        failed = true;
      }
      if(modelSettings->getGenerateSeismic() == false) {
        setupExtendedBackgroundSimbox(timeSimbox, correlationDirection, timeBGSimbox, outputFormat, outputFlag);
        error = timeBGSimbox->checkError(modelSettings->getLzLimit(),errText);
        if(error == 0)
        {
          LogKit::LogFormatted(LogKit::LOW,"\nTime interval used for background modelling:\n");
          double zmin, zmax;
          timeBGSimbox->getMinMaxZ(zmin,zmax);
          LogKit::LogFormatted(LogKit::LOW,"  Two-way-time          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
            zmin+timeBGSimbox->getlz()*timeBGSimbox->getAvgRelThick()*0.5,
            zmin,zmax); 
          LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
            timeBGSimbox->getlz()*timeBGSimbox->getAvgRelThick(),
            timeBGSimbox->getlz()*timeBGSimbox->getMinRelThick(),
            timeBGSimbox->getlz());
          LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
            timeBGSimbox->getdz()*timeBGSimbox->getAvgRelThick(),
            timeBGSimbox->getdz(),
            timeBGSimbox->getdz()*timeBGSimbox->getMinRelThick());
        }
        else
        {
          sprintf(errText,"%s Could not make the grid for background model.\n",errText);
          failed = true;
        }
      }
    }

    if(failed == false) {
      estimateXYPaddingSizes(timeSimbox, modelSettings);
      //
      // Check if CRAVA has enough memory to run calculation without buffering to disk
      //
      checkAvailableMemory(timeSimbox, modelSettings, seismicFile); 
    }

    //
    // Make time simbox with constant thicknesses (needed for log filtering and facies probabilities)
    //
    timeSimboxConstThick = new Simbox(timeSimbox);
    Surface * tsurf = new Surface(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface()));
    timeSimboxConstThick->setDepth(tsurf, 0, timeSimbox->getlz(), timeSimbox->getdz());

    if((outputFlag & ModelSettings::EXTRA_SURFACES) > 0 && (outputFlag & ModelSettings::NOTIME) == 0)
      timeSimboxConstThick->writeTopBotGrids("Surface_Top_Time_ConstThick", 
      "Surface_Base_Time_ConstThick", 
      outputFormat);
  }
  else
  {
    timeSimbox->externalFailure();
    failed = true;
  }
}

void 
Model::setSimboxSurfaces(Simbox                        *& simbox, 
                         const std::vector<std::string> & surfFile, 
                         bool                             parallelSurfaces, 
                         double                           dTop,
                         double                           lz, 
                         double                           dz, 
                         int                              nz,
                         int                              outputFormat,
                         int                              outputFlag,
                         char                           * errText,
                         int                            & error)
{
  const std::string & topName  = surfFile[0]; 
  const std::string & baseName = surfFile[1]; 

  Surface * z0Grid = NULL;
  Surface * z1Grid = NULL;
  try {
    if (NRLib2::IsNumber(topName)) {
      // Find the smallest surface that covers the simbox. For simplicity 
      // we use only four nodes (nx=ny=2).
      double xMin, xMax;
      double yMin, yMax;
      findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                  simbox->getlx(), simbox->getly(),
                                  simbox->getAngle(), 
                                  xMin,yMin,xMax,yMax);
      z0Grid = new Surface(xMin, yMin, xMax-xMin, yMax-yMin, 2, 2, atof(topName.c_str()));
    } 
    else {
      Surface tmpSurf = NRLib2::ReadStormSurf(topName);
      z0Grid = new Surface(tmpSurf);
    }
  }
  catch (NRLib2::Exception & e) {
    sprintf(errText,"%s%s",errText,e.what());
    error = 1;
  }

  if(error == 0) {
    if(parallelSurfaces) { //Only one reference surface
      simbox->setDepth(z0Grid, dTop, lz, dz);
    }
    else {
      try {
        if (NRLib2::IsNumber(baseName)) {
          // Find the smallest surface that covers the simbox. For simplicity 
          // we use only four nodes (nx=ny=2).
          double xMin, xMax;
          double yMin, yMax;
          findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                      simbox->getlx(), simbox->getly(),
                                      simbox->getAngle(), 
                                      xMin,yMin,xMax,yMax);
          z1Grid = new Surface(xMin, yMin, xMax-xMin, yMax-yMin, 2, 2, atof(baseName.c_str()));
        }
        else {
          Surface tmpSurf = NRLib2::ReadStormSurf(baseName);
          z1Grid = new Surface(tmpSurf);
        }
      }
      catch (NRLib2::Exception & e) {
        sprintf(errText,"%s%s",errText,e.what());
        error = 1;
      }
      if(error == 0) {
        try {
          simbox->setDepth(z0Grid, z1Grid, nz);
        }
        catch (NRLib2::Exception & e) {
          sprintf(errText,"%s%s",errText,e.what());
          std::string text(std::string("Seismic data"));
          writeAreas(modelSettings_->getAreaParameters(),simbox,text);
          error = 1;
        }
      }
    }
    if (error == 0) {
      if((outputFlag & ModelSettings::NOTIME) == 0)
        simbox->writeTopBotGrids("Surface_Top_Time", 
                                 "Surface_Base_Time", 
                                 outputFormat);
    }
  }
}

void
Model::setupExtendedTimeSimbox(Simbox   * timeSimbox, 
                               Surface  * corrSurf, 
                               Simbox  *& timeCutSimbox,
                               int        outputFormat,
                               int        outputFlag)
{
  timeCutSimbox = new Simbox(timeSimbox);
  double * corrPlanePars = findPlane(corrSurf);
  Surface * meanSurf = new Surface(*corrSurf);
  int i;
  for(i=0;i<meanSurf->GetN();i++)
    (*meanSurf)(i) = 0;

  meanSurf->Add(&(timeSimbox->GetTopSurface()));
  meanSurf->Add(&(timeSimbox->GetBotSurface()));
  meanSurf->Multiply(0.5);
  double * refPlanePars = findPlane(meanSurf);
  delete meanSurf;

  for(i=0;i<3;i++)
    refPlanePars[i] -= corrPlanePars[i];
  gradX_ = refPlanePars[1];
  gradY_ = refPlanePars[2];
  Surface * refPlane = createPlaneSurface(refPlanePars, corrSurf);

  writeSurfaceToFile(refPlane,"CorrelationRotationPlane",outputFormat);

  refPlane->Add(corrSurf);
  delete [] corrPlanePars;
  delete [] refPlanePars;

  Surface * topSurf = new Surface(*refPlane);
  topSurf->Subtract(&(timeSimbox->GetTopSurface()));
  double shiftTop = topSurf->Max();
  shiftTop *= -1.0;
  topSurf->Add(shiftTop);
  topSurf->Add(&(timeSimbox->GetTopSurface()));

  Surface * botSurf = new Surface(*refPlane);
  botSurf->Subtract(&(timeSimbox->GetBotSurface()));
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
  botSurf->Add(&(timeSimbox->GetBotSurface()));

  timeSimbox->setDepth(topSurf, botSurf, nz);
  
  if((outputFlag & ModelSettings::EXTRA_SURFACES) > 0 && (outputFlag & ModelSettings::NOTIME) == 0)
    timeSimbox->writeTopBotGrids("Surface_Top_Time_Extended", 
                                 "Surface_Base_Time_Extended", 
                                 outputFormat);

  delete refPlane;
}

void
Model::setupExtendedBackgroundSimbox(Simbox   * timeSimbox, 
                                     Surface  * corrSurf, 
                                     Simbox  *& timeBGSimbox,
                                     int        outputFormat,
                                     int        outputFlag)
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
  dTop->Subtract(&(timeSimbox->GetTopSurface()));
  dTop->Multiply(-1.0);
  double shiftTop = dTop->Min();
  delete dTop;
  Surface * topSurf = new Surface(*tmpSurf);
  topSurf->Add(shiftTop);

  //
  // Find base surface of background simbox
  //
  Surface * dBot = new Surface(*tmpSurf);
  dBot->Subtract(&(timeSimbox->GetBotSurface()));
  dBot->Multiply(-1.0);
  double shiftBot = dBot->Max();
  delete dBot;
  Surface * botSurf = new Surface(*tmpSurf);
  botSurf->Add(shiftBot);

  //
  // Calculate number of layers of background simbox 
  //
  tmpSurf->Assign(0.0);
  tmpSurf->Add(botSurf);
  tmpSurf->Subtract(topSurf);
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

  if((outputFlag & ModelSettings::EXTRA_SURFACES) > 0 && (outputFlag & ModelSettings::NOTIME) == 0)
    timeBGSimbox->writeTopBotGrids("Surface_Top_Time_BG", 
                                   "Surface_Base_Time_BG", 
                                   outputFormat);
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
  for(i=0;i<surf->GetN();i++) {
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
  for(i=0;i<result->GetN();i++) {
    result->GetXY(i,x,y);
    (*result)(i) = planeParams[0]+planeParams[1]*x+planeParams[2]*y;
  }
  return(result);
}


void
Model::estimateXYPaddingSizes(Simbox         * timeSimbox, 
                              ModelSettings *& modelSettings)
{
  bool newPaddings = false;
  float xPad    = 0.0f;
  float yPad    = 0.0f;
  float xPadFac = 0.0f;
  float yPadFac = 0.0f;
  if (modelSettings->getXpad() == 0.0f && modelSettings->getYpad() == 0.0f)
  {
    float range1  = modelSettings->getLateralCorr()->getRange();
    float range2  = modelSettings->getLateralCorr()->getSubRange();
    float angle   = modelSettings->getLateralCorr()->getAngle();
    float factor  = 1.54f;  // At dist = 1.54*range, an exponential variogram is reduced to 0.01

    xPad          = factor * MAXIM(fabs(range1*cos(angle)),fabs(range2*sin(angle)));
    yPad          = factor * MAXIM(fabs(range1*sin(angle)),fabs(range2*cos(angle)));

    xPadFac       = static_cast<float>(MINIM(1.0f, xPad / timeSimbox->getlx())); // A padding of more than 100% is insensible
    yPadFac       = static_cast<float>(MINIM(1.0f, yPad / timeSimbox->getly()));

    modelSettings->setXpad(xPadFac);
    modelSettings->setYpad(yPadFac);
    newPaddings = true;
  }
  int nxPad = setPaddingSize(timeSimbox->getnx(), modelSettings->getXpad());
  int nyPad = setPaddingSize(timeSimbox->getny(), modelSettings->getYpad());
  modelSettings->setNXpad(nxPad);
  modelSettings->setNYpad(nyPad);

  if (newPaddings)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nPadding sizes estimated from lateral correlation ranges:\n");
    LogKit::LogFormatted(LogKit::LOW,"  xPad, xPadFac, nxPad                     : %8.f, %6.2f, %6d\n", xPad, xPadFac, nxPad);
    LogKit::LogFormatted(LogKit::LOW,"  yPad, yPadFac, nyPad                     : %8.f, %6.2f, %6d\n", yPad, yPadFac, nyPad);
  }
}

void
Model::estimateZPaddingSize(Simbox         * timeSimbox,
                            ModelSettings *& modelSettings)
{
  bool newPadding = false;
  float zPad    = 0.0f;
  float zPadFac = 0.0f;
  if (modelSettings->getZpad() == 0.0f)
  {
    float wLength = 300.0f;           // Assume a wavelet is approx 300ms.
    zPad          = wLength / 2.0f;   // Use one wavelet as padding
    zPadFac       = static_cast<float>(MINIM(1.0f,zPad / (timeSimbox->getlz()*timeSimbox->getMinRelThick())));
    
    modelSettings->setZpad(zPadFac);
    newPadding = true;
  }
  int nzPad = setPaddingSize(timeSimbox->getnz(), modelSettings->getZpad());
  modelSettings->setNZpad(nzPad);

  if (newPadding)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nPadding sizes estimated from an assumed wavelet length:\n");
    LogKit::LogFormatted(LogKit::LOW,"  zPad, zPadFac, nzPad                     : %8.f, %6.2f, %6d\n", zPad, zPadFac, nzPad);
  }
}

void
Model::processSeismic(FFTGrid      **& seisCube,
                      Simbox        *& timeSimbox,
                      ModelSettings *& modelSettings, 
                      InputFiles     * inputFiles,
                      char           * errText,
                      bool           & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  int error = 0;

  if(inputFiles->getNumberOfSeismicFiles() > 0)
  {
    char tmpErrText[MAX_STRING];
    bool areaInModelFile = modelSettings->getAreaParameters() != NULL;
    int nAngles = modelSettings->getNumberOfAngles();
    const SegyGeometry ** geometry = new const SegyGeometry * [nAngles];
    seisCube = new FFTGrid * [nAngles];

    int readerror = 0;
    for (int i =0 ; i < nAngles ; i++) {
      geometry[i] = NULL;
      
      readerror = readSegyFile(inputFiles->getSeismicFile(i), seisCube[i],
                               timeSimbox, modelSettings, 
                               tmpErrText, geometry[i],
                               FFTGrid::DATA, i);
      error += readerror;
    }

    if(error == 0)
    {
      seisCube[0]->setOutputFormat(modelSettings->getFormatFlag()); //static, controls all grids.
      if(modelSettings->getDebugFlag() == 1)
      {
        char sName[100];
        for(int i=0 ; i<nAngles ; i++)
        {
          sprintf(sName, "origSeis%d",i);
          std::string sgriLabel("Original seismic data for incidence angle ");
          sgriLabel += NRLib2::ToString(modelSettings->getAngle(i));
          seisCube[i]->writeFile(sName, timeSimbox, sgriLabel);
        }
      }

      LogKit::LogFormatted(LogKit::LOW,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::LOW,"-------------------------------------------------------------------------------------------------\n");
      if (areaInModelFile) {
        double azimuth = (-1)*timeSimbox->getAngle()*(180.0/M_PI);
        if (azimuth < 0)
          azimuth += 360.0;
        LogKit::LogFormatted(LogKit::LOW,"Model file       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n", 
                             timeSimbox->getx0(), timeSimbox->gety0(), 
                             timeSimbox->getlx(), timeSimbox->getly(), azimuth, 
                             timeSimbox->getdx(), timeSimbox->getdy());
      }
      for (int i = 0 ; i < nAngles ; i++)
        if (geometry[i] != NULL) {
          double geoAngle = (-1)*timeSimbox->getAngle()*(180/M_PI);
          if (geoAngle < 0)
            geoAngle += 360.0f;
          LogKit::LogFormatted(LogKit::LOW,"Seismic data %d   %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",i,
                               geometry[i]->GetX0(), geometry[i]->GetY0(), 
                               geometry[i]->Getlx(), geometry[i]->Getly(), geoAngle,
                               geometry[i]->GetDx(), geometry[i]->GetDy());
          }
      LogKit::LogFormatted(LogKit::LOW,"\nTime simulation grids:\n");
      LogKit::LogFormatted(LogKit::LOW,"  Output grid         %4i * %4i * %4i   : %10i\n",
                           timeSimbox->getnx(),timeSimbox->getny(),timeSimbox->getnz(),
                           timeSimbox->getnx()*timeSimbox->getny()*timeSimbox->getnz()); 
      LogKit::LogFormatted(LogKit::LOW,"  FFT grid            %4i * %4i * %4i   : %10i\n",
                           modelSettings->getNXpad(),modelSettings->getNYpad(),modelSettings->getNZpad(),
                           modelSettings->getNXpad()*modelSettings->getNYpad()*modelSettings->getNZpad());

      for (int i =0 ; i < nAngles ; i++)
        if (geometry[i] != NULL)
          delete geometry[i];
      delete [] geometry;
    }
    else
    {
      sprintf(errText, "%sReading of seismic data failed:\n %s\n", errText, tmpErrText);
    }
  }
  failed = error > 0;
  Timings::setTimeSeismic(wall,cpu);
}

void 
Model::processWells(WellData     **& wells,
                    Simbox         * timeSimbox,
                    Simbox         * timeBGSimbox,
                    Simbox         * timeSimboxConstThick,
                    FFTGrid       ** seisCube,
                    RandomGen      * randomGen,
                    ModelSettings *& modelSettings, 
                    InputFiles     * inputFiles,
                    char           * errText,
                    bool           & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
  LogKit::LogFormatted(LogKit::LOW,"\nReading well data:\n");

  bool    faciesLogGiven = modelSettings->getFaciesLogGiven();
  int     nWells         = modelSettings->getNumberOfWells();
  int     nFacies        = 0;

  int error = 0;

  char tmpErrText[MAX_STRING];
  sprintf(tmpErrText,"%c",'\0');
  wells = new WellData *[nWells];
  for(int i=0 ; i<nWells ; i++) {
    wells[i] = new WellData(inputFiles->getWellFile(i), 
                            inputFiles->getLogNames(),
                            modelSettings, 
                            modelSettings->getIndicatorFacies(i),
                            modelSettings->getIndicatorWavelet(i),
                            modelSettings->getIndicatorBGTrend(i),
                            faciesLogGiven);
    if(wells[i]->checkError(tmpErrText) != 0) {
      sprintf(errText,"%s%s", errText, tmpErrText);
      error = 1;
    }
  }

  if (error == 0) {
    if(modelSettings->getFaciesLogGiven()) { 
      checkFaciesNames(wells, modelSettings, tmpErrText, error);
      nFacies = modelSettings->getNumberOfFacies(); // nFacies is set in checkFaciesNames()
    }

    Utils::writeHeader("Processing Wells");
    
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
      LogKit::LogFormatted(LogKit::LOW,"%s : \n",wells[i]->getWellname());
      if(wells[i]!=NULL) {
        if(wells[i]->checkSimbox(timeSimbox) == 1) {
          skip = true;
          nohit++;
        }
        if(wells[i]->getNd() == 0) {
          LogKit::LogFormatted(LogKit::LOW,"  IGNORED (no log entries found)\n");
          skip = true;
          empty++;
        }
        if(wells[i]->isFaciesOk()==0) {
          LogKit::LogFormatted(LogKit::LOW,"   IGNORED (facies log has wrong entries)\n");
          skip = true;
          facieslognotok++;
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
          int nAngles = modelSettings->getNumberOfAngles();
          for (int iAngle = 0 ; iAngle < nAngles ; iAngle++) 
            wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
          if (timeBGSimbox==NULL)
            wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeSimbox, randomGen) ); // Need a copy constructor?
          else
            wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeBGSimbox, randomGen) );
          if (nFacies > 0)
            wells[i]->countFacies(timeSimbox,faciesCount[i]);
          if((modelSettings->getOutputFlag() & ModelSettings::WELLS) > 0) 
            wells[i]->writeRMSWell();
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
                             wells[i]->getWellname(),
                             nMerges[i],
                             nInvalidAlpha[i], 
                             nInvalidBeta[i], 
                             nInvalidRho[i],
                             (rankCorr[i] > modelSettings->getMaxRankCorr() ? "yes" : " no"),
                             rankCorr[i],
                             (devAngle[i] > modelSettings->getMaxDevAngle() ? "yes" : " no"),
                             devAngle[i]);
      else  
        LogKit::LogFormatted(LogKit::LOW,"%-23s      -       -    -    -       - /     -       -  /    -\n",
                             wells[i]->getWellname());
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
          LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[i]->getWellname());
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
          LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[i]->getWellname());
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
          LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[i]->getWellname());
          for (int f = 0 ; f < nFacies ; f++) {
            LogKit::LogFormatted(LogKit::MEDIUM,"%12d ",faciesCount[i][f]);
          }
          LogKit::LogFormatted(LogKit::MEDIUM,"\n");
        } 
        else {
          LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[i]->getWellname());
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
    if (nWells==0) {
      LogKit::LogFormatted(LogKit::LOW,"\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
      LogKit::LogFormatted(LogKit::LOW,"\n       below is correct. If it is not, you probably have problems with coordinate scaling. To");
      LogKit::LogFormatted(LogKit::LOW,"\n       correct this problem the seismic data must be re-exported from your data base using");
      LogKit::LogFormatted(LogKit::LOW,"\n       the \'by-pass coordinate scaling (byte 71-72)\' toggle that matches your version of");
      LogKit::LogFormatted(LogKit::LOW,"\n       CRAVA (see beginning of CRAVA log file). Alternatively, CRAVA can be recompiled.\n");
      LogKit::LogFormatted(LogKit::LOW,"\n                                   X0          Y0        DeltaX      DeltaY      Angle");
      LogKit::LogFormatted(LogKit::LOW,"\n       -------------------------------------------------------------------------------");
      LogKit::LogFormatted(LogKit::LOW,"\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n", 
                           timeSimbox->getx0(), timeSimbox->gety0(), 
                           timeSimbox->getlx(), timeSimbox->getly(), 
                           (timeSimbox->getAngle()*180)/PI);
      error = 1;
    }
  }
  failed = error > 0;
  Timings::setTimeWells(wall,cpu);
}

void Model::checkFaciesNames(WellData      ** wells,
                             ModelSettings *& modelSettings,
                             char           * tmpErrText,
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
          sprintf(tmpErrText,"Problem with facies logs. Facies names and numbers are not uniquely defined.\n");
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
}

void 
Model::processBackground(Background   *& background,
                         WellData     ** wells,
                         Simbox        * timeSimbox,
                         Simbox        * timeBGSimbox,
                         ModelSettings * modelSettings, 
                         InputFiles    * inputFiles,
                         char          * errText,
                         bool          & failed)
{
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::BACKGROUND) > 0 || 
      (modelSettings->getOutputFlag() & ModelSettings::WAVELETS)   > 0 )
  {
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
        loadVelocity(velocity, timeSimbox, modelSettings_, 
                     backVelFile, errText, failed);
      if (!failed) 
      {
        if(modelSettings->getBackgroundVario() == NULL)
        {
          sprintf(errText,"%sThere is no variogram available for the background modelling\n",errText);
          failed = true;
        }
        for (int i=0 ; i<3 ; i++)
        {
          backModel[i] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);              
          backModel[i]->setType(FFTGrid::PARAMETER);
        }
        background = new Background(backModel, wells, velocity, timeSimbox, timeBGSimbox, modelSettings);
      }
      if(velocity != NULL)
        delete velocity;
    }
    else 
    {
      const char * parName[]={"Vp background","Vs background","Rho background"};
      for(int i=0 ; i<3 ; i++)
      {
        float constBackValue = modelSettings->getConstBackValue(i);

        const std::string & backFile = inputFiles->getBackFile(i);

        if(constBackValue < 0)
        {
          if(backFile.size() > 0)
          {
            char tmpErrText[MAX_STRING];
            sprintf(tmpErrText,"%c",'\0');
            int readerror = 0;
            if(findFileType(backFile) == SEGYFILE) {
              const SegyGeometry * geometry = NULL;
              readerror = readSegyFile(backFile, backModel[i], 
                                       timeSimbox, modelSettings, 
                                       tmpErrText, geometry, 
                                       FFTGrid::PARAMETER);
            }
            else
              readerror = readStormFile(backFile, backModel[i], parName[i], 
                                        timeSimbox, modelSettings,
                                        tmpErrText);
            backModel[i]->logTransf();
            
            if(readerror != 0)
            {
              sprintf(errText,"%sReading of file \'%s\' for parameter \'%s\' failed\n",
                      errText,backFile.c_str(),parName[i]);
              failed = true;
            }
          }
          else
          {
            sprintf(errText,"%sReading of file for parameter \'%s\' failed. No file name is given.\n",
                    errText,parName[i]);
            failed = true;
          }
        }
        else if(constBackValue > 0)
        {
          if(modelSettings->getFileGrid() == 1)
            backModel[i] = new FFTFileGrid(nx, ny, nz, nxPad, nyPad, nzPad);
          else
            backModel[i] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);              
          backModel[i]->setType(FFTGrid::PARAMETER);
          backModel[i]->fillInConstant(float( log( constBackValue )));
        }
        else
        {
          sprintf(errText,"%sTrying to set background model to 0 for parameter %s\n\n",
                  errText,parName[i]);
          failed = true;
        }
      }
      background = new Background(backModel);
    }
    if((modelSettings->getOutputFlag() & ModelSettings::BACKGROUND) > 0)
      background->writeBackgrounds(timeSimbox, timeDepthMapping_, timeCutMapping_); 
    Timings::setTimePriorExpectation(wall,cpu);
  }
}

void 
Model::processPriorCorrelations(Corr         *& correlations,
                                Background    * background,
                                WellData     ** wells,
                                Simbox        * timeSimbox,
                                ModelSettings * modelSettings, 
                                InputFiles    * inputFiles,
                                char          * errText,
                                bool          & failed)
{
  bool printResult = (modelSettings->getOutputFlag() & (ModelSettings::PRIORCORRELATIONS + ModelSettings::CORRELATION)) > 0;
  if (modelSettings->getDoInversion() || printResult)
  {
    Utils::writeHeader("Prior Covariance");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);
    //
    // Parameter correlation can be set in model file.
    // Default NULL, indicating that estimate will be used.
    //
    float ** paramCorr = NULL;
    const std::string & paramCorrFile = inputFiles->getParamCorrFile();

    if(paramCorrFile != "") 
    {
      char tmpErrText[MAX_STRING];
      sprintf(tmpErrText,"%c",'\0');
      paramCorr = readMatrix(paramCorrFile, 3, 3, "parameter correlation", tmpErrText);
      if(paramCorr == NULL) 
      {
        sprintf(errText,"%sReading of file \'%s\' for parameter correlation matrix failed\n%s\n",
                errText,paramCorrFile.c_str(),tmpErrText);
        failed = true;
      }
      LogKit::LogFormatted(LogKit::LOW,"Parameter correlation read from file.\n\n");
    }

    Surface * CorrXY = findCorrXYGrid();

    if(modelSettings->getLateralCorr()==NULL) { // NBNB-PAL: this will never be true (default lateral corr)
      estimateCorrXYFromSeismic(CorrXY,
                                seisCube_,
                                modelSettings->getNumberOfAngles());
    }

    Analyzelog * analyze = new Analyzelog(wells, 
                                          background,
                                          timeSimbox, 
                                          modelSettings);

    correlations = new Corr(analyze->getPointVar0(), 
                            analyze->getVar0(), 
                            analyze->getCorrT(), 
                            analyze->getNumberOfLags(),
                            static_cast<float>(timeSimbox->getdz()), 
                            CorrXY);
    delete analyze;

    if(correlations == NULL)
    {
      sprintf(errText,"%sCould not construct prior covariance. Unknown why...\n",errText);
      failed = true;
    }

    //
    // NBNB-PAL: 
    //
    // CorrT   (PriorCorrT)
    // CorrXY  (PriorCorrXY)
    // Var0    (PriorVar0)
    //
    // may be read from file.... currently only Var0 can be read.
    //
    //
    if(paramCorr != NULL)
      correlations->setPriorVar0(paramCorr);

    if(printResult)
      correlations->writeFilePriorVariances();
    
    correlations->printPriorVariances();

    Timings::setTimePriorCorrelation(wall,cpu);
  }
}

Surface * 
Model::findCorrXYGrid(void)
{
  int npix, nx, ny;
  int i,j,refi,refj;

  float dx  = static_cast<float>(timeSimbox_->getdx());
  float dy  = static_cast<float>(timeSimbox_->getdy());
  float snx = static_cast<float>(timeSimbox_->getnx());
  float sny = static_cast<float>(timeSimbox_->getny());

  nx = findClosestFactorableNumber(static_cast<int>(ceil(snx*(1.0f+modelSettings_->getXpad())))); //Use padded grid
  ny = findClosestFactorableNumber(static_cast<int>(ceil(sny*(1.0f+modelSettings_->getYpad()))));
  npix = nx*ny;
  Surface * grid = new Surface(0, 0, dx*nx, dy*ny, nx, ny, RMISSING);

  if(modelSettings_->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
  {
    for(j=0;j<ny;j++)
    {
      for(i=0;i<nx;i++)
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
        (*grid)(j*nx+i) = modelSettings_->getLateralCorr()->corr(refi*dx, refj*dy);
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
Model::processReflectionMatrix(float       **& reflectionMatrix,
                               Background    * background,
                               ModelSettings * modelSettings, 
                               InputFiles    * inputFiles,
                               char          * errText,
                               bool          & failed)
{
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model.
  //
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::WAVELETS)   > 0 )
  {
    const std::string & reflMatrFile = inputFiles->getReflMatrFile();

    if(reflMatrFile != "") {
      reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", errText);
      if(reflectionMatrix == NULL) {
        sprintf(errText,"%sReading of file \'%s\' for reflection matrix failed\n",errText,reflMatrFile.c_str());
        failed = true;
      }
      LogKit::LogFormatted(LogKit::LOW,"Reflection parameters read from file.\n\n");
    }
    else {
      if (background != NULL) {
        setupDefaultReflectionMatrix(reflectionMatrix,
                                     background,
                                     modelSettings);
      }
      else {
        sprintf(errText,"%s\nFailed to set up reflection matrix. Background model is empty.\n",errText);
        failed = true;
      }
    }
  }
}

void
Model::setupDefaultReflectionMatrix(float       **& reflectionMatrix,
                                    Background    * background,
                                    ModelSettings * modelSettings)
{
  int i;
  float ** A = new float * [modelSettings->getNumberOfAngles()];
  // For debugging
  //background->setClassicVsVp();
  double vsvp  = background->getMeanVsVp();
  double vsvp2 = vsvp*vsvp;
  for(i = 0; i < modelSettings->getNumberOfAngles(); i++)
  {
    double angle = static_cast<double>(modelSettings->getAngle(i));
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(modelSettings->getSeismicType(i) == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float((1.0 +tan2t )/2.0) ; 
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0);
    }
    else if(modelSettings->getSeismicType(i) == ModelSettings::PSSEIS) {
      double cost  = cos(angle);
      double cosp  = sqrt(1-vsvp2*sint2);
      double fac   = 0.5*sint/cosp;

      A[i][0] = 0;
      A[i][1] = float(4.0*fac*(vsvp2*sint2+vsvp*cost*cosp));
      A[i][2] = float(fac*(-1.0+2*vsvp2*sint2+2*vsvp*cost*cosp));
    }
  }
  reflectionMatrix = A;
  LogKit::LogFormatted(LogKit::LOW,"\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",1.0f/vsvp);
}

void
Model::processWavelets(Wavelet     **& wavelet,
                       FFTGrid      ** seisCube,
                       WellData     ** wells,
                       float        ** reflectionMatrix,
                       Simbox        * timeSimbox,
                       Surface      ** waveletEstimInterval,
                       Surface      ** shiftGrids,
                       Surface      ** gainGrids,
                       ModelSettings * modelSettings, 
                       InputFiles    * inputFiles,
                       char          * errText,
                       bool          & failed)
{
  int error = 0;
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::WAVELETS) > 0 )
  {
    Utils::writeHeader("Processing/generating wavelets");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    bool estimateStuff = false;
    for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {  
      estimateStuff = estimateStuff || modelSettings->getEstimateWavelet(i); 
      estimateStuff = estimateStuff || modelSettings->getEstimateSNRatio(i); 
    }
    if (estimateStuff) 
    {
      LogKit::LogFormatted(LogKit::HIGH,"\n\nWells that cannot be used in wavelet generation or noise estimation:");
      LogKit::LogFormatted(LogKit::HIGH,"\n  Deviated wells.");
      LogKit::LogFormatted(LogKit::HIGH,"\n  Wells with too little data.");
    }
    wavelet    = new Wavelet * [modelSettings->getNumberOfAngles()];
    shiftGrids = new Surface * [modelSettings->getNumberOfAngles()];
    gainGrids  = new Surface * [modelSettings->getNumberOfAngles()];

    for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {  
      float angle = modelSettings->getAngle(i);
      shiftGrids[i] = NULL;
      gainGrids[i]  = NULL;
      LogKit::LogFormatted(LogKit::LOW,"\nAngle stack : %.1f deg",angle*180.0/PI);
      if (modelSettings->getEstimateWavelet(i)) 
      {
        if (timeSimbox->getdz() > 4.01f) { // Require this density for wavelet estimation
          LogKit::LogFormatted(LogKit::LOW,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
          LogKit::LogFormatted(LogKit::LOW,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
          LogKit::LogFormatted(LogKit::LOW,"         the number of layers must be increased.                                    \n");
        }
        wavelet[i] = new Wavelet1D(timeSimbox, 
                                   seisCube[i], 
                                   wells, 
                                   waveletEstimInterval,                                   
                                   modelSettings, 
                                   reflectionMatrix[i],
                                   i);
      }
      else
      {
        const std::string & waveletFile = inputFiles->getWaveletFile(i);

        int fileFormat = getWaveletFileFormat(waveletFile,errText);
        if(fileFormat < 0)
        {
          sprintf(errText, "%s Unknown file format of file  %s.\n", errText, waveletFile.c_str());
          error += 1;
        }
        else {
          if (fileFormat == Wavelet::SGRI)
            wavelet[i] = new Wavelet3D(waveletFile, 
                                       modelSettings, 
                                       timeSimbox, 
                                       angle, 
                                       reflectionMatrix[i],
                                       error, 
                                       errText);
          else {
            wavelet[i] = new Wavelet1D(waveletFile, 
                                       fileFormat, 
                                       modelSettings, 
                                       reflectionMatrix[i],
                                       error, 
                                       errText);
            if (error == 0) {
              //wavelet[i]->write1DWLas3DWL(); //Frode: For debugging and testing
              //wavelet[i]->write3DWLfrom1DWL();
              wavelet[i]->resample(static_cast<float>(timeSimbox->getdz()), timeSimbox->getnz(), 
                                   modelSettings->getZpad(), angle);
            }
          }
        }
      }
      if (error == 0) {
        if (modelSettings->getMatchEnergies(i)) // If true we later scale wavelet to get EmpSN = TheoSN.
          wavelet[i]->scale(modelSettings->getWaveletScale(i));
        if ((wavelet[i]->getDim() == 3) && !timeSimbox->getIsConstantThick()) {
          sprintf(errText, "%s Simbox must have constant thickness when 3D wavelet.\n", errText);
          error += 1;
        }
        if (modelSettings->getEstimateSNRatio(i))
        {
          if (wavelet[i]->getDim() == 3) { //Not possible to estimate signal-to-noise ratio for 3D wavelets
            sprintf(errText, "%s Estimation of signal-to-noise ratio is not possible for 3D wavelets.\n", errText);
            sprintf(errText, "%s The s/n ratio must be specified in the model file\n", errText);
            error += 1;
          }
          else {
            float SNRatio = wavelet[i]->calculateSNRatio(timeSimbox, seisCube[i], wells, 
                                                         modelSettings->getNumberOfWells(), 
                                                         errText, error, shiftGrids[i], gainGrids[i]);
            modelSettings->setSNRatio(i,SNRatio);
          }
        }
        else
        {
          float SNRatio = modelSettings->getSNRatio(i);
          if (SNRatio <= 1.0f || SNRatio > 10.f)
          {
            sprintf(errText, "%s Illegal signal-to-noise ratio of %.3f for cube %d\n", errText, SNRatio,i);
            sprintf(errText, "%s Ratio must be in interval 1.0 < S/N ratio < 10.0\n", errText);
            error += 1;
          }
        }
        if (error == 0) {
          if((modelSettings->getOutputFlag() & ModelSettings::WAVELETS) > 0) 
          {
            char fileName[MAX_STRING];
            sprintf(fileName,"Wavelet_Scaled");
            wavelet[i]->writeWaveletToFile(fileName, 1.0, timeSimbox); // dt_max = 1.0;
          }
          if(shiftGrids != NULL && shiftGrids[i] != NULL)
            wavelet[i]->setShiftGrid(shiftGrids[i], timeSimbox);
          if(gainGrids != NULL && gainGrids[i] != NULL)
            wavelet[i]->setGainGrid(gainGrids[i], timeSimbox);
        }
      }
    }
    Timings::setTimeWavelets(wall,cpu);
  }
  failed = error > 0;
}

int
Model::getWaveletFileFormat(const std::string & fileName, char * errText)
{
  int fileformat = -1;
  char* dummyStr = new char[MAX_STRING];
  // test for old file format
  FILE* file = fopen(fileName.c_str(),"r");
  for(int i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName.c_str());
      return 0;
    } // endif
  }  // end for i
  fclose(file);
  char* targetString =new char[MAX_STRING];
  strcpy(targetString,"CMX");
  int  pos = findEnd(dummyStr, 0, targetString);
  if(pos>=0)
    fileformat= Wavelet::OLD;

  if(fileformat<0) // not old format
  {
    // Test for Sgri format
    file = fopen(fileName.c_str(), "r");
    if (fscanf(file, "%s", dummyStr) == EOF)
    {
      sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName.c_str());
      return 0;
    }
    strcpy(targetString, "NORSAR");
    readToEOL(file);
    pos = findEnd(dummyStr, 0, targetString);
    if (pos >= 0)
      fileformat = Wavelet::SGRI;
    fclose(file);

    if (fileformat != Wavelet::SGRI) {
      // test for jason file format
      file = fopen(fileName.c_str(),"r");
      bool lineIsComment = true; 
      while( lineIsComment ==true)
      {
        if(fscanf(file,"%s",dummyStr) == EOF)
        {
          readToEOL(file);
          sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName.c_str());
          return 0;
        } // endif
        else
        {
          readToEOL(file);
          if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
          {
            lineIsComment = false;
          }
        }
      }  // end while
      fclose(file);
      if(atof(dummyStr)!=0) // not convertable number
      {
        fileformat= Wavelet::JASON;
      }
    }
  }
  delete[] dummyStr;
  delete [] targetString;
  return fileformat;
}

void Model::processPriorFaciesProb(float         *& priorFacies,
                                   WellData      ** wells,
                                   RandomGen      * randomGen,
                                   int              nz,
                                   ModelSettings  * modelSettings)
{
  int outputFlag = modelSettings->getOutputFlag();
  int nFacies    = modelSettings->getNumberOfFacies();
  if(nFacies > 0 && (outputFlag & (ModelSettings::FACIESPROB + ModelSettings::FACIESPROBRELATIVE)) > 0)
  {
    Utils::writeHeader("Prior Facies Probabilities");
    //
    // NBNB-PAL: We should be able to read priorFacies from file. 
    //
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

    for (int w = 0 ; w < nWells ; w++)
    {
      if(wells[w]->getUseForFaciesProbabilities())
      { 
        //
        // Note that we use timeSimbox to calculate prior facies probabilities
        // instead of the simbox with parallel top and base surfaces. This
        // will make the prior probabilities slightly different, but that
        // should not be a problem.
        //
        BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
        bl->getVerticalTrend(bl->getAlpha(),vtAlpha);
        bl->getVerticalTrend(bl->getBeta(),vtBeta);
        bl->getVerticalTrend(bl->getRho(),vtRho);
        bl->getVerticalTrend(bl->getFacies(),vtFacies,randomGen);
        for(int i=0 ; i<nz ; i++)
        {
          if(vtAlpha[i] != RMISSING && vtBeta[i] != RMISSING && vtRho[i] != RMISSING) {
            if (vtFacies[i] != IMISSING)
              faciesCount[w][vtFacies[i]]++;
            faciesLog[w*nz + i] = vtFacies[i];
          }
          else
            faciesLog[w*nz + i] = IMISSING;
        }
      }
    }

    delete [] vtAlpha;
    delete [] vtBeta; 
    delete [] vtRho;  
    delete [] vtFacies; 

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
      if(wells[w]->getUseForFaciesProbabilities())
      { 
        float tot = 0.0;
        for (int i = 0 ; i < nFacies ; i++)
          tot += static_cast<float>(faciesCount[w][i]);
        LogKit::LogFormatted(LogKit::LOW,"%-23s ",wells[w]->getWellname());
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
        LogKit::LogFormatted(LogKit::MEDIUM,"%-23s ",wells[w]->getWellname());
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
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: No valid facies log entries have been found\n");
    }
    delete [] nData;
  }
}

void
Model::loadExtraSurfaces(Surface  **& waveletEstimInterval,
                         Surface  **& faciesEstimInterval,
                         Simbox     * timeSimbox,
                         InputFiles * inputFiles,
                         char       * errText,
                         bool       & failed)
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
    waveletEstimInterval = new Surface*[2];
    try {
      if (NRLib2::IsNumber(topWEI)) 
        waveletEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWEI.c_str()));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(topWEI);
        waveletEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s\n", errText,e.what());
      failed = true;
    }
    
    try {
      if (NRLib2::IsNumber(baseWEI)) 
        waveletEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWEI.c_str()));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(baseWEI);
        waveletEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s\n", errText,e.what());
      failed = true;
    }
  }
  //
  // Get facies estimation interval
  //
  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" && baseFEI != "") {  
    faciesEstimInterval = new Surface*[2];
    try {
      if (NRLib2::IsNumber(topFEI)) 
        faciesEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topFEI.c_str()));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(topFEI);
        faciesEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s\n", errText,e.what());
      failed = true;
    }

    try {
      if (NRLib2::IsNumber(baseFEI)) 
        faciesEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseFEI.c_str()));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(baseFEI);
        faciesEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s\n", errText,e.what());
      failed = true;
    }
  }
}

// NBNB The following routine is stupid, and assumes SEGY if file 
// does not start with 'storm_petro_binary'
int
Model::findFileType(const std::string & fileName)
{
  FILE * file = fopen(fileName.c_str(),"r");
  char header[19];
  fread(header, 1, 18, file);
  fclose(file);
  header[18] = 0;
  if(strcmp(header,"storm_petro_binary") == 0)
    return(STORMFILE);
  else
    return(SEGYFILE);
}

void
Model::printSettings(ModelSettings * modelSettings,
                     InputFiles    * inputFiles)
{
  LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings:\n");
  int logLevel = modelSettings->getLogLevel();
  std::string logText("*NONE*");
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
  LogKit::LogFormatted(LogKit::LOW,"  Log level                                : %10s\n",logText.c_str());
  if (modelSettings->getNumberOfSimulations() == 0)
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
  LogKit::LogFormatted(LogKit::LOW,"  Kriging                                  : %10s\n",(modelSettings->getKrigingParameters()==NULL ? "no" : "yes"));

  LogKit::LogFormatted(LogKit::HIGH,"\nUnit settings/assumptions:\n");
  LogKit::LogFormatted(LogKit::HIGH,"  Time                                     : %10s\n","ms TWT");
  LogKit::LogFormatted(LogKit::HIGH,"  Frequency                                : %10s\n","Hz");
  LogKit::LogFormatted(LogKit::HIGH,"  Length                                   : %10s\n","m");
  LogKit::LogFormatted(LogKit::HIGH,"  Velocities                               : %10s\n","m/s");
  LogKit::LogFormatted(LogKit::HIGH,"  Density                                  : %10s\n","kg/dm3");
  LogKit::LogFormatted(LogKit::HIGH,"  Angles                                   : %10s\n","   degrees (clockwise relative to north)");

  TraceHeaderFormat * thf = modelSettings_->getTraceHeaderFormat(); 
  if (thf != NULL) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nSegY trace header format:\n");
    LogKit::LogFormatted(LogKit::LOW,"  Format name                              : %10s\n",thf->GetFormatName().c_str());
    if (thf->GetBypassCoordScaling())
      LogKit::LogFormatted(LogKit::LOW,"  Bypass coordinate scaling                :        yes\n");
    if (!thf->GetStandardType()) 
    {
      LogKit::LogFormatted(LogKit::LOW,"  Start pos coordinate scaling             : %10d\n",thf->GetScalCoLoc());
      LogKit::LogFormatted(LogKit::LOW,"  Start pos trace x coordinate             : %10d\n",thf->GetUtmxLoc());
      LogKit::LogFormatted(LogKit::LOW,"  Start pos trace y coordinate             : %10d\n",thf->GetUtmyLoc());
      LogKit::LogFormatted(LogKit::LOW,"  Start pos inline index                   : %10d\n",thf->GetInlineLoc());
      LogKit::LogFormatted(LogKit::LOW,"  Start pos crossline index                : %10d\n",thf->GetCrosslineLoc());
      LogKit::LogFormatted(LogKit::LOW,"  Coordinate system                        : %10s\n",thf->GetCoordSys()==0 ? "UTM" : "ILXL" );
    }
  }
  //
  // WELL PROCESSING
  //
  LogKit::LogFormatted(LogKit::LOW,"\nSettings for well processing:\n");
  LogKit::LogFormatted(LogKit::LOW,"  Threshold for merging log entries        : %10.2f ms\n",modelSettings->getMaxMergeDist());
  LogKit::LogFormatted(LogKit::LOW,"  Threshold for Vp-Vs rank correlation     : %10.2f\n",modelSettings->getMaxRankCorr());
  LogKit::LogFormatted(LogKit::LOW,"  Threshold for deviation angle            : %10.1f (=%.2fm/ms TWT)\n",
                       modelSettings->getMaxDevAngle(),tan(modelSettings->getMaxDevAngle()*PI/180.0));
  LogKit::LogFormatted(LogKit::LOW,"  High cut for background modelling        : %10.1f\n",modelSettings->getMaxHzBackground());
  LogKit::LogFormatted(LogKit::LOW,"  High cut for seismic resolution          : %10.1f\n",modelSettings->getMaxHzSeismic());
  LogKit::LogFormatted(LogKit::LOW,"\nRange of allowed parameter values:\n");
  LogKit::LogFormatted(LogKit::LOW,"  Vp  - min                                : %10.0f\n",modelSettings->getAlphaMin());
  LogKit::LogFormatted(LogKit::LOW,"  Vp  - max                                : %10.0f\n",modelSettings->getAlphaMax());
  LogKit::LogFormatted(LogKit::LOW,"  Vs  - min                                : %10.0f\n",modelSettings->getBetaMin());
  LogKit::LogFormatted(LogKit::LOW,"  Vs  - max                                : %10.0f\n",modelSettings->getBetaMax());
  LogKit::LogFormatted(LogKit::LOW,"  Rho - min                                : %10.1f\n",modelSettings->getRhoMin());
  LogKit::LogFormatted(LogKit::LOW,"  Rho - max                                : %10.1f\n",modelSettings->getRhoMax());  

  //
  // WELL DATA
  //
  if (modelSettings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nWell logs:\n");
    const std::vector<std::string> & logNames = inputFiles->getLogNames();

    if (logNames.size() > 0)
    {
      LogKit::LogFormatted(LogKit::LOW,"  Time                                     : %10s\n",  logNames[0].c_str());
      if(NRLib2::Uppercase(logNames[1])=="VP" || 
         NRLib2::Uppercase(logNames[1])=="LFP_VP")
        LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity                          : %10s\n",logNames[1].c_str());
      else
        LogKit::LogFormatted(LogKit::LOW,"  Sonic                                    : %10s\n",logNames[1].c_str());
      if(NRLib2::Uppercase(logNames[3])=="VS" || 
         NRLib2::Uppercase(logNames[3])=="LFP_VS")
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
      LogKit::LogFormatted(LogKit::LOW,"\nUse well in estimation of:                  ");
      if (generateBackground) LogKit::LogFormatted(LogKit::LOW,"  BGTrend");
      if (estimateWavelet)    LogKit::LogFormatted(LogKit::LOW,"  Wavelet");
      if (estimateFaciesProb) LogKit::LogFormatted(LogKit::LOW,"  Facies");
      LogKit::LogFormatted(LogKit::LOW,"\n");
      for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++) 
      {
        LogKit::LogFormatted(LogKit::LOW,"  %-2d                                       :  ",i+1);

        if (generateBackground) LogKit::LogFormatted(LogKit::LOW,"  %-7s",(modelSettings->getIndicatorBGTrend(i) ? "yes" : " no"));
        if (estimateWavelet)    LogKit::LogFormatted(LogKit::LOW,"  %-7s",(modelSettings->getIndicatorWavelet(i) ? "yes" : " no"));
        if (estimateFaciesProb) LogKit::LogFormatted(LogKit::LOW,"  %-7s",(modelSettings->getIndicatorFacies(i)  ? "yes" : " no"));
        LogKit::LogFormatted(LogKit::LOW,"\n");
      }
    }
  }

  //
  // SURFACES
  // 
  LogKit::LogFormatted(LogKit::LOW,"\nTime surfaces:\n");
  if (modelSettings->getParallelTimeSurfaces())
  {
    LogKit::LogFormatted(LogKit::LOW,"  Surface                                  : %s\n",     inputFiles->getTimeSurfFile(1).c_str());
    LogKit::LogFormatted(LogKit::LOW,"  Shift to top surface                     : %10.1f\n", modelSettings->getTimeDTop());
    LogKit::LogFormatted(LogKit::LOW,"  Time slice                               : %10.1f\n", modelSettings->getTimeLz());
    LogKit::LogFormatted(LogKit::LOW,"  Sampling density                         : %10.1f\n", modelSettings->getTimeDz());
    LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n",   int(modelSettings->getTimeLz()/modelSettings->getTimeDz()+0.5));
  }
  else
  {
    const std::string & topName  = inputFiles->getTimeSurfFile(0); 
    const std::string & baseName = inputFiles->getTimeSurfFile(1); 

    if (NRLib2::IsNumber(topName))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topName.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Top surface                              : %s\n",    topName.c_str());

    if (NRLib2::IsNumber(baseName))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n", atof(baseName.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Base surface                             : %s\n",   baseName.c_str());
      LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n", modelSettings->getTimeNz());

    LogKit::LogFormatted(LogKit::LOW,"  Minimum allowed value for lmin/lmax      : %10.2f\n", modelSettings->getLzLimit());
  }
  if (inputFiles->getCorrDirFile() != "")
    LogKit::LogFormatted(LogKit::LOW,"\n  Correlation direction                    : %10s\n", inputFiles->getCorrDirFile().c_str());

  if (modelSettings->getDoDepthConversion())
  {
    LogKit::LogFormatted(LogKit::LOW,"\nDepth conversion:\n");
    if (inputFiles->getDepthSurfFile(0) != "")
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : %s\n", inputFiles->getDepthSurfFile(0).c_str());
    else
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : %s\n", "Not given");
    if (inputFiles->getDepthSurfFile(1) != "")
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : %s\n", inputFiles->getDepthSurfFile(1).c_str());
    else
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : %s\n", "Not given");
    LogKit::LogFormatted(LogKit::LOW,"  Velocity field                           : %s\n", inputFiles->getVelocityField().c_str());
  }

  const std::string & topWEI  = inputFiles->getWaveletEstIntFile(0);
  const std::string & baseWEI = inputFiles->getWaveletEstIntFile(1);

  if (topWEI != "" || baseWEI != "") {
    LogKit::LogFormatted(LogKit::LOW,"\nWavelet estimation interval:\n");
    if (NRLib2::IsNumber(topWEI))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10s\n",topWEI.c_str());
    
    if (NRLib2::IsNumber(baseWEI))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10s\n",baseWEI.c_str());
  }

  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" || baseFEI != "") {
    LogKit::LogFormatted(LogKit::LOW,"\nFacies estimation interval:\n");
    if (NRLib2::IsNumber(topFEI))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10s\n",topFEI.c_str());
    
    if (NRLib2::IsNumber(baseFEI))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10s\n",baseFEI.c_str());
  }

  //
  // BACKGROUND
  //
  if (modelSettings->getGenerateBackground()) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBackground model (estimated):\n");
    if (inputFiles->getBackVelFile() != "")
      LogKit::LogFormatted(LogKit::LOW,"  Trend for p-wave velocity                : %10s\n",inputFiles->getBackVelFile().c_str());
    Vario       * vario  = modelSettings->getBackgroundVario();
    GenExpVario * pVario = dynamic_cast<GenExpVario*>(vario);
    LogKit::LogFormatted(LogKit::LOW,"  Variogram\n");
    LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",vario->getType());
    if (pVario != NULL)
    LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pVario->getPower());
    LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",vario->getRange());
    if (vario->getAnisotropic()) 
    {
      LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",vario->getSubRange());
      LogKit::LogFormatted(LogKit::LOW,"    Angle                                  : %10.1f\n",vario->getAngle()*(-1)*(180/M_PI));
    }
    LogKit::LogFormatted(LogKit::LOW,"  High cut frequency for well logs         : %10.1f\n",modelSettings->getMaxHzBackground());
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBackground model:\n");
    if (modelSettings->getConstBackValue(0) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(0));
    else
      LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity read from file           : %10s\n",inputFiles->getBackFile(0).c_str());
    
    if (modelSettings->getConstBackValue(1) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(1));
    else
      LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity read from file           : %10s\n",inputFiles->getBackFile(1).c_str());
      
    if (modelSettings->getConstBackValue(2) > 0)
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10.1f\n",modelSettings->getConstBackValue(2));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Density read from file                   : %10s\n",inputFiles->getBackFile(2).c_str());
  }
  if (modelSettings->getGenerateSeismic())
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
      LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : %s\n",inputFiles->getWaveletFile(i).c_str());
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
      LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",corr->getType());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",corr->getRange());
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",corr->getSubRange());
        LogKit::LogFormatted(LogKit::LOW,"    Angle                                  : %10.1f\n",corr->getAngle()*(-1)*(180/M_PI));
      }
    }
    //
    // SEISMIC
    //
    LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for seismic:\n");
    LogKit::LogFormatted(LogKit::LOW,"  White noise component                    : %10.2f\n",modelSettings->getWNC());
    LogKit::LogFormatted(LogKit::LOW,"  Low cut for inversion                    : %10.1f\n",modelSettings->getLowCut());
    LogKit::LogFormatted(LogKit::LOW,"  High cut for inversion                   : %10.1f\n",modelSettings->getHighCut());
    corr  = modelSettings->getAngularCorr();
    GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
    LogKit::LogFormatted(LogKit::LOW,"  Angular correlation:\n");
    LogKit::LogFormatted(LogKit::LOW,"    Model                                  : %10s\n",corr->getType());
    if (pCorr != NULL)
      LogKit::LogFormatted(LogKit::LOW,"    Power                                  : %10.1f\n",pCorr->getPower());
    LogKit::LogFormatted(LogKit::LOW,"    Range                                  : %10.1f\n",corr->getRange()*180.0/PI);
    if (corr->getAnisotropic())
    {
      LogKit::LogFormatted(LogKit::LOW,"    Subrange                               : %10.1f\n",corr->getSubRange()*180.0/PI);
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
    
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {
      LogKit::LogFormatted(LogKit::LOW,"\nSettings for AVO stack %d:\n",i+1);
      LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10.1f\n",(modelSettings->getAngle(i)*180/PI));
      LogKit::LogFormatted(LogKit::LOW,"  Segy offset                              : %10.1f\n",modelSettings->getSegyOffset());
      LogKit::LogFormatted(LogKit::LOW,"  Data                                     : %s\n",inputFiles->getSeismicFile(i).c_str());
      if (modelSettings->getEstimateWavelet(i))
        LogKit::LogFormatted(LogKit::LOW,"  Estimate wavelet                         : %10s\n", "yes");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : %s\n",inputFiles->getWaveletFile(i).c_str());
      if (modelSettings->getUseLocalWavelet()) {
        bool estimateLocalWaveletData = true;
        if (estimateLocalWaveletData) {
          LogKit::LogFormatted(LogKit::LOW,"  Estimate local shift map                 : %10s\n", "yes");
          LogKit::LogFormatted(LogKit::LOW,"  Estimate local gain map                  : %10s\n", "yes");
        }
        else {
          LogKit::LogFormatted(LogKit::LOW,"  Local shift map                          : %10s\n", "filnavn1");
          LogKit::LogFormatted(LogKit::LOW,"  Local gain map                           : %10s\n", "filnavn2");
        }
      }
      if (modelSettings->getEstimateSNRatio(i)) 
        LogKit::LogFormatted(LogKit::LOW,"  Estimate signal-to-noise ratio           : %10s\n", "yes");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Signal-to-noise ratio                    : %10.1f\n",modelSettings->getSNRatio(i));
      if (modelSettings->getMatchEnergies(i)) 
        LogKit::LogFormatted(LogKit::LOW,"  Match empirical and theoretical energies : %10s\n", "yes");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Wavelet scale                            : %10.2e\n",modelSettings->getWaveletScale(i));
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
                              char          * errText, 
                              bool          & failed)
{
  if(modelSettings->getDoDepthConversion()) 
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
        timeDepthMapping_->calculateSurfaceFromVelocity(velocity, timeSimbox);
        timeDepthMapping_->setDepthSimbox(timeSimbox, timeSimbox->getnz(), 
                                          modelSettings->getFormatFlag(),
                                          failed, errText);            // NBNB-PAL: Er dettet riktig nz (timeCut vs time)? 
        timeDepthMapping_->makeTimeDepthMapping(velocity, timeSimbox);
      }
    }
    if(velocity != NULL)
      delete velocity;
  }
}

void 
Model::loadVelocity(FFTGrid          *& velocity,
                    Simbox            * timeSimbox,
                    ModelSettings     * modelSettings, 
                    const std::string & velocityField, 
                    char              * errText,
                    bool              & failed)
{
  if(velocityField == "CONSTANT")
    velocity = NULL;
  else if(velocityField == "FROM_INVERSION")
  {
    velocityFromInversion_ = true;
    velocity = NULL;
  }
  else
  {
    const char * parName = "velocity field";
    char tmpErrText[MAX_STRING];
    sprintf(tmpErrText,"%c",'\0');
    int readerror = 0;

    if(findFileType(velocityField) == SEGYFILE) {
      const SegyGeometry * geometry = NULL;
      readerror = readSegyFile(velocityField, velocity,
                               timeSimbox, modelSettings, 
                               tmpErrText, geometry, 
                               FFTGrid::VELOCITY);
    }
    else
      readerror = readStormFile(velocityField, velocity, parName, 
                                timeSimbox, modelSettings,
                                tmpErrText);
    if (readerror==0) {
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
        text += "\n  Minimum Vp = "+NRLib2::ToString(logMin,2)+"    Number of too low values  : "+NRLib2::ToString(tooLow);
        text += "\n  Maximum Vp = "+NRLib2::ToString(logMax,2)+"    Number of too high values : "+NRLib2::ToString(tooHigh);
        text += "\nThe range of allowed values can changed using the ALLOWED_PARAMETER_VALUES keyword\n";
        text += "\naborting...\n";
        sprintf(errText,"%sReading of file \'%s\' for parameter \'%s\' failed\n%s\n", 
                errText,velocityField.c_str(),parName,text.c_str());
        failed = true;
      } 
    }
    else {
      sprintf(errText,"%sReading of file \'%s\' for parameter \'%s\' failed\n%s\n", 
              errText,velocityField.c_str(),parName,tmpErrText);
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

  LogKit::LogFormatted(LogKit::LOW,"\nThe top and/or base time surfaces do not cover the area specified by the %s",text.c_str());
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
  const NRLib2::Surface<double> & top  = timeSimbox->GetTopSurface();
  const NRLib2::Surface<double> & base = timeSimbox->GetBotSurface();
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
Model::writeSurfaceToFile(Surface     * surface,
                          std::string   name,
                          int           format)
{
  std::string fileName = ModelSettings::makeFullFileName(name);
  
  if((format & FFTGrid::ASCIIFORMAT) == FFTGrid::ASCIIFORMAT)
    NRLib2::WriteIrapClassicAsciiSurf(*surface, fileName+".irap");
  else  
    NRLib2::WriteStormBinarySurf(*surface, fileName+".storm");
}
