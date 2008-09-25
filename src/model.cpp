
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

#include "lib/random.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/global_def.h"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"

Model::Model(char * fileName)
{
  modelSettings_          = NULL;
  timeSimbox_             = new Simbox();
  wells_                  = NULL;
  background_             = NULL;
  priorCorrelations_      = NULL;
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
  gradX_                  = 0;
  gradY_                  = 0;
 
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
  bool failedVelocity     = false;
  bool failedLoadingModel = false;

  ModelFile * modelFile   = new ModelFile(fileName);
  Simbox * timeCutSimbox = NULL;
  if (modelFile->getParsingFailed()) {
    failedModelFile = true;
  }
  else
  {
    modelSettings_ = modelFile->getModelSettings();

    hasSignalToNoiseRatio_ = modelFile->getHasSignalToNoiseRatio();     // NBNB-PAL: ==> UT

    LogKit::SetScreenLog(modelSettings_->getLogLevel());

    char * logFileName = ModelSettings::makeFullFileName("logFile.txt");
    LogKit::SetFileLog(logFileName,modelSettings_->getLogLevel());
    delete [] logFileName;

    if(modelSettings_->getDebugFlag() > 0)
    {
      char * fName = ModelSettings::makeFullFileName("debug",".txt");
      LogKit::SetFileLog(std::string(fName), LogKit::DEBUGHIGH+LogKit::DEBUGLOW);
      delete [] fName;
    }
    LogKit::EndBuffering();
    
    if(modelFile->getSeedFile()==NULL)
      randomGen_ = new RandomGen(modelFile->getSeed());
    else
      randomGen_ = new RandomGen(modelFile->getSeedFile());

    if(modelSettings_->getNumberOfSimulations() == 0)
      modelSettings_->setOutputFlag(ModelSettings::PREDICTION); //write predicted grids. 
    
    printSettings(modelSettings_, modelFile, hasSignalToNoiseRatio_);
    
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
    LogKit::LogFormatted(LogKit::LOW,"\n***                       Reading input data                        ***"); 
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");

    char errText[MAX_STRING];
    sprintf(errText,"%c",'\0');

    int format = 0;
    if((modelSettings_->getFormatFlag() & FFTGrid::STORMASCIIFORMAT) == FFTGrid::STORMASCIIFORMAT)
      format += 1;
    if((modelSettings_->getFormatFlag() & FFTGrid::STORMFORMAT) == FFTGrid::STORMFORMAT)
      format +=2;

    makeTimeSimbox(timeSimbox_, modelSettings_, modelFile, //Handles correlation direction too.
                   errText, failedSimbox, timeCutSimbox);

    if(!failedSimbox)
    { 
      //
      // FORWARD MODELLING
      //
      if (modelSettings_->getGenerateSeismic() == true)
      {
        processBackground(background_, wells_, timeSimbox_, 
                          modelSettings_, modelFile,
                          errText, failedBackground);
        if (!failedBackground)
        {
          processReflectionMatrix(reflectionMatrix_, background_, 
                                  modelSettings_, modelFile, 
                                  errText, failedReflMat);  
          if (!failedReflMat)
          {
            processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                            timeSimbox_, waveletEstimInterval_, shiftGrids_, gainGrids_,
                            modelSettings_, modelFile, hasSignalToNoiseRatio_,
                            errText, failedWavelet);
          }              
        }
        if(modelSettings_->getFormatFlag() == 0)
          modelSettings_->setFormatFlag(1);  //Default, but not initialized due to possible double output.
        background_->getAlpha()->setOutputFormat(modelSettings_->getFormatFlag()); //static, controls all grids.
      }
      else
      {
        //
        // INVERSION
        //
        processSeismic(seisCube_, timeSimbox_, 
                       modelSettings_, modelFile, 
                       errText, failedSeismic);
        
        if(failedSeismic==false)
        {
          completeTimeCutSimbox(timeCutSimbox, modelSettings_, errText,failedSimbox);   // Copies area to timeCutSimbox if needed. 
          if(timeCutSimbox!=NULL)
            timeCutMapping_ = new GridMapping(timeCutSimbox,modelFile, modelSettings_,0, failedSimbox, errText, format);
          if(failedSimbox==false)
          {
            FFTGrid * velocity = NULL;
            if(timeCutSimbox!=NULL)
              processVelocity(velocity, timeCutSimbox,
                              modelSettings_, modelFile, 
                              errText, failedVelocity);
            else
              processVelocity(velocity,timeSimbox_,
                              modelSettings_, modelFile, 
                              errText, failedVelocity);
            if(!failedVelocity && modelFile->getDoDepthConversion()==1)
              timeDepthMapping_ = new GridMapping(timeSimbox_, modelFile, modelSettings_, 1, failedSimbox, errText, format, velocity);
             if(velocity !=NULL)
              delete velocity;
           }
        }
        if(!(failedSeismic || failedSimbox))
        { 
          processWells(wells_, timeSimbox_, randomGen_, 
                       modelSettings_, modelFile, 
                       errText, failedWells);
          loadExtraSurfaces(waveletEstimInterval_, 
                            faciesEstimInterval_,
                            timeSimbox_, modelFile,
                            errText, failedExtraSurf);
          if (!failedWells)
          {
            processBackground(background_, wells_, timeSimbox_, 
                              modelSettings_, modelFile,
                              errText, failedBackground);

            if (!failedBackground)
            {
              processPriorCorrelations(priorCorrelations_, background_, wells_, timeSimbox_,
                                       modelSettings_,modelFile,
                                       errText, failedPriorCorr);
              processReflectionMatrix(reflectionMatrix_, background_, 
                                      modelSettings_, modelFile, 
                                      errText, failedReflMat);
              if (!failedReflMat && !failedExtraSurf)
              {
                processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                                timeSimbox_, waveletEstimInterval_, shiftGrids_, gainGrids_,
                                modelSettings_, modelFile, hasSignalToNoiseRatio_, 
                                errText, failedWavelet);
              }
            }
          }
        }
        if (!failedWells && !failedExtraSurf)
        {
          processPriorFaciesProb(priorFacies_,
                                 wells_,
                                 randomGen_,
                                 modelFile->getTimeNz(),
                                 modelSettings_);
        }
      }
    }

    failedLoadingModel = failedSimbox  || failedSeismic  || failedPriorCorr  ||
                         failedWells   || failedReflMat  || failedBackground ||
                         failedWavelet || failedVelocity || failedExtraSurf;

    if (failedLoadingModel) 
      LogKit::LogFormatted(LogKit::ERROR,"\nERROR(s) while loading model. \n %s", errText);
  }
  failed_ = failedModelFile || failedLoadingModel;
  
  delete modelFile;
  if(timeCutSimbox !=NULL)
    delete timeCutSimbox;
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

  if (priorCorrelations_ != NULL)
    delete priorCorrelations_;

  if(timeDepthMapping_!=NULL)
    delete timeDepthMapping_;

  if(timeCutMapping_!=NULL)
    delete timeCutMapping_;

  if(correlationDirection_ !=NULL)
    delete correlationDirection_;
  delete randomGen_;
  delete modelSettings_;
  delete timeSimbox_;
}

void
Model::releaseGrids(void)
{
  delete background_;
  seisCube_ = NULL;
}

float **
Model::readMatrix(char * fileName, int n1, int n2, const char * readReason, char * errText)
{
  float * tmpRes = new float[n1*n2+1];
  FILE * inFile = fopen(fileName,"r");
  LogKit::LogFormatted(LogKit::LOW,"Reading %s from file %s .... ",readReason, fileName);
  char storage[MAX_STRING];
  int index = 0;
  int error = 0;
  while(error == 0 && fscanf(inFile,"%s",storage) != EOF) {
    if(index < n1*n2) {
      if(isNumber(storage) == 1) {
        tmpRes[index] = float(atof(storage));
      }
      else {
        sprintf(errText,"Found '%s' in file %s, expected a number.\n", storage, fileName);
        error = 1;
      }
    }
    index++;
  }
  if(error == 0) {
    if(index != n1*n2) {
      error = 1;
      sprintf(errText,"Found %d elements in file %s, expected %d.\n",
        index, fileName,n1*n2);
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
                            ModelSettings * modelSettings)
{
  if (modelSettings->getFileGrid() > 0) // Disk buffering is turn on
    return;

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
  int workMem = 2500+int( 0.65*gridSize); //Size of memory used beyond grids.

  char   * memchunk0 = new char[workMem];
  float ** memchunk  = new float*[nGrids];

  float    megaBytes = static_cast<float>(4*nGrids)*static_cast<float>(gridSize)/(1024.f*1024.f);
  float    gigaBytes = megaBytes/1024.f;
  if (gigaBytes < 0.01f)
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f megaBytes\n",megaBytes);
  else
    LogKit::LogFormatted(LogKit::LOW,"\nMemory needed by CRAVA:  %.2f gigaBytes\n",gigaBytes);
  
  int i;
  for(i = 0 ; i < nGrids ; i++)
    memchunk[i]  = new float[gridSize];

  if(memchunk[nGrids-1] == NULL)  //Could not allocate memory
  {
    modelSettings->setFileGrid(1);
    LogKit::LogFormatted(LogKit::LOW,"Not enough memory to hold all grids. Using file storage.\n");
  }
  else
  {
    modelSettings->setFileGrid(0);
  }

  for(i = 0;i < nGrids ;i++)
    if(memchunk[i] != NULL) delete [] memchunk[i];
  if(memchunk != NULL) delete [] memchunk;

  if(memchunk0 != NULL) delete[] memchunk0;
}


// readSegyFiles: reads SegY files form fNames table and stores in target table.
int
Model::readSegyFiles(char          ** fNames, 
                     int              nFiles, 
                     FFTGrid       ** target, 
                     Simbox        *& timeSimbox, 
                     SegyGeometry **& geometries,
                     ModelSettings *& modelSettings, 
                     char           * errText, 
                     int              fileno)
{
  char tmpErr[MAX_STRING];
  int error = 0;
  int sbError = 0;
  strcpy(errText, "");
  strcpy(tmpErr,"");
  SegY * segy = NULL;
  // flag fileno used if only one of the files should be read from. Used in processBackground.
  for(int i=0 ; i<nFiles ; i++)
  {
    if(fileno==-1 || fileno==i)
    {
      target[i] = NULL;
      if(error == 0)
      {
        try
        {
          segy = new SegY(fNames[i], 
                          modelSettings->getSegyOffset(),
                          *(modelSettings->getTraceHeaderFormat()));
          bool onlyVolume = modelSettings->getAreaParameters() != NULL;;
          segy->readAllTraces(timeSimbox, 
                              modelSettings->getZpad(),
                              onlyVolume);
          segy->createRegularGrid();
        }
        catch (NRLib2::Exception & e)
        {
          sprintf(errText,"%s%s",errText,e.what());
          error++;
        }
        
        if (error == 0)
        {
            const SegyGeometry *geometry;
            geometry = segy->getGeometry();
            geometry->writeGeometry();
            geometries[i] = new SegyGeometry(geometry);
          if(timeSimbox->status() == Simbox::NOAREA)
          {
            modelSettings->setAreaParameters(geometry);
            sbError = timeSimbox->setArea(geometry, tmpErr);
            if(sbError==1)
            {
              sprintf(errText,"%s Could not define time simulation grid.",tmpErr);
              error++;
            }
            else
            {
              sbError = timeSimbox->checkError(modelSettings->getLzLimit(), tmpErr);

              if(sbError == 0)
              {
                estimateXYPaddingSizes(timeSimbox, modelSettings);
                //
                // Check if CRAVA has enough memory to run calculation without buffering to disk
                //
                checkAvailableMemory(timeSimbox, modelSettings);
              }
              else if(sbError == Simbox::INTERNALERROR)
              {
                sprintf(errText,"%s%s", errText, tmpErr);
                error++;
              }
            }
          }
          else
          {
            sbError = timeSimbox->insideRectangle(geometry);
            if(sbError == 1)
            {
              sprintf(errText,"%sSpecified area in command AREA is larger than the seismic data from volume %s\n",
                      errText,fNames[i]);
              error++;
            }
          }
          if(sbError == 0)
          {
            if(modelSettings->getFileGrid() == 1)
              target[i] = new FFTFileGrid(timeSimbox->getnx(),
              timeSimbox->getny(), 
              timeSimbox->getnz(),
              modelSettings->getNXpad(), 
              modelSettings->getNYpad(), 
              modelSettings->getNZpad());
            else
              target[i] = new FFTGrid(timeSimbox->getnx(), 
              timeSimbox->getny(), 
              timeSimbox->getnz(),
              modelSettings->getNXpad(), 
              modelSettings->getNYpad(), 
              modelSettings->getNZpad());

            target[i]->setType(FFTGrid::DATA);
            target[i]->fillInFromSegY(segy, timeSimbox);
            target[i]->setAngle(modelSettings->getAngle()[i]);
          }
        }
        if (segy != NULL)
          delete segy;
      }
    }
  }
  return(error);
}



//NBNB Following routine only to be used for parameters!
int
Model::readStormFile(char           * fName, 
                     FFTGrid       *& target, 
                     const char     * parName, 
                     Simbox         * timeSimbox, 
                     ModelSettings *& modelSettings, 
                     char           * errText)
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
Model::makeTimeSimbox(Simbox        *& timeSimbox,
                      ModelSettings *& modelSettings, 
                      ModelFile      * modelFile,
                      char           * errText,
                      bool           & failed,
                      Simbox         *& timeCutSimbox)
{
  int error = 0;
  setSimboxSurfaces(timeSimbox, 
                    modelFile->getTimeSurfFile(), 
                    modelFile->getParallelTimeSurfaces(), 
                    modelFile->getTimeDTop(), 
                    modelFile->getTimeLz(), 
                    modelFile->getTimeDz(), 
                    modelFile->getTimeNz(),
                    error);
  const SegyGeometry * areaParams = NULL;
  if(error == 0)
  {
    sprintf(errText,"%c",'\0');
    const char * topname = "toptime.storm";
    const char * botname = "bottime.storm";
    if(!((modelSettings->getOutputFlag() & ModelSettings::NOTIME) > 0))
      timeSimbox->writeTopBotGrids(topname, botname);

    if(modelFile->getNWaveletTransfArgs() > 0 && timeSimbox->getIsConstantThick() == true)
      LogKit::LogFormatted(LogKit::WARNING,"\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");

    areaParams = modelSettings->getAreaParameters(); 
    estimateZPaddingSize(timeSimbox, modelSettings);   
    if (areaParams != NULL)
    {
      error = timeSimbox->setArea(areaParams, errText);
      if(error==1)
      {
        sprintf(errText,"%s The AREA specified in the model file extends outside the surface(s).\n",errText);
        failed = true;
      }
      else
      {
        error = timeSimbox->checkError(modelSettings->getLzLimit(),errText);

        if(error == 0)
        {
          LogKit::LogFormatted(LogKit::LOW,"\nTime output interval:\n");
          LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                         timeSimbox->getlz()*timeSimbox->getAvgRelThick(),
                         timeSimbox->getlz()*timeSimbox->getMinRelThick(),
                         timeSimbox->getlz());
          LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                         timeSimbox->getdz()*timeSimbox->getAvgRelThick(),
                         timeSimbox->getdz(),
                         timeSimbox->getdz()*timeSimbox->getMinRelThick());
        }
        else
        {
          sprintf(errText,"%s. Could not make time simulation grid.\n",errText);
          failed = true;
        }
      }
    }
  }
  else
  {
    timeSimbox->externalFailure();
    failed = true;
  }

  if(modelFile->getCorrDirFile() != NULL) {
    //
    // Get correlation direction
    //
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getCorrDirFile());
      correlationDirection_ = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText,"%s%s",errText,e.what());
      failed = true;
    }
    if(failed == false && modelSettings->getGenerateSeismic() == false)
      setupExtendedTimeSimbox(timeSimbox, correlationDirection_, timeCutSimbox); 
      //Extends timeSimbox for correlation coverage. Original stored in timeCutSimbox_
    error = timeSimbox->checkError(modelSettings->getLzLimit(),errText);
    if(error == 0)
    {
      if(areaParams!=NULL)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nTime inversion interval:\n");
        LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                       timeSimbox->getlz()*timeSimbox->getAvgRelThick(),
                       timeSimbox->getlz()*timeSimbox->getMinRelThick(),
                       timeSimbox->getlz());
        LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                       timeSimbox->getdz()*timeSimbox->getAvgRelThick(),
                       timeSimbox->getdz(),
                       timeSimbox->getdz()*timeSimbox->getMinRelThick());
      }
    }
    else
    {
      sprintf(errText,"%s Could not make the time simulation grid.\n",errText);
      failed = true;
    }
  }

  if(failed == false && areaParams!=NULL) {
    estimateXYPaddingSizes(timeSimbox, modelSettings);
    //
    // Check if CRAVA has enough memory to run calculation without buffering to disk
    //
    checkAvailableMemory(timeSimbox, modelSettings); 
  }
}

void 
Model::setSimboxSurfaces(Simbox *& simbox, 
                         char   ** surfFile, 
                         bool      parallelSurfaces, 
                         double    dTop,
                         double    lz, 
                         double    dz, 
                         int       nz,
                         int     & error)
{
  Surface * z0Grid = NULL;
  try {
    Surface tmpSurf = NRLib2::ReadStormSurf(surfFile[0]);
    z0Grid = new Surface(tmpSurf);
  }
  catch (NRLib2::Exception & e) {
    LogKit::LogFormatted(LogKit::ERROR,e.what());
    error = 1;
  }

  if(error == 0)
  {
    if(parallelSurfaces) //Only one reference surface
    {
      simbox->setDepth(z0Grid, dTop, lz, dz);
    }
    else
    {
      Surface * z1Grid = NULL;
      try {
        Surface tmpSurf = NRLib2::ReadStormSurf(surfFile[1]);
        z1Grid = new Surface(tmpSurf);
      }
      catch (NRLib2::Exception & e) {
        LogKit::LogFormatted(LogKit::ERROR,e.what());
        error = 1;
      }
      if(error == 0)
        simbox->setDepth(z0Grid, z1Grid, nz);
    }
  }
}

void
Model::setupExtendedTimeSimbox(Simbox  * timeSimbox, 
                               Surface * corrSurf, Simbox *& timeCutSimbox)
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
  botSurf->Add(shiftBot);
  botSurf->Add(&(timeSimbox->GetBotSurface()));

  double thick = shiftBot-shiftTop;
  int nz = int(0.5+thick/timeCutSimbox->getdz());//NBNB Ragnar: Rounding - is it ok?

  timeSimbox->setDepth(topSurf, botSurf, nz);

  delete refPlane;
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
  for(i=0;i<surf->GetN();i++) {
    surf->GetXY(i, x, y);
    z = (*surf)(i);
    A[0][1] += x;
    A[0][2] += y;
    A[1][1] += x*x;
    A[1][2] += x*y;
    A[2][2] += y*y;
    b[0] += z;
    b[1] += x*z;
    b[2] += y*z;
  }

  A[0][0] = surf->GetN();
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

    xPadFac       = float(MINIM(1.0f, xPad / timeSimbox->getlx())); // A padding of more than 100% is insensible
    yPadFac       = float(MINIM(1.0f, yPad / timeSimbox->getly()));

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
    float factor  = 1.0f;
    float wLength = 300.0f;                  // Assume a wavelet is approx 300ms.
    zPad          = factor * wLength / 2.0f; // Use one wavelet as padding
    zPadFac       = float(MINIM(1.0f,zPad / (timeSimbox->getlz()*timeSimbox->getMinRelThick())));
    
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
Model::completeTimeCutSimbox(Simbox        *& timeCutSimbox,
                             ModelSettings  * modelSettings,
                             char            * errText,
                             bool            &failed)
{
  if(timeCutSimbox != NULL && timeCutSimbox->status() == Simbox::NOAREA)
  {
    const SegyGeometry * areaParams = modelSettings->getAreaParameters(); 
    if (areaParams != NULL)
    {
      int error = timeCutSimbox->setArea(areaParams, errText);
      if(error==1)
      {
        sprintf(errText,"%s Problem with definition of time cut simbox.",errText);
        failed = true;
      }

    }
  }
}

void
Model::processSeismic(FFTGrid      **& seisCube,
                      Simbox        *& timeSimbox,
                      ModelSettings *& modelSettings, 
                      ModelFile      * modelFile,
                      char           * errText,
                      bool           & failed)
{
  char ** seismicFile = modelFile->getSeismicFile();
  int error = 0;

  if(seismicFile != NULL)
  {
    char tmpErrText[MAX_STRING];
    bool areaInModelFile = modelSettings->getAreaParameters() != NULL;
    int nAngles = modelSettings->getNumberOfAngles();
    SegyGeometry ** geometry = new SegyGeometry * [nAngles];
    for (int i =0 ; i < nAngles ; i++)
      geometry[i] = NULL;

    seisCube = new FFTGrid * [nAngles];

    if(readSegyFiles(seismicFile, nAngles, seisCube, timeSimbox, 
                     geometry, modelSettings, tmpErrText) == 0)
    {
      int formatFlag = modelSettings->getFormatFlag();
      if(formatFlag == 0)
        formatFlag = 1;   //Default, but not initialized due to possible double output.
      seisCube[0]->setOutputFormat(formatFlag); //static, controls all grids.
      if(modelSettings->getDebugFlag() == 1)
      {
        char sName[100];
        for(int i=0 ; i<nAngles ; i++)
        {
          sprintf(sName, "origSeis%d",i);
          seisCube[i]->writeFile(sName, timeSimbox);
        }
      }

      LogKit::LogFormatted(LogKit::LOW,"\nArea/resolution           x0           y0         dx      dy     azimuth\n");
      LogKit::LogFormatted(LogKit::LOW,"------------------------------------------------------------------------\n");
      if (areaInModelFile) {
        double azimuth = (-1)*timeSimbox->getAngle()*(180.0/M_PI);
        if (azimuth < 0)
          azimuth += 360.0;
        LogKit::LogFormatted(LogKit::LOW,"Model file       %11.2f  %11.2f    %7.2f %7.2f    %8.3f\n", 
                             timeSimbox->getx0(), timeSimbox->gety0(), timeSimbox->getdx(), 
                             timeSimbox->getdy(), azimuth);
      }
      for (int i = 0 ; i < nAngles ; i++)
        if (geometry[i] != NULL) {
          double geoAngle = (-1)*timeSimbox->getAngle()*(180/M_PI);
          if (geoAngle < 0)
            geoAngle += 360.0f;
          LogKit::LogFormatted(LogKit::LOW,"Seismic data %d   %11.2f  %11.2f    %7.2f %7.2f    %8.3f\n",i,
                               geometry[i]->getX0(), geometry[i]->getY0(), geometry[i]->getDx(), 
                               geometry[i]->getDy(), geoAngle);
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
      error = 1;
    }
  }
  failed = error > 0;
}

void 
Model::processWells(WellData     **& wells,
                    Simbox         * timeSimbox,
                    RandomGen      * randomGen,
                    ModelSettings *& modelSettings, 
                    ModelFile      * modelFile,
                    char           * errText,
                    bool           & failed)
{
  LogKit::LogFormatted(LogKit::LOW,"\nReading well data:\n");

  char ** wellFile       = modelFile->getWellFile();
  char ** headerList     = modelFile->getHeaderList();
  bool    faciesLogGiven = modelFile->getFaciesLogGiven();
  int     nWells         = modelSettings->getNumberOfWells();
  int     nFacies        = 0;

  int error = 0;

  char tmpErrText[MAX_STRING];
  sprintf(tmpErrText,"%c",'\0');
  wells = new WellData *[nWells];
  for(int i=0 ; i<nWells ; i++) {
    wells[i] = new WellData(wellFile[i], modelSettings, headerList, faciesLogGiven, i);
    if(wells[i]->checkError(tmpErrText) != 0) {
      sprintf(errText,"%s%s", errText, tmpErrText);
      error = 1;
    }
  }

  if (error == 0) {
    if(modelFile->getFaciesLogGiven()) { 
      checkFaciesNames(wells, modelSettings, tmpErrText, error);
      nFacies = modelSettings->getNumberOfFacies(); // nFacies is set in checkFaciesNames()
    }
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
    LogKit::LogFormatted(LogKit::LOW,"\n***                       Processing Wells                          ***"); 
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");
    
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
          wells[i]->setBlockedLogsPropThick( new BlockedLogs(wells[i], timeSimbox, randomGen) );
          if (nFacies > 0)
            wells[i]->countFacies(timeSimbox,faciesCount[i]);
          if((modelSettings->getOutputFlag() & ModelSettings::WELLS) > 0) 
            wells[i]->writeRMSWell();
          if(modelSettings->getDebugFlag() > 0)
            wells[i]->getBlockedLogsPropThick()->writeToFile(static_cast<float>(timeSimbox->getdz()), 2, true); // 2 = BG logs
          
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
        LogKit::LogFormatted(LogKit::LOW,"%12s ",modelSettings->getFaciesName(i));
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
            float faciesProb = static_cast<float>(faciesCount[i][f])/tot;
            LogKit::LogFormatted(LogKit::LOW,"%12.4f ",faciesProb);
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
        LogKit::LogFormatted(LogKit::MEDIUM,"%12s ",modelSettings->getFaciesName(i));
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
      LogKit::LogFormatted(LogKit::LOW,"\nAborting\n");
      error = 1;
    }
  }
  failed = error > 0;
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

  int     nnames  = globalmax - globalmin + 1;
  char ** names   = new char * [nnames];

  for (int i =0 ; i < nnames ; i++) {
    names[i] = NULL;
  }
  
  for(int w=0 ; w<modelSettings->getNumberOfWells() ; w++)
  {
    if(wells[w]->isFaciesLogDefined())
    {
      for(int i=0 ; i < wells[w]->getNFacies() ; i++)
      {
        char * name = wells[w]->getFaciesName(i);
        int    fnr  = wells[w]->getFaciesNr(i) - globalmin;

        if(names[fnr] == NULL) {
          names[fnr]   = name;
        }
        else if(strcmp(names[fnr],name) != 0)
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
    if(names[i] != NULL) 
      LogKit::LogFormatted(LogKit::LOW,"    %2d           %-20s\n",i+globalmin,names[i]);

  int nFacies = 0;
  for(int i=0 ; i<nnames ; i++)
    if(names[i] != NULL) 
      nFacies++;

  char ** faciesNames = new char * [nFacies];

  int j = -1;
  for(int i=0 ; i<nnames ; i++)
  {
    if(names[i] != NULL)
    {
      j++;
      faciesNames[j]   = names[i];
    }
  }
  modelSettings->setNumberOfFacies(nFacies);
  modelSettings->setFaciesNames(faciesNames,nFacies);

  delete [] faciesNames;
  delete [] names;
}

void 
Model::processBackground(Background   *& background,
                         WellData     ** wells,
                         Simbox        * timeSimbox,
                         ModelSettings * modelSettings, 
                         ModelFile     * modelFile,
                         char          * errText,
                         bool          & failed)
{
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::BACKGROUND) > 0 || 
      (modelSettings->getOutputFlag() & ModelSettings::WAVELETS)   > 0 )
  {
    FFTGrid * backModel[3];
    const int nx    = timeSimbox->getnx();
    const int ny    = timeSimbox->getny();
    const int nz    = timeSimbox->getnz();
    const int nxPad = modelSettings->getNXpad();
    const int nyPad = modelSettings->getNYpad();
    const int nzPad = modelSettings->getNZpad();
    if (modelFile->getGenerateBackground()) 
    {
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
      LogKit::LogFormatted(LogKit::LOW,"\n***              Prior Expectations / Background Model              ***"); 
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");
      
      if(modelSettings->getBackgroundVario() == NULL)
      {
        sprintf(errText,"%sDid not manage to make variogram for background modelling\n",errText);
        failed = true;
      }
      for (int i=0 ; i<3 ; i++)
      {
        backModel[i] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);              
        backModel[i]->setType(FFTGrid::PARAMETER);
      }
      background = new Background(backModel, wells, timeSimbox, modelSettings);
    }
    else 
    {
      const char * parName[]={"Vp background","Vs background","Rho background"};
      for(int i=0 ; i<3 ; i++)
      {
        float  * constBack = modelFile->getConstBack();
        char  ** backFile  = modelFile->getBackFile();

        if(constBack[i] < 0)
        {
          if(backFile[i] != NULL)
          {
            char tmpErrText[MAX_STRING];
            sprintf(tmpErrText,"%c",'\0');
            int readerror = 0;
            if(findFileType(backFile[i]) == SEGYFILE) {
              int nAngles = modelSettings->getNumberOfAngles();
              SegyGeometry ** geometry = new SegyGeometry * [nAngles];
              for (int j =0 ; j < nAngles ; j++)
                geometry[i] = NULL;
              readerror = readSegyFiles(backFile, 3, backModel, timeSimbox, 
                                        geometry, modelSettings,
                                        tmpErrText, i);
              for (int j =0 ; j < nAngles ; j++)
                if (geometry[i] != NULL)
                  delete geometry[i];
              delete [] geometry;
            }
            else
              readerror = readStormFile(backFile[i], backModel[i], parName[i], 
                                        timeSimbox, modelSettings,
                                        tmpErrText);
            backModel[i]->logTransf();
            
            if(readerror != 0)
            {
              sprintf(errText,"%sReading of file \'%s\' for parameter \'%s\' failed\n",
                      errText,backFile[i],parName[i]);
              failed = true;
            }
          }
          else
          {
            sprintf(errText,"%sReading of file for parameter \'%s\' failed. File pointer is NULL\n",
                    errText,parName[i]);
            failed = true;
          }
        }
        else if(constBack[i] > 0)
        {
          if(modelSettings->getFileGrid() == 1)
            backModel[i] = new FFTFileGrid(nx, ny, nz, nxPad, nyPad, nzPad);
          else
            backModel[i] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);              
          backModel[i]->setType(FFTGrid::PARAMETER);
          backModel[i]->fillInConstant(float( log( constBack[i] )));
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
      background->writeBackgrounds(timeSimbox); 
  }
}

void 
Model::processPriorCorrelations(Corr         *& priorCorrelations,
                                Background    * background,
                                WellData     ** wells,
                                Simbox        * timeSimbox,
                                ModelSettings * modelSettings, 
                                ModelFile     * modelFile,
                                char          * errText,
                                bool          & failed)
{
  bool printResult = (modelSettings->getOutputFlag() & (ModelSettings::PRIORCORRELATIONS + ModelSettings::CORRELATION)) > 0;
  if (modelSettings->getDoInversion() || printResult)
  {
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
    LogKit::LogFormatted(LogKit::LOW,"\n***                        Prior Covariance                         ***"); 
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");
    time_t timestart, timeend;
    time(&timestart);

    //
    // Parameter correlation can be set in model file.
    // Default NULL, indicating that estimate will be used.
    //
    float ** paramCorr = NULL;
    char   * paramCorrFile = modelFile->getParamCorrFile();
    if(paramCorrFile != NULL) 
    {
      char tmpErrText[MAX_STRING];
      sprintf(tmpErrText,"%c",'\0');
      paramCorr = readMatrix(paramCorrFile, 3, 3, "parameter correlation", tmpErrText);
      if(paramCorrFile == NULL) 
      {
        sprintf(errText,"%sReading of file \'%s\' for parameter correlation matrix failed\n%s\n",
                errText,paramCorrFile,tmpErrText);
        failed = true;
      }
      LogKit::LogFormatted(LogKit::LOW,"Parameter correlation read from file.\n\n");
    }

    Surface * CorrXY = findCorrXYGrid();

    if(modelSettings->getLateralCorr()==NULL) { // NBNB-PAL: this will never be true (default lateral corr)
      estimateCorrXYFromSeismic(CorrXY,
                                seisCube_,
                                modelSettings_->getNumberOfAngles());
      time(&timeend);
      LogKit::LogFormatted(LogKit::LOW,"\nEstimate parameter lateral correlation from seismic in %d seconds.\n",
                           static_cast<int>(timeend-timestart));
    }

    Analyzelog * analyze = new Analyzelog(wells, 
                                          background,
                                          timeSimbox, 
                                          modelSettings);

    priorCorrelations = new Corr(analyze->getPointVar0(), 
                                 analyze->getVar0(), 
                                 analyze->getCorrT(), 
                                 analyze->getNumberOfLags(),
                                 static_cast<float>(timeSimbox->getdz()), 
                                 CorrXY);
    delete analyze;

    if(priorCorrelations == NULL)
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
      priorCorrelations->setVar0(paramCorr);

    if(printResult)
      priorCorrelations->dumpResult();
    
    priorCorrelations->printVariancesToScreen();

    time(&timeend);
    LogKit::LogFormatted(LogKit::DEBUGLOW,"\n\nTime elapsed :  %d\n",timeend-timestart);  
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
                               ModelFile     * modelFile,                  
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
    
    char * reflMatrFile = modelFile->getReflMatrFile();
    if(reflMatrFile != NULL) {
      reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", errText);
      if(reflectionMatrix == NULL) {
        sprintf(errText,"%sReading of file \'%s\' for reflection matrix failed\n",errText,reflMatrFile);
        failed = true;
      }
      LogKit::LogFormatted(LogKit::LOW,"Reflection parameters read from file.\n\n");
    }
    else {
      if (background != NULL) {
        setupDefaultReflectionMatrix(reflectionMatrix,
                                     background,
                                     modelSettings,
                                     modelFile);
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
                                    ModelSettings * modelSettings,
                                    ModelFile     * modelFile)
{
  int i;
  float ** A = new float * [modelSettings->getNumberOfAngles()];
  // For debugging
  //background->setClassicVsVp();
  double vsvp  = background->getMeanVsVp();
  double vsvp2 = vsvp*vsvp;
  for(i = 0; i < modelSettings->getNumberOfAngles(); i++)
  {
    double angle = static_cast<double>(modelSettings->getAngle()[i]);
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(modelFile->getSeisType(i) == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float((1.0 +tan2t )/2.0) ; 
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0);
    }
    else if(modelFile->getSeisType(i) == ModelSettings::PSSEIS) {
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
                       ModelFile     * modelFile,
                       bool          & hasSignalToNoiseRatio,
                       char          * errText,
                       bool          & failed)
{
  int error = 0;
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::WAVELETS) > 0 )
  {
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
    LogKit::LogFormatted(LogKit::LOW,"\n***                 Processing/generating wavelets                  ***"); 
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");
    bool estimateStuff = false;
    for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {  
      estimateStuff = estimateStuff || (modelFile->getWaveletFile()[i][0] == '*'); 
      estimateStuff = estimateStuff || (modelSettings->getNoiseEnergy()[i]==RMISSING); 
    }
    if (estimateStuff) 
    {
      LogKit::LogFormatted(LogKit::HIGH,"\nWells that cannot be used in wavelet generation or noise estimation:");
      LogKit::LogFormatted(LogKit::HIGH,"\n  Deviated wells.");
      LogKit::LogFormatted(LogKit::HIGH,"\n  Wells with too little data.\n");
    }
    wavelet = new Wavelet *[modelSettings->getNumberOfAngles()];
    
    char ** waveletFile = modelFile->getWaveletFile();
    float * waveScale   = modelFile->getWaveletScale();
    float * noiseEnergy = modelSettings->getNoiseEnergy();
    
    for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {  
      LogKit::LogFormatted(LogKit::LOW,"\nAngle stack : %.1f deg\n",modelSettings->getAngle()[i]*180.0/PI);
      if (waveletFile[i][0] == '*') 
      {
        if (timeSimbox->getdz() > 4.0f) { // Require this density for wavelet estimation
          LogKit::LogFormatted(LogKit::LOW,"\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
          LogKit::LogFormatted(LogKit::LOW,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
          LogKit::LogFormatted(LogKit::LOW,"         the number of layers must be increased.                                    \n");
        }
        wavelet[i] = new Wavelet1D(timeSimbox, 
                                   seisCube[i], 
                                   wells, 
                                   waveletEstimInterval,                                   
                                   modelSettings, 
                                   reflectionMatrix[i]);
      }
      else
      {
        int fileFormat = getWaveletFileFormat(waveletFile[i],errText);
        if(fileFormat < 0)
        {
          sprintf(errText, "%s Unknown file format of file  %s.\n", errText, waveletFile[i]);
          error += 1;
        }
        else {
          if (fileFormat == Wavelet::SGRI)
            wavelet[i] = new Wavelet3D(waveletFile[i], 
                                       modelSettings, 
                                       timeSimbox, 
                                       modelSettings->getAngle()[i], 
                                       reflectionMatrix[i],
                                       error, 
                                       errText);
          else {
            wavelet[i] = new Wavelet1D(waveletFile[i], 
                                       fileFormat, 
                                       modelSettings, 
                                       reflectionMatrix[i],
                                       error, 
                                       errText);
            if (error == 0) {
              //wavelet[i]->write1DWLas3DWL(); //Frode: For debugging and testing
              //wavelet[i]->write3DWLfrom1DWL();
              wavelet[i]->resample(static_cast<float>(timeSimbox->getdz()), timeSimbox->getnz(), 
                                   modelSettings->getZpad(), modelSettings->getAngle()[i]);
            }
          }
        }
      }
      if (error == 0) {
        if (waveScale[i] != RMISSING)        // If RMISSING we will later scale wavelet to get EmpSN = TheoSN.
          wavelet[i]->scale(waveScale[i]);
        if ((wavelet[i]->getDim() == 3) && !timeSimbox->getIsConstantThick()) {
          sprintf(errText, "%s Simbox must have constant thickness when 3D wavelet.\n", errText);
          error += 1;
        }
        if (noiseEnergy[i] == RMISSING)
        {
          if (wavelet[i]->getDim() == 3) { //Not possible to estimate noise variance when 3D wavelet
            sprintf(errText, "%s Estimation of noise st. dev. is not possible for 3D wavelet.\n", errText);
            sprintf(errText, "%s        Noise st.dev. or s/n ratio must be specified in modelfile\n", errText);
            error += 1;
          }
          else {
            noiseEnergy[i] = wavelet[i]->getNoiseStandardDeviation(timeSimbox, seisCube[i], wells, modelSettings->getNumberOfWells(), 
                                                                   errText, error);
            if (hasSignalToNoiseRatio) 
            {
              //LogKit::LogFormatted(LogKit::LOW,"\n  Since seismic noise is estimated directly, keyword GIVESIGNALTONOISERATIO has no effect.\n");
              hasSignalToNoiseRatio = false;
            }
          }
        }
        else
        {
          if (hasSignalToNoiseRatio && (noiseEnergy[i] <= 1.0 || noiseEnergy[i] > 10.0))
          {
            sprintf(errText, "%s Illegal signal-to-noise ratio of %.3f for cube %d\n", errText, noiseEnergy[i],i);
            sprintf(errText, "%s Ratio must be in interval 1.0 < SNratio < 10.0\n", errText);
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
  }
  failed = error > 0;
}

int
Model::getWaveletFileFormat(char * fileName, char * errText)
{
  int fileformat = -1;
  char* dummyStr = new char[MAX_STRING];
  // test for old file format
  FILE* file = fopen(fileName,"r");
  for(int i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName);
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
    file = fopen(fileName, "r");
    if (fscanf(file, "%s", dummyStr) == EOF)
    {
      sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName);
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
      file = fopen(fileName,"r");
      bool lineIsComment = true; 
      while( lineIsComment ==true)
      {
        if(fscanf(file,"%s",dummyStr) == EOF)
        {
          readToEOL(file);
          sprintf(errText,"%sEnd of wavelet file %s is premature\n",errText,fileName);
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
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
    LogKit::LogFormatted(LogKit::LOW,"\n***                     Prior Facies Probabilities                  ***"); 
    LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");

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
        BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
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
    // Print facies count for each well
    //
    LogKit::LogFormatted(LogKit::LOW,"BlockedWell             ");
    for (int i = 0 ; i < nFacies ; i++)
      LogKit::LogFormatted(LogKit::LOW,"%12s ",modelSettings->getFaciesName(i));
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
          LogKit::LogFormatted(LogKit::LOW," %d %12.4f ",faciesCount[w][i],faciesProb);
        }
        LogKit::LogFormatted(LogKit::LOW,"\n");
      }
    }
    LogKit::LogFormatted(LogKit::LOW,"\n");
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
        LogKit::LogFormatted(LogKit::LOW,"%-15s %10.4f\n",modelSettings->getFaciesName(i),priorFacies[i]);
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
                         ModelFile  * modelFile,
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
  if (modelFile->getWaveletEstIntFile() != NULL) {  
    waveletEstimInterval = new Surface*[2];

    char * topName  = modelFile->getWaveletEstIntFile()[0]; 
    char * baseName = modelFile->getWaveletEstIntFile()[1]; 

    try {
      if (isNumber(topName)) 
        waveletEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topName));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(topName);
        waveletEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s", errText,e.what());
      failed = true;
    }
    
    try {
      if (isNumber(baseName)) 
        waveletEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseName));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(baseName);
        waveletEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s", errText,e.what());
      failed = true;
    }
  }
  //
  // Get facies estimation interval
  //
  if (modelFile->getFaciesEstIntFile() != NULL) {  

    char * topName  = modelFile->getFaciesEstIntFile()[0]; 
    char * baseName = modelFile->getFaciesEstIntFile()[1]; 

    try {
      if (isNumber(topName)) 
        faciesEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topName));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(topName);
        faciesEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s", errText,e.what());
      failed = true;
    }

    try {
      if (isNumber(baseName)) 
        faciesEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseName));
      else { 
        Surface tmpSurf = NRLib2::ReadStormSurf(baseName);
        faciesEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText, "%s%s", errText,e.what());
      failed = true;
    }
  }
}

// NBNB The following routine is stupid, and assumes SEGY if file 
// does not start with 'storm_petro_binary'
int
Model::findFileType(char * fileName)
{
  FILE * file = fopen(fileName,"r");
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
                     ModelFile     * modelFile,
                     bool            hasSignalToNoiseRatio)
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
    if(modelFile->getSeedFile()==NULL) {
      if (modelFile->getSeed() == 0)
        LogKit::LogFormatted(LogKit::LOW,"  Seed                                     :          0 (default seed)\n");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Seed                                     : %10d\n",modelFile->getSeed());
    }
    else
      LogKit::LogFormatted(LogKit::LOW,"  Seed read from file                    : %10s\n",modelFile->getSeedFile());
    LogKit::LogFormatted(LogKit::LOW,"  Number of realisations                   : %10d\n",modelSettings->getNumberOfSimulations());
  }
  LogKit::LogFormatted(LogKit::LOW,"  Kriging                                  : %10s\n",(modelSettings->getKrigingParameters()==NULL ? "no" : "yes"));

  LogKit::LogFormatted(LogKit::LOW,"\nUnit settings/assumptions:\n");
  LogKit::LogFormatted(LogKit::LOW,"  Time                                     : %10s\n","ms TWT");
  LogKit::LogFormatted(LogKit::LOW,"  Frequency                                : %10s\n","Hz");
  LogKit::LogFormatted(LogKit::LOW,"  Length                                   : %10s\n","m");
  LogKit::LogFormatted(LogKit::LOW,"  Velocities                               : %10s\n","m/s");
  LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10s\n","kg/dm3");
  LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10s\n","degrees");

  TraceHeaderFormat * thf = modelSettings_->getTraceHeaderFormat(); 
  LogKit::LogFormatted(LogKit::LOW,"\nSegY trace header format:\n");
  LogKit::LogFormatted(LogKit::LOW,"  Format name                              : %10s\n",thf->getFormatName().c_str());
  if (thf->getBypassCoordScaling())
    LogKit::LogFormatted(LogKit::LOW,"  Bypass coordinate scaling                :        yes\n");
  if (!thf->getStandardType()) 
  {
    LogKit::LogFormatted(LogKit::LOW,"  Start pos coordinate scaling             : %10d\n",thf->getScalCoLoc());
    LogKit::LogFormatted(LogKit::LOW,"  Start pos trace x coordinate             : %10d\n",thf->getUtmxLoc());
    LogKit::LogFormatted(LogKit::LOW,"  Start pos trace y coordinate             : %10d\n",thf->getUtmyLoc());
    LogKit::LogFormatted(LogKit::LOW,"  Start pos inline index                   : %10d\n",thf->getInlineLoc());
    LogKit::LogFormatted(LogKit::LOW,"  Start pos crossline index                : %10d\n",thf->getCrosslineLoc());
    LogKit::LogFormatted(LogKit::LOW,"  Coordinate system                        : %10s\n",thf->getCoordSys()==0 ? "UTM" : "ILXL" );
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
    char ** headerList = modelFile->getHeaderList();
    if (headerList != NULL)
    {
      LogKit::LogFormatted(LogKit::LOW,"  Time                                     : %10s\n",  headerList[0]);
      if(strcmp(uppercase(headerList[1]),"VP"    )==0 ||
         strcmp(uppercase(headerList[1]),"LFP_VP")==0)
        LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity                          : %10s\n",headerList[1]);
      else
        LogKit::LogFormatted(LogKit::LOW,"  Sonic                                    : %10s\n",headerList[1]);
      if(strcmp(uppercase(headerList[3]),"VS"    )==0 || 
         strcmp(uppercase(headerList[3]),"LFP_VS")==0)
        LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity                          : %10s\n",headerList[3]);
      else
        LogKit::LogFormatted(LogKit::LOW,"  Shear sonic                              : %10s\n",headerList[3]);
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10s\n",  headerList[2]);
      if (modelFile->getFaciesLogGiven())
        LogKit::LogFormatted(LogKit::LOW,"  Facies                                   : %10s\n",headerList[4]);
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
      LogKit::LogFormatted(LogKit::LOW,"  %-2d                                       : %s\n",i+1,modelFile->getWellFile()[i]);
    }
    bool generateBackground = modelFile->getGenerateBackground();
    bool estimateFaciesProb = modelFile->getFaciesLogGiven();
    bool estimateWavelet    = false;
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
      estimateWavelet = estimateWavelet || modelFile->getWaveletFile()[i][0] == '*';
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
  if (modelFile->getParallelTimeSurfaces())
  {
    LogKit::LogFormatted(LogKit::LOW,"  Surface                                  : %s\n",     modelFile->getTimeSurfFile()[0]);
    LogKit::LogFormatted(LogKit::LOW,"  Shift to top surface                     : %10.1f\n", modelFile->getTimeDTop());
    LogKit::LogFormatted(LogKit::LOW,"  Time slice                               : %10.1f\n", modelFile->getTimeLz());
    LogKit::LogFormatted(LogKit::LOW,"  Sampling density                         : %10.1f\n", modelFile->getTimeDz());
    LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n",   int(modelFile->getTimeLz()/modelFile->getTimeDz()+0.5));
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"  Top surface                              : %s\n",     modelFile->getTimeSurfFile()[0]);
    LogKit::LogFormatted(LogKit::LOW,"  Base surface                             : %s\n",     modelFile->getTimeSurfFile()[1]);
    LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n",   modelFile->getTimeNz());
    LogKit::LogFormatted(LogKit::LOW,"  Minimum allowed value for lmin/lmax      : %10.2f\n", modelSettings->getLzLimit());
  }
  if (modelFile->getCorrDirFile() != NULL)
    LogKit::LogFormatted(LogKit::LOW,"\n  Correlation direction                    : %10s\n",   modelFile->getCorrDirFile());

  if (modelFile->getDoDepthConversion())
  {
    LogKit::LogFormatted(LogKit::LOW,"\nDepth conversion:\n");
    if (modelFile->getDepthSurfFile()[0] != NULL)
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : %s\n", modelFile->getDepthSurfFile()[0]);
    else
      LogKit::LogFormatted(LogKit::LOW,"  Top depth surface                        : %s\n", "Not given");
    if (modelFile->getDepthSurfFile()[1] != NULL)
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : %s\n", modelFile->getDepthSurfFile()[1]);
    else
      LogKit::LogFormatted(LogKit::LOW,"  Base depth surface                       : %s\n", "Not given");
    LogKit::LogFormatted(LogKit::LOW,"  Velocity field                           : %s\n", modelFile->getVelocityField());
  }

  if (modelFile->getWaveletEstIntFile() != NULL) {
    char * topName  = modelFile->getWaveletEstIntFile()[0]; 
    char * baseName = modelFile->getWaveletEstIntFile()[1]; 
    LogKit::LogFormatted(LogKit::LOW,"\nWavelet estimation interval:\n");
    if (isNumber(topName))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topName));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10s\n",topName);
    
    if (isNumber(baseName))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseName));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10s\n",baseName);
  }

  if (modelFile->getFaciesEstIntFile() != NULL) {
    char * topName  = modelFile->getFaciesEstIntFile()[0]; 
    char * baseName = modelFile->getFaciesEstIntFile()[1]; 
    LogKit::LogFormatted(LogKit::LOW,"\nFacies estimation interval:\n");
    if (isNumber(topName))
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10.2f\n",atof(topName));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Start time                               : %10s\n",topName);
    
    if (isNumber(baseName))
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10.2f\n",atof(baseName));
    else
      LogKit::LogFormatted(LogKit::LOW,"  Stop time                                : %10s\n",baseName);
  }

  //
  // BACKGROUND
  //
  if (modelFile->getGenerateBackground()) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBackground model (estimated):\n");
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
      LogKit::LogFormatted(LogKit::LOW,"    Angle                                  : %10.1f\n",vario->getAngle());
    }
    LogKit::LogFormatted(LogKit::LOW,"  High cut frequency for well logs         : %10.1f\n",modelSettings->getMaxHzBackground());
  }
  else
  {
    float * constBack = modelFile->getConstBack();
    LogKit::LogFormatted(LogKit::LOW,"\nBackground model:\n");
    if (constBack[0] > 0)
      LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity                          : %10.1f\n",constBack[0]);
    else
      LogKit::LogFormatted(LogKit::LOW,"  p-wave velocity read from file           : %10s\n",modelFile->getBackFile()[0]);
    
    if (constBack[1] > 0)
      LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity                          : %10.1f\n",constBack[1]);
    else
      LogKit::LogFormatted(LogKit::LOW,"  s-wave velocity read from file           : %10s\n",modelFile->getBackFile()[1]);
      
    if (constBack[2] > 0)
      LogKit::LogFormatted(LogKit::LOW,"  Density                                  : %10.1f\n",constBack[2]);
    else
      LogKit::LogFormatted(LogKit::LOW,"  Density read from file                   : %10s\n",modelFile->getBackFile()[2]);
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
      LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10.1f\n",(modelSettings->getAngle()[i]*180/PI));
      LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : %s\n",modelFile->getWaveletFile()[i]);
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
        LogKit::LogFormatted(LogKit::LOW,"    Angle                                  : %10.1f\n",corr->getAngle());
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
      estimateNoise = estimateNoise || (modelSettings->getNoiseEnergy()[i]==RMISSING); 
    }
    LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings for wavelet:\n");
    if (estimateNoise)
      LogKit::LogFormatted(LogKit::LOW,"  Maximum shift in noise estimation        : %10.1f\n",modelSettings->getMaxWaveletShift());
    LogKit::LogFormatted(LogKit::LOW,"  Minimum relative amplitude               : %10.3f\n",modelSettings->getMinRelWaveletAmp());
    LogKit::LogFormatted(LogKit::LOW,"  Wavelet tapering length                  : %10.1f\n",modelSettings->getWaveletTaperingL());
    
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {
      LogKit::LogFormatted(LogKit::LOW,"\nSettings for AVO stack %d:\n",i+1);
      LogKit::LogFormatted(LogKit::LOW,"  Angle                                    : %10.1f\n",(modelSettings->getAngle()[i]*180/PI));
      LogKit::LogFormatted(LogKit::LOW,"  Segy offset                              : %10.1f\n",modelSettings->getSegyOffset());
      LogKit::LogFormatted(LogKit::LOW,"  Data                                     : %s\n",modelFile->getSeismicFile()[i]);
      if (modelFile->getWaveletFile()[i][0] == '*')
        LogKit::LogFormatted(LogKit::LOW,"  Estimate wavelet                         : %10s\n", "yes");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Read wavelet from file                   : %s\n",modelFile->getWaveletFile()[i]);
      
      float * noiseEnergy   = modelSettings->getNoiseEnergy();
      bool  * matchEnergies = modelSettings->getMatchEnergies();
      if (noiseEnergy[i] == RMISSING) 
        LogKit::LogFormatted(LogKit::LOW,"  Estimate signal-to-noise ratio           : %10s\n", "yes");
      else
        if (hasSignalToNoiseRatio)
          LogKit::LogFormatted(LogKit::LOW,"  Signal-to-noise ratio                    : %10.1f\n",noiseEnergy[i]);
        else
          LogKit::LogFormatted(LogKit::LOW,"  Error std.dev. in seismic                : %10.2e\n",noiseEnergy[i]);
      if (matchEnergies[i]) 
        LogKit::LogFormatted(LogKit::LOW,"  Match empirical and theoretical energies : %10s\n", "yes");
      else
        LogKit::LogFormatted(LogKit::LOW,"  Wavelet scale                            : %10.2e\n",modelFile->getWaveletScale()[i]);
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

  double cI = dx*cosrot*gradX_+dy*sinrot*gradY_;
  double cJ = -dx*sinrot*gradX_+dy*cosrot*gradY_;

  corrGradI = float(cI/timeSimbox_->getdz());
  corrGradJ = float(cJ/timeSimbox_->getdz());
}


void 
Model::processVelocity(FFTGrid      *& velocity,
                       Simbox        * timeSimbox,
                       ModelSettings * modelSettings, 
                       ModelFile     * modelFile, 
                       char          * errText,
                       bool          & failed)
{
  if(modelFile->getDoDepthConversion() == true)
  {
    char * velocityField = modelFile->getVelocityField();

    if(strcmp(velocityField,"CONSTANT") == 0)
      velocity = NULL;
    else if(strcmp(velocityField,"FROM_INVERSION")==0)
    {
      velocityFromInversion_ = true;
      velocity = NULL;
    }
    else
    {
      const char * parName = "Velocity";
      char tmpErrText[MAX_STRING];
      sprintf(tmpErrText,"%c",'\0');
      int readerror = 0;
      if(findFileType(velocityField) == SEGYFILE) {
        int nAngles = modelSettings->getNumberOfAngles();
        SegyGeometry ** geometry = new SegyGeometry * [nAngles];
        for (int i =0 ; i < nAngles ; i++)
          geometry[i] = NULL;
        readerror = readSegyFiles((&velocityField), 1,(&velocity),
                                  timeSimbox, geometry, modelSettings,
                                  tmpErrText,0);
        for (int i =0 ; i < nAngles ; i++)
          if (geometry[i] != NULL)
            delete geometry[i];
        delete [] geometry;
      }
      else
        readerror = readStormFile(velocityField, velocity, parName, 
                                  timeSimbox, modelSettings,
                                  tmpErrText);
      if(readerror != 0)
      {
        sprintf(errText,"%sReading of file \'%s\' for parameter \'%s\' failed\n%s\n", 
                errText,velocityField,parName,tmpErrText);
        failed = true;
      }
    }
  }
}

