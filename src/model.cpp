#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

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

#include "lib/random.h"
#include "lib/lib_misc.h"
#include "lib/global_def.h"
#include "lib/segy.h"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"

Model::Model(char * fileName)
{
  modelSettings_         = NULL;
  timeSimbox_            = new Simbox();
  depthSimbox_           = NULL;
  wells_                 = NULL;
  background_            = NULL;
  priorCorrelations_     = NULL;
  priorFacies_           = NULL;
  seisCube_              = NULL;
  wavelet_               = NULL;
  shiftGrids_            = NULL;
  gainGrids_             = NULL;
  waveletEstimInterval_  = NULL;
  faciesEstimInterval_   = NULL; 
  correlationDirection_  = NULL;     
  reflectionMatrix_      = NULL;
  randomGen_             = NULL;
  failed_                = false;

  ModelFile * modelFile = new ModelFile(fileName);
  if (modelFile->getParsingFailed()) {
    failed_ = true;
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
    sprintf(errText, "");

    makeTimeSimbox(timeSimbox_, modelSettings_, modelFile, 
                   errText, failed_);

    if(failed_ == false)
    { 
      if (modelSettings_->getGenerateSeismic() == true)
      {
        processBackground(background_, wells_, timeSimbox_, 
                          modelSettings_, modelFile,
                          errText);
        processReflectionMatrix(reflectionMatrix_, background_, 
                                modelSettings_, modelFile, 
                                errText);  
        if (processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                          timeSimbox_, shiftGrids_, gainGrids_,
                          modelSettings_, modelFile,
                          hasSignalToNoiseRatio_, errText) != 0) {
            failed_ = true;
            LogKit::LogFormatted(LogKit::ERROR,"\nERROR(S) while processing wavelets: \n %s", errText);
            sprintf(errText, "");
        }
        if(modelSettings_->getFormatFlag() == 0)
          modelSettings_->setFormatFlag(1);  //Default, but not initialized due to possible double output.
        background_->getAlpha()->setOutputFormat(modelSettings_->getFormatFlag()); //static, controls all grids.
      }
      else
      {
        bool failedSeismic = false;
        if (processSeismic(seisCube_, timeSimbox_, modelSettings_, modelFile, errText) != 0) {
          failedSeismic = true;
          LogKit::LogFormatted(LogKit::ERROR,"\nERROR(s) while processing seismic data:\n %s", errText);
          sprintf(errText, "");
        }
        makeDepthSimbox(depthSimbox_, modelSettings_, modelFile, 
                        errText, failed_);
        if(failed_ == false)
        { 
          processWells(wells_, timeSimbox_, randomGen_, 
                       modelSettings_, modelFile, 
                       errText);
          processBackground(background_, wells_, timeSimbox_, 
                            modelSettings_, modelFile,
                            errText);
          processPriorCorrelations(priorCorrelations_, background_, wells_, timeSimbox_,
                                   modelSettings_,modelFile,
                                   errText);
          processReflectionMatrix(reflectionMatrix_, background_, 
                                  modelSettings_, modelFile, 
                                  errText);
          if (processWavelets(wavelet_, seisCube_, wells_, reflectionMatrix_,
                          timeSimbox_, shiftGrids_, gainGrids_,
                          modelSettings_, modelFile,
                          hasSignalToNoiseRatio_, errText) != 0) {
              failed_ = true;
              LogKit::LogFormatted(LogKit::ERROR,"\nERROR(s) while processing wavelets: \n%s", errText);
              sprintf(errText, "");
          }
          processPriorFaciesProb(priorFacies_,
                                 wells_,
                                 randomGen_,
                                 modelFile->getTimeNz(),
                                 modelSettings_);
          loadExtraSurfaces(waveletEstimInterval_,
                            faciesEstimInterval_,
                            correlationDirection_,
                            modelFile);
        }
        if (failedSeismic)
          failed_ = true;
      }
    }
    if (failed_) 
      LogKit::LogFormatted(LogKit::ERROR,"\nERROR(s) in input from modelfile. \n %s", errText);
  }
  delete modelFile;
}

Model::~Model(void)
{
  delete randomGen_;
  delete timeSimbox_;

  if(!modelSettings_->getGenerateSeismic()) 
  {
    releaseWells();
  }

  if (depthSimbox_ != NULL)
    delete depthSimbox_;

  if (priorCorrelations_ != NULL)
    delete priorCorrelations_;

  for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
    if(wavelet_[i] != NULL)
      delete wavelet_[i];
  delete [] wavelet_;

  if(shiftGrids_ != NULL)
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(shiftGrids_[i] != NULL)
      {
        delete shiftGrids_[i];
      }
  }

  if(gainGrids_ != NULL)
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(gainGrids_[i] != NULL)
      {
        delete gainGrids_[i];
      }
  }
  delete modelSettings_;
}

void
Model::releaseWells(void)
{
  int i;
  for(i=0;i<modelSettings_->getNumberOfWells();i++)
  {
    if(wells_[i] != NULL)
      delete wells_[i];
  }
  delete [] wells_;
}

void
Model::releaseGrids(void)
{
  delete background_;
  seisCube_ = NULL;
}


// checkFileOpen: checks if the files in fNames table can be opened. Error message is
//                given in errText if not, and error is binary coded numbering of failures.  
int 
Model::checkFileOpen(char ** fNames, int nFiles, const char * command, char * errText, int start,
                     bool details)
{
  if(details == true)
    assert(nFiles<32); //Otherwise, error overflows.
  int i, error = 0;
  int flag = 1;
  int nErr = 0;
  char errFiles[MAX_STRING];
  strcpy(errFiles," ");
  for(i=start;i<nFiles+start;i++)
  {
    //    LogKit::LogFormatted(LogKit::LOW,"CFO: %s\n",fNames[i]);
    if(fNames[i][0] != '*' && fNames[i][0] != '?')
    {
      FILE * file = fopen(fNames[i],"r");
      if(file == 0)
      {
        error+= flag;
        nErr++;
        sprintf(errFiles,"%s%s ", errFiles, fNames[i]);
      }
      else
        fclose(file);
      if(details == true)
        flag *= 2;
    }
  }
  if(error > 0)
  {
    if(nErr == 1)
      sprintf(errText,"Could not open file %s (command %s).\n",errFiles, command);
    else
      sprintf(errText,"Could not open these files:%s (command %s).\n",errFiles, command);
  }
  return(error);
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
        error = 1;
        sprintf(errText,"Found '%s' in file %s, expected a number.\n",
          storage, fileName);
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

  float    megaBytes = static_cast<float>(4*nGrids*gridSize)/(1024.f*1024.f);
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
                     ModelSettings *& modelSettings, 
                     char           * errText)
{
  //  long int timestart, timeend;
  //  time(&timestart);
  char tmpErr[MAX_STRING];
  int error = 0;
  int sbError = 0;
  int okFiles = checkFileOpen(fNames, nFiles, "DUMMY", tmpErr, 0);
  strcpy(errText, "");
  SegY * segy;
  int i, flag = 1;
  for(i=0 ; i<nFiles ; i++)
  {
    target[i] = NULL;
    if((okFiles & flag) == 0)
    {
      segy = new SegY(fNames[i], timeSimbox, modelSettings->getZpad(), modelSettings->getSegyOffset());
      if(segy->checkError(tmpErr) != 0)
      {
        error++;
        sprintf(errText,"%s%s", errText, tmpErr);
      }
      else
      {
        if(timeSimbox->status() == Simbox::NOAREA)
        {
          //
          // Needed in modelSettings (at least) for depth simbox 
          //
          double * ap = new double[7]; 
          segy->getAreaParameters(ap);
          modelSettings->setAreaParameters(ap);
          timeSimbox->setArea(ap);
          delete [] ap;

          //
          // Get information about seismic lines
          //
          int * sl = new int[4]; 
          segy->getSeisLines(sl);
          //modelSettings->setSeisLines(sl);
          timeSimbox->setSeisLines(sl);
          delete [] sl;

          sbError = timeSimbox->checkError(modelSettings->getLzLimit(), errText);

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
            error++;
            sprintf(errText,"%s%s", errText, tmpErr);
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
          target[i]->fillInFromSegY(segy);
          target[i]->setAngle(modelSettings->getAngle()[i]);
        }
      }
      delete segy;
    }
  }
  //  time(&timeend);
  //  LogKit::LogFormatted(LogKit::LOW,"SEGY read and interpolate %i grids in %ld seconds \n",nFiles,timeend-timestart);
  return(error);
}

int
Model::readSegyFile(char          * fName, 
                    FFTGrid       * target, 
                    Simbox        * timeSimbox, 
                    ModelSettings * modelSettings, 
                    char          * errText) 
{
  char tmpErr[MAX_STRING];
  int error = 0;

  SegY * segy;
  target = NULL;
  segy = new SegY(fName, timeSimbox, modelSettings->getZpad(), modelSettings->getSegyOffset());
  if(segy->checkError(tmpErr) != 0)
  {
    sprintf(errText,"%s%s", errText, tmpErr);
    error = 1;
  }
  else
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
    target->fillInFromSegY(segy);
    target->logTransf();
  }
  delete segy;

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
  FILE * file = fopen(fName,"rb");
  float x0, y0, z0, lx, ly, lz, rot;
  int nx, ny, nz, i;
  char temp[MAX_STRING];
  for(i=0;i<5;i++)
    fscanf(file,"%s",temp);
  fscanf(file,"%f %f %f %f",&x0, &lx, &y0, &ly);
  fscanf(file,"%s",temp);
  if(!isNumber(temp))
  {
    z0 = RMISSING; // Dummy setting to avoid warning of uninitialized variable
    error = 1;
    sprintf(errText,"Top surface in %s should be constant value\n", fName);
  }
  else
    z0 = float(atof(temp));
  fscanf(file,"%s",temp);
  if(!isNumber(temp))
  {
    if(error == 1)
      sprintf(errText,"%sThe same goes for the bottom surface.\n",
      errText);
    else
    {
      error = 1;
      sprintf(errText,"Bottom surface in %s should be constant value\n",
        fName);
    }
  }
  fscanf(file,"%s",temp);
  if(!isNumber(temp) || atof(temp) != 0)
  {
    sprintf(errText,"%sErosion should be 0 in %s.\n",errText,fName);
    error = 1;
    fscanf(file,"%s",temp);
  }
  else
  {
    fscanf(file,"%s",temp);
    if(!isNumber(temp) || atof(temp) != 0)
    {
      sprintf(errText,"%sErosion should be 0 in %s.\n",errText,fName);
      error = 1;
    }
  }
  if(error == 0)
  {
    fscanf(file,"%f",&lz);
    fscanf(file,"%f",&rot);
    rot = rot*float(PI/180.0);
    fscanf(file,"%d %d %d",&nx, &ny, &nz);
    float lmax = lx;
    if(ly > lmax)
      lmax = ly;
    lmax *= 1.1f;
    Surface * tmpTop = new Surface(x0-1.5*lmax,y0-1.5*lmax,3*lmax,3*lmax,2,2,z0);
    Simbox * tmpSimbox = new Simbox(x0, y0, tmpTop, lx, ly, lz, rot,
      lx/float(nx), ly/float(ny), lz/float(nz));
    readToEOL(file);
    int nData = nx*ny*nz;
    float * grid = new float[nData];
    fread(grid, 4, nData, file);
#ifndef BIGENDIAN
    swap4Bytes((char * ) grid, nData);
#endif
    int baseInd, i, j, k, k1, kInd, index;
    for(i=0; i<nx; i++) 
      for(j=0;j<ny;j++)
      {
        baseInd = i+j*nx;
        k1 = nz;
        kInd = -1;
        for(k=0;k<nz;k++)
        {
          index = baseInd + k*nx*ny;
          if(grid[index] != WELLMISSING) 
          {
            if(k < k1)
              k1 = k;
            else if(index > kInd)
              kInd = index;
          }
          else if(kInd >= 0)
            grid[index] = grid[kInd];
        }
        if(k1 < nz)
        {
          kInd = baseInd + k1*nx*ny;
          for(k=0;k<k1;k++)
            grid[baseInd+k*nx*ny] = grid[baseInd+k1*nx*ny];
        }
      }

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
      target->fillInFromStorm(tmpSimbox, timeSimbox, grid, parName);
      target->logTransf();
      delete [] grid;
      delete tmpSimbox;
  }    
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
                      bool           & failed)
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
  if(error == 0)
  {
    const char * topname = "toptime.storm";
    const char * botname = "bottime.storm";
    if(!((modelSettings->getOutputFlag() & ModelSettings::NOTIME) > 0))
      timeSimbox->writeTopBotGrids(topname, botname);

    if(modelFile->getNWaveletTransfArgs() > 0 && timeSimbox->getIsConstantThick() == true)
      LogKit::LogFormatted(LogKit::LOW,"\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");

    double * areaParams = modelSettings->getAreaParameters(); 
    estimateZPaddingSize(timeSimbox, modelSettings);   
    if (areaParams != NULL)
    {
      timeSimbox->setArea(areaParams);
      
      error = timeSimbox->checkError(modelSettings->getLzLimit(),errText);

      if(error == 0)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nTime surfaces:\n");
        LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                         timeSimbox->getlz()*timeSimbox->getAvgRelThick(),
                         timeSimbox->getlz()*timeSimbox->getMinRelThick(),
                         timeSimbox->getlz());
        LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                         timeSimbox->getdz()*timeSimbox->getAvgRelThick(),
                         timeSimbox->getdz(),
                         timeSimbox->getdz()*timeSimbox->getMinRelThick());
        
        estimateXYPaddingSizes(timeSimbox, modelSettings);

        //
        // Check if CRAVA has enough memory to run calculation without buffering to disk
        //
        checkAvailableMemory(timeSimbox, modelSettings); 
      }
      else
      {
        if(error == Simbox::INTERNALERROR)
        {
          LogKit::LogFormatted(LogKit::LOW,"ERROR: A problems was encountered for simulation grid\n");
          LogKit::LogFormatted(LogKit::LOW,"       %s\n",errText);
        }
        failed = true;
      }
    }
  }
  else
  {
    timeSimbox->externalFailure();
    failed = true;
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
Model::makeDepthSimbox(Simbox       *& depthSimbox,
                       ModelSettings * modelSettings, 
                       ModelFile     * modelFile,
                       char          * errText,
                       bool          & failed)
{
  if (modelFile->getHasDepthSurfaces())
  {
    depthSimbox = new Simbox();
    int error = 0;
    setSimboxSurfaces(depthSimbox, 
                      modelFile->getDepthSurfFile(), 
                      modelFile->getParallelDepthSurfaces(), 
                      modelFile->getDepthDTop(), 
                      modelFile->getDepthLz(), 
                      modelFile->getDepthDz(), 
                      modelFile->getDepthNz(),
                      error);
    if(error == 0)
    {
      const char * topname = "topdepth.storm";
      const char * botname = "botdepth.storm";
      depthSimbox->writeTopBotGrids(topname, botname);

      double * areaParams = modelSettings->getAreaParameters(); 
      if (areaParams != NULL)
      {
        depthSimbox->setArea(areaParams);
        error = depthSimbox->checkError(modelSettings->getLzLimit(),errText);

        if(error == Simbox::INTERNALERROR)
        {
          LogKit::LogFormatted(LogKit::LOW,"ERROR: A problems was encountered for simulation grid\n");
          LogKit::LogFormatted(LogKit::LOW,"       %s\n",errText);
          failed = true;
        }
      }
      else
      {
        LogKit::LogFormatted(LogKit::ERROR,"ERROR: No area available for depth simbox\n");
        failed = true;
      }
      if (depthSimbox->status() == Simbox::BOXOK)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nDepth surfaces:\n");
        LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                         depthSimbox->getlz()*depthSimbox->getAvgRelThick(),
                         depthSimbox->getlz()*depthSimbox->getMinRelThick(),
                         depthSimbox->getlz());
        LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                         depthSimbox->getdz()*depthSimbox->getAvgRelThick(),
                         depthSimbox->getdz(),
                         depthSimbox->getdz()*depthSimbox->getMinRelThick());
      }
    }
    else
    {
      depthSimbox->externalFailure();
      failed = true;
    }
  }
}

int 
Model::processSeismic(FFTGrid      **& seisCube,
                      Simbox        *& timeSimbox,
                      ModelSettings *& modelSettings, 
                      ModelFile      * modelFile,
                      char           * errText)
{
  char ** seismicFile = modelFile->getSeismicFile();
  int error = 0;

  if(seismicFile != NULL)
  {
    seisCube = new FFTGrid * [modelSettings->getNumberOfAngles()];

    char tmpErrText[MAX_STRING];
    if(readSegyFiles(seismicFile, modelSettings->getNumberOfAngles(), seisCube, 
                     timeSimbox, modelSettings, tmpErrText) == 0)
    {
      int formatFlag = modelSettings->getFormatFlag();
      if(formatFlag == 0)
        formatFlag = 1;   //Default, but not initialized due to possible double output.
      seisCube[0]->setOutputFormat(formatFlag); //static, controls all grids.
      if(modelSettings->getDebugFlag() == 1)
      {
        char sName[100];
        for(int i=0 ; i<modelSettings->getNumberOfAngles() ; i++)
        {
          sprintf(sName, "origSeis%d",i);
          seisCube[i]->writeFile(sName, timeSimbox);
        }
      }
      LogKit::LogFormatted(LogKit::LOW,"\nTime simulation grids:\n");
      LogKit::LogFormatted(LogKit::LOW,"  Output grid         %4i * %4i * %4i   : %10i\n",
                           timeSimbox->getnx(),timeSimbox->getny(),timeSimbox->getnz(),
                           timeSimbox->getnx()*timeSimbox->getny()*timeSimbox->getnz()); 
      LogKit::LogFormatted(LogKit::LOW,"  FFT grid            %4i * %4i * %4i   : %10i\n",
                           modelSettings->getNXpad(),modelSettings->getNYpad(),modelSettings->getNZpad(),
                           modelSettings->getNXpad()*modelSettings->getNYpad()*modelSettings->getNZpad());
    }
    else
    {
      sprintf(errText, "%sERROR: Reading of seismic data files failed:\n %s\n", errText, tmpErrText);
      error = 1;
    }
  }
  return (error);
}

void 
Model::processWells(WellData     **& wells,
                    Simbox         * timeSimbox,
                    RandomGen      * randomGen,
                    ModelSettings *& modelSettings, 
                    ModelFile      * modelFile,
                    char           * errText)
{
  LogKit::LogFormatted(LogKit::LOW,"\nReading well data:\n");

  char ** wellFile       = modelFile->getWellFile();
  char ** headerList     = modelFile->getHeaderList();
  bool    faciesLogGiven = modelFile->getFaciesLogGiven();
  int     nWells         = modelSettings->getNumberOfWells();
  int     nFacies        = 0;

  int error = checkFileOpen(wellFile, nWells, "DUMMY", errText, 0);

  if(error > 0)
  {
    printf("ERROR: %s\n",errText);
    exit(1);
  }
  else
  {
    wells = new WellData *[nWells];
    for(int i=0 ; i<nWells ; i++)
    {
      wells[i] = new WellData(wellFile[i], modelSettings, headerList, faciesLogGiven, i);

      if(wells[i]->checkError(errText) != 0)
      {
        printf("%s\n",errText);
        exit(1);
      }
    }
  }

  if(modelFile->getFaciesLogGiven()) { 
    checkFaciesNames(wells, modelSettings);
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
    if(wells[i]!=NULL)
    {
      if(wells[i]->checkSimbox(timeSimbox) == 1) 
      {
        skip = true;
        nohit++;
      }
      if(wells[i]->getNd() == 0) 
      {
        LogKit::LogFormatted(LogKit::LOW,"  IGNORED (no log entries found)\n");
        skip = true;
        empty++;
      }
      if(wells[i]->isFaciesOk()==0)
      {
        LogKit::LogFormatted(LogKit::LOW,"   IGNORED (facies log has wrong entries)\n");
        skip = true;
        facieslognotok++;
      }
      if(skip)
        validIndex[i] = false;
      else
      {
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
  for(int i=0 ; i<nWells ; i++)
  {
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
          LogKit::LogFormatted(LogKit::LOW,"%12.4f ",faciesCount[i][f],faciesProb);
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
  if (nWells==0)
  {
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
    exit(1);
  }
}

void Model::checkFaciesNames(WellData      ** wells,
                             ModelSettings *& modelSettings)
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
    names[i]   = NULL;
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
          LogKit::LogFormatted(LogKit::LOW,"Error in facies logs. Facies names and numbers are not uniquely defined.\n");
          exit(1);
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
                         char          * errText)
{
  if (modelSettings->getDoInversion() || 
      modelSettings->getGenerateSeismic() || 
      (modelSettings->getOutputFlag() & ModelSettings::BACKGROUND) > 0 )
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
        printf("ERROR: Did not manage to make variogram for background modelling\n");
        exit(1);
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
            sprintf(errText,"%c",'\0');
            int okFiles = checkFileOpen(backFile, 1, "DUMMY", errText, i);
            if(okFiles == 0)
            {
              int readerror = 0;
              if(constBack[i] == ModelFile::SEGYFILE)
                readerror = readSegyFile(backFile[i], backModel[i],
                                         timeSimbox, modelSettings,
                                         errText);
              else
                readerror = readStormFile(backFile[i], backModel[i], parName[i], 
                                          timeSimbox, modelSettings,
                                          errText);
              if(readerror != 0)
              {
                LogKit::LogFormatted(LogKit::LOW,"ERROR: Reading of file \'%s\' for parameter \'%s\' failed\n",
                                     backFile[i],parName[i]);
                printf("       %s\n",errText);
                exit(1);
              }
            }
            else 
            {
              LogKit::LogFormatted(LogKit::LOW,errText);
              printf("       %s\n",errText);
              exit(1);
            }
          }
          else
          {
            LogKit::LogFormatted(LogKit::LOW,"ERROR: Reading of file for parameter \'%s\' failed. File pointer is NULL\n",
                                 parName[i]);
            exit(1);
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
          LogKit::LogFormatted(LogKit::LOW,"ERROR: Could not set background model constBack[%d] == NULL\n",i);
          exit(1);
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
                                char          * errText)
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
      paramCorr = readMatrix(paramCorrFile, 3, 3, "parameter correlation", errText);
      if(paramCorrFile == NULL) 
      {
        LogKit::LogFormatted(LogKit::LOW,"ERROR: Reading of file \'%s\' for parameter correlation matrix failed\n",paramCorrFile);
        LogKit::LogFormatted(LogKit::LOW,"       %s\n",errText);
        exit(1);
      }
      LogKit::LogFormatted(LogKit::LOW,"Parameter correlation read from file.\n\n");
    }

    Surface * CorrXY = getCorrXYGrid();

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
      LogKit::LogFormatted(LogKit::LOW,"\nErrors detected when estimating prior covariance.\nAborting\n");
      exit(1);
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
Model::getCorrXYGrid(void)
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
      transf = new FFTFileGrid((FFTFileGrid *) seisCube[i]); //move new out of loop? Copy grid instead
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
                               char          * errText)
{
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model.
  //
  char * reflMatrFile = modelFile->getReflMatrFile();
  if(reflMatrFile != NULL) {
    reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", errText);
    if(reflectionMatrix == NULL) {
      LogKit::LogFormatted(LogKit::LOW,"ERROR: Reading of file \'%s\' for reflection matrix failed\n",reflMatrFile);
      LogKit::LogFormatted(LogKit::LOW,"       %s\n",errText);
      exit(1);
    }
    LogKit::LogFormatted(LogKit::LOW,"Reflection parameters read from file.\n\n");
  }
  else
  {
    setupDefaultReflectionMatrix(reflectionMatrix,
                                 background,
                                 modelSettings,
                                 modelFile);
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

int 
Model::processWavelets(Wavelet     **& wavelet,
                       FFTGrid      ** seisCube,
                       WellData     ** wells,
                       float        ** reflectionMatrix,
                       Simbox        * timeSimbox,
                       Surface      ** shiftGrids,
                       Surface      ** gainGrids,
                       ModelSettings * modelSettings, 
                       ModelFile     * modelFile,
                       bool          & hasSignalToNoiseRatio,
                       char          * errText)
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
      LogKit::LogFormatted(LogKit::LOW,"\nWells that cannot be used in wavelet generation or noise estimation:");
      LogKit::LogFormatted(LogKit::LOW,"\n  Deviated wells.");
      LogKit::LogFormatted(LogKit::LOW,"\n  Wells with too little data.\n");
    }
    wavelet = new Wavelet *[modelSettings->getNumberOfAngles()];
    
    char ** waveletFile     = modelFile->getWaveletFile();
    float * waveScale       = modelFile->getWaveletScale();
    float * noiseEnergy     = modelSettings->getNoiseEnergy();
    
    for(int i=0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {  
      LogKit::LogFormatted(LogKit::LOW,"\nAngle stack : %.1f deg",modelSettings->getAngle()[i]*180.0/PI);
      if (waveletFile[i][0] == '*') 
      {
        if (timeSimbox->getdz() > 4.0f) { // Require this density for wavelet estimation
          LogKit::LogFormatted(LogKit::LOW,"\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
          LogKit::LogFormatted(LogKit::LOW,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
          LogKit::LogFormatted(LogKit::LOW,"         the number of layers must be increased.                                    \n");
        }
        wavelet[i] = new Wavelet1D(timeSimbox, seisCube[i], wells, modelSettings, reflectionMatrix[i]);
      }
      else
      {
        int fileFormat = getWaveletFileFormat(waveletFile[i]);
        if(fileFormat < 0)
        {
          sprintf(errText, "%sERROR: Unknown file format of file  %s.\n", errText, waveletFile[i]);
          error += 1;
        }
        else {
          if (fileFormat == Wavelet::SGRI)
            wavelet[i] = new Wavelet3D(waveletFile[i], modelSettings, timeSimbox, modelSettings->getAngle()[i]);
          else {
            wavelet[i] = new Wavelet1D(waveletFile[i], modelSettings, fileFormat);
            wavelet[i]->resample(static_cast<float>(timeSimbox->getdz()), timeSimbox->getnz(), 
              modelSettings->getZpad(), modelSettings->getAngle()[i]);
          }
          wavelet[i]->setReflCoeff(reflectionMatrix[i]);
        }
      }
      if (error == 0) {
        if (waveScale[i] != RMISSING)        // If RMISSING we will later scale wavelet to get EmpSN = TheoSN.
          wavelet[i]->scale(waveScale[i]);
        if ((wavelet[i]->getDim() == 3) && !timeSimbox->getIsConstantThick()) {
          sprintf(errText, "%s ERROR: Simbox must have constant thickness when 3D wavelet.\n", errText);
          error += 1;
        }
        if (noiseEnergy[i] == RMISSING)
        {
          if (wavelet[i]->getDim() == 3) { //Not possible to estimate noise variance when 3D wavelet
            sprintf(errText, "%s ERROR: Estimation of noise st. dev. is not possible for 3D wavelet.\n", errText);
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
            sprintf(errText, "%s ERROR: Illegal signal-to-noise ratio of %.3f for cube %d\n", errText, noiseEnergy[i],i);
            sprintf(errText, "%s       Ratio must be in interval 1.0 < SNratio < 10.0\n", errText);
            error += 1;
          }
        }
        if (error == 0) {
          if((modelSettings->getOutputFlag() & ModelSettings::WAVELETS) > 0) 
          {
            char fileName[MAX_STRING];
            sprintf(fileName,"Wavelet_Scaled");
            wavelet[i]->writeWaveletToFile(fileName, 1.0); // dt_max = 1.0;
          }
          if(shiftGrids != NULL && shiftGrids[i] != NULL)
            wavelet[i]->setShiftGrid(shiftGrids[i], timeSimbox);
          if(gainGrids != NULL && gainGrids[i] != NULL)
            wavelet[i]->setGainGrid(gainGrids[i], timeSimbox);
        }
      }
    }
  }
  return (error);
}

int
Model::getWaveletFileFormat(char * fileName)
{
  int fileformat=-1;
  char* dummyStr = new char[MAX_STRING];
  // test for old file format
  FILE* file = fopen(fileName,"r");
  for(int i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      LogKit::LogFormatted(LogKit::LOW,"ERROR: End of file %s is premature\n",fileName);
      exit(1);
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
      LogKit::LogFormatted(LogKit::LOW,"ERROR: End of file %s is premature\n",fileName);
      exit(1);
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
          LogKit::LogFormatted(LogKit::LOW,"ERROR: End of file %s is premature\n",fileName);
          exit(1);  
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
          LogKit::LogFormatted(LogKit::LOW,"%12.4f ",faciesCount[w][i],faciesProb);
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
                         Surface   *& correlationDirection,
                         ModelFile  * modelFile)
{
  bool error = false;
  //
  // Get correlation direction
  //
  if (modelFile->getCorrDirFile() != NULL) {  
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getCorrDirFile());
      correlationDirection = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      LogKit::LogFormatted(LogKit::ERROR,e.what());
      error = true;
    }
  }
  //
  // Get wavelet estimation interval
  //
  if (modelFile->getWaveletEstIntFile() != NULL) {  
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getWaveletEstIntFile()[0]);
      waveletEstimInterval[0] = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      LogKit::LogFormatted(LogKit::ERROR,e.what());
      error = true;
    }
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getWaveletEstIntFile()[1]);
      waveletEstimInterval[1] = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      LogKit::LogFormatted(LogKit::ERROR,e.what());
      error = true;
    }
  }
  //
  // Get facies estimation interval
  //
  if (modelFile->getFaciesEstIntFile() != NULL) {  
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getFaciesEstIntFile()[0]);
      faciesEstimInterval[0] = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      LogKit::LogFormatted(LogKit::ERROR,e.what());
      error = true;
    }
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(modelFile->getFaciesEstIntFile()[1]);
      faciesEstimInterval[1] = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      LogKit::LogFormatted(LogKit::ERROR,e.what());
      error = true;
    }
  }
  if (error)
    LogKit::LogFormatted(LogKit::ERROR,"ERROR loading on or more surfaces\nAborting\n");
}

void
Model::printSettings(ModelSettings * modelSettings,
                     ModelFile     * modelFile,
                     bool            hasSignalToNoiseRatio)
{
  if (modelFile->getCorrDirFile() != NULL)
    LogKit::LogFormatted(LogKit::LOW,"  Correlation direction file: %10s\n",modelFile->getCorrDirFile());
  if (modelFile->getWaveletEstIntFile() != NULL) {
    LogKit::LogFormatted(LogKit::LOW,"  Wavelet estimation file 1:  %10s\n",modelFile->getWaveletEstIntFile()[0]);
    LogKit::LogFormatted(LogKit::LOW,"  Wavelet estimation file 2:  %10s\n",modelFile->getWaveletEstIntFile()[1]);
  }
  if (modelFile->getFaciesEstIntFile() != NULL) {
    LogKit::LogFormatted(LogKit::LOW,"  Facies estimation file 1:   %10s\n",modelFile->getFaciesEstIntFile()[0]);
    LogKit::LogFormatted(LogKit::LOW,"  Facies estimation file 2:   %10s\n",modelFile->getFaciesEstIntFile()[1]);
  }

  LogKit::LogFormatted(LogKit::LOW,"\nGeneral settings:\n");
  int logLevel = modelSettings->getLogLevel();
  std::string logText("*NONE*");
  if (logLevel == LogKit::ERROR)
    logText = "ERROR";
  else if (logLevel == LogKit::WARNING)
    logText = "WARNING";
  else if (logLevel == LogKit::LOW)
    logText = "LOW";
  else if (logLevel == LogKit::MEDIUM)
    logText = "MEDIUM";
  else if (logLevel == LogKit::HIGH)
    logText = "HIGH";
  else if (logLevel == LogKit::DEBUGLOW)
     logText = "DEBUGLOW";
  else if (logLevel == LogKit::DEBUGHIGH)
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
        LogKit::LogFormatted(LogKit::LOW,"  %-2d                                       :   ",i+1);

        if (generateBackground) LogKit::LogFormatted(LogKit::LOW,"  %-6d",modelSettings->getIndicatorBGTrend(i));
        if (estimateWavelet)    LogKit::LogFormatted(LogKit::LOW,"  %-6d",modelSettings->getIndicatorWavelet(i));
        if (estimateFaciesProb) LogKit::LogFormatted(LogKit::LOW,"  %-6d",modelSettings->getIndicatorFacies(i));
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
  if (modelFile->getDepthSurfFile() != NULL)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nDepth surfaces:\n");
    LogKit::LogFormatted(LogKit::LOW,"  Top surface                              : %s\n",   modelFile->getDepthSurfFile()[0]);
    LogKit::LogFormatted(LogKit::LOW,"  Base surface                             : %s\n",   modelFile->getDepthSurfFile()[1]);
    LogKit::LogFormatted(LogKit::LOW,"  Number of layers                         : %10d\n", modelFile->getDepthNz());
    LogKit::LogFormatted(LogKit::LOW,"\nVelocity model:\n");
    LogKit::LogFormatted(LogKit::LOW,"  Constant\n");
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
  //
  // SEISMIC
  //
  if (modelSettings->getGenerateSeismic())
  {
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
    // PRIOR CORRELATION
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

