#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

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
#include "lib/irapgrid.h"
#include "lib/global_def.h"
#include "lib/segy.h"
#include "lib/log.h"

Model::Model(char * fileName)
{
  modelFile_             = new ModelFile(fileName);
  timeSimbox_            = new Simbox();
  depthSimbox_           = NULL;
  wells_                 = NULL;
  background_            = NULL;
  priorCorrelations_     = NULL;
  seisCube_              = NULL;
  wavelet_               = NULL;
  shiftGrids_            = NULL;
  gainGrids_             = NULL;
  reflectionMatrix_      = NULL;
  nxPad_                 =  0;
  nyPad_                 =  0;
  nzPad_                 =  0;
  failed_                = false;

  if (modelFile_->getParsingFailed())
  {
    failed_ = true;
  }
  else
  {
    modelSettings_ = modelFile_->getModelSettings();

    hasSignalToNoiseRatio_ = modelFile_->getHasSignalToNoiseRatio();     // NBNB-PAL: ==> UT

    if(modelSettings_->getDebugFlag() > 0)
    {
      char * fName = LogKit::makeFullFileName("debug",".txt");
      LogKit::setDebug(modelSettings_->getDebugFlag(), fName);
      delete [] fName;
    }

    if(modelFile_->getSeedFile()==NULL)
      randomGen_ = new RandomGen(modelFile_->getSeed());
    else
      randomGen_ = new RandomGen(modelFile_->getSeedFile());

    if(modelSettings_->getNumberOfSimulations() == 0)
      modelSettings_->setOutputFlag(ModelSettings::PREDICTION); //write predicted grids. 
    
    printSettings();
    
    LogKit::writeLog("\n***********************************************************************");
    LogKit::writeLog("\n***                       Reading input data                        ***"); 
    LogKit::writeLog("\n***********************************************************************\n");

    char errText[MAX_STRING];

    makeTimeSimbox(errText);

    if(failed_ == false)
    { 
      if (modelSettings_->getGenerateSeismic() == true)
      {
        processBackground(errText);        
        processReflectionMatrix(errText);  
        processWavelets();        
        if(modelSettings_->getFormatFlag() == 0)
          modelSettings_->setFormatFlag(1);  //Default, but not initialized due to possible double output.
        background_->getAlpha()->setOutputFormat(modelSettings_->getFormatFlag()); //static, controls all grids.
      }
      else
      {
        processSeismic(errText);
        makeDepthSimbox(errText);
        if(failed_ == false)
        { 
          processWells(errText);
          processBackground(errText);        // Needs filtered wells if background are to be generated
          processPriorCorrelations(errText);
          processReflectionMatrix(errText);  // Needs background
          processWavelets();
          processPriorFaciesProb();
        }
      }
    }
  }
  delete modelFile_;
}

Model::~Model()
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
        freeIrapgrid(shiftGrids_[i]);
        free(shiftGrids_[i]);
      }
  }

  if(gainGrids_ != NULL)
  {
    for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
      if(gainGrids_[i] != NULL)
      {
        freeIrapgrid(gainGrids_[i]);
        free(gainGrids_[i]);
      }
  }
  delete modelSettings_;
}

void
Model::releaseWells()
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
Model::releaseGrids()
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
    //    LogKit::writeLog("CFO: %s\n",fNames[i]);
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
  LogKit::writeLog("Reading %s from file %s .... ",readReason, fileName);
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
    LogKit::writeLog("ok.\n");
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
    LogKit::writeLog("failed.\n");
  delete [] tmpRes;
  return(result);
}

void
Model::checkAvailableMemory()
{
  if (modelSettings_->getFileGrid() > 0) // Disk buffering is turn on
    return;

  FFTGrid * dummyGrid = new FFTGrid(timeSimbox_->getnx(), 
                                    timeSimbox_->getny(), 
                                    timeSimbox_->getnz(),
                                    nxPad_, 
                                    nyPad_, 
                                    nzPad_);

  int nGrids, gridSize = dummyGrid->getrsize();
  delete dummyGrid;
  if((modelSettings_->getOutputFlag() & ModelSettings::PREDICTION) == 1)
  {
    nGrids = 10 + modelSettings_->getNumberOfAngles();
    if(modelSettings_->getNumberOfSimulations() > 0 && nGrids < 13)
      nGrids = 13;
  }
  else
    nGrids = 12;
  int workMem = 2500+int( 0.65*gridSize); //Size of memory used beyond grids.

  char * memchunk0 = new char[workMem];
  float **memchunk  = new float*[nGrids];

  int i;
  for(i = 0;i < nGrids ;i++)
    memchunk[i]  = new float[gridSize];

  if(memchunk[nGrids-1] == NULL)  //Could not allocate memory
  {
    modelSettings_->setFileGrid(1);
    LogKit::writeLog("Not enough memory to hold all grids. Using file storage.\n");
  }
  else
  {
    modelSettings_->setFileGrid(0);
  }

  for(i = 0;i < nGrids ;i++)
    if(memchunk[i] != NULL) delete [] memchunk[i];
  if(memchunk != NULL) delete [] memchunk;

  if(memchunk0 != NULL) delete[] memchunk0;
}



// readSegyFiles: reads SegY files form fNames table and stores in target table.
// NB: Also deallocates file names.
int
Model::readSegyFiles(char ** fNames, int nFiles, FFTGrid ** target, char * errText, 
                     int gridType, int start)
{
  //  long int timestart, timeend;
  //  time(&timestart);
  char tmpErr[MAX_STRING];
  int error = 0;
  int sbError = 0;
  int okFiles = checkFileOpen(fNames, nFiles, "DUMMY", tmpErr, start);
  strcpy(errText, "");
  SegY * segy;
  int i, flag = 1;
  for(i=start;i<nFiles+start;i++)
  {
    target[i] = NULL;
    if((okFiles & flag) == 0)
    {
      segy = new SegY(fNames[i], timeSimbox_, modelSettings_->getZpad(), modelSettings_->getSegyOffset());
      if(segy->checkError(tmpErr) != 0)
      {
        error++;
        sprintf(errText,"%s%s", errText, tmpErr);
      }
      else
      {
        if(timeSimbox_->status() == Simbox::NOAREA)
        {
          sbError = segy->completeTimeSimbox(timeSimbox_, modelSettings_->getLzLimit(), tmpErr);
          if(sbError == 0)
          {
            estimateXYPaddingSizes();
            //
            // Check if CRAVA has enough memory to run calculation without buffering to disk
            //
            checkAvailableMemory();
          }
          else if(sbError == Simbox::INTERNALERROR)
          {
            error++;
            sprintf(errText,"%s%s", errText, tmpErr);
          }
        }
        if(sbError == 0)
        {
          if(depthSimbox_ != NULL && depthSimbox_->status() == Simbox::NOAREA)
            segy->completeDepthSimbox(depthSimbox_);

          if(modelSettings_->getFileGrid() == 1)
            target[i] = new FFTFileGrid(timeSimbox_->getnx(), timeSimbox_->getny(), timeSimbox_->getnz(),
                                        nxPad_, nyPad_, nzPad_);
          else
            target[i] = new FFTGrid(timeSimbox_->getnx(), timeSimbox_->getny(), timeSimbox_->getnz(),
                                    nxPad_, nyPad_, nzPad_);

          target[i]->setType(gridType);
          target[i]->fillInFromSegY(segy);
          if(target[i]->getType() == FFTGrid::DATA)
          {
            target[i]->setAngle(modelSettings_->getAngle()[i]);
            //	    char fName[20];
            //	    sprintf(fName, "paddedSeis%i",i);
            //          target[i]->writeStormFile(fName, simbox_,true,true);
          }
          else if(target[i]->getType() == FFTGrid::PARAMETER)
          {
            target[i]->logTransf();
          }
        }
      }
      delete segy;
    }
  }
  //  time(&timeend);
  //  LogKit::writeLog("SEGY read and interpolate %i grids in %ld seconds \n",nFiles,timeend-timestart);
  return(error);
}


//NBNB Following routine only to be used for parameters!
int
Model::readStormFile(char *fName, FFTGrid * & target, const char * parName, char * errText)
{
  int error = 0;
  sprintf(errText,"%c",'\0');
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
    int tmpErr = 0;
    float lmax = lx;
    if(ly > lmax)
      lmax = ly;
    lmax *= 1.1f;
    irapgrid * tmpTop = irapgridConstGrid(x0-1.5*lmax,y0-1.5*lmax,3*lmax,3*lmax,2,2,z0,
      &tmpErr);
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

      if(modelSettings_->getFileGrid() == 1)
        target = new FFTFileGrid(timeSimbox_->getnx(), timeSimbox_->getny(), 
                                 timeSimbox_->getnz(), nxPad_, nyPad_, nzPad_);
      else
        target = new FFTGrid(timeSimbox_->getnx(), timeSimbox_->getny(), 
                             timeSimbox_->getnz(), nxPad_, nyPad_, nzPad_);
      target->setType(FFTGrid::PARAMETER);
      target->fillInFromStorm(tmpSimbox, timeSimbox_, grid, parName);
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
Model::setupDefaultReflectionMatrix()
{
  int i;
  float ** A = new float * [modelSettings_->getNumberOfAngles()];
  double vsvp  = background_->getMeanVsVp();
  double vsvp2 = vsvp*vsvp;
  for(i = 0; i < modelSettings_->getNumberOfAngles(); i++)
  {
    double angle = static_cast<double>(modelSettings_->getAngle()[i]);
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(modelFile_->getSeisType(i) == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float((1.0 +tan2t )/2.0) ; 
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0);
    }
    else if(modelFile_->getSeisType(i) == ModelSettings::PSSEIS) {
      double cost  = cos(angle);
      double cosp  = sqrt(1-vsvp2*sint2);
      double fac   = 0.5*sint/cosp;

      A[i][0] = 0;
      A[i][1] = float(4.0*fac*(vsvp2*sint2+vsvp*cost*cosp));
      A[i][2] = float(fac*(-1.0+2*vsvp2*sint2+2*vsvp*cost*cosp));
    }
  }
  reflectionMatrix_ = A;
  LogKit::writeLog("\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",1.0f/vsvp);
}

float **
Model::getClassicAMatrix()
{
  background_->setClassicVsVp();
  setupDefaultReflectionMatrix();
  return(reflectionMatrix_);
}


void 
Model::makeTimeSimbox(char * errText)
{
  int error = 0;
  setSimboxSurfaces(timeSimbox_, 
                    modelFile_->getTimeSurfFile(), 
                    modelFile_->getParallelTimeSurfaces(), 
                    modelFile_->getTimeDTop(), 
                    modelFile_->getTimeLz(), 
                    modelFile_->getTimeDz(), 
                    modelFile_->getTimeNz(),
                    error);
  if(error == 0)
  {
    const char * topname = "toptime.irap";
    const char * botname = "bottime.irap";
    if(!((modelSettings_->getOutputFlag() & ModelSettings::NOTIME) > 0))
      timeSimbox_->writeTopBotGrids(topname, botname);

    if(modelFile_->getNWaveletTransfArgs() > 0 && timeSimbox_->getIsConstantThick() == true)
      LogKit::writeLog("\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");

    double * areaParams = modelFile_->getAreaParameters(); 
    estimateZPaddingSize();   
    if (areaParams != NULL)
    {
      timeSimbox_->setArea(areaParams[0], areaParams[1], areaParams[2], 
                           areaParams[3], areaParams[4], areaParams[5], 
                           areaParams[6]);
      
      error = timeSimbox_->checkError(modelSettings_->getLzLimit(),errText);

      if(error == 0)
      {
        LogKit::writeLog("\nTime surfaces:\n");
        LogKit::writeLog("  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                         timeSimbox_->getlz()*timeSimbox_->getAvgRelThick(),
                         timeSimbox_->getlz()*timeSimbox_->getMinRelThick(),
                         timeSimbox_->getlz());
        LogKit::writeLog("  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                         timeSimbox_->getdz()*timeSimbox_->getAvgRelThick(),
                         timeSimbox_->getdz(),
                         timeSimbox_->getdz()*timeSimbox_->getMinRelThick());
        
        if (timeSimbox_->getdz() > 4.0f) { // Require this density for wavelet estimation or facies probabilities
          bool estimateWavelet = false;
          for (int i = 0 ; i < modelSettings_->getNumberOfAngles() ; i++)
            if (modelFile_->getWaveletFile()[i][0] == '*')
              estimateWavelet = true;
          if (estimateWavelet) {
            LogKit::writeLog("\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
            LogKit::writeLog("         CRAVA are not reliable and the output results should be treated accordingly.\n");
          }
          if((modelSettings_->getOutputFlag() & ModelSettings::FACIESPROB) > 0) {
            LogKit::writeLog("\nWARNING: The minimum sampling density is lower than 4.0. The FACIES PROBABILITIES\n");
            LogKit::writeLog("         generated by CRAVA are not reliable. To get more reliable probabilities    \n");
            LogKit::writeLog("         the number of layers must be increased.                                    \n");
          }
        }
        estimateXYPaddingSizes();
        //
        // Check if CRAVA has enough memory to run calculation without buffering to disk
        //
        checkAvailableMemory(); 
      }
      else
      {
        if(error == Simbox::INTERNALERROR)
        {
          LogKit::writeLog("ERROR: A problems was encountered for simulation grid\n");
          LogKit::writeLog("       %s\n",errText);
        }
        failed_ = true;
      }
    }
  }
  else
  {
    timeSimbox_->externalFailure();
    failed_ = true;
  }
}

void 
Model::setSimboxSurfaces(Simbox * simbox, char ** surfFile, bool parallelSurfaces, 
                         double dTop, double lz, double dz, int nz,
                         int & error)
{
  irapgrid * z0Grid = irapgridRead(surfFile[0], &error);
  if(error != 0)
  {
    if(error == -3)
      LogKit::writeLog("ERROR: Failed to allocate memory for top surface %s.\n",surfFile[0]);
    else 
      LogKit::writeLog("ERROR: Reading of top surface from file %s failed.\n", surfFile[0]);
    error = 1;
    free(z0Grid);
  }
  else
  {
    if(parallelSurfaces) //Only one reference surface
    {
      simbox->setDepth(z0Grid, dTop, lz, dz);
    }
    else
    {
      irapgrid * z1Grid = irapgridRead(surfFile[1], &error);
      if(error != 0)
      {
        if(error == -3)
          LogKit::writeLog("ERROR: Failed to allocate memory for base surface %s.\n",surfFile[1]);
        else 
          LogKit::writeLog("ERROR: Reading of base surface from file %s failed.\n", surfFile[1]);
        error = 1;

        freeIrapgrid(z0Grid); //Note that the use of z0 and z1 here is correct.
        free(z0Grid);
        free(z1Grid);
      }
      else
      {
        simbox->setDepth(z0Grid, z1Grid, nz);
      }
    }
  }
}

void
Model::estimateXYPaddingSizes()
{
  bool newPaddings = false;
  float xPad    = 0.0f;
  float yPad    = 0.0f;
  float xPadFac = 0.0f;
  float yPadFac = 0.0f;
  if (modelSettings_->getXpad() == 0.0f && modelSettings_->getYpad() == 0.0f)
  {
    float range1  = modelSettings_->getLateralCorr()->getRange();
    float range2  = modelSettings_->getLateralCorr()->getSubRange();
    float angle   = modelSettings_->getLateralCorr()->getAngle();
    float factor  = 1.54f;  // At dist = 1.54*range, an exponential variogram is reduced to 0.01

    xPad          = factor * MAXIM(fabs(range1*cos(angle)),fabs(range2*sin(angle)));
    yPad          = factor * MAXIM(fabs(range1*sin(angle)),fabs(range2*cos(angle)));

    xPadFac       = float(MINIM(1.0f, xPad / timeSimbox_->getlx())); // A padding of more than 100% is insensible
    yPadFac       = float(MINIM(1.0f, yPad / timeSimbox_->getly()));

    modelSettings_->setXpad(xPadFac);
    modelSettings_->setYpad(yPadFac);
    newPaddings = true;
  }
  nxPad_ = setPaddingSize(timeSimbox_->getnx(), modelSettings_->getXpad());
  nyPad_ = setPaddingSize(timeSimbox_->getny(), modelSettings_->getYpad());

  if (newPaddings)
  {
    LogKit::writeLog("\nPadding sizes estimated from lateral correlation ranges:\n");
    LogKit::writeLog("  xPad, xPadFac, nxPad                     : %8.f, %6.2f, %6d\n", xPad, xPadFac, nxPad_);
    LogKit::writeLog("  yPad, yPadFac, nyPad                     : %8.f, %6.2f, %6d\n", yPad, yPadFac, nyPad_);
  }
}

void
Model::estimateZPaddingSize()
{
  bool newPadding = false;
  float zPad    = 0.0f;
  float zPadFac = 0.0f;
  if (modelSettings_->getZpad() == 0.0f)
  {
    float factor  = 1.0f;
    float wLength = 300.0f;                  // Assume a wavelet is approx 300ms.
    zPad          = factor * wLength / 2.0f; // Use one wavelet as padding
    zPadFac       = float(MINIM(1.0f,zPad / (timeSimbox_->getlz()*timeSimbox_->getMinRelThick())));
    
    modelSettings_->setZpad(zPadFac);
    newPadding = true;
  }
  nzPad_ = setPaddingSize(timeSimbox_->getnz(), modelSettings_->getZpad());

  if (newPadding)
  {
    LogKit::writeLog("\nPadding sizes estimated from an assumed wavelet length:\n");
    LogKit::writeLog("  zPad, zPadFac, nzPad                     : %8.f, %6.2f, %6d\n", zPad, zPadFac, nzPad_);
  }
}

void 
Model::makeDepthSimbox(char * errText)
{
  if (modelFile_->getHasDepthSurfaces())
  {
    depthSimbox_ = new Simbox();
    int error = 0;
    setSimboxSurfaces(depthSimbox_, 
                      modelFile_->getDepthSurfFile(), 
                      modelFile_->getParallelDepthSurfaces(), 
                      modelFile_->getDepthDTop(), 
                      modelFile_->getDepthLz(), 
                      modelFile_->getDepthDz(), 
                      modelFile_->getDepthNz(),
                      error);
    if(error == 0)
    {
      const char * topname = "topdepth.irap";
      const char * botname = "botdepth.irap";
      depthSimbox_->writeTopBotGrids(topname, botname);

      double * areaParams = modelFile_->getAreaParameters(); 
      if (areaParams != NULL)
      {
        depthSimbox_->setArea(areaParams[0], areaParams[1], areaParams[2], 
                              areaParams[3], areaParams[4], areaParams[5], 
                              areaParams[6]);
        
        error = depthSimbox_->checkError(modelSettings_->getLzLimit(),errText);

        if(error == Simbox::INTERNALERROR)
        {
          LogKit::writeLog("ERROR: A problems was encountered for simulation grid\n");
          LogKit::writeLog("       %s\n",errText);
          failed_ = true;
        }
      }
      if (depthSimbox_->status() == Simbox::BOXOK)
      {
        LogKit::writeLog("\nDepth surfaces:\n");
        LogKit::writeLog("  Interval thickness    avg / min / max    : %6.1f /%6.1f /%6.1f\n", 
                         depthSimbox_->getlz()*depthSimbox_->getAvgRelThick(),
                         depthSimbox_->getlz()*depthSimbox_->getMinRelThick(),
                         depthSimbox_->getlz());
        LogKit::writeLog("  Sampling density      avg / min / max    : %6.2f /%6.2f /%6.2f\n", 
                         depthSimbox_->getdz()*depthSimbox_->getAvgRelThick(),
                         depthSimbox_->getdz(),
                         depthSimbox_->getdz()*depthSimbox_->getMinRelThick());
      }
    }
    else
    {
      depthSimbox_->externalFailure();
      failed_ = true;
    }
  }
}

void 
Model::processSeismic(char * errText)
{
  char ** seismicFile = modelFile_->getSeismicFile();

  if(seismicFile != NULL)
  {
    seisCube_ = new FFTGrid * [modelSettings_->getNumberOfAngles()];

    if(readSegyFiles(seismicFile, modelSettings_->getNumberOfAngles(), seisCube_, errText, FFTGrid::DATA) == 0)
    {
      int formatFlag = modelSettings_->getFormatFlag();
      if(formatFlag == 0)
        formatFlag = 1;   //Default, but not initialized due to possible double output.
      seisCube_[0]->setOutputFormat(formatFlag); //static, controls all grids.
      if(modelSettings_->getDebugFlag() == 1)
      {
        char sName[100];
        for(int i=0 ; i<modelSettings_->getNumberOfAngles() ; i++)
        {
          sprintf(sName, "origSeis%d",i);
          seisCube_[i]->writeFile(sName, timeSimbox_);
        }
      }
      LogKit::writeLog("\nTime simulation grids:\n");
      LogKit::writeLog("  Output grid         %4i * %4i * %4i   : %10i\n",
                       timeSimbox_->getnx(),timeSimbox_->getny(),timeSimbox_->getnz(),
                       timeSimbox_->getnx()*timeSimbox_->getny()*timeSimbox_->getnz()); 
      LogKit::writeLog("  FFT grid            %4i * %4i * %4i   : %10i\n",
                       nxPad_,nyPad_,nzPad_,nxPad_*nyPad_*nzPad_);
    }
    else
    {
      printf("ERROR: Reading of seismic data files failed:\n");
      printf("       %s\n",errText);
      exit(1);
    }
  }
}

void 
Model::processWells(char * errText)
{
  LogKit::writeLog("\nReading well data:\n");

  char ** wellFile   = modelFile_->getWellFile();
  char ** headerList = modelFile_->getHeaderList();
  int     nWells     = modelSettings_->getNumberOfWells();

  int error = checkFileOpen(wellFile, nWells, "DUMMY", errText, 0);

  if(error > 0)
  {
    printf("ERROR: %s\n",errText);
    exit(1);
  }
  else
  {
    wells_ = new WellData *[nWells];
    for(int i=0 ; i<nWells ; i++)
    {
      wells_[i] = new WellData(wellFile[i], modelSettings_, headerList, modelFile_->getFaciesLogGiven());

      if(wells_[i]->checkError(errText) != 0)
      {
        printf("%s\n",errText);
        exit(1);
      }
    }
  }

  if(modelFile_->getFaciesLogGiven()) 
    checkFaciesNames();

  LogKit::writeLog("\n***********************************************************************");
  LogKit::writeLog("\n***                       Processing Wells                          ***"); 
  LogKit::writeLog("\n***********************************************************************\n\n");
  
  int   * validWells    = new int[nWells];
  bool  * validIndex    = new bool[nWells];
  int   * nMerges       = new int[nWells];
  int   * nInvalidAlpha = new int[nWells];
  int   * nInvalidBeta  = new int[nWells];
  int   * nInvalidRho   = new int[nWells];
  float * rankCorr      = new float[nWells];
  float * devAngle      = new float[nWells];
  
  int count = 0;
  int nohit=0;
  int empty=0;
  int facieslognotok = 0;
  for (int i=0 ; i<nWells ; i++)
  {
    bool skip = false;
    LogKit::writeLog("%s : \n",wells_[i]->getWellname());
    if(wells_[i]!=NULL)
    {
      if(wells_[i]->checkSimbox(timeSimbox_) == 1) 
      {
        skip = true;
        nohit++;
      }
      if(wells_[i]->getNd() == 0) 
      {
        LogKit::writeLog("  IGNORED (no log entries found)\n");
        skip = true;
        empty++;
      }
      if(wells_[i]->isFaciesOk()==0)
      {
        LogKit::writeLog("   IGNORED (facies log has wrong entries)\n");
        skip = true;
        facieslognotok++;
      }
      if(skip)
        validIndex[i] = false;
      else
      {
        validIndex[i] = true;
        wells_[i]->removeDuplicateLogEntries(nMerges[i]);
        wells_[i]->setWrongLogEntriesUndefined(nInvalidAlpha[i], nInvalidBeta[i], nInvalidRho[i]);
        wells_[i]->filterLogs();
        wells_[i]->lookForSyntheticVsLog(rankCorr[i]);
        wells_[i]->calculateDeviation(devAngle[i]);
        wells_[i]->setBlockedLogsPropThick( new BlockedLogs(wells_[i], 
                                                            timeSimbox_,
                                                            randomGen_) );
        if((modelSettings_->getOutputFlag() & ModelSettings::WELLS) > 0) 
          wells_[i]->writeRMSWell();
        if(modelSettings_->getDebugFlag() > 0)
          wells_[i]->getBlockedLogsPropThick()->writeToFile(static_cast<float>(timeSimbox_->getdz()), 2, true); // 2 = BG logs
        
        validWells[count] = i;
        count++;      
      }
    }
  }
  //
  // Write summary.
  //
  LogKit::writeLog("\n");
  LogKit::writeLog("                                      Invalid                                    \n");
  LogKit::writeLog("Well                    Merges      Vp   Vs  Rho  synthVs/Corr    Deviated/Angle \n");
  LogKit::writeLog("---------------------------------------------------------------------------------\n");
  for(int i=0 ; i<nWells ; i++)
  {
    if (validIndex[i]) 
      LogKit::writeLog("%-23s %6d    %4d %4d %4d     %3s / %5.3f      %3s / %4.1f\n",
                       wells_[i]->getWellname(),
                       nMerges[i],
                       nInvalidAlpha[i], 
                       nInvalidBeta[i], 
                       nInvalidRho[i],
                       (rankCorr[i] > modelSettings_->getMaxRankCorr() ? "yes" : " no"),
                       rankCorr[i],
                       (devAngle[i] > modelSettings_->getMaxDevAngle() ? "yes" : " no"),
                       devAngle[i]);
     else  
       LogKit::writeLog("%-23s      -       -    -    -       - /     -       -  /    -\n",wells_[i]->getWellname());
  }

  //
  // Remove invalid wells
  //
  for(int i=0 ; i<nWells ; i++)
    if (!validIndex[i]) 
      delete wells_[i];
  for(int i=0 ; i<count ; i++)
    wells_[i] = wells_[validWells[i]];
  for(int i=count ; i<nWells ; i++)
    wells_[i] = NULL;
  nWells = count;
  modelSettings_->setNumberOfWells(nWells);

  delete [] validWells;
  delete [] validIndex;
  delete [] nMerges;
  delete [] nInvalidAlpha;
  delete [] nInvalidBeta;
  delete [] nInvalidRho;
  delete [] rankCorr;
  delete [] devAngle;
  
  if (nohit>0)
    LogKit::writeLog("\nWARNING: %d well(s) do not hit the inversion volume and will be ignored.\n",nohit);
  if (empty>0)
    LogKit::writeLog("\nWARNING: %d well(s) contain no log entries and will be ignored.\n",empty);
  if(facieslognotok>0)
    LogKit::writeLog("\nWARNING: %d well(s) have wrong facies logs and will be ignored.\n",facieslognotok);
  if (nWells==0)
  {
    LogKit::writeLog("\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
    LogKit::writeLog("\n       below is correct. If it is not, you probably have problems with coordinate scaling. To");
    LogKit::writeLog("\n       correct this problem the seismic data must be re-exported from your data base using");
    LogKit::writeLog("\n       the \'by-pass coordinate scaling (byte 71-72)\' toggle that matches your version of");
    LogKit::writeLog("\n       CRAVA (see beginning of CRAVA log file). Alternatively, CRAVA can be recompiled.\n");
    LogKit::writeLog("\n                                   X0          Y0        DeltaX      DeltaY      Angle");
    LogKit::writeLog("\n       -------------------------------------------------------------------------------");
    LogKit::writeLog("\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n", 
                     timeSimbox_->getx0(), timeSimbox_->gety0(), 
                     timeSimbox_->getlx(), timeSimbox_->getly(), 
                     (timeSimbox_->getAngle()*180)/PI);
    LogKit::writeLog("\nAborting\n");
    exit(1);
  }
}

void Model::checkFaciesNames(void)
{
  int i, j, w;
  int min,max, globalmin = 0,globalmax = 0;
  bool first = true;
  for (w = 0; w < modelSettings_->getNumberOfWells(); w++) {
    if(wells_[w]->isFaciesLogDefined())
    {
      wells_[w]->getMinMaxFnr(min,max);
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

  int     nnames = globalmax-globalmin+1;
  char ** names  = new char*[nnames];

  for (i =0 ; i < nnames ; i++)
    names[i] = NULL;
  
  for(w=0 ; w<modelSettings_->getNumberOfWells() ; w++)
  {
    if(wells_[w]->isFaciesLogDefined())
    {
      for(i=0 ; i < wells_[w]->getNFacies() ; i++)
      {
        char * name = wells_[w]->getFaciesName(i);
        int    fnr  = wells_[w]->getFaciesNr(i) - globalmin;

        if(names[fnr] == NULL)
          names[fnr] = name;
        else if(strcmp(names[fnr],name) != 0)
        {
          LogKit::writeLog("Error in facies logs. Facies names and numbers are not uniquely defined.\n");
          exit(1);
        }
      }
    }
  }

  LogKit::writeLog("\nFaciesNumber      FaciesName            ");
  LogKit::writeLog("\n--------------------------------------\n");
  for(i=0 ; i<nnames ; i++)
    if(names[i] != NULL) 
      LogKit::writeLog("     %2d           %-20s\n",i+globalmin,names[i]);


  int nFacies = 0;
  for(i=0 ; i<nnames ; i++)
    if(names[i] != NULL) 
      nFacies++;

  char ** faciesNames = new char*[nFacies];

  j = -1;
  for(i=0 ; i<nnames ; i++)
  {
    if(names[i] != NULL)
    {
      j++;
      faciesNames[j] = names[i];
    }
  }
  modelSettings_->setNumberOfFacies(nFacies);
  modelSettings_->setFaciesNames(faciesNames,nFacies);

  delete [] faciesNames;
  delete [] names;
}

void 
Model::processBackground(char * errText)
{
  if (modelSettings_->getDoInversion() || 
      modelSettings_->getGenerateSeismic() || 
      (modelSettings_->getOutputFlag() & ModelSettings::BACKGROUND) > 0 )
  {
    FFTGrid * backModel[3];
    if (modelFile_->getGenerateBackground()) 
    {
      LogKit::writeLog("\n***********************************************************************");
      LogKit::writeLog("\n***              Prior Expectations / Background Model              ***"); 
      LogKit::writeLog("\n***********************************************************************\n");
      
      if(modelSettings_->getBackgroundVario() == NULL)
      {
        printf("ERROR: Did not manage to make variogram for background modelling\n");
        exit(1);
      }
      
      for (int i=0 ; i<3 ; i++)
      {
        backModel[i] = new FFTGrid(timeSimbox_->getnx(), timeSimbox_->getny(), timeSimbox_->getnz(), 
                                   nxPad_, nyPad_, nzPad_);              
        backModel[i]->setType(FFTGrid::PARAMETER);
      }
      
      background_ = new Background(backModel, wells_, timeSimbox_, modelSettings_);
    }
    else 
    {
      const char * parName[]={"Vp background","Vs background","Rho background"};
      for(int i=0 ; i<3 ; i++)
      {
        float  * constBack = modelFile_->getConstBack();
        char  ** backFile  = modelFile_->getBackFile();
        
        if(constBack[i] < 0)
        {
          if(backFile[i] != NULL)
          {
            int readerror = 0;
            if(constBack[i] == ModelFile::SEGYFILE)
              readerror = readSegyFiles(backFile, 1, backModel, errText, FFTGrid::PARAMETER, i);
            else
              readerror = readStormFile(backFile[i], backModel[i], parName[i], errText);
            if(readerror != 0)
            {
              LogKit::writeLog("ERROR: Reading of file \'%s\' for parameter \'%s\' failed\n",
                               backFile[i],parName[i]);
              printf("       %s\n",errText);
              exit(1);
            }
          }
          else
          {
            LogKit::writeLog("ERROR: Reading of file for parameter \'%s\' failed. File pointer is NULL\n",
                             parName[i]);
            exit(1);
          }
        }
        else if(constBack[i] > 0)
        {
          if(modelSettings_->getFileGrid() == 1)
            backModel[i] = new FFTFileGrid(timeSimbox_->getnx(), timeSimbox_->getny(), timeSimbox_->getnz(), 
                                           nxPad_, nyPad_, nzPad_);
          else
            backModel[i] = new FFTGrid(timeSimbox_->getnx(), timeSimbox_->getny(), timeSimbox_->getnz(), 
                                       nxPad_, nyPad_, nzPad_);              
          backModel[i]->setType(FFTGrid::PARAMETER);
          backModel[i]->fillInConstant(float( log( constBack[i] )));
        }
        else
        {
          LogKit::writeLog("ERROR: Could not set background model constBack[%d] == NULL\n",i);
          exit(1);
        }
      }
      background_ = new Background(backModel);
    }
    if((modelSettings_->getOutputFlag() & ModelSettings::BACKGROUND) > 0)
      background_->writeBackgrounds(timeSimbox_); 
  }
}

void 
Model::processPriorCorrelations(char * errText)
{
  bool printResult = (modelSettings_->getOutputFlag() & (ModelSettings::PRIORCORRELATIONS + ModelSettings::CORRELATION)) > 0;
  if (modelSettings_->getDoInversion() || printResult)
  {
    LogKit::writeLog("\n***********************************************************************");
    LogKit::writeLog("\n***                        Prior Covariance                         ***"); 
    LogKit::writeLog("\n***********************************************************************\n");
    time_t timestart, timeend;
    time(&timestart);

    //
    // Parameter correlation can be set in model file.
    // Default NULL, indicating that estimate will be used.
    //
    float ** paramCorr = NULL;
    char   * paramCorrFile = modelFile_->getParamCorrFile();
    if(paramCorrFile != NULL) 
    {
      paramCorr = readMatrix(paramCorrFile, 3, 3, "parameter correlation", errText);
      if(paramCorrFile == NULL) 
      {
        LogKit::writeLog("ERROR: Reading of file \'%s\' for parameter correlation matrix failed\n",paramCorrFile);
        LogKit::writeLog("       %s\n",errText);
        exit(1);
      }
      LogKit::writeLog("Parameter correlation read from file.\n\n");
    }

    irapgrid * CorrXY = getCorrXYGrid();

    if(CorrXY->grid == NULL) { // NBNB-PAL: this will never be true (default lateral corr)
      estimateCorrXYFromSeismic(CorrXY);
      time(&timeend);
      LogKit::writeLog("\nEstimate parameter lateral correlation from seismic in %d seconds.\n",
                       static_cast<int>(timeend-timestart));
    }

    Analyzelog * analyze = new Analyzelog(wells_, 
                                          background_,
                                          timeSimbox_, 
                                          modelSettings_);

    priorCorrelations_ = new Corr(analyze->getPointVar0(), 
                                  analyze->getVar0(), 
                                  analyze->getCorrT(), 
                                  analyze->getNumberOfLags(),
                                  static_cast<float>(timeSimbox_->getdz()), 
                                  CorrXY);
    delete analyze;



    if(priorCorrelations_ == NULL)
    {
      LogKit::writeLog("\nErrors detected when estimating prior covariance.\nAborting\n");
      exit(1);
    }

    //
    // NBNB-PAL: 
    //
    // CorrT   (PriorCorrT)
    // CorrXY  (PriorCorrXY)
    // Var0    (PriorVar0)
    //
    // may be read from file.... currently only Var0_ may be read.
    //
    //
    if(paramCorr != NULL)
      priorCorrelations_->setVar0(paramCorr);

    if(printResult)
      priorCorrelations_->dumpResult();
    
    priorCorrelations_->printVariancesToScreen();

    time(&timeend);
    LogKit::writeDebugLog("\n\nTime elapsed :  %d\n",timeend-timestart);  
  }
}

irapgrid * 
Model::getCorrXYGrid()
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
  struct irapgrid *rIrap;
  rIrap = (struct irapgrid *) calloc(1, sizeof(struct irapgrid));
  rIrap->xmin = 0;
  rIrap->xmax = dx*nx;
  rIrap->ymin = 0;
  rIrap->ymax = dy*ny;
  rIrap->xinc = dx;
  rIrap->yinc = dy;
  rIrap->nx = nx;
  rIrap->ny = ny;
  rIrap->bin = 0;
  rIrap->missingcode = IRAPMISSING;
  rIrap->constValue = IRAPMISSING;
  rIrap->filename = 0;
  if(modelSettings_->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
  {
    rIrap->grid = (double *) malloc(((unsigned) npix)*sizeof(double));
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
        rIrap->grid[j*nx+i] = modelSettings_->getLateralCorr()->corr(refi*dx, refj*dy);
      }
    }
  }
  else
  {
    rIrap->grid = NULL;
  }
  return(rIrap);
}

void 
Model::estimateCorrXYFromSeismic(irapgrid * CorrXY)
{
  FFTGrid * transf;
  float   * grid;

  int n = CorrXY->nx*CorrXY->ny;
  grid = new float[n];

  for(int i=0 ; i<n ; i++)
    grid[i] = 0.0;

  for(int i=0;i<modelSettings_->getNumberOfAngles();i++)
  {
    if(seisCube_[i]->isFile())
      transf = new FFTFileGrid((FFTFileGrid *) seisCube_[i]); //move new out of loop? Copy grid instead
    else
      transf = new FFTGrid(seisCube_[i]); //move new out of loop? Copy grid instead

    transf->setAccessMode(FFTGrid::RANDOMACCESS);
    transf->fftInPlace();
    transf->square();
    transf->invFFTInPlace();
    transf->collapseAndAdd( grid ); //the result of the collapse (the result for z=0) is is added to grid
    transf->endAccess();
    delete transf;
  }
  float sill = grid[0];
  CorrXY->grid = (double *) malloc(((unsigned) n)*sizeof(double));
  for(int i=0;i<n;i++)
    CorrXY->grid[i] = grid[i]/sill;
  delete [] grid;

  char * fName= LogKit::makeFullFileName("PriorCorrXY_EstimatedFromSeismic",".irap");
  FILE * file = fopen(fName, "w");
  irapgridWritept(file,CorrXY);
  fclose(file);
  delete [] fName;
}

void 
Model::processReflectionMatrix(char * errText)
{
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model.
  //
  char * reflMatrFile = modelFile_->getReflMatrFile();
  if(reflMatrFile != NULL) {
    reflectionMatrix_ = readMatrix(reflMatrFile, modelSettings_->getNumberOfAngles(), 3, "reflection matrix", errText);
    if(reflectionMatrix_ == NULL) {
      LogKit::writeLog("ERROR: Reading of file \'%s\' for reflection matrix failed\n",reflMatrFile);
      LogKit::writeLog("       %s\n",errText);
      exit(1);
    }
    LogKit::writeLog("Reflection parameters read from file.\n\n");
  }
  else
  {
    setupDefaultReflectionMatrix();
  }
}

void 
Model::processWavelets(void)
{
  if (modelSettings_->getDoInversion() || 
      modelSettings_->getGenerateSeismic() || 
      (modelSettings_->getOutputFlag() & ModelSettings::WAVELETS) > 0 )
  {
    LogKit::writeLog("\n***********************************************************************");
    LogKit::writeLog("\n***                 Processing/generating wavelets                  ***"); 
    LogKit::writeLog("\n***********************************************************************\n");
    bool estimateStuff = false;
    for(int i=0 ; i < modelSettings_->getNumberOfAngles() ; i++)
    {  
      estimateStuff = estimateStuff || (modelFile_->getWaveletFile()[i][0] == '*'); 
      estimateStuff = estimateStuff || (modelSettings_->getNoiseEnergy()[i]==RMISSING); 
    }
    if (estimateStuff) 
    {
      LogKit::writeLog("\nWells that cannot be used in wavelet generation or noise estimation:");
      LogKit::writeLog("\n  Deviated wells.");
      LogKit::writeLog("\n  Wells with too little data.\n");
    }
    wavelet_ = new Wavelet *[modelSettings_->getNumberOfAngles()];
    
    char ** waveletFile     = modelFile_->getWaveletFile();
    float * waveScale       = modelFile_->getWaveletScale();
    float * noiseEnergy     = modelSettings_->getNoiseEnergy();
    
    for(int i=0 ; i < modelSettings_->getNumberOfAngles() ; i++)
    {  
      LogKit::writeLog("\nAngle stack : %.1f deg",modelSettings_->getAngle()[i]*180.0/PI);
      if (waveletFile[i][0] == '*') 
        wavelet_[i] = new Wavelet1D(timeSimbox_, seisCube_[i], wells_, modelSettings_, reflectionMatrix_[i]);
      else
      {
        int fileFormat = getWaveletFileFormat(waveletFile[i]);
        if(fileFormat < 0)
        {
          LogKit::writeLog("ERROR: Unknown file format of file  %s.\n",waveletFile[i]);
          exit(1);
        }
        if (fileFormat == Wavelet::SGRI)
          wavelet_[i] = new Wavelet3D(waveletFile[i], modelSettings_, timeSimbox_, modelSettings_->getAngle()[i]);
        else {
          wavelet_[i] = new Wavelet1D(waveletFile[i], modelSettings_, fileFormat);
          wavelet_[i]->resample(static_cast<float>(timeSimbox_->getdz()), timeSimbox_->getnz(), 
                                modelSettings_->getZpad(), modelSettings_->getAngle()[i]);
        }
        wavelet_[i]->setReflCoeff(reflectionMatrix_[i]);
      }
      if (waveScale[i] != RMISSING)        // If RMISSING we will later scale wavelet to get EmpSN = TheoSN.
        wavelet_[i]->scale(waveScale[i]);
      if (noiseEnergy[i] == RMISSING)
      {
        noiseEnergy[i] = wavelet_[i]->getNoiseStandardDeviation(timeSimbox_, seisCube_[i], wells_, modelSettings_->getNumberOfWells());
        if (hasSignalToNoiseRatio_) 
        {
          //LogKit::writeLog("\n  Since seismic noise is estimated directly, keyword GIVESIGNALTONOISERATIO has no effect.\n");
          hasSignalToNoiseRatio_ = false;
        }
      }
      else
      {
        if (hasSignalToNoiseRatio_ && (noiseEnergy[i] <= 1.0 || noiseEnergy[i] > 10.0))
        {
          LogKit::writeLog("ERROR: Illegal signal-to-noise ratio of %.3f for cube %d\n",noiseEnergy[i],i);
          LogKit::writeLog("       Ratio must be in interval 1.0 < SNratio < 10.0\n");
          exit(1);
        }
      }
      if((modelSettings_->getOutputFlag() & ModelSettings::WAVELETS) > 0) 
      {
        char fileName[MAX_STRING];
        sprintf(fileName,"Wavelet_Scaled");
        wavelet_[i]->writeWaveletToFile(fileName, 1.0); // dt_max = 1.0;
      }
      if(shiftGrids_ != NULL && shiftGrids_[i] != NULL)
        wavelet_[i]->setShiftGrid(shiftGrids_[i], timeSimbox_);
      if(gainGrids_ != NULL && gainGrids_[i] != NULL)
        wavelet_[i]->setGainGrid(gainGrids_[i], timeSimbox_);
    }
  }
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
      LogKit::writeLog("ERROR: End of file %s is premature\n",fileName);
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
      LogKit::writeLog("ERROR: End of file %s is premature\n",fileName);
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
          LogKit::writeLog("ERROR: End of file %s is premature\n",fileName);
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

void
Model::printSettings(void)
{
  LogKit::writeLog("\nGeneral settings:\n");
  if (modelSettings_->getNumberOfSimulations() == 0)
  {
    LogKit::writeLog("  Modelling mode                           : prediction\n");
  }
  else
  {
    LogKit::writeLog("  Modelling mode                           : simulation\n");
    if(modelFile_->getSeedFile()==NULL) {
      if (modelFile_->getSeed() == 0)
        LogKit::writeLog("  Seed                                     :          0 (default seed)\n");
      else
        LogKit::writeLog("  Seed                                     : %10d\n",modelFile_->getSeed());
    }
    else
      LogKit::writeLog("  Seed read from file                    : %10s\n",modelFile_->getSeedFile());
    LogKit::writeLog("  Number of realisations                   : %10d\n",modelSettings_->getNumberOfSimulations());
  }
  LogKit::writeLog("  Kriging                                  : %10s\n",(modelSettings_->getKrigingParameters()==NULL ? "no" : "yes"));

  LogKit::writeLog("\nUnit setting/assumptions:\n");
  LogKit::writeLog("  Time                                     : %10s\n","ms TWT");
  LogKit::writeLog("  Frequency                                : %10s\n","Hz");
  LogKit::writeLog("  Length                                   : %10s\n","m");
  LogKit::writeLog("  Velocities                               : %10s\n","m/s");
  LogKit::writeLog("  Density                                  : %10s\n","kg/m3");
  LogKit::writeLog("  Angle                                    : %10s\n","degrees");

  LogKit::writeLog("\nSettings for well processing:\n");
  LogKit::writeLog("  Threshold for merging log entries        : %10.2f ms\n",modelSettings_->getMaxMergeDist());
  LogKit::writeLog("  Threshold for Vp-Vs rank correlation     : %10.2f\n",modelSettings_->getMaxRankCorr());
  LogKit::writeLog("  Threshold for deviation angle            : %10.1f (=%.2fm/ms TWT)\n",
                   modelSettings_->getMaxDevAngle(),tan(modelSettings_->getMaxDevAngle()*PI/180.0));
  LogKit::writeLog("  High cut for background modelling        : %10.1f\n",modelSettings_->getMaxHzBackground());
  LogKit::writeLog("  High cut for seismic resolution          : %10.1f\n",modelSettings_->getMaxHzSeismic());
  LogKit::writeLog("\nRange of allowed parameter values:\n");
  LogKit::writeLog("  Vp  - min                                : %10.0f\n",modelSettings_->getAlphaMin());
  LogKit::writeLog("  Vp  - max                                : %10.0f\n",modelSettings_->getAlphaMax());
  LogKit::writeLog("  Vs  - min                                : %10.0f\n",modelSettings_->getBetaMin());
  LogKit::writeLog("  Vs  - max                                : %10.0f\n",modelSettings_->getBetaMax());
  LogKit::writeLog("  Rho - min                                : %10.1f\n",modelSettings_->getRhoMin());
  LogKit::writeLog("  Rho - max                                : %10.1f\n",modelSettings_->getRhoMax());  

  //
  // WELL DATA
  //
  if (modelSettings_->getNumberOfWells() > 0)
  {
    LogKit::writeLog("\nWell logs:\n");
    char ** headerList = modelFile_->getHeaderList();
    if (headerList != NULL)
    {
      LogKit::writeLog("  Time                                     : %10s\n",  headerList[0]);
      if(strcmp(uppercase(headerList[1]),"VP"    )==0 ||
         strcmp(uppercase(headerList[1]),"LFP_VP")==0)
        LogKit::writeLog("  p-wave velocity                          : %10s\n",headerList[1]);
      else
        LogKit::writeLog("  Sonic                                    : %10s\n",headerList[1]);
      if(strcmp(uppercase(headerList[3]),"VS"    )==0 || 
         strcmp(uppercase(headerList[3]),"LFP_VS")==0)
        LogKit::writeLog("  s-wave velocity                          : %10s\n",headerList[3]);
      else
        LogKit::writeLog("  Shear sonic                              : %10s\n",headerList[3]);
      LogKit::writeLog("  Density                                  : %10s\n",  headerList[2]);
      if (modelFile_->getFaciesLogGiven())
        LogKit::writeLog("  Facies                                   : %10s\n",headerList[4]);
    }
    else
    {
      LogKit::writeLog("  Time                                     : %10s\n","TWT");
      LogKit::writeLog("  Sonic                                    : %10s\n","DT");
      LogKit::writeLog("  Shear sonic                              : %10s\n","DTS");
      LogKit::writeLog("  Density                                  : %10s\n","RHOB");
      LogKit::writeLog("  Facies                                   : %10s\n","FACIES");
    }
    LogKit::writeLog("\nWell files:\n");
    for (int i = 0 ; i < modelSettings_->getNumberOfWells() ; i++) 
    {
      LogKit::writeLog("  %-2d                                       : %s\n",i+1,modelFile_->getWellFile()[i]);
    }
  }

  //
  // SURFACES
  // 
  LogKit::writeLog("\nTime surfaces:\n");
  if (modelFile_->getParallelTimeSurfaces())
  {
    LogKit::writeLog("  Surface                                  : %s\n",     modelFile_->getTimeSurfFile()[0]);
    LogKit::writeLog("  Shift to top surface                     : %10.1f\n", modelFile_->getTimeDTop());
    LogKit::writeLog("  Time slice                               : %10.1f\n", modelFile_->getTimeLz());
    LogKit::writeLog("  Sampling density                         : %10.1f\n", modelFile_->getTimeDz());
    LogKit::writeLog("  Number of layers                         : %10d\n",   int(modelFile_->getTimeLz()/modelFile_->getTimeDz()+0.5));
  }
  else
  {
    LogKit::writeLog("  Top surface                              : %s\n",     modelFile_->getTimeSurfFile()[0]);
    LogKit::writeLog("  Base surface                             : %s\n",     modelFile_->getTimeSurfFile()[1]);
    LogKit::writeLog("  Number of layers                         : %10d\n",   modelFile_->getTimeNz());
    LogKit::writeLog("  Minimum allowed value for lmin/lmax      : %10.2f\n", modelSettings_->getLzLimit());
  }
  if (modelFile_->getDepthSurfFile() != NULL)
  {
    LogKit::writeLog("\nDepth surfaces:\n");
    LogKit::writeLog("  Top surface                              : %s\n",   modelFile_->getDepthSurfFile()[0]);
    LogKit::writeLog("  Base surface                             : %s\n",   modelFile_->getDepthSurfFile()[1]);
    LogKit::writeLog("  Number of layers                         : %10d\n", modelFile_->getDepthNz());
  }
  //
  // BACKGROUND
  //
  if (modelFile_->getGenerateBackground()) 
  {
    LogKit::writeLog("\nBackground model (estimated):\n");
    Vario       * vario  = modelSettings_->getBackgroundVario();
    GenExpVario * pVario = dynamic_cast<GenExpVario*>(vario);
    LogKit::writeLog("  Variogram\n");
    LogKit::writeLog("    Model                                  : %10s\n",vario->getType());
    if (pVario != NULL)
    LogKit::writeLog("    Power                                  : %10.1f\n",pVario->getPower());
    LogKit::writeLog("    Range                                  : %10.1f\n",vario->getRange());
    if (vario->getAnisotropic()) 
    {
      LogKit::writeLog("    Subrange                               : %10.1f\n",vario->getSubRange());
      LogKit::writeLog("    Angle                                  : %10.1f\n",vario->getAngle());
    }
    LogKit::writeLog("  High cut frequency for well logs         : %10.1f\n",modelSettings_->getMaxHzBackground());
  }
  else
  {
    float * constBack = modelFile_->getConstBack();
    LogKit::writeLog("\nBackground model:\n");
    if (constBack[0] > 0)
      LogKit::writeLog("  p-wave velocity                          : %10.1f\n",constBack[0]);
    else
      LogKit::writeLog("  p-wave velocity read from file           : %10s\n",modelFile_->getBackFile()[0]);
    
    if (constBack[1] > 0)
      LogKit::writeLog("  s-wave velocity                          : %10.1f\n",constBack[1]);
    else
      LogKit::writeLog("  s-wave velocity read from file           : %10s\n",modelFile_->getBackFile()[1]);
      
    if (constBack[2] > 0)
      LogKit::writeLog("  Density                                  : %10.1f\n",constBack[2]);
    else
      LogKit::writeLog("  Density read from file                   : %10s\n",modelFile_->getBackFile()[2]);
  }
  //
  // SEISMIC
  //
  if (modelSettings_->getGenerateSeismic())
  {
    LogKit::writeLog("\nGeneral settings for seismic:\n");
    LogKit::writeLog("  Generating seismic                       : %10s\n","yes");
    for (int i = 0 ; i < modelSettings_->getNumberOfAngles() ; i++)
    {
      LogKit::writeLog("\nSettings for AVO stack %d:\n",i+1);
      LogKit::writeLog("  Angle                                    : %10.1f\n",(modelSettings_->getAngle()[i]*180/PI));
      LogKit::writeLog("  Read wavelet from file                   : %s\n",modelFile_->getWaveletFile()[i]);
    }
  }
  else
  {
    // PRIOR CORRELATION
    Vario * corr = modelSettings_->getLateralCorr();
    if (corr != NULL) {
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::writeLog("\nPrior correlation:\n");
      LogKit::writeLog("  Range of allowed parameter values:\n");
      LogKit::writeLog("    Var{Vp}  - min                         : %10.1e\n",modelSettings_->getVarAlphaMin());
      LogKit::writeLog("    Var{Vp}  - max                         : %10.1e\n",modelSettings_->getVarAlphaMax());
      LogKit::writeLog("    Var{Vs}  - min                         : %10.1e\n",modelSettings_->getVarBetaMin());
      LogKit::writeLog("    Var{Vs}  - max                         : %10.1e\n",modelSettings_->getVarBetaMax());
      LogKit::writeLog("    Var{Rho} - min                         : %10.1e\n",modelSettings_->getVarRhoMin());
      LogKit::writeLog("    Var{Rho} - max                         : %10.1e\n",modelSettings_->getVarRhoMax());  
      LogKit::writeLog("  Lateral correlation:\n");
      LogKit::writeLog("    Model                                  : %10s\n",corr->getType());
      if (pCorr != NULL)
        LogKit::writeLog("    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::writeLog("    Range                                  : %10.1f\n",corr->getRange());
      if (corr->getAnisotropic())
      {
        LogKit::writeLog("    Subrange                               : %10.1f\n",corr->getSubRange());
        LogKit::writeLog("    Angle                                  : %10.1f\n",corr->getAngle());
      }
    }
    LogKit::writeLog("\nGeneral settings for seismic:\n");
    LogKit::writeLog("  White noise component                    : %10.2f\n",modelSettings_->getWNC());
    LogKit::writeLog("  Low cut for inversion                    : %10.1f\n",modelSettings_->getLowCut());
    LogKit::writeLog("  High cut for inversion                   : %10.1f\n",modelSettings_->getHighCut());
    corr  = modelSettings_->getAngularCorr();
    GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
    LogKit::writeLog("  Angular correlation:\n");
    LogKit::writeLog("    Model                                  : %10s\n",corr->getType());
    if (pCorr != NULL)
      LogKit::writeLog("    Power                                  : %10.1f\n",pCorr->getPower());
    LogKit::writeLog("    Range                                  : %10.1f\n",corr->getRange()*180.0/PI);
    if (corr->getAnisotropic())
    {
      LogKit::writeLog("    Subrange                               : %10.1f\n",corr->getSubRange()*180.0/PI);
      LogKit::writeLog("    Angle                                  : %10.1f\n",corr->getAngle());
    }
    bool estimateNoise = false;
    for (int i = 0 ; i < modelSettings_->getNumberOfAngles() ; i++) {
      estimateNoise = estimateNoise || (modelSettings_->getNoiseEnergy()[i]==RMISSING); 
    }
    LogKit::writeLog("\nGeneral settings for wavelet:\n");
    if (estimateNoise)
      LogKit::writeLog("  Maximum shift in noise estimation        : %10.1f\n",modelSettings_->getMaxWaveletShift());
    LogKit::writeLog("  Minimum relative amplitude               : %10.3f\n",modelSettings_->getMinRelWaveletAmp());
    LogKit::writeLog("  Wavelet tapering length                  : %10.1f\n",modelSettings_->getWaveletTaperingL());
    
    for (int i = 0 ; i < modelSettings_->getNumberOfAngles() ; i++)
    {
      LogKit::writeLog("\nSettings for AVO stack %d:\n",i+1);
      LogKit::writeLog("  Angle                                    : %10.1f\n",(modelSettings_->getAngle()[i]*180/PI));
      LogKit::writeLog("  Segy offset                              : %10.1f\n",modelSettings_->getSegyOffset());
      LogKit::writeLog("  Data                                     : %s\n",modelFile_->getSeismicFile()[i]);
      if (modelFile_->getWaveletFile()[i][0] == '*')
        LogKit::writeLog("  Estimate wavelet                         : %10s\n", "yes");
      else
        LogKit::writeLog("  Read wavelet from file                   : %s\n",modelFile_->getWaveletFile()[i]);
      
      float * noiseEnergy   = modelSettings_->getNoiseEnergy();
      bool  * matchEnergies = modelSettings_->getMatchEnergies();
      if (noiseEnergy[i] == RMISSING) 
        LogKit::writeLog("  Estimate signal-to-noise ratio           : %10s\n", "yes");
      else
        if (hasSignalToNoiseRatio_)
          LogKit::writeLog("  Signal-to-noise ratio                    : %10.1f\n",noiseEnergy[i]);
        else
          LogKit::writeLog("  Error std.dev. in seismic                : %10.2e\n",noiseEnergy[i]);
      if (matchEnergies[i]) 
        LogKit::writeLog("  Match empirical and theoretical energies : %10s\n", "yes");
      else
        LogKit::writeLog("  Wavelet scale                            : %10.2e\n",modelFile_->getWaveletScale()[i]);
    }
  }
}
void Model::processPriorFaciesProb()
{
  int outputFlag = modelFile_->getModelSettings()->getOutputFlag();
  if((outputFlag & ModelSettings::FACIESPROB) >0 || (outputFlag & ModelSettings::FACIESPROBRELATIVE)>0)
  {
    int w1, i;
    int nz = modelFile_->getTimeNz();
    float * vtAlpha = new float[nz];            // vt = vertical trend
    float * vtBeta  = new float[nz];
    float * vtRho   = new float[nz];
    
    int nWells = modelFile_->getModelSettings()->getNumberOfWells();
    int nFacies = modelFile_->getModelSettings()->getNumberOfFacies();
    int ndata = nWells*nz;
    int * facieslog = new int[ndata];
    for (w1 = 0 ; w1 < nWells ; w1++)
    {
      if(!(wells_[w1]->isDeviated()))
      { 
        BlockedLogs * bw = new BlockedLogs(wells_[w1], timeSimbox_, randomGen_) ;
        bw->getVerticalTrend(bw->getAlpha(),vtAlpha);
        bw->getVerticalTrend(bw->getBeta(),vtBeta);
        bw->getVerticalTrend(bw->getRho(),vtRho);
        int * facies = new int[nz];
        bw->getVerticalTrendDiscrete(bw->getFacies(),facies,randomGen_);
        for(i=0;i<nz;i++)
        {
          if(vtAlpha[i]!=RMISSING && vtBeta[i]!=RMISSING && vtRho[i]!=RMISSING)
            facieslog[i+w1*nz] = facies[i];
          else
            facieslog[i+w1*nz] = IMISSING;
        }
      }
    }

    int f;
    float sum = 0.0f;
    int * nData = new int[nFacies];
    priorFacies_ = new float[nFacies];
    for(i=0;i<nFacies;i++)
      nData[i] = 0;

    for(i=0;i<ndata;i++)
    {
      if(facieslog[i]!=IMISSING)
      {
        f = facieslog[i];
        nData[f]++;
      }
    }
    for(i=0;i<nFacies;i++)
      sum+=nData[i];

    for(i=0;i<nFacies;i++)
      priorFacies_[i] = float(nData[i])/sum;
  }
  else
    priorFacies_ = NULL;
}
