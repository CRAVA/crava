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
#include "src/modelavodynamic.h"
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
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


ModelAVODynamic::ModelAVODynamic(ModelSettings       *& modelSettings,
                                 const InputFiles     * inputFiles,
                                 std::vector<bool>      failedGeneralDetails,
                                 std::vector<bool>      failedStaticDetails,
                                 Simbox               * timeSimbox,
                                 Simbox              *& timeBGSimbox,
                                 RandomGen            * randomGen,
                                 GridMapping          * timeDepthMapping,
                                 GridMapping          * timeCutMapping,
                                 std::vector<Surface *> waveletEstimInterval,
                                 std::vector<Surface *> wellMoveInterval,
                                 std::vector<Surface *> faciesEstimInterval,
                                 ModelAVOStatic       * modelAVOstatic)
{
  numberOfAngles_         = modelSettings->getNumberOfAngles();
  background_             = NULL;
  correlations_           = NULL;
  seisCube_               = NULL;
  wavelet_                = NULL;
  reflectionMatrix_       = NULL;
  failed_                 = false;

  bool failedSimbox       = failedGeneralDetails[0];;
  bool failedDepthConv    = failedGeneralDetails[1];
  bool failedWells        = failedStaticDetails[0];
  bool failedExtraSurf    = failedStaticDetails[1];
  // ModelAVOStatic's failedPriorFacies is not used here.

  bool failedWavelet      = false;
  bool failedSeismic      = false;
  bool failedReflMat      = false;
  bool failedBackground   = false;
  bool failedPriorCorr    = false;

  bool failedLoadingModel = false;

  std::string errText("");

  if(!failedSimbox)
  {
    //
    // FORWARD MODELLING
    //
    if (modelSettings->getForwardModeling() == true)
    {
      processBackground(background_, modelAVOstatic->getWells(), timeSimbox, timeBGSimbox,
                        timeDepthMapping, timeCutMapping,
                        modelSettings, inputFiles,
                        errText, failedBackground);
      if (!failedBackground)
      {
        processReflectionMatrix(reflectionMatrix_, background_, modelAVOstatic->getWells(), modelSettings,
                                inputFiles, errText, failedReflMat);
        if (!failedReflMat)
        {
          processWavelets(wavelet_, seisCube_, modelAVOstatic->getWells(), reflectionMatrix_,
                          timeSimbox, waveletEstimInterval,
                          modelSettings, inputFiles, errText, failedWavelet);
        }
      }
    }
    else {
      //
      // INVERSION/ESTIMATION
      //
      bool estimationMode = modelSettings->getEstimationMode();

      if (!failedWells && !failedDepthConv) {
        bool backgroundDone = false;

        if(!(modelSettings->getOptimizeWellLocation() == true &&
             modelSettings->getGenerateBackground() == true))
        {
          if(estimationMode == false ||
             modelSettings->getEstimateBackground() == true ||
             modelSettings->getEstimateCorrelations() == true ||
             modelSettings->getOptimizeWellLocation() == true)
          {
            processBackground(background_, modelAVOstatic->getWells(), timeSimbox, timeBGSimbox,
                              timeDepthMapping, timeCutMapping,
                              modelSettings, inputFiles,
                              errText, failedBackground);
            backgroundDone = true;
          }
          if(failedBackground == false && backgroundDone == true &&
            (estimationMode == false || modelSettings->getEstimateWaveletNoise() ||
             modelSettings->getOptimizeWellLocation() == true))
          {
            processReflectionMatrix(reflectionMatrix_, background_, modelAVOstatic->getWells(), modelSettings,
                                      inputFiles, errText, failedReflMat);
          }
          else if(estimationMode == true &&
                  modelSettings->getEstimateWaveletNoise() == true &&
                  modelSettings->getEstimateBackground() == false &&
                  modelSettings->getEstimateCorrelations() == false)
          {
            processReflectionMatrix(reflectionMatrix_, background_, modelAVOstatic->getWells(), modelSettings,
                                      inputFiles, errText, failedReflMat);
            backgroundDone = true; //Not really, but do not need it in this case.
          }
          if(failedBackground == false && backgroundDone == true &&
             failedReflMat == false && failedExtraSurf == false &&
             (estimationMode == false || modelSettings->getEstimateWaveletNoise() ||
              modelSettings->getOptimizeWellLocation() == true))
          {
            processSeismic(seisCube_, timeSimbox, timeDepthMapping, timeCutMapping,
                           modelSettings, inputFiles, errText, failedSeismic);
            if(failedSeismic == false && modelSettings->getOptimizeWellLocation() == true)
            {
              for(int i=0;i<modelSettings->getNumberOfAngles();i++)
                seisCube_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
              modelAVOstatic->processWellLocation(seisCube_, modelAVOstatic->getWells(), reflectionMatrix_,
                                                  timeSimbox, modelSettings, wellMoveInterval);
              for(int i=0;i<modelSettings->getNumberOfAngles();i++)
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
          processReflectionMatrix(reflectionMatrix_, background_, modelAVOstatic->getWells(), modelSettings,
                                    inputFiles, errText, failedReflMat);
          if (failedReflMat == false && failedExtraSurf == false)
          {
            processSeismic(seisCube_, timeSimbox, timeDepthMapping, timeCutMapping,
                           modelSettings, inputFiles, errText, failedSeismic);
            if(failedSeismic == false)
            {
              modelAVOstatic->processWellLocation(seisCube_, modelAVOstatic->getWells(), reflectionMatrix_,
                                                  timeSimbox, modelSettings, wellMoveInterval);
            }
          }
          if(estimationMode == false ||
             modelSettings->getEstimateBackground() == true ||
             modelSettings->getEstimateCorrelations() == true ||
             modelSettings->getEstimateWaveletNoise() == true)
          {
            processBackground(background_, modelAVOstatic->getWells(), timeSimbox, timeBGSimbox,
                              timeDepthMapping, timeCutMapping,
                              modelSettings, inputFiles, errText, failedBackground);
            backgroundDone = true;
          }
        }

        if(failedBackground == false && backgroundDone == true &&
           failedSeismic == false && failedReflMat == false &&
           (estimationMode == false || modelSettings->getEstimateCorrelations() == true))
        {
          processPriorCorrelations(correlations_, background_, modelAVOstatic->getWells(), timeSimbox, modelSettings,
                                   inputFiles, errText, failedPriorCorr);
        }

        if(failedSeismic == false && failedBackground == false &&
          (estimationMode == false || modelSettings->getEstimateWaveletNoise() || modelSettings->getOptimizeWellLocation()))
        {
          modelAVOstatic->addSeismicLogs(modelAVOstatic->getWells(), seisCube_, modelSettings);
          if(failedReflMat == false && failedExtraSurf == false)
          {
            processWavelets(wavelet_, seisCube_, modelAVOstatic->getWells(), reflectionMatrix_,
                            timeSimbox, waveletEstimInterval,
                            modelSettings, inputFiles, errText, failedWavelet);
          }
        }
        if(failedSeismic == false && failedWavelet == false && failedReflMat == false && failedBackground == false &&
          (modelSettings->getOptimizeWellLocation() || modelSettings->getEstimateWaveletNoise() ))
        {
          modelAVOstatic->generateSyntheticSeismic(wavelet_, modelAVOstatic->getWells(), reflectionMatrix_, timeSimbox, modelSettings);
        }
      }

      if (!failedWells) {
        if(estimationMode || (modelSettings->getWellOutputFlag() & IO::WELLS) > 0)
          modelAVOstatic->writeWells(modelAVOstatic->getWells(), modelSettings);
        if(estimationMode)
          modelAVOstatic->writeBlockedWells(modelAVOstatic->getWells(), modelSettings);
      }
    }
  }
  failedLoadingModel = failedSeismic    || failedPriorCorr  ||  failedReflMat   ||
                       failedBackground || failedWavelet;

  if (failedLoadingModel) {
    LogKit::WriteHeader("Error(s) while loading data");
    LogKit::LogMessage(LogKit::Error,"\n"+errText);
    LogKit::LogMessage(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel;
  failed_details_.push_back(failedSeismic);
  failed_details_.push_back(failedPriorCorr);
  failed_details_.push_back(failedReflMat);
  failed_details_.push_back(failedBackground);
  failed_details_.push_back(failedWavelet);
}


ModelAVODynamic::~ModelAVODynamic(void)
{
  if (wavelet_ != NULL)
  {
    for(int i=0;i<numberOfAngles_;i++)
      if(wavelet_[i] != NULL)
        delete wavelet_[i];
    delete [] wavelet_;
  }

  if(reflectionMatrix_ != NULL) {
    for(int i = 0;i<numberOfAngles_;i++)
      delete [] reflectionMatrix_[i] ;
    delete [] reflectionMatrix_ ;
  }

  for (int i = 0 ; i < static_cast<int>(localNoiseScale_.size()) ; i++) {
    if (localNoiseScale_[i] != NULL)
      delete localNoiseScale_[i];
  }

  if (correlations_ != NULL)
    delete correlations_;

  if(background_ != NULL)
    delete background_;

  if(seisCube_ != NULL) {
    for(int i=0;i<numberOfAngles_;i++)
      delete seisCube_[i];
    delete [] seisCube_;
  }
}

void
ModelAVODynamic::releaseGrids(void)
{
  background_->releaseGrids();
  delete background_;
  background_ = NULL;
  seisCube_ = NULL;
}

float **
ModelAVODynamic::readMatrix(const std::string & fileName, int n1, int n2,
                            const std::string & readReason,
                            std::string       & errText)
{
  float * tmpRes = new float[n1*n2+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,fileName);
  std::string text = "Reading "+readReason+" from file "+fileName+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
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
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
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
    LogKit::LogFormatted(LogKit::Low,"failed.\n");
  delete [] tmpRes;
  return(result);
}

void
ModelAVODynamic::processSeismic(FFTGrid        **& seisCube,
                                Simbox          *& timeSimbox,
                                GridMapping     *& timeDepthMapping,
                                GridMapping     *& timeCutMapping,
                                ModelSettings   *& modelSettings,
                                const InputFiles * inputFiles,
                                std::string      & errText,
                                bool             & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  if(inputFiles->getNumberOfSeismicFiles() > 0)
  {
    LogKit::WriteHeader("Reading seismic data");

    int nAngles = modelSettings->getNumberOfAngles();
    const SegyGeometry ** geometry = new const SegyGeometry * [nAngles];
    seisCube = new FFTGrid * [nAngles];

    bool outsideWarning = false;
    for (int i = 0 ; i < nAngles ; i++) {
      geometry[i] = NULL;
      std::string tmpErrText("");
      std::string fileName = inputFiles->getSeismicFile(i);
      std::string angle    = NRLib::ToString(modelSettings->getAngle(i)*(180/M_PI), 1);
      std::string dataName = "Seismic data angle stack "+angle;
      float       offset = modelSettings->getLocalSegyOffset(i);
      if(offset < 0)
        offset = modelSettings->getSegyOffset();

      int outsideTraces = 0;
      ModelGeneral::readGridFromFile(fileName,
                                     dataName,
                                     offset,
                                     seisCube[i],
                                     geometry[i],
                                     modelSettings->getTraceHeaderFormat(i),
                                     FFTGrid::DATA,
                                     timeSimbox,
                                     modelSettings,
                                     outsideTraces,
                                     tmpErrText);
      if(tmpErrText != "")
      {
        tmpErrText += "\nReading of file \'"+fileName+"\' for "+dataName+" failed.\n";
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
            LogKit::LogMessage(LogKit::Warning, "Warning: "+NRLib::ToString(outsideTraces)+" traces in the grid were outside the data area in file "
              +fileName+". Note that this includes traces in the padding.\n");
            outsideWarning = true;
          }
        }
      }
    }
    LogKit::LogFormatted(LogKit::Low,"\n");
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
        LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
        LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
        for (int i = 0 ; i < nAngles ; i++)
        {
          if (geometry[i] != NULL) {
            double geoAngle = (-1)*timeSimbox->getAngle()*(180/M_PI);
            if (geoAngle < 0)
              geoAngle += 360.0;
            LogKit::LogFormatted(LogKit::Low,"Seismic data %d   %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",i,
                                 geometry[i]->GetX0(), geometry[i]->GetY0(),
                                 geometry[i]->Getlx(), geometry[i]->Getly(), geoAngle,
                                 geometry[i]->GetDx(), geometry[i]->GetDy());
          }
        }
      }

      if((modelSettings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0) {
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
                                 timeDepthMapping,
                                 timeCutMapping, *modelSettings->getTraceHeaderFormatOutput());
        }
      }
      if((modelSettings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0) {
        for(int i=0;i<nAngles;i++) {
          std::string angle    = NRLib::ToString(i);
          std::string baseName = IO::FileTemporarySeismic() + angle;
          seisCube[i]->writeCravaFile(baseName, timeSimbox);
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
ModelAVODynamic::processBackground(Background         *& background,
                                   WellData           ** wells,
                                   Simbox              * timeSimbox,
                                   Simbox              * timeBGSimbox,
                                   GridMapping        *& timeDepthMapping,
                                   GridMapping        *& timeCutMapping,
                                   ModelSettings       * modelSettings,
                                   const InputFiles    * inputFiles,
                                   std::string         & errText,
                                   bool                & failed)
{
  if (modelSettings->getForwardModeling())
    LogKit::WriteHeader("Earth Model");
  else
    LogKit::WriteHeader("Prior Expectations / Background Model");

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
    if (backVelFile != ""){
      bool dummy;
      ModelGeneral::loadVelocity(velocity, timeSimbox, modelSettings, backVelFile, dummy, errText, failed);
    }
    if (!failed)
    {
      if(modelSettings->getBackgroundVario() == NULL)
      {
        errText += "There is no variogram available for the background modelling.\n";
        failed = true;
      }
      for (int i=0 ; i<3 ; i++)
      {
        backModel[i] = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());
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
    if (modelSettings->getUseAIBackground())
      parName.push_back("AI "+modelSettings->getBackgroundType());
    else
      parName.push_back("Vp "+modelSettings->getBackgroundType());
    if (modelSettings->getUseSIBackground())
      parName.push_back("SI "+modelSettings->getBackgroundType());
    else if (modelSettings->getUseVpVsBackground())
      parName.push_back("Vp/Vs "+modelSettings->getBackgroundType());
    else
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
          int outsideTraces = 0;
          ModelGeneral::readGridFromFile(backFile,
                                         parName[i],
                                         offset,
                                         backModel[i],
                                         dummy1,
                                         dummy2,
                                         FFTGrid::PARAMETER,
                                         timeSimbox,
                                         modelSettings,
                                         outsideTraces,
                                         errorText);
          if(errorText != "")
          {
            errorText += "Reading of file '"+backFile+"' for parameter '"+parName[i]+"' failed\n";
            errText += errorText;
            failed = true;
          }
          else {
            backModel[i]->calculateStatistics();
            backModel[i]->setUndefinedCellsToGlobalAverage();
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
        backModel[i] = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());
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

      LogKit::LogFormatted(LogKit::Low, "\nSummary                Average   Minimum   Maximum\n");
      LogKit::LogFormatted(LogKit::Low, "--------------------------------------------------\n");
      for(int i=0 ; i<3 ; i++) {
        LogKit::LogFormatted(LogKit::Low, "%-20s %9.2f %9.2f %9.2f\n",
                             parName[i].c_str(),
                             backModel[i]->getAvgReal(),
                             backModel[i]->getMinReal(),
                             backModel[i]->getMaxReal());
      }
      if (modelSettings->getUseAIBackground())   { // Vp = AI/Rho     ==> lnVp = lnAI - lnRho
        LogKit::LogMessage(LogKit::Low, "\nMaking Vp background from AI and Rho\n");
        backModel[0]->subtract(backModel[2]);
      }
      if (modelSettings->getUseVpVsBackground()) { // Vs = SI/Rho     ==> lnVs = lnSI - lnRho
        LogKit::LogMessage(LogKit::Low, "\nMaking Vs background from SI and Rho\n");
        backModel[1]->subtract(backModel[2]);
      }
      else if (modelSettings->getUseVpVsBackground()) { // Vs = Vp/(Vp/Vs) ==> lnVs = lnVp - ln(Vp/Vs)
        LogKit::LogMessage(LogKit::Low, "\nMaking Vs background from Vp and Vp/Vs\n");
        backModel[1]->subtract(backModel[0]);
        backModel[1]->changeSign();
      }
      background = new Background(backModel);
    }
  }

  if (failed == false) {
    if((modelSettings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
      background->writeBackgrounds(timeSimbox,
                                   timeDepthMapping,
                                   timeCutMapping,
                                   modelSettings->getFileGrid(),
                                   *modelSettings->getTraceHeaderFormatOutput());
    }
  }

  Timings::setTimePriorExpectation(wall,cpu);
}


void
ModelAVODynamic::processPriorCorrelations(Corr                   *& correlations,
                                          Background              * background,
                                          WellData               ** wells,
                                          Simbox                  * timeSimbox,
                                          ModelSettings           * modelSettings,
                                          const InputFiles        * inputFiles,
                                          std::string             & errText,
                                          bool                    & failed)
{
  bool printResult = ((modelSettings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0 ||
                      modelSettings->getEstimationMode() == true);
  if (modelSettings->getDoInversion() || printResult)
  {
    LogKit::WriteHeader("Prior Covariance");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    const std::string & paramCovFile = inputFiles->getParamCorrFile();
    const std::string & corrTFile    = inputFiles->getTempCorrFile();

    bool estimateParamCov = paramCovFile == "";
    bool estimateTempCorr = corrTFile    == "";

    //
    // Read parameter covariance (Var0) from file
    //
    float ** paramCorr = NULL;
    bool failedParamCorr = false;
    if(!estimateParamCov)
    {
      std::string tmpErrText("");
      paramCorr = readMatrix(paramCovFile, 3, 3, "parameter covariance", tmpErrText);
      if(paramCorr == NULL)
      {
        errText += "Reading of file "+paramCovFile+" for parameter covariance matrix failed\n";
        errText += tmpErrText;
        failedParamCorr = true;
      }
    }

    //
    // Estimate lateral correlation from seismic data
    //
    Surface * CorrXY = findCorrXYGrid(timeSimbox, modelSettings);

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
    if (estimateParamCov || estimateTempCorr) //Need well estimation
    {
      std::string tmpErrTxt;
      Analyzelog * analyze = new Analyzelog(wells,
                                            background,
                                            timeSimbox,
                                            modelSettings,
                                            tmpErrTxt);
      if (tmpErrTxt != "") {
        errText += tmpErrTxt;
        failedParamCorr = true;
      }

      if(estimateParamCov)
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
          LogKit::LogFormatted(LogKit::High,
            "\nOnly able to estimate %d of %d lags needed in temporal correlation. The rest are set to 0.\n", nEst, nCorrT);
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
ModelAVODynamic::findCorrXYGrid(Simbox * timeSimbox, ModelSettings * modelSettings)
{
  float dx  = static_cast<float>(timeSimbox->getdx());
  float dy  = static_cast<float>(timeSimbox->getdy());

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
ModelAVODynamic::estimateCorrXYFromSeismic(Surface *& corrXY,
                                           FFTGrid ** seisCube,
                                           int        nAngles)
{
  FFTGrid * transf;
  float   * grid;

  int n = static_cast<int>(corrXY->GetNI()*corrXY->GetNJ());
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
ModelAVODynamic::processReflectionMatrix(float            **& reflectionMatrix,
                                         Background         * background,
                                         WellData          ** wells,
                                         ModelSettings      * modelSettings,
                                         const InputFiles   * inputFiles,
                                         std::string        & errText,
                                         bool               & failed)
{
  LogKit::WriteHeader("Reflection matrix");
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model or wells.
  //
  const std::string & reflMatrFile = inputFiles->getReflMatrFile();
  const double        vpvs         = modelSettings->getVpVsRatio();

  if (reflMatrFile != "") {
    std::string tmpErrText("");
    reflectionMatrix = readMatrix(reflMatrFile, modelSettings->getNumberOfAngles(), 3, "reflection matrix", tmpErrText);
    if(reflectionMatrix == NULL) {
      errText += "Reading of file "+reflMatrFile+ " for reflection matrix failed\n";
      errText += tmpErrText;
      failed = true;
    }
    LogKit::LogFormatted(LogKit::Low,"\nReflection parameters read from file.\n\n");
  }
  else if (vpvs != RMISSING) {
    LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs ratio specified in model file.\n");
    double vsvp = 1.0/vpvs;
    setupDefaultReflectionMatrix(reflectionMatrix, vsvp, modelSettings);
  }
  else if (background == NULL || modelSettings->getVpVsRatioFromWells()) {
    LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp and Vs from wells\n");
    double vsvp = vsvpFromWells(wells, modelSettings->getNumberOfWells());
    setupDefaultReflectionMatrix(reflectionMatrix, vsvp, modelSettings);
  }
  else {
    if (modelSettings->getForwardModeling())
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp and Vs from earth model\n");
    else
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp and Vs from background model\n");

    double vsvp = background->getMeanVsVp();
    setupDefaultReflectionMatrix(reflectionMatrix, vsvp, modelSettings);
  }
}

void
ModelAVODynamic::setupDefaultReflectionMatrix(float       **& reflectionMatrix,
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
  double vpvs = 1.0f/vsvp;
  LogKit::LogFormatted(LogKit::Low,"\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",vpvs);
  std::string text;
  if (vpvs < modelSettings->getVpVsRatioMin()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very small Vp/Vs-ratio has been detected. Values below %.2f are regarded unlikely. \n",modelSettings->getVpVsRatioMin());
    text  = "Check the Vp/Vs-ratio. A small value has been found. If the value is acceptable,\n";
    text += "   you can remove this task using the <minimim-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
  else if (vpvs > modelSettings->getVpVsRatioMax()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very large Vp/Vs-ratio has been detected. Values above %.2f are regarded unlikely. \n",modelSettings->getVpVsRatioMax());
    text  = "Check the Vp/Vs-ratio. A large value has been found. If the value is acceptable,\n";
    text += "   you can remove this task using the <maximum-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
}

double ModelAVODynamic::vsvpFromWells(WellData ** wells,
                                      int         nWells)
{
  int   N      = 0;
  float VsVp   = 0.0f;

  for(int i=0 ; i < nWells ; i++) {
    N    += wells[i]->getNumberOfVsVpSamples();
    VsVp += wells[i]->getMeanVsVp()*N;
  }
  VsVp /= N;

  return static_cast<double>(VsVp);
}

void
ModelAVODynamic::processWavelets(Wavelet                    **& wavelet,
                                 FFTGrid                     ** seisCube,
                                 WellData                    ** wells,
                                 float                       ** reflectionMatrix,
                                 Simbox                       * timeSimbox,
                                 const std::vector<Surface *> & waveletEstimInterval,
                                 ModelSettings                * modelSettings,
                                 const InputFiles             * inputFiles,
                                 std::string                  & errText,
                                 bool                         & failed)
{
  int error = 0;
  LogKit::WriteHeader("Processing/generating wavelets");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  unsigned int nAngles = modelSettings->getNumberOfAngles();

  wavelet = new Wavelet * [nAngles];
  localNoiseScale_.resize(nAngles);

  bool has3Dwavelet = false;
  for(unsigned int i=0 ; i < nAngles ; i++) {
    localNoiseScale_[i] = NULL;
    if (modelSettings->getWaveletDim(i) == Wavelet::THREE_D)
      has3Dwavelet = true;
    if(modelSettings->getEstimateWavelet(i) == true)
      modelSettings->setWaveletScale(i,1.0);
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
    float distance, sigma_m;
    modelSettings->getTimeGradientSettings(distance, sigma_m);
    std::vector<std::vector<double> > SigmaXY;
    for (unsigned int w=0; w<nWells; w++) {
      BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick();
      bl->setTimeGradientSettings(distance, sigma_m);
      bl->findSeismicGradient(seisCube, timeSimbox, nAngles,tGradX[w], tGradY[w],SigmaXY);
    }
  }

  if (timeSimbox->getdz() > 4.01f && modelSettings->getEstimateNumberOfWavelets() > 0)
  { // Require this density for wavelet estimation
    LogKit::LogFormatted(LogKit::Low,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
    LogKit::LogFormatted(LogKit::Low,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
    LogKit::LogFormatted(LogKit::Low,"         The number of layers must be increased.                                  \n");
    std::string text("");
    text += "Increase the number of layers to improve the quality of the wavelet estimation.\n";
    text += "   The minimum sampling density is "+NRLib::ToString(timeSimbox->getdz())+", and it should be";
    text += "lower than 4.0.\n   To obtain the desired density, the number of layers should be at least ";
    text += NRLib::ToString(static_cast<int>(timeSimbox->GetLZ()/4.0))+"\n";
    TaskList::addTask(text);
  }

  // check if local noise is set for some angles.
  bool localNoiseSet = false;
  for(unsigned int i=0 ; i < nAngles ; i++) {
    float angle = float(modelSettings->getAngle(i)*180.0/M_PI);
    LogKit::LogFormatted(LogKit::Low,"\nAngle stack : %.1f deg",angle);
    if(modelSettings->getForwardModeling()==false)
      seisCube[i]->setAccessMode(FFTGrid::RANDOMACCESS);
    if (modelSettings->getWaveletDim(i) == Wavelet::ONE_D)
      error += process1DWavelet(modelSettings,
                                inputFiles,
                                timeSimbox,
                                seisCube,
                                wells,
                                waveletEstimInterval,
                                reflectionMatrix[i],
                                errText,
                                wavelet[i],
                                i);
    else
      error += process3DWavelet(modelSettings,
                                inputFiles,
                                timeSimbox,
                                seisCube,
                                wells,
                                waveletEstimInterval,
                                reflectionMatrix[i],
                                errText,
                                wavelet[i],
                                i,
                                refTimeGradX,
                                refTimeGradY,
                                tGradX,
                                tGradY);

    if(localNoiseScale_[i]!=NULL)
      localNoiseSet = true;
    if(modelSettings->getForwardModeling()==false) // else, no seismic data
      seisCube[i]->endAccess();
  } // end i (angles)

  if(localNoiseSet==true) {
    for(unsigned int i=0;i<nAngles;i++)
      if(localNoiseScale_[i]==NULL)
        localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(),
                                         timeSimbox->getny(),
                                         1.0);
  }

  Timings::setTimeWavelets(wall,cpu);
  failed = error > 0;
}

int
ModelAVODynamic::process1DWavelet(ModelSettings                * modelSettings,
                                  const InputFiles             * inputFiles,
                                  Simbox                       * timeSimbox,
                                  FFTGrid                     ** seisCube,
                                  WellData                    ** wells,
                                  const std::vector<Surface *> & waveletEstimInterval,
                                  float                        * reflectionMatrix,
                                  std::string                  & errText,
                                  Wavelet                     *& wavelet,
                                  unsigned int                   i)
{
  int error = 0;
  Grid2D * shiftGrid(NULL);
  Grid2D * gainGrid(NULL);
  if(modelSettings->getUseLocalWavelet() && inputFiles->getScaleFile(i) != "") {
      Surface help(inputFiles->getScaleFile(i));
      gainGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
      resampleSurfaceToGrid2D(timeSimbox, &help, gainGrid);
  }
  if (modelSettings->getUseLocalWavelet() && inputFiles->getShiftFile(i) != ""){
    Surface helpShift(inputFiles->getShiftFile(i));
    shiftGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
    resampleSurfaceToGrid2D(timeSimbox, &helpShift, shiftGrid);
  }
  if (modelSettings->getUseLocalNoise() && inputFiles->getLocalNoiseFile(i) != ""){
    Surface helpNoise(inputFiles->getLocalNoiseFile(i));
    localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(), timeSimbox->getny(), 0.0);
    resampleSurfaceToGrid2D(timeSimbox, &helpNoise, localNoiseScale_[i]);
  }

  if (modelSettings->getEstimateWavelet(i))
    wavelet = new Wavelet1D(timeSimbox,
                            seisCube[i],
                            wells,
                            waveletEstimInterval,
                            modelSettings,
                            reflectionMatrix,
                            i,
                            error,
                            errText);

  else { //Not estimation modus
    if(modelSettings->getUseRickerWavelet(i))
        wavelet = new Wavelet1D(modelSettings,
                                reflectionMatrix,
                                modelSettings->getAngle(i),
                                modelSettings->getRickerPeakFrequency(i),
                                error);
    else {
      const std::string & waveletFile = inputFiles->getWaveletFile(i);
      int fileFormat = getWaveletFileFormat(waveletFile,errText);
      if(fileFormat < 0) {
        errText += "Unknown file format of file '"+waveletFile+"'.\n";
        error++;
      }
      else
        wavelet = new Wavelet1D(waveletFile,
                                fileFormat,
                                modelSettings,
                                reflectionMatrix,
                                modelSettings->getAngle(i),
                                error,
                                errText);
    }
      // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
      if(seisCube!=NULL) {// If forward modeling, we have no seismic, can not prescale wavelet.
        float       prescale = wavelet->findGlobalScaleForGivenWavelet(modelSettings, timeSimbox, seisCube[i], wells);
        const float limHigh  = 3.0f;
        const float limLow   = 0.33f;

        if(modelSettings->getEstimateGlobalWaveletScale(i)) // prescale, then we have correct size order, and later scale estimation will be ok.
           wavelet->multiplyRAmpByConstant(prescale);
        else {
          if(modelSettings->getWaveletScale(i)!= 1.0f && (prescale>limHigh || prescale<limLow)) {
             std::string text = "The wavelet given for angle no "+NRLib::ToString(i)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
            if(modelSettings->getEstimateLocalScale(i)) {
              errText += text;
              error++;
            }
            else {
              LogKit::LogFormatted(LogKit::Warning,"\nWARNING: "+text);
              TaskList::addTask("The wavelet is badly scaled. Consider having CRAVA estimate global wavelet scale");
            }
          }
        }
      }
      if (error == 0)
        wavelet->resample(static_cast<float>(timeSimbox->getdz()),
                          timeSimbox->getnz(),
                          modelSettings->getNZpad());
  }

  if (error == 0) {
    wavelet->scale(modelSettings->getWaveletScale(i));

    if (modelSettings->getForwardModeling() == false) {
      float SNRatio = wavelet->calculateSNRatioAndLocalWavelet(timeSimbox,
                                                               seisCube[i],
                                                               wells,
                                                               modelSettings,
                                                               errText,
                                                               error,
                                                               i,
                                                               localNoiseScale_[i],
                                                               shiftGrid,
                                                               gainGrid);

      if(modelSettings->getEstimateSNRatio(i))
        modelSettings->setSNRatio(i,SNRatio);
    }

    if (error == 0) {
      if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (modelSettings->getEstimationMode() && modelSettings->getEstimateWavelet(i))) {
        std::string type;
        if (modelSettings->getEstimateWavelet(i))
          type = "Estimated_";
        else if (modelSettings->getWaveletScale(i) == 1.00)
          type = "";
        else
          type = "Scaled_";
        wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0); // dt_max = 1.0;
      }

      float SNRatio = modelSettings->getSNRatio(i);
      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((SNRatio <=SNLow  || SNRatio > SNHigh) && modelSettings->getForwardModeling()==false) {
        errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio)+" for cube "+NRLib::ToString(i+1)+".\n";
        errText += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }

      bool useLocalNoise = modelSettings->getEstimateLocalNoise(i) || inputFiles->getLocalNoiseFile(i) != "";
      bool useLocalShift = modelSettings->getEstimateLocalShift(i) || inputFiles->getShiftFile(i) != "";
      bool useLocalGain  = modelSettings->getEstimateLocalScale(i) || inputFiles->getScaleFile(i) != "";

      if (useLocalNoise)
        readAndWriteLocalGridsToFile(inputFiles->getLocalNoiseFile(i),
                                     IO::PrefixLocalNoise(),
                                     1.0,  // Scale map with this factor before writing to disk
                                     modelSettings,
                                     i,
                                     timeSimbox,
                                     localNoiseScale_[i]);

      if (useLocalShift) {
        readAndWriteLocalGridsToFile(inputFiles->getShiftFile(i),
                                     IO::PrefixLocalWaveletShift(),
                                     1.0,
                                     modelSettings,
                                     i,
                                     timeSimbox,
                                     shiftGrid);
        wavelet->setShiftGrid(shiftGrid);
      }

      if (useLocalGain) {
        readAndWriteLocalGridsToFile(inputFiles->getScaleFile(i),
                                     IO::PrefixLocalWaveletGain(),
                                     1.0,
                                     modelSettings,
                                     i,
                                     timeSimbox,
                                     gainGrid);
        wavelet->setGainGrid(gainGrid);
      }
    }
  }
  return error;
}

int
ModelAVODynamic::process3DWavelet(ModelSettings                           * modelSettings,
                                  const InputFiles                        * inputFiles,
                                  Simbox                                  * timeSimbox,
                                  FFTGrid                                ** seisCube,
                                  WellData                               ** wells,
                                  const std::vector<Surface *>            & waveletEstimInterval,
                                  float                                   * reflectionMatrix,
                                  std::string                             & errText,
                                  Wavelet                                *& wavelet,
                                  unsigned int                              i,
                                  const NRLib::Grid2D<float>              & refTimeGradX,
                                  const NRLib::Grid2D<float>              & refTimeGradY,
                                  const std::vector<std::vector<double> > & tGradX,
                                  const std::vector<std::vector<double> > & tGradY)
{
  int error = 0;
  if (modelSettings->getEstimateWavelet(i)) {
    wavelet = new Wavelet3D(inputFiles->getWaveletFilterFile(i),
                            waveletEstimInterval,
                            refTimeGradX,
                            refTimeGradY,
                            tGradX,
                            tGradY,
                            seisCube[i],
                            modelSettings,
                            wells,
                            timeSimbox,
                            reflectionMatrix,
                            i,
                            error,
                            errText);
  }
  else { //Not estimation modus
    const std::string & waveletFile = inputFiles->getWaveletFile(i);
    int fileFormat = getWaveletFileFormat(waveletFile,errText);
    if(fileFormat < 0) {
      errText += "Unknown file format of file '"+waveletFile+"'.\n";
      error++;
    }
    else {
      wavelet = new Wavelet3D(waveletFile,
                              fileFormat,
                              modelSettings,
                              reflectionMatrix,
                              modelSettings->getAngle(i),
                              error,
                              errText,
                              inputFiles->getWaveletFilterFile(i));
      if (error == 0)
        wavelet->resample(static_cast<float>(timeSimbox->getdz()),
                          timeSimbox->getnz(),
                          modelSettings->getNZpad());
    }
  }
  if ((modelSettings->getEstimationMode() == false) && !timeSimbox->getIsConstantThick()) {
    errText += "Simbox must have constant thickness if forward modelling or inversion when 3D wavelet.\n";
    error++;
  }
  if (error == 0) {
    wavelet->scale(modelSettings->getWaveletScale(i));
    bool localEst = (modelSettings->getEstimateLocalScale(i) || modelSettings->getEstimateLocalShift(i) ||
                     modelSettings->getEstimateLocalNoise(i) || modelSettings->getEstimateGlobalWaveletScale(i) ||
                     modelSettings->getEstimateSNRatio(i));

    if (localEst && modelSettings->getForwardModeling() == false) {
      float SNRatio = wavelet->calculateSNRatio(timeSimbox,
                                                seisCube[i],
                                                wells,
                                                modelSettings,
                                                errText,
                                                error,
                                                refTimeGradX,
                                                refTimeGradY,
                                                tGradX,
                                                tGradY,
                                                i);
      if(modelSettings->getEstimateSNRatio(i))
        modelSettings->setSNRatio(i,SNRatio);

    }
    if (error == 0) {
      if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (modelSettings->getEstimationMode() && modelSettings->getEstimateWavelet(i))) {
        std::string type;
        if (modelSettings->getEstimateWavelet(i))
          type = "Estimated_";
        else if (modelSettings->getWaveletScale(i) == 1.00)
          type = "";
        else
          type = "Scaled_";
        wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0); // dt_max = 1.0;
      }

      float SNRatio = modelSettings->getSNRatio(i);
      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((SNRatio <=SNLow  || SNRatio > SNHigh) && modelSettings->getForwardModeling()==false) {
        errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio)+" for cube "+NRLib::ToString(i+1)+".\n";
        errText += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }
    }
  }

  return error;
}

int
ModelAVODynamic::getWaveletFileFormat(const std::string & fileName, std::string & errText)
{
  int fileformat = -1;
  int line       = 0;
  int pos;
  std::string dummyStr;
  std::string targetString;

  std::ifstream file;
  NRLib::OpenRead(file,fileName);

  std::getline(file,dummyStr);
  line++;
  targetString = "pulse file-3";
  pos = Utils::findEnd(dummyStr, 0, targetString);
  if (pos >= 0)
    fileformat = Wavelet::NORSAR;
  file.close();
  file.clear();

  if(fileformat<0) { // not norsar format
      // test for jason file format
    NRLib::OpenRead(file,fileName);
    line         = 0;
    int thisLine = 0;
    bool lineIsComment = true;
    while (lineIsComment == true) {
      NRLib::ReadNextToken(file,dummyStr,line);
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
    if (NRLib::IsNumber(dummyStr)) // not convertable number
      fileformat= Wavelet::JASON;
  }
  return fileformat;
}

void
ModelAVODynamic::readAndWriteLocalGridsToFile(const std::string   & fileName,
                                              const std::string   & type,
                                              const float           scaleFactor,
                                              const ModelSettings * modelSettings,
                                              const unsigned int    i,
                                              Simbox              * timeSimbox,
                                              Grid2D             *& grid)
{
  bool   estimationMode   = modelSettings->getEstimationMode();
  int    outputFormat     = modelSettings->getOutputGridFormat();
  double angle            = modelSettings->getAngle(i)*180.0/M_PI;

  Surface * help = NULL;

  if(fileName != "") {
    std::string toPath = NRLib::RemovePath(fileName);

    if (type == IO::PrefixLocalNoise())
      toPath = NRLib::PrependDir(IO::PathToNoise(), toPath);
    else
      toPath = NRLib::PrependDir(IO::PathToWavelets(), toPath);

    NRLib::CreateDirIfNotExists(toPath);
    NRLib::CopyFile(fileName, toPath, true);
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
ModelAVODynamic::resampleSurfaceToGrid2D(const Simbox  * simbox,
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
ModelAVODynamic::resampleGrid2DToSurface(const Simbox   * simbox,
                                         const Grid2D   * grid,
                                         Surface       *& surface)
{
  double xmin,xmax,ymin,ymax;
  simbox->getMinAndMaxXY(xmin,xmax,ymin,ymax);
  int nx,ny;
  double angle = simbox->getAngle()*180.0/M_PI;
  if(angle > -45 || angle < 45)
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::cos(simbox->getAngle())+0.5)) * 2;
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::cos(simbox->getAngle())+0.5)) * 2;
  }
  else
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::sin(simbox->getAngle())+0.5)) * 2;
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::sin(simbox->getAngle())+0.5)) * 2;
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
ModelAVODynamic::findTimeGradientSurface(const std::string    & refTimeFile,
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
    NRLib::RegularSurfaceRotated<float> t0surface(refTimeFile);

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
      for (int i = nx-1; i >= 0; i--) {
        for (int j = ny-1; j >= 0; j--) {
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

