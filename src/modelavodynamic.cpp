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
#include "src/modelavodynamic.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
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
#include "src/seismicparametersholder.h"

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
#include "nrlib/volume/volume.hpp"


ModelAVODynamic::ModelAVODynamic(ModelSettings               *& modelSettings,
                                 const InputFiles             * inputFiles,
                                 std::vector<bool>              failedGeneralDetails,
                                 std::vector<bool>              failedStaticDetails,
                                 const Simbox                 * timeSimbox,
                                 const Simbox                 * timeBGSimbox,
                                 const Surface                * correlationDirection,
                                 RandomGen                    * /*randomGen*/,
                                 GridMapping                  * timeDepthMapping,
                                 const GridMapping            * timeCutMapping,
                                 const std::vector<Surface *> & waveletEstimInterval,
                                 const std::vector<Surface *> & /*wellMoveInterval*/,
                                 const std::vector<Surface *> & /*faciesEstimInterval*/,
                                 ModelAVOStatic               * modelAVOstatic,
                                 ModelGeneral                 * modelGeneral,
                                 int                            t,
                                 SeismicParametersHolder      & seismicParameters)
{
  seisCube_               = NULL;
  wavelet_                = NULL;
  reflectionMatrix_       = NULL;
  failed_                 = false;
  thisTimeLapse_          = t;
  angularCorr_            = modelSettings->getAngularCorr(thisTimeLapse_);
  SNRatio_                = modelSettings->getSNRatio(thisTimeLapse_);
  useLocalNoise_          = modelSettings->getUseLocalNoise(thisTimeLapse_);
  angle_                  = modelSettings->getAngle(thisTimeLapse_);
  estimateWavelet_        = modelSettings->getEstimateWavelet(thisTimeLapse_);
  matchEnergies_          = modelSettings->getMatchEnergies(thisTimeLapse_);
  numberOfAngles_         = modelSettings->getNumberOfAngles(thisTimeLapse_);

  bool failedSimbox       = failedGeneralDetails[0];
  bool failedDepthConv    = failedGeneralDetails[1];
  bool failedWells        = failedGeneralDetails[2];
  bool failedExtraSurf    = failedStaticDetails[0];
  // ModelAVOStatic's failedPriorFacies is not used here.

  bool failedWavelet      = false;
  bool failedSeismic      = false;
  bool failedReflMat      = false;
  bool failedBackground   = false;
  bool failedPriorCorr    = false;

  bool failedLoadingModel = false;

  Background * background = NULL;

  std::string errText("");

  if(!failedSimbox)
  {
    //
    // FORWARD MODELLING
    //
    if (modelSettings->getForwardModeling() == true)
    {
      processBackground(background,
                        modelGeneral->getWells(),
                        timeSimbox,
                        timeBGSimbox,
                        timeDepthMapping,
                        timeCutMapping,
                        modelSettings,
                        modelGeneral,
                        inputFiles,
                        thisTimeLapse_,
                        errText,
                        failedBackground);

      if (!failedBackground)
      {
        processReflectionMatrix(reflectionMatrix_, background, modelGeneral->getWells(), modelSettings,
                                inputFiles, errText, failedReflMat);
        if (!failedReflMat)
        {
          processWavelets(wavelet_, seisCube_, modelGeneral->getWells(), reflectionMatrix_,
                          timeSimbox, correlationDirection, waveletEstimInterval,
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
            processBackground(background,
                              modelGeneral->getWells(),
                              timeSimbox,
                              timeBGSimbox,
                              timeDepthMapping,
                              timeCutMapping,
                              modelSettings,
                              modelGeneral,
                              inputFiles,
                              thisTimeLapse_,
                              errText,
                              failedBackground);

            backgroundDone = true;
          }
          if(failedBackground == false && backgroundDone == true &&
            (estimationMode == false || modelSettings->getEstimateWaveletNoise() ||
             modelSettings->getOptimizeWellLocation() == true))
          {
            processReflectionMatrix(reflectionMatrix_, background, modelGeneral->getWells(), modelSettings,
                                      inputFiles, errText, failedReflMat);
          }
          else if(estimationMode == true &&
                  modelSettings->getEstimateWaveletNoise() == true &&
                  modelSettings->getEstimateBackground() == false &&
                  modelSettings->getEstimateCorrelations() == false)
          {
            processReflectionMatrix(reflectionMatrix_, background, modelGeneral->getWells(), modelSettings,
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
              for(int i=0;i<numberOfAngles_;i++)
                seisCube_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
              modelGeneral->processWellLocation(seisCube_, reflectionMatrix_,
                                                modelSettings, modelAVOstatic->getWellMoveInterval());
              for(int i=0;i<numberOfAngles_;i++)
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
          processReflectionMatrix(reflectionMatrix_, background, modelGeneral->getWells(), modelSettings,
                                    inputFiles, errText, failedReflMat);
          if (failedReflMat == false && failedExtraSurf == false)
          {
            processSeismic(seisCube_, timeSimbox, timeDepthMapping, timeCutMapping,
                           modelSettings, inputFiles, errText, failedSeismic);
            if(failedSeismic == false)
            {
              modelGeneral->processWellLocation(seisCube_, reflectionMatrix_,
                                                modelSettings, modelAVOstatic->getWellMoveInterval());
            }
          }
          if(estimationMode == false ||
             modelSettings->getEstimateBackground() == true ||
             modelSettings->getEstimateCorrelations() == true ||
             modelSettings->getEstimateWaveletNoise() == true)
          {
            processBackground(background,
                              modelGeneral->getWells(),
                              timeSimbox,
                              timeBGSimbox,
                              timeDepthMapping,
                              timeCutMapping,
                              modelSettings,
                              modelGeneral,
                              inputFiles,
                              thisTimeLapse_,
                              errText,
                              failedBackground);

            backgroundDone = true;
          }
        }

        if(failedBackground == false && backgroundDone == true &&
           failedSeismic == false && failedReflMat == false &&
           (estimationMode == false || modelSettings->getEstimateCorrelations() == true))
        {
          modelGeneral->processPriorCorrelations(background,
                                                 modelGeneral->getWells(),
                                                 timeSimbox,
                                                 modelSettings,
                                                 modelGeneral->getPriorFacies(),
                                                 seisCube_,
                                                 inputFiles,
                                                 seismicParameters,
                                                 errText,
                                                 failedPriorCorr);
        }

        if(failedSeismic == false && failedBackground == false &&
          (estimationMode == false || modelSettings->getEstimateWaveletNoise() || modelSettings->getOptimizeWellLocation()))
        {
          modelAVOstatic->addSeismicLogs(modelGeneral->getWells(), seisCube_, modelSettings, numberOfAngles_);
          if(failedReflMat == false && failedExtraSurf == false)
          {
            processWavelets(wavelet_, seisCube_, modelGeneral->getWells(), reflectionMatrix_,
                            timeSimbox, correlationDirection, waveletEstimInterval,
                            modelSettings, inputFiles, errText, failedWavelet);
          }
        }
        if(failedSeismic == false && failedWavelet == false && failedReflMat == false && failedBackground == false &&
          (modelSettings->getOptimizeWellLocation() || modelSettings->getEstimateWaveletNoise() ))
        {
          modelAVOstatic->generateSyntheticSeismic(wavelet_, modelGeneral->getWells(), reflectionMatrix_, timeSimbox, modelSettings, numberOfAngles_);
        }
      }

      if (!failedWells) {
        if(estimationMode || (modelSettings->getWellOutputFlag() & IO::WELLS) > 0)
          modelAVOstatic->writeWells(modelGeneral->getWells(), modelSettings);
        if(estimationMode)
          modelAVOstatic->writeBlockedWells(modelGeneral->getWells(), modelSettings, modelGeneral->getFaciesNames(), modelGeneral->getFaciesLabel());
      }
    }
  }

  if(failedBackground == false) {
    seismicParameters.setBackgroundParameters(background->getAlpha(), background->getBeta(), background->getRho());
    background->releaseGrids();
    delete background;
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

ModelAVODynamic::ModelAVODynamic(ModelSettings          *& modelSettings,
                                 const InputFiles        * inputFiles,
                                 ModelAVOStatic          * modelAVOstatic,
                                 ModelGeneral            * modelGeneral,
                                 SeismicParametersHolder & seismicParameters,
                                 const Simbox            * timeSimbox,
                                 const Surface           * correlationDirection,
                                 const GridMapping       * timeDepthMapping,
                                 const GridMapping       * timeCutMapping,
                                 int                       t)
{ //Time lapse constructor
  seisCube_               = NULL;
  wavelet_                = NULL;
  reflectionMatrix_       = NULL;
  failed_                 = false;
  thisTimeLapse_          = t;
  angularCorr_            = modelSettings->getAngularCorr(thisTimeLapse_);
  SNRatio_                = modelSettings->getSNRatio(thisTimeLapse_);
  useLocalNoise_          = modelSettings->getUseLocalNoise(thisTimeLapse_);
  angle_                  = modelSettings->getAngle(thisTimeLapse_);
  estimateWavelet_        = modelSettings->getEstimateWavelet(thisTimeLapse_);
  matchEnergies_          = modelSettings->getMatchEnergies(thisTimeLapse_);
  numberOfAngles_         = modelSettings->getNumberOfAngles(thisTimeLapse_);

  bool failedWavelet      = false;
  bool failedSeismic      = false;
  bool failedReflMat      = false;
  bool failedPriorCorr    = false;
  bool failedLoadingModel = false;

  Background * background = NULL;

  std::string errText("");

  bool estimationMode = modelSettings->getEstimationMode();

  if(estimationMode == false) {
    FFTGrid * backModel[3];
    backModel[0] = seismicParameters.GetMuAlpha();
    backModel[1] = seismicParameters.GetMuBeta();
    backModel[2] = seismicParameters.GetMuRho();

    background = new Background(backModel);
  }
  if(estimationMode == false || modelSettings->getEstimateWaveletNoise() ){
    processReflectionMatrix(reflectionMatrix_, background, modelGeneral->getWells(), modelSettings,
                              inputFiles, errText, failedReflMat);
  }

  if(failedReflMat == false && (estimationMode == false || modelSettings->getEstimateWaveletNoise() ))
    processSeismic(seisCube_, timeSimbox, timeDepthMapping, timeCutMapping,
                   modelSettings, inputFiles, errText, failedSeismic);

  if(failedSeismic == false && (estimationMode == false || modelSettings->getEstimateWaveletNoise() )){
    modelAVOstatic->addSeismicLogs(modelGeneral->getWells(), seisCube_, modelSettings, numberOfAngles_);
    if(failedReflMat == false)
      processWavelets(wavelet_, seisCube_, modelGeneral->getWells(), reflectionMatrix_,
                      timeSimbox, correlationDirection, modelAVOstatic->getWaveletEstimInterval(),
                      modelSettings, inputFiles, errText, failedWavelet);
  }
  if(failedSeismic == false && failedWavelet == false && failedReflMat == false && modelSettings->getEstimateWaveletNoise() )
    modelAVOstatic->generateSyntheticSeismic(wavelet_, modelGeneral->getWells(), reflectionMatrix_, timeSimbox, modelSettings, numberOfAngles_);

  if(estimationMode || (modelSettings->getWellOutputFlag() & IO::WELLS) > 0)
    modelAVOstatic->writeWells(modelGeneral->getWells(), modelSettings);
  if(estimationMode)
    modelAVOstatic->writeBlockedWells(modelGeneral->getWells(), modelSettings, modelGeneral->getFaciesNames(), modelGeneral->getFaciesLabel());

  if(background != NULL) {
    // New background model has not been made, and background is owned by seismicParameters
    background->releaseGrids();
    delete background;
  }

  failedLoadingModel = failedSeismic || failedPriorCorr || failedReflMat || failedWavelet;

  if (failedLoadingModel) {
    LogKit::WriteHeader("Error(s) while loading data");
    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel;
  failed_details_.push_back(failedSeismic);
  failed_details_.push_back(failedPriorCorr);
  failed_details_.push_back(failedReflMat);
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

  if(seisCube_ != NULL) {
    for(int i=0;i<numberOfAngles_;i++)
      delete seisCube_[i];
    delete [] seisCube_;
  }
}

void
ModelAVODynamic::releaseGrids(void)
{
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
ModelAVODynamic::processSeismic(FFTGrid             **& seisCube,
                                const Simbox          * timeSimbox,
                                const GridMapping     * timeDepthMapping,
                                const GridMapping     * timeCutMapping,
                                const ModelSettings   * modelSettings,
                                const InputFiles      * inputFiles,
                                std::string           & errText,
                                bool                  & failed)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  if(inputFiles->getNumberOfSeismicFiles(thisTimeLapse_) > 0)
  {
    LogKit::WriteHeader("Reading seismic data");

    std::vector<float> offset            = modelSettings->getLocalSegyOffset(thisTimeLapse_);
    const SegyGeometry ** geometry       = new const SegyGeometry * [numberOfAngles_];
    seisCube                             = new FFTGrid * [numberOfAngles_];

    const Simbox * timeCutSimbox = NULL;
    if (timeCutMapping != NULL)
      timeCutSimbox = timeCutMapping->getSimbox(); // For the got-enough-data test
    else
      timeCutSimbox = timeSimbox;

    for (int i = 0 ; i < numberOfAngles_ ; i++) {
      geometry[i] = NULL;
      std::string tmpErrText("");
      std::string fileName = inputFiles->getSeismicFile(thisTimeLapse_,i);
      std::string angle    = NRLib::ToString(angle_[i]*(180/M_PI), 1);
      std::string dataName = "Seismic data angle stack "+angle;
      if(offset[i] < 0)
        offset[i] = modelSettings->getSegyOffset(thisTimeLapse_);

      ModelGeneral::readGridFromFile(fileName,
                                     dataName,
                                     offset[i],
                                     seisCube[i],
                                     geometry[i],
                                     modelSettings->getTraceHeaderFormat(thisTimeLapse_,i),
                                     FFTGrid::DATA,
                                     timeSimbox,
                                     timeCutSimbox,
                                     modelSettings,
                                     tmpErrText);
      if(tmpErrText != "")
      {
        tmpErrText += "\nReading of file \'"+fileName+"\' for "+dataName+" failed.\n";
        errText += tmpErrText;
        failed = true;
      }
      else {
        seisCube[i]->setAngle(angle_[i]);
      }
    }

    LogKit::LogFormatted(LogKit::Low,"\n");

    if(failed == false)
    {
      bool segyVolumesRead = false;
      for (int i = 0 ; i < numberOfAngles_ ; i++)
      {
        if (geometry[i] != NULL)
          segyVolumesRead = true;
      }
      if (segyVolumesRead)
      {
        LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
        LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
        for (int i = 0 ; i < numberOfAngles_ ; i++)
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
        for(int i=0;i<numberOfAngles_;i++) {
          std::string angle    = NRLib::ToString(angle_[i]*(180/M_PI), 1);
          std::string baseName = IO::PrefixOriginalSeismicData() + angle;
          std::string sgriLabel = std::string("Original seismic data for angle stack ") + angle;
          if(offset[i] < 0)
            offset[i] = modelSettings->getSegyOffset(thisTimeLapse_);
          seisCube[i]->writeFile(baseName,
                                 IO::PathToSeismicData(),
                                 timeSimbox,
                                 sgriLabel,
                                 offset[i],
                                 timeDepthMapping,
                                 timeCutMapping,
                                 *modelSettings->getTraceHeaderFormatOutput());
        }
      }
      if((modelSettings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0) {
        for(int i=0;i<numberOfAngles_;i++) {
          std::string angle    = NRLib::ToString(i);
          std::string fileName = IO::makeFullFileName(IO::PathToSeismicData(), IO::FileTemporarySeismic()+angle);
          seisCube[i]->writeCravaFile(fileName, timeSimbox);
        }
      }

      for (int i =0 ; i < numberOfAngles_ ; i++)
        if (geometry[i] != NULL)
          delete geometry[i];
      delete [] geometry;
    }
  }
  Timings::setTimeSeismic(wall,cpu);
}


void
ModelAVODynamic::processBackground(Background                    *& background,
                                   const std::vector<WellData *>  & wells,
                                   const Simbox                   * timeSimbox,
                                   const Simbox                   * timeBGSimbox,
                                   GridMapping                   *& timeDepthMapping,
                                   const GridMapping              * timeCutMapping,
                                   const ModelSettings            * modelSettings,
                                   ModelGeneral                   * modelGeneral,
                                   const InputFiles               * inputFiles,
                                   const int                      & thisTimeLapse,
                                   std::string                    & errText,
                                   bool                           & failed)
{
  if (modelSettings->getForwardModeling())
    LogKit::WriteHeader("Earth Model");
  else
    LogKit::WriteHeader("Prior Expectations / Background Model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  const Simbox * timeCutSimbox = NULL;
  if (timeCutMapping != NULL)
    timeCutSimbox = timeCutMapping->getSimbox(); // For the got-enough-data test
  else
    timeCutSimbox = timeSimbox;

  FFTGrid * backModel[3];
  const int nx    = timeSimbox->getnx();
  const int ny    = timeSimbox->getny();
  const int nz    = timeSimbox->getnz();
  const int nxPad = modelSettings->getNXpad();
  const int nyPad = modelSettings->getNYpad();
  const int nzPad = modelSettings->getNZpad();
  if (modelSettings->getGenerateBackground()) {

    if(modelSettings->getGenerateBackgroundFromRockPhysics() == false) {

      FFTGrid * velocity = NULL;
      std::string backVelFile = inputFiles->getBackVelFile();
      if (backVelFile != ""){
        bool dummy;
        ModelGeneral::loadVelocity(velocity,
                                   timeSimbox,
                                   timeCutSimbox,
                                   modelSettings,
                                   backVelFile,
                                   dummy,
                                   errText,
                                   failed);
      }
      if (!failed) {

        if(modelSettings->getBackgroundVario() == NULL) {
          errText += "There is no variogram available for the background modelling.\n";
          failed = true;
        }
        for (int i=0 ; i<3 ; i++)
        {
          backModel[i] = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());
          backModel[i]->setType(FFTGrid::PARAMETER);
        }

        //if(modelSettings->getMultizoneBackground() == false) //Background() is changed because of CommonData
        //  background = new Background(backModel, wells, velocity, timeSimbox, timeBGSimbox, modelSettings);
        //else
        //  background = new Background(backModel, wells, timeSimbox, modelSettings, inputFiles->getMultizoneSurfaceFiles());

        if(velocity != NULL)
          delete velocity;
      }
    }
    else {

      for (int i=0 ; i<3 ; i++) {
        backModel[i] = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, modelSettings->getFileGrid());
        backModel[i]->createRealGrid();
        backModel[i]->setType(FFTGrid::PARAMETER);
      }


      // Get prior probabilities for the facies in a vector
      std::vector<std::string> facies_names = modelGeneral->getFaciesNames();
      int                      n_facies     = static_cast<int>(facies_names.size());

      std::vector<float> priorProbability = modelGeneral->getPriorFacies();

      std::vector<DistributionsRock *> rock_distribution(n_facies);
      typedef std::map<std::string, DistributionsRock *> rfMapType;
      rfMapType rfMap = modelGeneral->getRockDistributionTime0();

      for(int i=0; i<n_facies; i++) {
        rfMapType::iterator iter = rfMap.find(facies_names[i]);
        if(iter != rfMap.end())
          rock_distribution[i] = iter->second;
      }

      // filling in the backModel in Background
      modelGeneral->generateRockPhysics3DBackground(rock_distribution,
                                                    priorProbability,
                                                    *backModel[0],
                                                    *backModel[1],
                                                    *backModel[2]);

      background = new Background(backModel);

    }
  }
  else {
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
          const float               offset = modelSettings->getSegyOffset(thisTimeLapse);
          std::string errorText("");
          ModelGeneral::readGridFromFile(backFile,
                                         parName[i],
                                         offset,
                                         backModel[i],
                                         dummy1,
                                         dummy2,
                                         FFTGrid::PARAMETER,
                                         timeSimbox,
                                         timeCutSimbox,
                                         modelSettings,
                                         errorText);
          if(errorText != "")
          {
            errorText += "Reading of file '"+backFile+"' for parameter '"+parName[i]+"' failed\n\n";
            errText += errorText;
            failed = true;
          }
          else {
            backModel[i]->calculateStatistics();
            backModel[i]->setUndefinedCellsToGlobalAverage();
            backModel[i]->logTransf();
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
        backModel[i]->calculateStatistics();
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
      if (modelSettings->getUseSIBackground()) { // Vs = SI/Rho     ==> lnVs = lnSI - lnRho
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
ModelAVODynamic::processReflectionMatrix(float               **& reflectionMatrix,
                                         const Background      * background,
                                         const std::vector<WellData *> & wells,
                                         const ModelSettings   * modelSettings,
                                         const InputFiles      * inputFiles,
                                         std::string           & errText,
                                         bool                  & failed)
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
    reflectionMatrix = readMatrix(reflMatrFile, numberOfAngles_, 3, "reflection matrix", tmpErrText);
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
ModelAVODynamic::setupDefaultReflectionMatrix(float             **& reflectionMatrix,
                                              double                vsvp,
                                              const ModelSettings * modelSettings)
{
  int      i;
  float ** A      = new float * [numberOfAngles_];

  // For debugging
  //background->setClassicVsVp();

  double           vsvp2       = vsvp*vsvp;
  std::vector<int> seismicType = modelSettings->getSeismicType(thisTimeLapse_);

  for(i = 0; i < numberOfAngles_; i++)
  {
    double angle = static_cast<double>(angle_[i]);
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(seismicType[i] == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float( (1.0 +tan2t )/2.0 ) ;
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0 );
    }
    else if(seismicType[i] == ModelSettings::PSSEIS) {
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

double ModelAVODynamic::vsvpFromWells(const std::vector<WellData *> & wells,
                                      int                             nWells)
{
  int   N      = 0;
  float VsVp   = 0.0f;

  //Change: Vp/Vs per interval. Only interval with wells.
  //Get interval-names (modelSettings->getIntervalNames()), set setVpVsRatioInterval() per interval.

  for(int i=0 ; i < nWells ; i++) {
    int n = wells[i]->getNumberOfVsVpSamples();
    N    += n;
    VsVp += wells[i]->getMeanVsVp()*n;
  }
  VsVp /= N;

  return static_cast<double>(VsVp);
}

void
ModelAVODynamic::processWavelets(Wavelet                    **& wavelet,
                                 FFTGrid                     ** seisCube,
                                 std::vector<WellData *>        wells,
                                 const float          * const * reflectionMatrix,
                                 const Simbox                 * timeSimbox,
                                 const Surface                * correlationDirection,
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

  wavelet = new Wavelet * [numberOfAngles_];
  localNoiseScale_.resize(numberOfAngles_);

  bool has3Dwavelet = false;
  for(int i=0 ; i < numberOfAngles_ ; i++) {
    localNoiseScale_[i] = NULL;
    if (modelSettings->getWaveletDim(i) == Wavelet::THREE_D)
      has3Dwavelet = true;
    if(estimateWavelet_[i] == true)
      modelSettings->setWaveletScale(thisTimeLapse_,i,1.0);
  }

  unsigned int                      nWells = modelSettings->getNumberOfWells();

  std::vector<std::vector<double> > tGradX(nWells);
  std::vector<std::vector<double> > tGradY(nWells);

  NRLib::Grid2D<float>              refTimeGradX;         ///< Time gradient in x-direction for reference time surface (t0)
  NRLib::Grid2D<float>              refTimeGradY;         ///< Time gradient in x-direction for reference time surface (t0)
  NRLib::Grid2D<float>              structureDepthGradX;  ///< Depth gradient in x-direction for structure ( correlationDirection-t0)*v0/2
  NRLib::Grid2D<float>              structureDepthGradY;  ///< Depth gradient in y-direction for structure ( correlationDirection-t0)*v0/2


  if (has3Dwavelet) {
    if (inputFiles->getRefSurfaceFile() != "") {
      Surface  t0Surf;
      try {
        t0Surf =Surface(inputFiles->getRefSurfaceFile());
      }
      catch (NRLib::Exception & e) {
        errText += e.what();
        failed = true;
      }
       if(!failed)
      {
        double v0=modelSettings->getAverageVelocity();
        computeStructureDepthGradient(v0,
                                      modelSettings->getGradientSmoothingRange(),
                                      timeSimbox,
                                      &t0Surf,
                                      correlationDirection,
                                      structureDepthGradX,
                                      structureDepthGradY);
        Wavelet3D::setGradientMaps(structureDepthGradX,
                                   structureDepthGradY);
        computeReferenceTimeGradient(timeSimbox,
                                     &t0Surf,
                                     refTimeGradX,
                                     refTimeGradY);
      }
      else{
        errText += "Problems reading reference time surface in (x,y).\n";
        error = 1;
      }
    }
    bool estimateWellGradient = modelSettings->getEstimateWellGradientFromSeismic();
    float distance, sigma_m;
    modelSettings->getTimeGradientSettings(distance, sigma_m, thisTimeLapse_);
    std::vector<std::vector<double> > SigmaXY;
    for (unsigned int w=0; w<nWells; w++) {
      BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick();
      if(!estimateWellGradient & ((structureDepthGradX.GetN()> 0) & (structureDepthGradY.GetN()>0))){ // then we allready have the
        double v0=modelSettings->getAverageVelocity();
        bl->setSeismicGradient(v0,structureDepthGradX,structureDepthGradY,refTimeGradX, refTimeGradY,tGradX[w], tGradY[w]);
      }else{
        bl->setTimeGradientSettings(distance, sigma_m);
        bl->findSeismicGradient(seisCube, timeSimbox, numberOfAngles_,tGradX[w], tGradY[w],SigmaXY);
      }
    }
  }

  if (timeSimbox->getdz() > 4.01f && modelSettings->getEstimateNumberOfWavelets(thisTimeLapse_) > 0)
  { // Require this density for wavelet estimation
    LogKit::LogFormatted(LogKit::Low,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
    LogKit::LogFormatted(LogKit::Low,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
    LogKit::LogFormatted(LogKit::Low,"         The number of layers must be increased.                                  \n");
    std::string text("");
    text += "Increase the number of layers to improve the quality of the wavelet estimation.\n";
    text += "   The minimum sampling density is "+NRLib::ToString(timeSimbox->getdz())+", and it should be ";
    text += "lower than 4.0.\n   To obtain the desired density, the number of layers should be at least ";
    text += NRLib::ToString(static_cast<int>(timeSimbox->GetLZ()/4.0))+"\n";
    TaskList::addTask(text);
  }

  // check if local noise is set for some angles.
  bool localNoiseSet = false;
  std::vector<bool> useRickerWavelet = modelSettings->getUseRickerWavelet(thisTimeLapse_);
  for(int i=0 ; i < numberOfAngles_ ; i++) {
    float angle = float(angle_[i]*180.0/M_PI);
    LogKit::LogFormatted(LogKit::Low,"\nAngle stack : %.1f deg",angle);
    if(modelSettings->getForwardModeling()==false)
      seisCube[i]->setAccessMode(FFTGrid::RANDOMACCESS);
    if (modelSettings->getWaveletDim(i) == Wavelet::ONE_D)
      error += process1DWavelet(modelSettings,
                                inputFiles,
                                timeSimbox,
                                //seisCube,
                                wells,
                                //waveletEstimInterval,
                                reflectionMatrix[i],
                                errText,
                                wavelet[i],
                                i,
                                useRickerWavelet[i]);
    else
      error += process3DWavelet(modelSettings,
                                inputFiles,
                                timeSimbox,
                                //seisCube,
                                //wells,
                                //waveletEstimInterval,
                                reflectionMatrix[i],
                                errText,
                                wavelet[i],
                                i);
                                //refTimeGradX,
                                //refTimeGradY,
                                //tGradX,
                                //tGradY

    if(localNoiseScale_[i]!=NULL)
      localNoiseSet = true;
    if(modelSettings->getForwardModeling()==false) // else, no seismic data
      seisCube[i]->endAccess();
  } // end i (angles)

  if(localNoiseSet==true) {
    for(int i=0;i<numberOfAngles_;i++)
      if(localNoiseScale_[i]==NULL)
        localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(),
                                         timeSimbox->getny(),
                                         1.0);
  }

  Timings::setTimeWavelets(wall,cpu);
  failed = error > 0;
}

int
ModelAVODynamic::process1DWavelet(const ModelSettings          * modelSettings,
                                  const InputFiles             * inputFiles,
                                  const Simbox                 * timeSimbox,
                                  //const FFTGrid        * const * seisCube,
                                  std::vector<WellData *>        wells,
                                  //const std::vector<Surface *> & waveletEstimInterval,
                                  const float                  * reflectionMatrix,
                                  std::string                  & errText,
                                  Wavelet                     *& wavelet,
                                  unsigned int                   i,
                                  bool                           useRickerWavelet)
{
  int error = 0;
  Grid2D * shiftGrid(NULL);
  Grid2D * gainGrid(NULL);
  if(modelSettings->getUseLocalWavelet() && inputFiles->getScaleFile(thisTimeLapse_,i) != "") {
      Surface help(inputFiles->getScaleFile(thisTimeLapse_,i));
      gainGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
      resampleSurfaceToGrid2D(timeSimbox, &help, gainGrid);
  }
  if (modelSettings->getUseLocalWavelet() && inputFiles->getShiftFile(thisTimeLapse_,i) != ""){
    Surface helpShift(inputFiles->getShiftFile(thisTimeLapse_,i));
    shiftGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
    resampleSurfaceToGrid2D(timeSimbox, &helpShift, shiftGrid);
  }
  if (useLocalNoise_ && inputFiles->getLocalNoiseFile(thisTimeLapse_,i) != ""){
    Surface helpNoise(inputFiles->getLocalNoiseFile(thisTimeLapse_,i));
    localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(), timeSimbox->getny(), 0.0);
    resampleSurfaceToGrid2D(timeSimbox, &helpNoise, localNoiseScale_[i]);
  }

  //if (estimateWavelet_[i])
  //  wavelet = new Wavelet1D(timeSimbox,
  //                          seisCube[i],
  //                          wells,
  //                          waveletEstimInterval,
  //                          modelSettings,
  //                          reflectionMatrix,
  //                          i,
  //                          error,
  //                          errText);

  else { //Not estimation modus
    if(useRickerWavelet)
        wavelet = new Wavelet1D(modelSettings,
                                reflectionMatrix,
                                angle_[i],
                                modelSettings->getRickerPeakFrequency(thisTimeLapse_,i),
                                error);
    else {
      const std::string & waveletFile = inputFiles->getWaveletFile(thisTimeLapse_,i);
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
                                angle_[i],
                                error,
                                errText);
    }
      // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
      //if(seisCube!=NULL) {// If forward modeling, we have no seismic, can not prescale wavelet.
      //  float       prescale = wavelet->findGlobalScaleForGivenWavelet(modelSettings, timeSimbox, seisCube[i], wells);
      //  const float limHigh  = 3.0f;
      //  const float limLow   = 0.33f;

      //  if(modelSettings->getEstimateGlobalWaveletScale(thisTimeLapse_,i)) // prescale, then we have correct size order, and later scale estimation will be ok.
      //     wavelet->multiplyRAmpByConstant(prescale);
      //  else {
      //    if(modelSettings->getWaveletScale(thisTimeLapse_,i)!= 1.0f && (prescale>limHigh || prescale<limLow)) {
      //       std::string text = "The wavelet given for angle no "+NRLib::ToString(i)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
      //      if(modelSettings->getEstimateLocalScale(thisTimeLapse_,i)) {
      //        errText += text;
      //        error++;
      //      }
      //      else {
      //        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: "+text);
      //        TaskList::addTask("The wavelet is badly scaled. Consider having CRAVA estimate global wavelet scale");
      //      }
      //    }
      //  }
      //}
      if (error == 0)
        wavelet->resample(static_cast<float>(timeSimbox->getdz()),
                          timeSimbox->getnz(),
                          modelSettings->getNZpad());
  }

  if (error == 0) {
    wavelet->scale(modelSettings->getWaveletScale(thisTimeLapse_,i));

    if (modelSettings->getForwardModeling() == false && modelSettings->getNumberOfWells() > 0) {
      //float SNRatio = wavelet->calculateSNRatioAndLocalWavelet(timeSimbox,
      //                                                         seisCube[i],
      //                                                         wells,
      //                                                         modelSettings,
      //                                                         errText,
      //                                                         error,
      //                                                         i,
      //                                                         localNoiseScale_[i],
      //                                                         shiftGrid,
      //                                                         gainGrid,
      //                                                         SNRatio_[i],
      //                                                         modelSettings->getWaveletScale(thisTimeLapse_,i),
      //                                                         modelSettings->getEstimateSNRatio(thisTimeLapse_,i),
      //                                                         modelSettings->getEstimateGlobalWaveletScale(thisTimeLapse_,i),
      //                                                         modelSettings->getEstimateLocalNoise(thisTimeLapse_,i),
      //                                                         modelSettings->getEstimateLocalShift(thisTimeLapse_,i),
      //                                                         modelSettings->getEstimateLocalScale(thisTimeLapse_,i),
      //                                                         estimateWavelet_[i]);

      //if(modelSettings->getEstimateSNRatio(thisTimeLapse_,i))
      //  SNRatio_[i] = SNRatio;
    }

    if (error == 0) {
      if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (modelSettings->getEstimationMode() && estimateWavelet_[i])) {
        std::string type;
        if (estimateWavelet_[i]) {
          type = "Estimated_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,true); // dt_max = 1.0;
        }
        else if (modelSettings->getWaveletScale(thisTimeLapse_,i) == 1.00) {
          type = "";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else {
          type = "Scaled_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
      }
      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((SNRatio_[i] <=SNLow  || SNRatio_[i] > SNHigh) && modelSettings->getForwardModeling()==false) {
        errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio_[i])+" for cube "+NRLib::ToString(i+1)+".\n";
        errText += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }

      bool useLocalNoise = modelSettings->getEstimateLocalNoise(thisTimeLapse_,i) || inputFiles->getLocalNoiseFile(thisTimeLapse_,i) != "";
      bool useLocalShift = modelSettings->getEstimateLocalShift(thisTimeLapse_,i) || inputFiles->getShiftFile(thisTimeLapse_,i)      != "";
      bool useLocalGain  = modelSettings->getEstimateLocalScale(thisTimeLapse_,i) || inputFiles->getScaleFile(thisTimeLapse_,i)      != "";

      if (useLocalNoise)
        readAndWriteLocalGridsToFile(inputFiles->getLocalNoiseFile(thisTimeLapse_,i),
                                     IO::PrefixLocalNoise(),
                                     1.0,  // Scale map with this factor before writing to disk
                                     modelSettings,
                                     i,
                                     timeSimbox,
                                     localNoiseScale_[i]);

      if (useLocalShift) {
        readAndWriteLocalGridsToFile(inputFiles->getShiftFile(thisTimeLapse_,i),
                                     IO::PrefixLocalWaveletShift(),
                                     1.0,
                                     modelSettings,
                                     i,
                                     timeSimbox,
                                     shiftGrid);
        wavelet->setShiftGrid(shiftGrid);
      }

      if (useLocalGain) {
        readAndWriteLocalGridsToFile(inputFiles->getScaleFile(thisTimeLapse_,i),
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
ModelAVODynamic::process3DWavelet(const ModelSettings                     * modelSettings,
                                  const InputFiles                        * inputFiles,
                                  const Simbox                            * timeSimbox,
                                  //const FFTGrid                   * const * seisCube,
                                  //const std::vector<WellData *>           & wells,
                                  //const std::vector<Surface *>            & waveletEstimInterval,
                                  const float                             * reflectionMatrix,
                                  std::string                             & errText,
                                  Wavelet                                *& wavelet,
                                  unsigned int                              i)
{
  int error = 0;
  if (estimateWavelet_[i]) {
    //wavelet = new Wavelet3D(inputFiles->getWaveletFilterFile(i),
    //                        waveletEstimInterval,
    //                        refTimeGradX,
    //                        refTimeGradY,
    //                        tGradX,
    //                        tGradY,
    //                        seisCube[i],
    //                        modelSettings,
    //                        wells,
    //                        timeSimbox,
    //                        reflectionMatrix,
    //                        i,
    //                        error,
    //                        errText);
  }
  else { //Not estimation modus
    const std::string & waveletFile = inputFiles->getWaveletFile(thisTimeLapse_,i);
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
                              angle_[i],
                              error,
                              errText,
                              inputFiles->getWaveletFilterFile(i));
      if (error == 0)
        wavelet->resample(static_cast<float>(timeSimbox->getdz()),
                          timeSimbox->getnz(),
                          modelSettings->getNZpad());
    }
  }
  if ((modelSettings->getEstimationMode() == false) && timeSimbox->getIsConstantThick()) {
    errText += "Simbox with constant thicknessis not implemented for modelling or inversion when 3D wavelet.\n";
    error++;
  }
  if (error == 0) {
    wavelet->scale(modelSettings->getWaveletScale(thisTimeLapse_,i));
/*    bool localEst = (modelSettings->getEstimateLocalScale(thisTimeLapse_,i) ||
                     modelSettings->getEstimateLocalShift(thisTimeLapse_,i) ||
                     modelSettings->getEstimateLocalNoise(thisTimeLapse_,i) ||
                     modelSettings->getEstimateGlobalWaveletScale(thisTimeLapse_,i) ||
                     modelSettings->getEstimateSNRatio(thisTimeLapse_,i)); */

    //if (localEst && modelSettings->getForwardModeling() == false) {
    //  float SNRatio = wavelet->calculateSNRatio(timeSimbox,
    //                                          seisCube[i],
    //                                          wells,
    //                                          modelSettings,
    //                                          errText,
    //                                          error,
    //                                          refTimeGradX,
    //                                          refTimeGradY,
    //                                          tGradX,
    //                                          tGradY,
    //                                          i,
    //                                          SNRatio_[i],
    //                                          modelSettings->getEstimateSNRatio(thisTimeLapse_,i),
    //                                          estimateWavelet_[i]);
    //  if(modelSettings->getEstimateSNRatio(thisTimeLapse_,i))
    //    SNRatio_[i] = SNRatio;

    //}
    if (error == 0) {
      if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (modelSettings->getEstimationMode() && estimateWavelet_[i])) {
        std::string type;
        if (estimateWavelet_[i]) {
          type = "Estimated_";
          if(wavelet->getDim()==1)
            wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,true); // dt_max = 1.0;
          else
            wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else if (modelSettings->getWaveletScale(thisTimeLapse_,i) == 1.00) {
          type = "";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else {
          type = "Scaled_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
      }

      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((SNRatio_[i] <=SNLow  || SNRatio_[i] > SNHigh) && modelSettings->getForwardModeling()==false) {
        errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio_[i])+" for cube "+NRLib::ToString(i+1)+".\n";
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
                                              const Simbox        * timeSimbox,
                                              const Grid2D        * grid)
{
  bool   estimationMode   = modelSettings->getEstimationMode();
  int    outputFormat     = modelSettings->getOutputGridFormat();
  double angle            = angle_[i]*180.0/M_PI;

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
                                         const Simbox         * simbox,
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

void
ModelAVODynamic::computeStructureDepthGradient(double                 v0,
                                               double                 radius,
                                               const Simbox         * timeSimbox,
                                               const Surface        * t0Surf,
                                               const Surface        * correlationDirection_,
                                               NRLib::Grid2D<float> & structureDepthGradX,
                                               NRLib::Grid2D<float> & structureDepthGradY)
 {
   double ds = 12.5;

   int nx = timeSimbox->getnx();
   int ny = timeSimbox->getny();
   structureDepthGradX.Resize(nx,ny);
   structureDepthGradY.Resize(nx,ny);
   double mp=v0*0.001*0.5; // 0.001 is due to s vs ms convension

   for(int i=0;i<nx;i++)
     for(int j=0;j<ny;j++)
     {
       double x,y;
       double gx,gy,gxTmp,gyTmp;
       gx=0.0;
       gy=0.0;
       timeSimbox->getXYCoord(i,j,x,y);
       calculateSmoothGrad( t0Surf, x, y, radius, ds,gxTmp, gyTmp);
       gx=-gxTmp;
       gy=-gyTmp;
       if( correlationDirection_!=NULL){
         calculateSmoothGrad( correlationDirection_, x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp;
         gy+=gyTmp;
       }else{
         calculateSmoothGrad( &(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface())), x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp*0.5;
         gy+=gyTmp*0.5;
         calculateSmoothGrad( &(dynamic_cast<const Surface &> (timeSimbox->GetBotSurface())), x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp*0.5;
         gy+=gyTmp*0.5;
       }

       gx*=mp;
       gy*=mp;
       structureDepthGradX(i,j) =float(gx);
       structureDepthGradY(i,j) =float(gy);
     }
 }

void
ModelAVODynamic::computeReferenceTimeGradient(const Simbox     *timeSimbox,
                              const Surface * t0Surf,
                              NRLib::Grid2D<float> &refTimeGradX,
                              NRLib::Grid2D<float> &refTimeGradY)
 {
   double radius = 50.0;
   double ds = 12.5;
   int nx = timeSimbox->getnx();
   int ny = timeSimbox->getny();
   refTimeGradX.Resize(nx,ny);
   refTimeGradY.Resize(nx,ny);
   for(int i=0;i<nx;i++)
     for(int j=0;j<ny;j++)
     {
       double x,y;
       double gx,gy;
       timeSimbox->getXYCoord(i,j,x,y);
       calculateSmoothGrad( t0Surf, x, y, radius, ds,gx, gy);
       refTimeGradX(i,j) =float(gx);
       refTimeGradY(i,j) =float(gy);
     }
 }


void
ModelAVODynamic::calculateSmoothGrad(const Surface * surf, double x, double y, double radius, double ds,  double& gx, double& gy)
{
  /// Return smoothed Gradient. Computes the gradient as a regression
  /// among points within a given distance from the central point.
  //  Returns missing if central point is outside the grid.
  /// Returns  otherwise it is ok,
  int i,j,k,l;
  int discRadius = int(floor(radius/ds));
  int nPoints    = (2*discRadius+1)*(2*discRadius+1);
  std::vector<double> Z(3*nPoints,0.0);
  std::vector<double> Y(nPoints,0.0);
  std::vector<double> cov(9,0.0);
  std::vector<double> invCov(9,0.0);
  std::vector<double> proj(3,0.0);
  double z0;
  int cy=0;
  int cz=0;

  double baseDepth =surf->GetZ(x,y);


  for(i=- discRadius;i<=discRadius;i++)
    for( j=- discRadius;j<=discRadius;j++)
    {
      double dx= i*ds;
      double dy= j*ds;
      double foo=surf->GetZ(x+dx,y+dy);
      Y[cy] = foo;
      Z[cz] = 1.0;
      Z[cz + 1] = dx;
      Z[cz + 2] = dy;
      cy++;
      cz += 3;
    }

  int nData=0;
  for(i=0; i < nPoints; i++)
  {
    if(!surf->IsMissing(Y[i]))
    {
      nData++;
      for(k=0;k<3;k++)
      {
        for(l=0;l<3;l++)
          cov[k+l*3]+=Z[k + 3*i] * Z[l + 3*i];

        proj[k]+=Z[k + 3*i]*(Y[i]-baseDepth);
      }
    }
  }

  double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  +   cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

  if(det != 0)
  {
      invCov[0] = (cov[4]*cov[8] - cov[5]*cov[7]) / det;
      invCov[1] = (cov[2]*cov[7] - cov[1]*cov[8]) / det;
      invCov[2] = (cov[1]*cov[5] - cov[2]*cov[4]) / det;
      invCov[3] = (cov[5]*cov[6] - cov[3]*cov[8]) / det;
      invCov[4] = (cov[0]*cov[8] - cov[2]*cov[6]) / det;
      invCov[5] = (cov[2]*cov[3] - cov[0]*cov[5]) / det;
      invCov[6] = (cov[3]*cov[7] - cov[4]*cov[6]) / det;
      invCov[7] = (cov[1]*cov[6] - cov[0]*cov[7]) / det;
      invCov[8] = (cov[0]*cov[4] - cov[1]*cov[3]) / det;

      z0 = baseDepth;
      gx = 0.0;
      gy = 0.0;
      for(k=0;k<3;k++)
      {
        z0 += invCov[k]*proj[k]; //NBNB check
        gx += invCov[3+k]*proj[k];
        gy += invCov[6+k]*proj[k];
      }
  }
  else
  {
   gx = RMISSING;
   gy = RMISSING;
  }
}
