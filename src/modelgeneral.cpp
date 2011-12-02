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


ModelGeneral::ModelGeneral(ModelSettings *& modelSettings, const InputFiles * inputFiles, Simbox *& timeBGSimbox)
{
  timeSimbox_             = new Simbox();
  timeSimboxConstThick_   = NULL;

  correlationDirection_   = NULL;
  randomGen_              = NULL;
  failed_                 = false;
  gradX_                  = 0.0;
  gradY_                  = 0.0;

  timeDepthMapping_       = NULL;
  timeCutMapping_         = NULL;
  velocityFromInversion_  = false;

  bool failedSimbox       = false;
  bool failedDepthConv    = false;

  bool failedLoadingModel = false;

  Simbox * timeCutSimbox  = NULL;

  {
    int debugLevel = modelSettings->getLogLevel();
    if(modelSettings->getDebugLevel() == 1)
      debugLevel = LogKit::L_DebugLow;
    else if(modelSettings->getDebugLevel() == 2)
      debugLevel = LogKit::L_DebugHigh;

    LogKit::SetScreenLog(debugLevel);

    std::string logFileName = IO::makeFullFileName("",IO::FileLog()+IO::SuffixTextFiles());
    LogKit::SetFileLog(logFileName,modelSettings->getLogLevel());

    if(modelSettings->getDebugFlag() > 0)
    {
      std::string fName = IO::makeFullFileName("",IO::FileDebug()+IO::SuffixTextFiles());
      LogKit::SetFileLog(fName, debugLevel);
    }

    if(modelSettings->getErrorFileFlag() == true)
    {
      std::string fName = IO::makeFullFileName("",IO::FileError()+IO::SuffixTextFiles());
      LogKit::SetFileLog(fName, LogKit::Error);
    }
    LogKit::EndBuffering();

    if(inputFiles->getSeedFile() == "")
      randomGen_ = new RandomGen(modelSettings->getSeed());
    else
      randomGen_ = new RandomGen(inputFiles->getSeedFile().c_str());

    if(modelSettings->getNumberOfSimulations() == 0)
      modelSettings->setWritePrediction(true); //write predicted grids.

    printSettings(modelSettings, inputFiles);

    //Set output for all FFTGrids.
    FFTGrid::setOutputFlags(modelSettings->getOutputGridFormat(),
                            modelSettings->getOutputGridDomain());

    std::string errText("");

    LogKit::WriteHeader("Defining modelling grid");
    makeTimeSimboxes(timeSimbox_, timeCutSimbox, timeBGSimbox, timeSimboxConstThick_,  //Handles correlation direction too.
                     correlationDirection_, modelSettings, inputFiles,
                     errText, failedSimbox);

    if(!failedSimbox)
    {
      //
      // FORWARD MODELLING
      //
      if (modelSettings->getForwardModeling() == true)
      {
        checkAvailableMemory(timeSimbox_, modelSettings, inputFiles);
      }
      else {
        //
        // INVERSION/ESTIMATION
        //
        if(timeCutSimbox!=NULL)  {
          timeCutMapping_ = new GridMapping();
          timeCutMapping_->makeTimeTimeMapping(timeCutSimbox);
        }

        checkAvailableMemory(timeSimbox_, modelSettings, inputFiles);

        bool estimationMode = modelSettings->getEstimationMode();

        if(estimationMode == false && modelSettings->getDoDepthConversion() == true)
        {
          processDepthConversion(timeCutSimbox, timeSimbox_, modelSettings,
                                 inputFiles, errText, failedDepthConv);
        }

      }
    }
    failedLoadingModel = failedSimbox  || failedDepthConv;

    if (failedLoadingModel) {
      LogKit::WriteHeader("Error(s) while loading data");
      LogKit::LogFormatted(LogKit::Error,"\n"+errText);
      LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
    }
  }

  failed_ = failedLoadingModel;
  failed_details_.push_back(failedSimbox);
  failed_details_.push_back(failedDepthConv);

  if(timeCutSimbox != NULL)
    delete timeCutSimbox;
}


ModelGeneral::~ModelGeneral(void)
{
  if(timeDepthMapping_!=NULL)
    delete timeDepthMapping_;

  if(timeCutMapping_!=NULL)
    delete timeCutMapping_;

  if(correlationDirection_ !=NULL)
    delete correlationDirection_;

  delete randomGen_;
  delete timeSimbox_;
  delete timeSimboxConstThick_;

}

void
ModelGeneral::checkAvailableMemory(Simbox        * timeSimbox,
                                   ModelSettings * modelSettings,
                                   const InputFiles    * inputFiles)
{
  LogKit::WriteHeader("Estimating amount of memory needed");
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
      if(modelSettings->getUseLocalNoise() == true || (modelSettings->getEstimateFaciesProb() && modelSettings->getFaciesProbRelative()))
        baseP += nGridBackground;
      int baseU = 0;
      if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        baseU += modelSettings->getNumberOfFacies();

      //First peak: At inversion
      int peak1P = baseP + nGridSeismicData; //Need seismic data as well here.
      int peak1U = baseU;

      long long int peakGridMem = peak1P*gridSizePad + peak1U*gridSizeBase; //First peak must be currently largest.
      int peakNGrid   = peak1P;                                             //Also in number of padded grids

      if(modelSettings->getNumberOfSimulations() > 0) { //Second possible peak when simulating.
        int peak2P = baseP + 3; //Three extra parameter grids for simulated parameters.
        if(modelSettings->getUseLocalNoise() == true &&
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
  float mem2        = static_cast<float>(modelSettings->getNumberOfAngles())*gridSizePad + memOneSeis; //Peak memory when reading seismic, overestimated.

  float neededMem   = mem0 + std::max(mem1, mem2);

  float megaBytes   = neededMem/(1024.f*1024.f);
  float gigaBytes   = megaBytes/1024.f;

  LogKit::LogFormatted(LogKit::High,"\nMemory needed for reading seismic data       : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding internal grids (%2d): %10.2f MB\n",nGrids, mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding other entities     : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (gigaBytes < 0.01f)
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.2f megaBytes\n",megaBytes);
  else
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.2f gigaBytes\n",gigaBytes);

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

int
ModelGeneral::readSegyFile(const std::string       & fileName,
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
ModelGeneral::readStormFile(const std::string  & fName,
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
ModelGeneral::setPaddingSize(int nx, double px)
{
  int leastint    = static_cast<int>(ceil(nx*(1.0f+px)));
  int closestprod = FFTGrid::findClosestFactorableNumber(leastint);
  return(closestprod);
}

void
ModelGeneral::makeTimeSimboxes(Simbox   *& timeSimbox,
                        Simbox          *& timeCutSimbox,
                        Simbox          *& timeBGSimbox,
                        Simbox          *& timeSimboxConstThick,
                        Surface         *& correlationDirection,
                        ModelSettings   *& modelSettings,
                        const InputFiles * inputFiles,
                        std::string      & errText,
                        bool             & failed)
{
  std::string gridFile("");

  int  areaSpecification      = modelSettings->getAreaSpecification();

  bool estimationModeNeedILXL = modelSettings->getEstimationMode() &&
                                (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA ||
                                (modelSettings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0 ||
                                (modelSettings->getOutputGridFormat() & IO::SEGY) > 0);

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
    LogKit::LogFormatted(LogKit::High,"\nArea information has been taken from model file\n");
    areaType = "Model file";
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_SURFACE)
  {
    LogKit::LogFormatted(LogKit::High,"\nFinding area information from surface \'"+inputFiles->getAreaSurfaceFile()+"\'\n");
    areaType = "Surface";
    Surface surf(inputFiles->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    modelSettings->setAreaParameters(&geometry);
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_GRID_DATA)
  {
    LogKit::LogFormatted(LogKit::High,"\nFinding inversion area from grid data in file \'"+gridFile+"\'\n");
    areaType = "Grid data";
    std::string tmpErrText;
    SegyGeometry * geometry;
    getGeometryFromGridOnFile(gridFile,
                              modelSettings->getTraceHeaderFormat(0),
                              geometry,
                              tmpErrText);

    if(geometry != NULL) {
      geometry->WriteGeometry(true);
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
      LogKit::LogFormatted(LogKit::Low,"\nResolution                x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
      double azimuth = (-1)*timeSimbox->getAngle()*(180.0/M_PI);
      if (azimuth < 0)
        azimuth += 360.0;
      LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                           areaType.c_str(),
                           timeSimbox->getx0(), timeSimbox->gety0(),
                           timeSimbox->getlx(), timeSimbox->getly(), azimuth,
                           timeSimbox->getdx(), timeSimbox->getdy());
    }

    float minHorRes = modelSettings->getMinHorizontalRes();
    if (timeSimbox->getdx() < minHorRes || timeSimbox->getdy() < minHorRes){
      failed = true;
      errText += "The horizontal resolution in dx and dy should be above "+NRLib::ToString(minHorRes)+" m. \n";
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
            LogKit::LogFormatted(LogKit::High,"\nFinding IL/XL information from grid data file \'"+gridFile+"\'\n");
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
            LogKit::LogFormatted(LogKit::High,"\nCannot extract IL/XL information from non-SEGY grid data file \'"+gridFile+"\'\n");
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

      setSimboxSurfaces(timeSimbox,
                        inputFiles->getTimeSurfFiles(),
                        modelSettings,
                        errText,
                        failed);

      if(!failed)
      {
        if(modelSettings->getUseLocalWavelet() && timeSimbox->getIsConstantThick())
        {
          LogKit::LogFormatted(LogKit::Warning,"\nWarning: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");
          TaskList::addTask("If local wavelet is to be used, constant thickness in depth should be removed.");
        }


        int status = timeSimbox->calculateDz(modelSettings->getLzLimit(),errText);
        estimateZPaddingSize(timeSimbox, modelSettings);

        float minSampDens = modelSettings->getMinSamplingDensity();
        if (timeSimbox->getdz()*timeSimbox->getMinRelThick() < minSampDens){
          failed   = true;
          errText += "The sampling density should be above "+NRLib::ToString(minSampDens)+" ms.\nPlease decrease the number of layers. \n";
        }

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
              if(timeSimbox->CheckSurface(tmpSurf) == true)
                correlationDirection = new Surface(tmpSurf);
              else {
                errText += "Error: Correlation surface does not cover volume.\n";
                failed = true;
              }
            }
            catch (NRLib::Exception & e) {
              errText += e.what();
              failed = true;
            }

            if(failed == false && modelSettings->getForwardModeling() == false) {
              //Extends timeSimbox for correlation coverage. Original stored in timeCutSimbox
              setupExtendedTimeSimbox(timeSimbox, correlationDirection,
                                      timeCutSimbox,
                                      modelSettings->getOutputGridFormat(),
                                      modelSettings->getOutputGridDomain(),
                                      modelSettings->getOtherOutputFlag());
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

            if(modelSettings->getForwardModeling() == false && failed == false) {
              setupExtendedBackgroundSimbox(timeSimbox, correlationDirection, timeBGSimbox,
                                            modelSettings->getOutputGridFormat(),
                                            modelSettings->getOutputGridDomain(),
                                            modelSettings->getOtherOutputFlag());
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
            unsigned long long int gridsize = (timeSimbox->getnx()+modelSettings->getNXpad())*(timeSimbox->getny()+modelSettings->getNYpad())*(timeSimbox->getnz()+modelSettings->getNZpad());

            if(gridsize > std::numeric_limits<unsigned int>::max())
            {
              failed = true;
              errText+= "Grid size is too large.  Reduce resolution or extend the grid.\n";
            }

            LogKit::LogFormatted(LogKit::Low,"\nTime simulation grids:\n");
            LogKit::LogFormatted(LogKit::Low,"  Output grid         %4i * %4i * %4i   : %10i\n",
                                 timeSimbox->getnx(),timeSimbox->getny(),timeSimbox->getnz(),
                                 timeSimbox->getnx()*timeSimbox->getny()*timeSimbox->getnz());
            LogKit::LogFormatted(LogKit::Low,"  FFT grid            %4i * %4i * %4i   : %10i\n",
                                 modelSettings->getNXpad(),modelSettings->getNYpad(),modelSettings->getNZpad(),
                                 modelSettings->getNXpad()*modelSettings->getNYpad()*modelSettings->getNZpad());
          }

          //
          // Make time simbox with constant thicknesses (needed for log filtering and facies probabilities)
          //
          timeSimboxConstThick = new Simbox(timeSimbox);
          Surface tsurf(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface()));
          timeSimboxConstThick->setDepth(tsurf, 0, timeSimbox->getlz(), timeSimbox->getdz());

          if((modelSettings->getOtherOutputFlag() & IO::EXTRA_SURFACES) > 0 && (modelSettings->getOutputGridDomain() & IO::TIMEDOMAIN) > 0) {
            std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_ConstThick";
            std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_ConstThick";
            timeSimboxConstThick->writeTopBotGrids(topSurf,
                                                   baseSurf,
                                                   IO::PathToInversionResults(),
                                                   modelSettings->getOutputGridFormat());
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
ModelGeneral::logIntervalInformation(const Simbox      * simbox,
                                     const std::string & header_text1,
                                     const std::string & header_text2)
{
  LogKit::LogFormatted(LogKit::Low,"\n"+header_text1+"\n");
  double zmin, zmax;
  simbox->getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::Low," %13s          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       header_text2.c_str(),
                       zmin+simbox->getlz()*simbox->getAvgRelThick()*0.5,
                       zmin,zmax);
  LogKit::LogFormatted(LogKit::Low,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       simbox->getlz()*simbox->getAvgRelThick(),
                       simbox->getlz()*simbox->getMinRelThick(),
                       simbox->getlz());
  LogKit::LogFormatted(LogKit::Low,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n",
                       simbox->getdz()*simbox->getAvgRelThick(),
                       simbox->getdz(),
                       simbox->getdz()*simbox->getMinRelThick());
}

void
ModelGeneral::setSimboxSurfaces(Simbox                        *& simbox,
                                const std::vector<std::string> & surfFile,
                                ModelSettings                  * modelSettings,
                                std::string                    & errText,
                                bool                           & failed)
{
  const std::string & topName = surfFile[0];

  bool   generateSeismic    = modelSettings->getForwardModeling();
  bool   estimationMode     = modelSettings->getEstimationMode();
  bool   generateBackground = modelSettings->getGenerateBackground();
  bool   parallelSurfaces   = modelSettings->getParallelTimeSurfaces();
  int    nz                 = modelSettings->getTimeNz();
  int    outputFormat       = modelSettings->getOutputGridFormat();
  int    outputDomain       = modelSettings->getOutputGridDomain();
  int    outputGridsElastic = modelSettings->getOutputGridsElastic();
  int    outputGridsOther   = modelSettings->getOutputGridsOther();
  int    outputGridsSeismic = modelSettings->getOutputGridsSeismic();
  double dTop               = modelSettings->getTimeDTop();
  double lz                 = modelSettings->getTimeLz();
  double dz                 = modelSettings->getTimeDz();

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
      simbox->setDepth(*z0Grid, dTop, lz, dz, modelSettings->getRunFromPanel());
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
          simbox->setDepth(*z0Grid, *z1Grid, nz, modelSettings->getRunFromPanel());
        }
        catch (NRLib::Exception & e) {
          errText += e.what();
          std::string text(std::string("Seismic data"));
          writeAreas(modelSettings->getAreaParameters(),simbox,text);
          failed = true;
        }
      }
    }
    if (!failed) {
      if((outputDomain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        simbox->setTopBotName(topSurf,baseSurf,outputFormat);
        if (generateSeismic) {
          simbox->writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   outputFormat);
        }
        else if (!estimationMode){
          if (outputGridsElastic > 0 || outputGridsOther > 0 || outputGridsSeismic > 0)
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToInversionResults(),
                                     outputFormat);
        }
        if((outputFormat & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((outputGridsElastic & IO::BACKGROUND) > 0 ||
              (outputGridsElastic & IO::BACKGROUND_TREND) > 0 ||
              (estimationMode && generateBackground)) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToBackground(),
                                     outputFormat);
          }
          if ((outputGridsOther & IO::CORRELATION) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     outputFormat);
          }
          if ((outputGridsSeismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     outputFormat);
          }
          if ((outputGridsOther & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     outputFormat);
          }
        }
      }
    }
  }
  delete z0Grid;
  delete z1Grid;
}

void
ModelGeneral::setupExtendedTimeSimbox(Simbox   * timeSimbox,
                                      Surface  * corrSurf,
                                      Simbox  *& timeCutSimbox,
                                      int        outputFormat,
                                      int        outputDomain,
                                      int        otherOutput)
{
  timeCutSimbox = new Simbox(timeSimbox);
  double * corrPlanePars = findPlane(corrSurf);
  Surface * meanSurf;
  if(corrSurf->GetNI() > 2)
    meanSurf = new Surface(*corrSurf);
  else {
    meanSurf = new Surface(dynamic_cast<const Surface &>(timeSimbox->GetTopSurface()));
    if(meanSurf->GetNI() == 2) { //Extend corrSurf to cover other surfaces.
      double minX = meanSurf->GetXMin();
      double maxX = meanSurf->GetXMax();
      double minY = meanSurf->GetYMin();
      double maxY = meanSurf->GetYMax();
      if(minX > corrSurf->GetXMin())
        minX = corrSurf->GetXMin();
      if(maxX < corrSurf->GetXMax())
        maxX = corrSurf->GetXMax();
      if(minY > corrSurf->GetYMin())
        minY = corrSurf->GetYMin();
      if(maxY < corrSurf->GetYMax())
        maxY = corrSurf->GetYMax();
      corrSurf->SetDimensions(minX, minY, maxX-minX, maxY-minY);
    }
  }
  int i;
  for(i=0;i<static_cast<int>(meanSurf->GetN());i++)
    (*meanSurf)(i) = 0;

  meanSurf->AddNonConform(&(timeSimbox->GetTopSurface()));
  meanSurf->AddNonConform(&(timeSimbox->GetBotSurface()));
  meanSurf->Multiply(0.5);
  double * refPlanePars = findPlane(meanSurf);

  for(i=0;i<3;i++)
    refPlanePars[i] -= corrPlanePars[i];
  gradX_ = refPlanePars[1];
  gradY_ = refPlanePars[2];

  Surface * refPlane = createPlaneSurface(refPlanePars, meanSurf);
  delete meanSurf;
  meanSurf = NULL;

  std::string fileName = "Correlation_Rotation_Plane";
  IO::writeSurfaceToFile(*refPlane, fileName, IO::PathToCorrelations(), outputFormat);

  refPlane->AddNonConform(corrSurf);
  delete [] corrPlanePars;
  delete [] refPlanePars;

  Surface topSurf(*refPlane);
  topSurf.SubtractNonConform(&(timeSimbox->GetTopSurface()));
  double shiftTop = topSurf.Max();
  shiftTop *= -1.0;
  topSurf.Add(shiftTop);
  topSurf.AddNonConform(&(timeSimbox->GetTopSurface()));

  Surface botSurf(*refPlane);
  botSurf.SubtractNonConform(&(timeSimbox->GetBotSurface()));
  double shiftBot = botSurf.Min();
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
    LogKit::LogFormatted(LogKit::High,"\nNumber of layers in inversion increased from %d",timeCutSimbox->getnz());
    LogKit::LogFormatted(LogKit::High," to %d in grid created using correlation direction.\n",nz);
  }
  botSurf.Add(shiftBot);
  botSurf.AddNonConform(&(timeSimbox->GetBotSurface()));

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
ModelGeneral::setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
                                            Surface  * corrSurf,
                                            Simbox  *& timeBGSimbox,
                                            int        outputFormat,
                                            int        outputDomain,
                                            int        otherOutput)
{
  //
  // Move correlation surface for easier handling.
  //
  Surface tmpSurf(*corrSurf);
  double avg = tmpSurf.Avg();
  if (avg > 0)
    tmpSurf.Subtract(avg);
  else
    tmpSurf.Add(avg); // This situation is not very likely, but ...

  //
  // Find top surface of background simbox.
  //
  // The funny/strange dTop->Multiply(-1.0) is due to NRLIB's current
  // inability to set dTop equal to Simbox top surface.
  //
  Surface dTop(tmpSurf);
  dTop.SubtractNonConform(&(timeSimbox->GetTopSurface()));
  dTop.Multiply(-1.0);
  double shiftTop = dTop.Min();
  Surface topSurf(tmpSurf);
  topSurf.Add(shiftTop);

  //
  // Find base surface of background simbox
  //
  Surface dBot(tmpSurf);
  dBot.SubtractNonConform(&(timeSimbox->GetBotSurface()));
  dBot.Multiply(-1.0);
  double shiftBot = dBot.Max();
  Surface botSurf(tmpSurf);
  botSurf.Add(shiftBot);

  //
  // Calculate number of layers of background simbox
  //
  tmpSurf.Assign(0.0);
  tmpSurf.AddNonConform(&botSurf);
  tmpSurf.SubtractNonConform(&topSurf);
  double dMax = tmpSurf.Max();
  double dt = timeSimbox->getdz();
  int nz;
  //
  // NBNB-PAL: I think it is a good idea to use a maximum dt of 10ms.
  //
  //if (dt < 10.0) {
  //  LogKit::LogFormatted(LogKit::High,"\nReducing sampling density for background",dt);
  //  LogKit::LogFormatted(LogKit::High," modelling from %.2fms to 10.0ms\n");
  //  dt = 10.0;  // A sampling density of 10.0ms is good enough for BG model
  // }
  nz = static_cast<int>(ceil(dMax/dt));

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
ModelGeneral::findPlane(Surface * surf)
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
ModelGeneral::createPlaneSurface(double * planeParams, Surface * templateSurf)
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
ModelGeneral::estimateXYPaddingSizes(Simbox         * timeSimbox,
                                     ModelSettings *& modelSettings)
{
  double dx      = timeSimbox->getdx();
  double dy      = timeSimbox->getdy();
  double lx      = timeSimbox->getlx();
  double ly      = timeSimbox->getly();
  int    nx      = timeSimbox->getnx();
  int    ny      = timeSimbox->getny();
  int    nz      = timeSimbox->getnz();

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
  int nzPad = modelSettings->getNZpad();

  double true_xPadFac = static_cast<double>(nxPad - nx)/static_cast<double>(nx);
  double true_yPadFac = static_cast<double>(nyPad - ny)/static_cast<double>(ny);
  double true_zPadFac = modelSettings->getZPadFac();
  double true_xPad    = true_xPadFac*lx;
  double true_yPad    = true_yPadFac*ly;
  double true_zPad    = true_zPadFac*(timeSimbox->getlz()*timeSimbox->getMinRelThick());

  modelSettings->setNXpad(nxPad);
  modelSettings->setNYpad(nyPad);
  modelSettings->setXPadFac(true_xPadFac);
  modelSettings->setYPadFac(true_yPadFac);

  std::string text1;
  std::string text2;
  int logLevel = LogKit::Medium;
  if (modelSettings->getEstimateXYPadding()) {
    text1 = " estimated from lateral correlation ranges in internal grid";
    logLevel = LogKit::Low;
  }
  if (modelSettings->getEstimateZPadding()) {
    text2 = " estimated from an assumed wavelet length";
    logLevel = LogKit::Low;
  }

  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text1+":\n");
  LogKit::LogFormatted(logLevel,"  xPad, xPadFac, nx, nxPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_xPad, true_xPadFac, nx, nxPad);
  LogKit::LogFormatted(logLevel,"  yPad, yPadFac, ny, nyPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_yPad, true_yPadFac, ny, nyPad);
  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text2+":\n");
  LogKit::LogFormatted(logLevel,"  zPad, zPadFac, nz, nzPad                 : %5.fms, %5.3f, %5d, %4d\n",
                       true_zPad, true_zPadFac, nz, nzPad);
}

void
ModelGeneral::estimateZPaddingSize(Simbox         * timeSimbox,
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
  int nzPad        = setPaddingSize(nz, zPadFac);
  zPadFac          = static_cast<double>(nzPad - nz)/static_cast<double>(nz);

  modelSettings->setNZpad(nzPad);
  modelSettings->setZPadFac(zPadFac);
}

void
ModelGeneral::readGridFromFile(const std::string       & fileName,
                               const std::string       & parName,
                               const float               offset,
                               FFTGrid                *& grid,
                               const SegyGeometry     *& geometry,
                               const TraceHeaderFormat * format,
                               int                       gridType,
                               Simbox                  * timeSimbox,
                               ModelSettings           * modelSettings,
                               int                     & outsideTraces,
                               std::string             & errText,
                               bool                      nopadding)
{
  int fileType = IO::findGridType(fileName);

  outsideTraces = 0;
  if(fileType == IO::CRAVA)
  {
    int nxPad, nyPad, nzPad;
    if(nopadding)
    {
      nxPad = timeSimbox->getnx();
      nyPad = timeSimbox->getny();
      nzPad = timeSimbox->getnz();
    }
    else
    {
      nxPad = modelSettings->getNXpad();
      nyPad = modelSettings->getNYpad();
      nzPad = modelSettings->getNZpad();
    }
    LogKit::LogFormatted(LogKit::Low,"\nReading grid \'"+parName+"\' from file "+fileName);
    grid = createFFTGrid(timeSimbox->getnx(),
                         timeSimbox->getny(),
                         timeSimbox->getnz(),
                         nxPad,
                         nyPad,
                         nzPad,
                         modelSettings->getFileGrid());

    grid->setType(gridType);
    grid->readCravaFile(fileName, errText, nopadding);
  }
  else if(fileType == IO::SEGY)
    outsideTraces = readSegyFile(fileName, grid, timeSimbox, modelSettings, geometry,
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
}

void
ModelGeneral::printSettings(ModelSettings     * modelSettings,
                            const InputFiles  * inputFiles)
{
  LogKit::WriteHeader("Model settings");

  LogKit::LogFormatted(LogKit::Low,"\nGeneral settings:\n");
  if(modelSettings->getForwardModeling()==true)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : forward\n");
  else if (modelSettings->getEstimationMode()==true)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : estimation\n");

  else if (modelSettings->getNumberOfSimulations() == 0)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : prediction\n");
  else
  {
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : simulation\n");
    if(inputFiles->getSeedFile()=="") {
      if (modelSettings->getSeed() == 0)
        LogKit::LogFormatted(LogKit::Low,"  Seed                                     :          0 (default seed)\n");
      else
        LogKit::LogFormatted(LogKit::Low,"  Seed                                     : %10d\n",modelSettings->getSeed());
    }
    else
      LogKit::LogFormatted(LogKit::Low,"  Seed read from file                      : %10s\n",inputFiles->getSeedFile().c_str());


    LogKit::LogFormatted(LogKit::Low,"  Number of realisations                   : %10d\n",modelSettings->getNumberOfSimulations());
  }
  if(modelSettings->getForwardModeling()==false)
  {
    LogKit::LogFormatted(LogKit::Low,"  Kriging                                  : %10s\n",(modelSettings->getKrigingParameter()>0 ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Facies probabilities                     : %10s\n",(modelSettings->getEstimateFaciesProb() ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Synthetic seismic                        : %10s\n",(modelSettings->getGenerateSeismicAfterInv() ? "yes" : "no" ));
  }

  if (modelSettings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::Low,"\nSettings for facies probability estimation:\n");
    LogKit::LogFormatted(LogKit::Low,"  Use elastic parameters relative to trend : %10s\n",(modelSettings->getFaciesProbRelative()     ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Include Vs information in estimation     : %10s\n",(modelSettings->getNoVsFaciesProb()         ? "no"  : "yes"));
    LogKit::LogFormatted(LogKit::Low,"  Use filtered well logs for estimation    : %10s\n",(modelSettings->getUseFilterForFaciesProb() ? "yes" : "no"));
  }

  LogKit::LogFormatted(LogKit::Low,"\nInput/Output settings:\n");
  std::string logText("*NONE*");
  int logLevel = modelSettings->getLogLevel();
  if (logLevel == LogKit::L_Error)
    logText = "ERROR";
  else if (logLevel == LogKit::L_Warning)
    logText = "WARNING";
  else if (logLevel == LogKit::L_Low)
    logText = "LOW";
  else if (logLevel == LogKit::L_Medium)
    logText = "MEDIUM";
  else if (logLevel == LogKit::L_High)
    logText = "HIGH";
  else if (logLevel == LogKit::L_DebugLow)
     logText = "DEBUGLOW";
  else if (logLevel == LogKit::L_DebugHigh)
    logText = "DEBUGHIGH";
  LogKit::LogFormatted(LogKit::Low, "  Log level                                : %10s\n",logText.c_str());
  if (inputFiles->getInputDirectory() != "")
    LogKit::LogFormatted(LogKit::High,"  Input directory                          : %10s\n",inputFiles->getInputDirectory().c_str());
  if (IO::getOutputPath() != "")
    LogKit::LogFormatted(LogKit::High,"  Output directory                         : %10s\n",IO::getOutputPath().c_str());

  int gridFormat         = modelSettings->getOutputGridFormat();
  int gridDomain         = modelSettings->getOutputGridDomain();
  int outputGridsOther   = modelSettings->getOutputGridsOther();
  int outputGridsElastic = modelSettings->getOutputGridsElastic();
  int outputGridsSeismic = modelSettings->getOutputGridsSeismic();

  if (outputGridsElastic > 0  || outputGridsSeismic > 0  || outputGridsOther > 0) {
    LogKit::LogFormatted(LogKit::Medium,"\nGrid output formats:\n");
    if (gridFormat & IO::SEGY) {
      const std::string & formatName = modelSettings->getTraceHeaderFormatOutput()->GetFormatName();
      LogKit::LogFormatted(LogKit::Medium,"  Segy - %-10s                        :        yes\n",formatName.c_str());
    }
    if (gridFormat & IO::STORM)
      LogKit::LogFormatted(LogKit::Medium,"  Storm                                    :        yes\n");
    if (gridFormat & IO::ASCII)
      LogKit::LogFormatted(LogKit::Medium,"  ASCII                                    :        yes\n");
    if (gridFormat & IO::SGRI)
      LogKit::LogFormatted(LogKit::Medium,"  Norsar                                   :        yes\n");
    if (gridFormat & IO::CRAVA)
      LogKit::LogFormatted(LogKit::Medium,"  Crava                                    :        yes\n");

    LogKit::LogFormatted(LogKit::Medium,"\nGrid output domains:\n");
    if (gridDomain & IO::TIMEDOMAIN)
      LogKit::LogFormatted(LogKit::Medium,"  Time                                     :        yes\n");
    if (gridDomain & IO::DEPTHDOMAIN)
      LogKit::LogFormatted(LogKit::Medium,"  Depth                                    :        yes\n");
  }

  if (outputGridsElastic > 0 &&
      modelSettings->getForwardModeling() == false) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of elastic parameters:\n");
    if ((outputGridsElastic & IO::VP) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Pressure-wave velocity  (Vp)             :        yes\n");
    if ((outputGridsElastic & IO::VS) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Shear-wave velocity  (Vs)                :        yes\n");
    if ((outputGridsElastic & IO::RHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Density  (Rho)                           :        yes\n");
    if ((outputGridsElastic & IO::AI) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Acoustic impedance  (AI)                 :        yes\n");
    if ((outputGridsElastic & IO::VPVSRATIO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Vp/Vs ratio                              :        yes\n");
    if ((outputGridsElastic & IO::SI) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Shear impedance  (SI)                    :        yes\n");
    if ((outputGridsElastic & IO::MURHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  MuRho  (SI*SI)                           :        yes\n");
    if ((outputGridsElastic & IO::LAMBDARHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  LambdaRho  (AI*AI - 2*SI*SI)             :        yes\n");
    if ((outputGridsElastic & IO::LAMELAMBDA) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Lame's first parameter                   :        yes\n");
    if ((outputGridsElastic & IO::LAMEMU) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Lame's second parameter (shear modulus)  :        yes\n");
    if ((outputGridsElastic & IO::POISSONRATIO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Poisson ratio  (X-1)/2(X-2), X=(Vp/Vs)^2 :        yes\n");
    if ((outputGridsElastic & IO::BACKGROUND) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Background (Vp, Vs, Rho)                 :        yes\n");
    if ((outputGridsElastic & IO::BACKGROUND_TREND) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Background trend (Vp, Vs, Rho)           :        yes\n");
  }

  if (modelSettings->getForwardModeling() ||
      outputGridsSeismic > 0) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of seismic data:\n");
    if ((outputGridsSeismic & IO::SYNTHETIC_SEISMIC_DATA) > 0 || modelSettings->getForwardModeling())
      LogKit::LogFormatted(LogKit::Medium,"  Synthetic seismic data (forward modelled):        yes\n");
    if ((outputGridsSeismic & IO::ORIGINAL_SEISMIC_DATA) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Original seismic data (in output grid)   :        yes\n");
    if ((outputGridsSeismic & IO::RESIDUAL) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Seismic data residuals                   :        yes\n");
  }

  if (modelSettings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of facies probability volumes:\n");
    if ((outputGridsOther & IO::FACIESPROB) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Facies probabilities                     :        yes\n");
    if ((outputGridsOther & IO::FACIESPROB_WITH_UNDEF) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Facies probabilities with undefined value:        yes\n");
  }

  if ((outputGridsOther & IO::CORRELATION)>0 ||
      (outputGridsOther & IO::EXTRA_GRIDS)  >0 ||
      (outputGridsOther & IO::TIME_TO_DEPTH_VELOCITY)>0) {
    LogKit::LogFormatted(LogKit::Medium,"\nOther grid output:\n");
    if ((outputGridsOther & IO::CORRELATION) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Posterior correlations                   :        yes\n");
    if ((outputGridsOther & IO::EXTRA_GRIDS) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Help grids (see use manual)              :        yes\n");
    if ((outputGridsOther & IO::TIME_TO_DEPTH_VELOCITY) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Time-to-depth velocity                   :        yes\n");
  }

  // NBNB-PAL: Vi får utvide testen nedenfor etter hvert...
  if (modelSettings->getFileGrid()) {
    LogKit::LogFormatted(LogKit::Medium,"\nAdvanced settings:\n");
    LogKit::LogFormatted(LogKit::Medium, "  Use intermediate disk storage for grids  :        yes\n");
  }

  LogKit::LogFormatted(LogKit::High,"\nUnit settings/assumptions:\n");
  LogKit::LogFormatted(LogKit::High,"  Time                                     : %10s\n","ms TWT");
  LogKit::LogFormatted(LogKit::High,"  Frequency                                : %10s\n","Hz");
  LogKit::LogFormatted(LogKit::High,"  Length                                   : %10s\n","m");
  LogKit::LogFormatted(LogKit::High,"  Velocities                               : %10s\n","m/s");
  LogKit::LogFormatted(LogKit::High,"  Density                                  : %10s\n","g/cm3");
  LogKit::LogFormatted(LogKit::High,"  Angles                                   : %10s\n","   degrees (clockwise relative to north when applicable)");

  //
  // WELL PROCESSING
  //
  if (modelSettings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::High,"\nSettings for well processing:\n");
    LogKit::LogFormatted(LogKit::High,"  Threshold for merging log entries        : %10.2f ms\n",modelSettings->getMaxMergeDist());
    LogKit::LogFormatted(LogKit::High,"  Threshold for Vp-Vs rank correlation     : %10.2f\n",modelSettings->getMaxRankCorr());
    LogKit::LogFormatted(LogKit::High,"  Threshold for deviation angle            : %10.1f (=%.2fm/ms TWT)\n",
                         modelSettings->getMaxDevAngle(),tan(modelSettings->getMaxDevAngle()*M_PI/180.0));
    LogKit::LogFormatted(LogKit::High,"  High cut for background modelling        : %10.1f\n",modelSettings->getMaxHzBackground());
    LogKit::LogFormatted(LogKit::High,"  High cut for seismic resolution          : %10.1f\n",modelSettings->getMaxHzSeismic());
  }
  LogKit::LogFormatted(LogKit::High,"\nRange of allowed parameter values:\n");
  LogKit::LogFormatted(LogKit::High,"  Vp  - min                                : %10.0f\n",modelSettings->getAlphaMin());
  LogKit::LogFormatted(LogKit::High,"  Vp  - max                                : %10.0f\n",modelSettings->getAlphaMax());
  LogKit::LogFormatted(LogKit::High,"  Vs  - min                                : %10.0f\n",modelSettings->getBetaMin());
  LogKit::LogFormatted(LogKit::High,"  Vs  - max                                : %10.0f\n",modelSettings->getBetaMax());
  LogKit::LogFormatted(LogKit::High,"  Rho - min                                : %10.1f\n",modelSettings->getRhoMin());
  LogKit::LogFormatted(LogKit::High,"  Rho - max                                : %10.1f\n",modelSettings->getRhoMax());

  //
  // WELL DATA
  //
  if (modelSettings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::Low,"\nWell logs:\n");
    const std::vector<std::string> & logNames = modelSettings->getLogNames();

    if (logNames.size() > 0)
    {
      LogKit::LogFormatted(LogKit::Low,"  Time                                     : %10s\n",  logNames[0].c_str());
      if(NRLib::Uppercase(logNames[1])=="VP" ||
         NRLib::Uppercase(logNames[1])=="LFP_VP")
        LogKit::LogFormatted(LogKit::Low,"  p-wave velocity                          : %10s\n",logNames[1].c_str());
      else
        LogKit::LogFormatted(LogKit::Low,"  Sonic                                    : %10s\n",logNames[1].c_str());
      if(NRLib::Uppercase(logNames[3])=="VS" ||
         NRLib::Uppercase(logNames[3])=="LFP_VS")
        LogKit::LogFormatted(LogKit::Low,"  s-wave velocity                          : %10s\n",logNames[3].c_str());
      else
        LogKit::LogFormatted(LogKit::Low,"  Shear sonic                              : %10s\n",logNames[3].c_str());
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10s\n",  logNames[2].c_str());
      if (modelSettings->getFaciesLogGiven())
        LogKit::LogFormatted(LogKit::Low,"  Facies                                   : %10s\n",logNames[4].c_str());
    }
    else
    {
      LogKit::LogFormatted(LogKit::Low,"  Time                                     : %10s\n","TWT");
      LogKit::LogFormatted(LogKit::Low,"  Sonic                                    : %10s\n","DT");
      LogKit::LogFormatted(LogKit::Low,"  Shear sonic                              : %10s\n","DTS");
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10s\n","RHOB");
      LogKit::LogFormatted(LogKit::Low,"  Facies                                   : %10s\n","FACIES");
    }
    LogKit::LogFormatted(LogKit::Low,"\nWell files:\n");
    for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++)
    {
      LogKit::LogFormatted(LogKit::Low,"  %-2d                                       : %s\n",i+1,inputFiles->getWellFile(i).c_str());
    }
    bool generateBackground = modelSettings->getGenerateBackground();
    bool estimateFaciesProb = modelSettings->getFaciesLogGiven();
    bool estimateWavelet    = false;
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
      estimateWavelet = estimateWavelet || modelSettings->getEstimateWavelet(i);
    if (generateBackground || estimateFaciesProb || estimateWavelet)
    {
      LogKit::LogFormatted(LogKit::Low,"\nUse well in estimation of:                   ");
      if (generateBackground) LogKit::LogFormatted(LogKit::Low,"BackgroundTrend  ");
      if (estimateWavelet)    LogKit::LogFormatted(LogKit::Low,"WaveletEstimation  ");
      if (estimateFaciesProb) LogKit::LogFormatted(LogKit::Low,"FaciesProbabilities");
      LogKit::LogFormatted(LogKit::Low,"\n");
      for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++)
      {
        LogKit::LogFormatted(LogKit::Low,"  %-2d                                       : ",i+1);
        if (generateBackground) {
          if (modelSettings->getIndicatorBGTrend(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","yes");
          else if (modelSettings->getIndicatorBGTrend(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","yes");
        }
        if (estimateWavelet) {
          if (modelSettings->getIndicatorWavelet(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-13s  ","yes");
          else if (modelSettings->getIndicatorWavelet(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-13s  ","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-13s  ","if possible");
        }
        if (estimateFaciesProb) {
          if (modelSettings->getIndicatorFacies(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-12s","yes");
          else if (modelSettings->getIndicatorFacies(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-12s","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-12s","if possible");
        }
        LogKit::LogFormatted(LogKit::Low,"\n");
      }
    }
    if ( modelSettings->getOptimizeWellLocation() )
    {
      LogKit::LogFormatted(LogKit::Low,"\nFor well, optimize position for            : Angle with Weight\n");
      for (int i = 0 ; i < modelSettings->getNumberOfWells() ; i++)
      {
        int nMoveAngles = modelSettings->getNumberOfWellAngles(i);
        if( nMoveAngles > 0 )
        {
          LogKit::LogFormatted(LogKit::Low," %2d %46.1f %10.1f\n",i+1,(modelSettings->getWellMoveAngle(i,0)*180/M_PI),modelSettings->getWellMoveWeight(i,0));
          for (int j=1; j<nMoveAngles; j++)
            LogKit::LogFormatted(LogKit::Low," %49.1f %10.1f\n",(modelSettings->getWellMoveAngle(i,j)*180/M_PI),modelSettings->getWellMoveWeight(i,j));
        }
        LogKit::LogFormatted(LogKit::Low,"\n");
      }
    }
  }

  //
  // AREA
  //
  std::string gridFile;
  int areaSpecification = modelSettings->getAreaSpecification();
  if(modelSettings->getForwardModeling()) {
    LogKit::LogFormatted(LogKit::Low,"\nSeismic area:\n");
    gridFile = inputFiles->getBackFile(0);    // Get geometry from earth model (Vp)
  }
  else {
    LogKit::LogFormatted(LogKit::Low,"\nInversion area");
    if(areaSpecification == ModelSettings::AREA_FROM_GRID_DATA)
      gridFile = inputFiles->getSeismicFile(0); // Get area from first seismic data volume
  }
  if (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA) {
    const std::vector<int> & areaILXL = modelSettings->getAreaILXL();
    LogKit::LogFormatted(LogKit::Low," taken from grid\n");
    LogKit::LogFormatted(LogKit::Low,"  Grid                                     : "+gridFile+"\n");
    if(areaILXL.size()>0)
    {
    if (areaILXL[0] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  In-line start                            : %10d\n", areaILXL[0]);
    if (areaILXL[1] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  In-line end                              : %10d\n", areaILXL[1]);
    if (areaILXL[4] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  In-line step                             : %10d\n", areaILXL[4]);
    if (areaILXL[2] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  Cross-line start                         : %10d\n", areaILXL[2]);
    if (areaILXL[3] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  Cross-line end                           : %10d\n", areaILXL[3]);
    if (areaILXL[5] != IMISSING)
      LogKit::LogFormatted(LogKit::Low,"  Cross-line step                          : %10d\n", areaILXL[5]);
    }
  }
  else if (areaSpecification == ModelSettings::AREA_FROM_UTM) {
    LogKit::LogFormatted(LogKit::Low," given as UTM coordinates\n");
    const SegyGeometry * geometry = modelSettings->getAreaParameters();
    LogKit::LogFormatted(LogKit::Low,"  Reference point x                        : %10.1f\n", geometry->GetX0());
    LogKit::LogFormatted(LogKit::Low,"  Reference point y                        : %10.1f\n", geometry->GetY0());
    LogKit::LogFormatted(LogKit::Low,"  Length x                                 : %10.1f\n", geometry->Getlx());
    LogKit::LogFormatted(LogKit::Low,"  Length y                                 : %10.1f\n", geometry->Getly());
    LogKit::LogFormatted(LogKit::Low,"  Sample density x                         : %10.1f\n", geometry->GetDx());
    LogKit::LogFormatted(LogKit::Low,"  Sample density y                         : %10.1f\n", geometry->GetDy());
    LogKit::LogFormatted(LogKit::Low,"  Rotation                                 : %10.4f\n", geometry->GetAngle()*(180.0/NRLib::Pi)*(-1));
  }
  else if (areaSpecification == ModelSettings::AREA_FROM_SURFACE) {
    LogKit::LogFormatted(LogKit::Low," taken from surface\n");
    LogKit::LogFormatted(LogKit::Low,"  Reference surface                        : "+inputFiles->getAreaSurfaceFile()+"\n");
  }

  //
  // SURFACES
  //
  LogKit::LogFormatted(LogKit::Low,"\nTime surfaces:\n");
  if (modelSettings->getParallelTimeSurfaces())
  {
    LogKit::LogFormatted(LogKit::Low,"  Surface                                  : "+inputFiles->getTimeSurfFile(0)+"\n");
    LogKit::LogFormatted(LogKit::Low,"  Shift to top surface                     : %10.1f\n", modelSettings->getTimeDTop());
    LogKit::LogFormatted(LogKit::Low,"  Time slice                               : %10.1f\n", modelSettings->getTimeLz());
    LogKit::LogFormatted(LogKit::Low,"  Sampling density                         : %10.1f\n", modelSettings->getTimeDz());
    LogKit::LogFormatted(LogKit::Low,"  Number of layers                         : %10d\n",   int(modelSettings->getTimeLz()/modelSettings->getTimeDz()+0.5));
  }
  else
  {
    const std::string & topName  = inputFiles->getTimeSurfFile(0);
    const std::string & baseName = inputFiles->getTimeSurfFile(1);

    if (NRLib::IsNumber(topName))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(topName.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Top surface                              : "+topName+"\n");

    if (NRLib::IsNumber(baseName))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n", atof(baseName.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Base surface                             : "+baseName+"\n");
      LogKit::LogFormatted(LogKit::Low,"  Number of layers                         : %10d\n", modelSettings->getTimeNz());

    LogKit::LogFormatted(LogKit::Low,"  Minimum allowed value for lmin/lmax      : %10.2f\n", modelSettings->getLzLimit());
  }
  if (inputFiles->getCorrDirFile() != "")
    LogKit::LogFormatted(LogKit::Low,"\n  Correlation direction                    : "+inputFiles->getCorrDirFile()+"\n");

  if (modelSettings->getDoDepthConversion())
  {
    LogKit::LogFormatted(LogKit::Low,"\nDepth conversion:\n");
    if (inputFiles->getDepthSurfFile(0) != "")
      LogKit::LogFormatted(LogKit::Low,"  Top depth surface                        : "+inputFiles->getDepthSurfFile(0)+"\n");
    else
      LogKit::LogFormatted(LogKit::Low,"  Top depth surface                        : %s\n", "Made from base depth surface and velocity field");
    if (inputFiles->getDepthSurfFile(1) != "")
      LogKit::LogFormatted(LogKit::Low,"  Base depth surface                       : "+inputFiles->getDepthSurfFile(1)+"\n");
    else
      LogKit::LogFormatted(LogKit::Low,"  Base depth surface                       : %s\n", "Made from top depth surface and velocity field");
    std::string velocityField = inputFiles->getVelocityField();
    if (modelSettings->getVelocityFromInversion()) {
      velocityField = "Use Vp from inversion";
    }
     LogKit::LogFormatted(LogKit::Low,"  Velocity field                           : "+velocityField+"\n");
  }

  const std::string & topWEI  = inputFiles->getWaveletEstIntFile(0);
  const std::string & baseWEI = inputFiles->getWaveletEstIntFile(1);

  if (topWEI != "" || baseWEI != "") {
    LogKit::LogFormatted(LogKit::Low,"\nWavelet estimation interval:\n");
    if (NRLib::IsNumber(topWEI))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(topWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : "+topWEI+"\n");

    if (NRLib::IsNumber(baseWEI))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n",atof(baseWEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : "+baseWEI+"\n");
  }

  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" || baseFEI != "") {
    LogKit::LogFormatted(LogKit::Low,"\nFacies estimation interval:\n");
    if (NRLib::IsNumber(topFEI))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(topFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : "+topFEI+"\n");

    if (NRLib::IsNumber(baseFEI))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n",atof(baseFEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : "+baseFEI+"\n");
  }

  //
  // BACKGROUND
  //
  if (modelSettings->getGenerateBackground())
  {
    LogKit::LogFormatted(LogKit::Low,"\nBackground model (estimated):\n");
    if (inputFiles->getBackVelFile() != "")
      LogKit::LogFormatted(LogKit::Low,"  Trend for p-wave velocity                : "+inputFiles->getBackVelFile()+"\n");
    Vario       * vario  = modelSettings->getBackgroundVario();
    GenExpVario * pVario = dynamic_cast<GenExpVario*>(vario);
    LogKit::LogFormatted(LogKit::Low,"  Variogram\n");
    LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(vario->getType()).c_str());
    if (pVario != NULL)
    LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pVario->getPower());
    LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",vario->getRange());
    if (vario->getAnisotropic())
    {
      LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",vario->getSubRange());
      LogKit::LogFormatted(LogKit::Low,"    Azimuth                                : %10.1f\n",90.0 - vario->getAngle()*(180/M_PI));
    }
    LogKit::LogFormatted(LogKit::Low,"  High cut frequency for well logs         : %10.1f\n",modelSettings->getMaxHzBackground());
  }
  else
  {
    if(modelSettings->getForwardModeling()==true)
      LogKit::LogFormatted(LogKit::Low,"\nEarth model:\n");
    else
      LogKit::LogFormatted(LogKit::Low,"\nBackground model:\n");
    if (modelSettings->getConstBackValue(0) > 0)
      LogKit::LogFormatted(LogKit::Low,"  P-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(0));
    else
      LogKit::LogFormatted(LogKit::Low,"  P-wave velocity read from file           : "+inputFiles->getBackFile(0)+"\n");

    if (modelSettings->getConstBackValue(1) > 0)
      LogKit::LogFormatted(LogKit::Low,"  S-wave velocity                          : %10.1f\n",modelSettings->getConstBackValue(1));
    else
      LogKit::LogFormatted(LogKit::Low,"  S-wave velocity read from file           : "+inputFiles->getBackFile(1)+"\n");

    if (modelSettings->getConstBackValue(2) > 0)
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10.1f\n",modelSettings->getConstBackValue(2));
    else
      LogKit::LogFormatted(LogKit::Low,"  Density read from file                   : "+inputFiles->getBackFile(2)+"\n");
  }

  TraceHeaderFormat * thf_old = modelSettings->getTraceHeaderFormat();
  if (thf_old != NULL)
  {
    LogKit::LogFormatted(LogKit::Low,"\nAdditional SegY trace header format:\n");
    if (thf_old != NULL) {
      LogKit::LogFormatted(LogKit::Low,"  Format name                              : "+thf_old->GetFormatName()+"\n");
      if (thf_old->GetBypassCoordScaling())
        LogKit::LogFormatted(LogKit::Low,"  Bypass coordinate scaling                :        yes\n");
      if (!thf_old->GetStandardType())
      {
        LogKit::LogFormatted(LogKit::Low,"  Start pos coordinate scaling             : %10d\n",thf_old->GetScalCoLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos trace x coordinate             : %10d\n",thf_old->GetUtmxLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos trace y coordinate             : %10d\n",thf_old->GetUtmyLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos inline index                   : %10d\n",thf_old->GetInlineLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos crossline index                : %10d\n",thf_old->GetCrosslineLoc());
        LogKit::LogFormatted(LogKit::Low,"  Coordinate system                        : %10s\n",thf_old->GetCoordSys()==0 ? "UTM" : "ILXL" );
      }
    }
  }

  if (modelSettings->getForwardModeling())
  {
    //
    // SEISMIC
    //
    LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for seismic:\n");
    LogKit::LogFormatted(LogKit::Low,"  Generating seismic                       : %10s\n","yes");
    for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
    {
      LogKit::LogFormatted(LogKit::Low,"\nSettings for AVO stack %d:\n",i+1);
      LogKit::LogFormatted(LogKit::Low,"  Angle                                    : %10.1f\n",(modelSettings->getAngle(i)*180/M_PI));
      LogKit::LogFormatted(LogKit::Low,"  Read wavelet from file                   : "+inputFiles->getWaveletFile(i)+"\n");
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
      LogKit::LogFormatted(LogKit::Low,"\nPrior correlation (of residuals):\n");
      LogKit::LogFormatted(LogKit::Low,"  Range of allowed parameter values:\n");
      LogKit::LogFormatted(LogKit::Low,"    Var{Vp}  - min                         : %10.1e\n",modelSettings->getVarAlphaMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vp}  - max                         : %10.1e\n",modelSettings->getVarAlphaMax());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vs}  - min                         : %10.1e\n",modelSettings->getVarBetaMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vs}  - max                         : %10.1e\n",modelSettings->getVarBetaMax());
      LogKit::LogFormatted(LogKit::Low,"    Var{Rho} - min                         : %10.1e\n",modelSettings->getVarRhoMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Rho} - max                         : %10.1e\n",modelSettings->getVarRhoMax());
      LogKit::LogFormatted(LogKit::Low,"  Lateral correlation:\n");
      LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",corr->getRange());
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",corr->getSubRange());
        LogKit::LogFormatted(LogKit::Low,"    Azimuth                                : %10.1f\n",90.0 - corr->getAngle()*(180/M_PI));
      }
    }
    //
    // PRIOR FACIES
    //
    if (modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE ||
        modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        // Can not be written when FACIES_FROM_WELLS as this information not is extracted yet
    {
      LogKit::LogFormatted(LogKit::Low,"\nPrior facies probabilities:\n");
      if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
      {
        typedef std::map<std::string,float> mapType;
        mapType myMap = modelSettings->getPriorFaciesProb();

        for(mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::Low,"   %-12s                            : %10.2f\n",(i->first).c_str(),i->second);
      }
      else if (modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
      {
        typedef std::map<std::string,std::string> mapType;
        mapType myMap = inputFiles->getPriorFaciesProbFile();

        for(mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::Low,"   %-12s                            : %10s\n",(i->first).c_str(),(i->second).c_str());
      }
    }
    //
    // SEISMIC
    //
    if (modelSettings->getNoSeismicNeeded()==false)
    {
      LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for seismic:\n");
      LogKit::LogFormatted(LogKit::Low,"  White noise component                    : %10.2f\n",modelSettings->getWNC());
      LogKit::LogFormatted(LogKit::Low,"  Low cut for inversion                    : %10.1f\n",modelSettings->getLowCut());
      LogKit::LogFormatted(LogKit::Low,"  High cut for inversion                   : %10.1f\n",modelSettings->getHighCut());
      corr  = modelSettings->getAngularCorr();
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::LogFormatted(LogKit::Low,"  Angular correlation:\n");
      LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",corr->getRange()*180.0/M_PI);
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",corr->getSubRange()*180.0/M_PI);
        LogKit::LogFormatted(LogKit::Low,"    Angle                                  : %10.1f\n",corr->getAngle());
      }
      bool estimateNoise = false;
      for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++) {
        estimateNoise = estimateNoise || modelSettings->getEstimateSNRatio(i);
      }
      LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for wavelet:\n");
      if (estimateNoise)
        LogKit::LogFormatted(LogKit::Low,"  Maximum shift in noise estimation        : %10.1f\n",modelSettings->getMaxWaveletShift());
      LogKit::LogFormatted(LogKit::Low,"  Minimum relative amplitude               : %10.3f\n",modelSettings->getMinRelWaveletAmp());
      LogKit::LogFormatted(LogKit::Low,"  Wavelet tapering length                  : %10.1f\n",modelSettings->getWaveletTaperingL());
      if (modelSettings->getOptimizeWellLocation()) {
        LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for well locations:\n");
        LogKit::LogFormatted(LogKit::Low,"  Maximum offset                           : %10.1f\n",modelSettings->getMaxWellOffset());
        LogKit::LogFormatted(LogKit::Low,"  Maximum vertical shift                   : %10.1f\n",modelSettings->getMaxWellShift());
      }
      for (int i = 0 ; i < modelSettings->getNumberOfAngles() ; i++)
      {
        LogKit::LogFormatted(LogKit::Low,"\nSettings for AVO stack %d:\n",i+1);
        LogKit::LogFormatted(LogKit::Low,"  Angle                                    : %10.1f\n",(modelSettings->getAngle(i)*180/M_PI));
        LogKit::LogFormatted(LogKit::Low,"  SegY start time                          : %10.1f\n",modelSettings->getSegyOffset());
        TraceHeaderFormat * thf = modelSettings->getTraceHeaderFormat(i);
        if (thf != NULL)
        {
          LogKit::LogFormatted(LogKit::Low,"  SegY trace header format:\n");
          LogKit::LogFormatted(LogKit::Low,"    Format name                            : "+thf->GetFormatName()+"\n");
          if (thf->GetBypassCoordScaling())
            LogKit::LogFormatted(LogKit::Low,"    Bypass coordinate scaling              :        yes\n");
          if (!thf->GetStandardType())
          {
            LogKit::LogFormatted(LogKit::Low,"    Start pos coordinate scaling           : %10d\n",thf->GetScalCoLoc());
            LogKit::LogFormatted(LogKit::Low,"    Start pos trace x coordinate           : %10d\n",thf->GetUtmxLoc());
            LogKit::LogFormatted(LogKit::Low,"    Start pos trace y coordinate           : %10d\n",thf->GetUtmyLoc());
            LogKit::LogFormatted(LogKit::Low,"    Start pos inline index                 : %10d\n",thf->GetInlineLoc());
            LogKit::LogFormatted(LogKit::Low,"    Start pos crossline index              : %10d\n",thf->GetCrosslineLoc());
            LogKit::LogFormatted(LogKit::Low,"    Coordinate system                      : %10s\n",thf->GetCoordSys()==0 ? "UTM" : "ILXL" );
          }
        }
        LogKit::LogFormatted(LogKit::Low,"  Data                                     : "+inputFiles->getSeismicFile(i)+"\n");
        if (modelSettings->getEstimateWavelet(i))
          LogKit::LogFormatted(LogKit::Low,"  Estimate wavelet                         : %10s\n", "yes");
        else
          LogKit::LogFormatted(LogKit::Low,"  Read wavelet from file                   : "+inputFiles->getWaveletFile(i)+"\n");
        if (modelSettings->getEstimateLocalShift(i))
          LogKit::LogFormatted(LogKit::Low,"  Estimate local shift map                 : %10s\n", "yes");
        else if (inputFiles->getShiftFile(i) != "")
          LogKit::LogFormatted(LogKit::Low,"  Local shift map                          : "+inputFiles->getShiftFile(i)+"\n");
        if (modelSettings->getEstimateLocalScale(i))
          LogKit::LogFormatted(LogKit::Low,"  Estimate local scale map                 : %10s\n", "yes");
        else if (inputFiles->getScaleFile(i) != "")
          LogKit::LogFormatted(LogKit::Low,"  Local scale map                          : "+inputFiles->getScaleFile(i)+"\n");
        if (modelSettings->getMatchEnergies(i))
          LogKit::LogFormatted(LogKit::Low,"  Match empirical and theoretical energies : %10s\n", "yes");
        if (!modelSettings->getEstimateWavelet(i) && !modelSettings->getMatchEnergies(i)) {
          if (modelSettings->getEstimateGlobalWaveletScale(i))
            LogKit::LogFormatted(LogKit::Low,"  Estimate global wavelet scale            : %10s\n","yes");
          else
            LogKit::LogFormatted(LogKit::Low,"  Global wavelet scale                     : %10.2f\n",modelSettings->getWaveletScale(i));
        }
        if (modelSettings->getEstimateSNRatio(i))
          LogKit::LogFormatted(LogKit::Low,"  Estimate signal-to-noise ratio           : %10s\n", "yes");
        else
          LogKit::LogFormatted(LogKit::Low,"  Signal-to-noise ratio                    : %10.1f\n",modelSettings->getSNRatio(i));
        if (modelSettings->getEstimateLocalNoise(i)) {
          if (inputFiles->getLocalNoiseFile(i) == "")
            LogKit::LogFormatted(LogKit::Low,"  Estimate local signal-to-noise ratio map : %10s\n", "yes");
          else
            LogKit::LogFormatted(LogKit::Low,"  Local signal-to-noise ratio map          : "+inputFiles->getLocalNoiseFile(i)+"\n");
        }
        if (modelSettings->getEstimateLocalNoise(i))
          LogKit::LogFormatted(LogKit::Low,"  Estimate local noise                     : %10s\n", "yes");
        if (inputFiles->getLocalNoiseFile(i) != "")
          LogKit::LogFormatted(LogKit::Low,"  Local noise                              : "+inputFiles->getLocalNoiseFile(i)+"\n");
      }
    }
  }
}

void
ModelGeneral::getCorrGradIJ(float & corrGradI, float &corrGradJ) const
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
ModelGeneral::processDepthConversion(Simbox            * timeCutSimbox,
                                     Simbox            * timeSimbox,
                                     ModelSettings     * modelSettings,
                                     const InputFiles  * inputFiles,
                                     std::string       & errText,
                                     bool              & failed)
{
  FFTGrid * velocity = NULL;
  if(timeCutSimbox != NULL)
    loadVelocity(velocity, timeCutSimbox, modelSettings,
                 inputFiles->getVelocityField(), velocityFromInversion_,
                 errText, failed);
  else
    loadVelocity(velocity, timeSimbox, modelSettings,
                 inputFiles->getVelocityField(), velocityFromInversion_,
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
                                        modelSettings->getOutputGridFormat(),
                                        failed, errText);            // NBNB-PAL: Er dettet riktig nz (timeCut vs time)?
      timeDepthMapping_->makeTimeDepthMapping(velocity, timeSimbox);
      velocity->endAccess();

      if((modelSettings->getOutputGridsOther() & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
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
                                        modelSettings->getOutputGridFormat(),
                                        failed,
                                        errText);

    }
  }
  if(velocity != NULL)
    delete velocity;
}

void
ModelGeneral::loadVelocity(FFTGrid          *& velocity,
                           Simbox            * timeSimbox,
                           ModelSettings     * modelSettings,
                           const std::string & velocityField,
                           bool              & velocityFromInversion,
                           std::string       & errText,
                           bool              & failed)
{
  LogKit::WriteHeader("Setup time-to-depth relationship");

  if(modelSettings->getVelocityFromInversion() == true)
  {
    velocityFromInversion = true;
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
    int outsideTraces = 0;
    readGridFromFile(velocityField,
                     "velocity field",
                     offset,
                     velocity,
                     dummy1,
                     dummy2,
                     FFTGrid::PARAMETER,
                     timeSimbox,
                     modelSettings,
                     outsideTraces,
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
ModelGeneral::writeAreas(const SegyGeometry * areaParams,
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

  LogKit::LogFormatted(LogKit::Low,"\nThe top and/or base time surfaces do not cover the area specified by the "+text);
  LogKit::LogFormatted(LogKit::Low,"\nPlease extrapolate surfaces or specify a smaller AREA in the model file.\n");
  LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx        ly     azimuth          dx      dy\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
  double azimuth = (-1)*areaRot*(180.0/M_PI);
  if (azimuth < 0)
    azimuth += 360.0;
  LogKit::LogFormatted(LogKit::Low,"Model area       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n\n",
                       areaX0, areaY0, areaLx, areaLy, areaDx, areaDy, azimuth);

  LogKit::LogFormatted(LogKit::Low,"Area                    xmin         xmax           ymin        ymax\n");
  LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %11.2f %11.2f\n",
                       text.c_str(),areaXmin, areaXmax, areaYmin, areaYmax);
  const NRLib::Surface<double> & top  = timeSimbox->GetTopSurface();
  const NRLib::Surface<double> & base = timeSimbox->GetBotSurface();
  LogKit::LogFormatted(LogKit::Low,"Top surface      %11.2f  %11.2f    %11.2f %11.2f\n",
                       top.GetXMin(), top.GetXMax(), top.GetYMin(), top.GetYMax());
  LogKit::LogFormatted(LogKit::Low,"Base surface     %11.2f  %11.2f    %11.2f %11.2f\n",
                       base.GetXMin(), base.GetXMax(), base.GetYMin(), base.GetYMax());
}

void
ModelGeneral::findSmallestSurfaceGeometry(const double   x0,
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
ModelGeneral::getGeometryFromGridOnFile(const std::string         gridFile,
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
ModelGeneral::geometryFromCravaFile(const std::string & fileName)
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
ModelGeneral::geometryFromStormFile(const std::string & fileName,
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
    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
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
    int    nx      = static_cast<int>(stormgrid->GetNI());
    int    ny      = static_cast<int>(stormgrid->GetNJ());
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
ModelGeneral::createFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid)
{
  FFTGrid* fftGrid;

  if(fileGrid)
    fftGrid =  new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fftGrid =  new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

  return(fftGrid);
}
