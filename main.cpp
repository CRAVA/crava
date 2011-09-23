#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#include "lib/timekit.hpp"
#include "lib/utils.h"

#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/program.h"
#include "src/definitions.h"
#include "src/corr.h"
#include "src/wavelet.h"
#include "src/crava.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/welldata.h"
#include "src/filterwelllogs.h"
#include "src/timings.h"
#include "src/spatialwellfilter.h"
#include "src/tasklist.h"

#include "src/xmlmodelfile.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/modelgeneral.h"

int main(int argc, char** argv)
{
  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::SetScreenLog(LogKit::L_Low);
  LogKit::StartBuffering();

  Program program( 1,                     // Major version
                   1,                     // Minor version
                   0,                     // Patch number
                   //"",                    // Use empty string "" for release versions
                   " beta",               // Use empty string "" for release versions
                   -1,                    // Validity of licence in days (-1 = infinite)
                  "Norsk Regnesentral");  // Who this copy of CRAVA is licensed to

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  try
  {
    {
    InputFiles * inputFiles;
    XmlModelFile modelFile(argv[1]);
    inputFiles     = modelFile.getInputFiles();
    ModelSettings * modelSettings = modelFile.getModelSettings();

    bool failedModelFile    = false;
    bool failedInputFiles   = false;
    bool failedLoadingModel = false;

    ModelGeneral    * modelGeneral    = NULL;
    ModelAVOStatic  * modelAVOstatic  = NULL;
    ModelAVODynamic * modelAVOdynamic = NULL;

    if (modelFile.getParsingFailed()) {
      failedModelFile = true;
    }

    std::string errTxt = inputFiles->addInputPathAndCheckFiles();
    if(errTxt != "") {
      LogKit::WriteHeader("Error opening files");
      LogKit::LogMessage(LogKit::Error, "\n"+errTxt);
      LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
      failedInputFiles = true;
    }

    if(!failedModelFile && !failedInputFiles){
      // Construct ModelGeneral object first.
      // For each data type, construct the static model class before the dynamic.
      Simbox * timeBGSimbox   = NULL;
      modelGeneral    = new ModelGeneral(modelSettings, inputFiles, timeBGSimbox);
      modelAVOstatic  = new ModelAVOStatic(modelSettings,
                                           inputFiles,
                                           modelGeneral->getFailedDetails(),
                                           modelGeneral->getTimeSimbox(),
                                           timeBGSimbox,
                                           modelGeneral->getTimeSimboxConstThick(),
                                           modelGeneral->getRandomGen());
      // Wells are adjusted by ModelAVODynamic constructor.
      modelAVOdynamic = new ModelAVODynamic(modelSettings, inputFiles,
                                            modelGeneral->getFailedDetails(),
                                            modelAVOstatic->getFailedDetails(),
                                            modelGeneral->getTimeSimbox(),
                                            timeBGSimbox,
                                            modelGeneral->getRandomGen(),
                                            modelGeneral->getTimeDepthMapping(),
                                            modelGeneral->getTimeCutMapping(),
                                            modelAVOstatic->getWaveletEstimInterval(),
                                            modelAVOstatic->getWellMoveInterval(),
                                            modelAVOstatic->getFaciesEstimInterval(),
                                            modelAVOstatic);
      if(timeBGSimbox != NULL)
        delete timeBGSimbox;
    }
    else {
      LogKit::SetFileLog(IO::FileLog()+IO::SuffixTextFiles(), modelSettings->getLogLevel());
      LogKit::EndBuffering();
    }
    delete inputFiles;
    failedLoadingModel =  modelGeneral    == NULL || modelGeneral->getFailed()   ||
                          modelAVOstatic  == NULL || modelAVOstatic->getFailed() ||
                          modelAVOdynamic == NULL || modelAVOdynamic->getFailed();

    if(failedModelFile || failedInputFiles || failedLoadingModel)
        return(1);

    Crava         * crava         = NULL;

    if(!modelSettings->getForwardModeling())
    {
      if (modelSettings->getDoInversion())
      {
        time_t timestart, timeend;
        time(&timestart);

        int nwells = modelSettings->getNumberOfWells();
        SpatialWellFilter *spatwellfilter = new SpatialWellFilter(nwells);
        crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter);

        std::string warningText("");

        if(crava->getWarning( warningText ) != 0)
         {
           LogKit::LogFormatted(LogKit::Low,"\nWarning  !!!\n");
           LogKit::LogFormatted(LogKit::Low,"%s",warningText.c_str());
           LogKit::LogFormatted(LogKit::Low,"\n");
         }
        crava->printEnergyToScreen();

        time(&timeend);
        LogKit::LogFormatted(LogKit::DebugLow,"\nTime elapsed :  %d\n",timeend-timestart);
        crava->computePostMeanResidAndFFTCov();
        time(&timeend);
        LogKit::LogFormatted(LogKit::DebugLow,"\nTime elapsed :  %d\n",timeend-timestart);

        if(modelSettings->getNumberOfSimulations() > 0)
        {
          crava->simulate(modelGeneral->getRandomGen());
        }

        Corr * corr = modelAVOdynamic->getCorrelations();
        corr->invFFT();
        if (!modelSettings->getUseLocalNoise()) // Already done in crava.cpp if local noise
          corr->createPostVariances();
        corr->printPostVariances();
        if((modelSettings->getOutputGridsOther() & IO::CORRELATION) > 0)
        {
          corr->writeFilePostVariances();
          corr->writeFilePostCovGrids(modelGeneral->getTimeSimbox());
        }

        int activeAngles = 0; //How many dimensions for local noise interpolation? Turn off for now.
        if(modelSettings->getUseLocalNoise()==true)
          activeAngles = modelSettings->getNumberOfAngles();
        spatwellfilter->doFiltering(corr,
                                    modelAVOstatic->getWells(),
                                    modelSettings->getNumberOfWells(),
                                    modelSettings->getNoVsFaciesProb(),
                                    modelSettings->getIndicatorFilter(),
                                    activeAngles,
                                    crava,
                                    modelAVOdynamic->getLocalNoiseScales());

        // FilterWellLogs * filteredlogs = NULL;
        //crava->filterLogs(modelGeneral->getTimeSimboxConstThick(),filteredlogs);

        if (modelSettings->getEstimateFaciesProb()) {
          bool useFilter = modelSettings->getUseFilterForFaciesProb();
          crava->computeFaciesProb(spatwellfilter, useFilter);
        }
        delete spatwellfilter;

        if(modelSettings->getKrigingParameter() > 0)
          crava->doPredictionKriging();

        if(modelSettings->getGenerateSeismicAfterInv())
          crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());

        //
        // Temporary placement.  crava.cpp needs a proper restructuring.
        //
        if((modelSettings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
          modelAVOstatic->writeBlockedWells(modelAVOstatic->getWells(),modelSettings);
        }
        if((modelSettings->getWellOutputFlag() & IO::BLOCKED_LOGS) > 0) {
          LogKit::LogFormatted(LogKit::Low,"\nWARNING: Writing of BLOCKED_LOGS is not implemented yet.\n");
        }

        delete crava;
      } //end doinversion
    }
    else // do forward modeling
    {
      LogKit::LogFormatted(LogKit::Low,"\nBuilding model ...\n");
      crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, 0);
      LogKit::LogFormatted(LogKit::Low,"\n               ... model built\n");

      crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());
      delete crava;
    }

    if (FFTGrid::getMaxAllowedGrids() != FFTGrid::getMaxAllocatedGrids() && modelSettings->getDoInversion()) {
      LogKit::LogFormatted(LogKit::DebugLow,"\nWARNING: A memory requirement inconsistency has been detected:");
      LogKit::LogFormatted(LogKit::DebugLow,"\n            Maximum number of grids requested  :  %2d",FFTGrid::getMaxAllowedGrids());
      LogKit::LogFormatted(LogKit::DebugLow,"\n            Maximum number of grids allocated  :  %2d",FFTGrid::getMaxAllocatedGrids());
      TaskList::addTask("The memory estimate in CRAVA failed. The developers would be interested to know about this, so if you inform support about this, and provide your xml-file, it would be appreciated.");
    }

    Timings::setTimeTotal(wall,cpu);
    Timings::reportAll(LogKit::Medium);

    TaskList::viewAllTasks(modelSettings->getTaskFileFlag());

    delete modelAVOdynamic;
    delete modelAVOstatic;
    delete modelGeneral;
    delete modelSettings;

    Timings::reportTotal();
    LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA closing  ***\n");
    LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA finished ***\n");
    LogKit::EndLog();

    }
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Out of memory: " << ba.what() << std::endl;
  }
  return(0);
}
