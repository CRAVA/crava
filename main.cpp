#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#include "lib/global_def.h"
#include "lib/timekit.hpp"
#include "lib/utils.h"

#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/program.h"
#include "src/definitions.h"
#include "src/corr.h"
#include "src/model.h"
#include "src/wavelet.h"
#include "src/crava.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/welldata.h"
#include "src/filterwelllogs.h"
#include "src/timings.h"
#include "src/spatialwellfilter.h"
#include "src/tasklist.h"

int main(int argc, char** argv)
{  
  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::SetScreenLog(LogKit::L_LOW);
  LogKit::StartBuffering();

  Program program( 0,                     // Major version
                   9,                     // Minor version 
                   7,                     // Patch number 
                   "beta",                // Use empty string "" for release versions 
                   -1,                    // Validity of licence in days (-1 = infinite)
                  "Norsk Regnesentral");  // Who this copy of CRAVA is licensed to

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  // Parsing modelfile and reading files
  Model * model = new Model(argv[1]);

  if(model->getFailed())
    return(1);
  
  ModelSettings * modelSettings = model->getModelSettings();
  Crava         * crava         = NULL;

  if(!modelSettings->getForwardModeling())
  {
    if (modelSettings->getDoInversion())
    {
      time_t timestart, timeend;
      time(&timestart);

      int nwells = modelSettings->getNumberOfWells();
      SpatialWellFilter *spatwellfilter = new SpatialWellFilter(nwells);
      crava = new Crava(model, spatwellfilter);
      
      std::string warningText("");
      
      if(crava->getWarning( warningText ) != 0)
       {
         LogKit::LogFormatted(LogKit::LOW,"\nWarning  !!!\n");
         LogKit::LogFormatted(LogKit::LOW,"%s",warningText.c_str());
         LogKit::LogFormatted(LogKit::LOW,"\n");
       }
      crava->printEnergyToScreen();
      
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  
      crava->computePostMeanResidAndFFTCov();
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  

      if(modelSettings->getNumberOfSimulations() > 0)
      {
        crava->simulate(model->getRandomGen());
      }

      Corr * corr = model->getCorrelations();
      corr->invFFT();
      if (!modelSettings->getUseLocalNoise()) // Already done in crava.cpp if local noise
        corr->createPostVariances();
      corr->printPostVariances();
      if((modelSettings->getGridOutputFlag() & IO::CORRELATION) > 0)
      {
        corr->writeFilePostVariances();
        corr->writeFilePostCovGrids(model->getTimeSimbox());
      }       

      int activeAngles = 0; //How many dimensions for local noise interpolation? Turn off for now.
      if(modelSettings->getUseLocalNoise()==true)
        activeAngles = modelSettings->getNumberOfAngles();
      spatwellfilter->doFiltering(corr,model->getWells(), 
                                  modelSettings->getNumberOfWells(), 
                                  modelSettings->getNoVsFaciesProb(),
                                  activeAngles, crava, model->getLocalNoiseScales());
      
      // FilterWellLogs * filteredlogs = NULL;
      //crava->filterLogs(model->getTimeSimboxConstThick(),filteredlogs);

      if (modelSettings->getEstimateFaciesProb())
        crava->computeFaciesProb(spatwellfilter);
      delete spatwellfilter;

      if(modelSettings->getKrigingParameter() > 0)
        crava->doPredictionKriging();

      if(modelSettings->getGenerateSeismicAfterInv())
        crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());

      //
      // Temporary placement.  crava.cpp needs a proper restructuring.
      //
      if((modelSettings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
        model->writeBlockedWells(model->getWells(),modelSettings);
      }
      if((modelSettings->getWellOutputFlag() & IO::BLOCKED_LOGS) > 0) {
        LogKit::LogFormatted(LogKit::LOW,"\nWARNING: Writing of BLOCKED_LOGS is not implemented yet.\n");
      }

      delete crava;
    } //end doinversion 
  }
  else // do forward modeling
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBuilding model ...\n");
    crava = new Crava(model, 0);
    LogKit::LogFormatted(LogKit::LOW,"\n               ... model built\n");
  
    crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());
    delete crava;
  } 
  delete model;

  if (FFTGrid::getMaxAllowedGrids() > FFTGrid::getMaxAllocatedGrids()) {
    LogKit::LogFormatted(LogKit::DEBUGLOW,"\nWARNING: A memory requirement inconsistency has been detected:"); 
    LogKit::LogFormatted(LogKit::DEBUGLOW,"\n            Maximum number of grids requested  :  %2d",FFTGrid::getMaxAllowedGrids()); 
    LogKit::LogFormatted(LogKit::DEBUGLOW,"\n            Maximum number of grids allocated  :  %2d",FFTGrid::getMaxAllocatedGrids()); 
    LogKit::LogFormatted(LogKit::DEBUGLOW,"\n         Consult method Model::checkAvailableMemory().\n"); 
    TaskList::addTask("CRAVA did not use as much memory as estimated. NR would be interested to konow about this, so if you could send your .xml-file to us, we would aprreciate it.");
  }
  
  Timings::setTimeTotal(wall,cpu);
  Timings::reportAll(LogKit::MEDIUM);

  TaskList::viewAllTasks();  

  Timings::reportTotal(LogKit::LOW);

  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA closing  ***\n"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA finished ***\n");
  LogKit::EndLog();
  return(0);
}
