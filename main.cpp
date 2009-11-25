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
                   5,                     // Patch number 
                   -1,                    // Validity of licence in days (-1 = infinite)
                  "Norsk Regnesentral");  // Who this copy of CRAVA is licensed to

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  // Parsing modelfile and reading files
  Model * model = new Model(argv[1]);

  if(model->getFailed())
    return(1);
  
  Crava * crava = NULL;

  if(model->getModelSettings()->getForwardModeling() == false)
  {
    if (model->getModelSettings()->getDoInversion())
    {
      time_t timestart, timeend;
      time(&timestart);

      int nwells = model->getModelSettings()->getNumberOfWells();
      SpatialWellFilter *spatwellfilter = new SpatialWellFilter(nwells);
      crava = new Crava(model, spatwellfilter);
      
      char * warningText = new char[12*MAX_STRING*model->getModelSettings()->getNumberOfAngles()];
      
      if(crava->getWarning( warningText ) != 0)
       {
         LogKit::LogFormatted(LogKit::LOW,"\nWarning  !!!\n");
         LogKit::LogFormatted(LogKit::LOW,"%s",warningText);
         LogKit::LogFormatted(LogKit::LOW,"\n");
       }
      crava->printEnergyToScreen();
      
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  
      crava->computePostMeanResidAndFFTCov();
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  

      if(model->getModelSettings()->getNumberOfSimulations() > 0)
      {
        crava->simulate(model->getRandomGen());
      }
      if(model->getModelSettings()->getGenerateSeismicAfterInversion() == true)
      {
        crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());
      }

      Corr * corr = model->getCorrelations();
      corr->invFFT();
      corr->getPostVariances();
      corr->printPostVariances();
      if((model->getModelSettings()->getGridOutputFlag() & IO::CORRELATION) > 0)
      {
        corr->writeFilePostVariances();
        corr->writeFilePostCovGrids(model->getTimeSimbox());
      }       

      spatwellfilter->doFiltering(corr,model->getWells(), 
                                  model->getModelSettings()->getNumberOfWells(), 
                                  model->getModelSettings()->getNoVsFaciesProb());
      
      // FilterWellLogs * filteredlogs = NULL;
      //crava->filterLogs(model->getTimeSimboxConstThick(),filteredlogs);

      if (model->getModelSettings()->getEstimateFaciesProb())
        crava->computeFaciesProb(spatwellfilter);
      delete spatwellfilter;
      //
      // Temporary placement.  crava.cpp needs a proper restructuring.
      //
      if((model->getModelSettings()->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
        WellData ** wells = model->getWells();
        for (int i=0 ; i<model->getModelSettings()->getNumberOfWells() ; i++)
          wells[i]->getBlockedLogsOrigThick()->writeWell(model->getModelSettings());
      }
      if((model->getModelSettings()->getWellOutputFlag() & IO::BLOCKED_LOGS) > 0) {
        LogKit::LogFormatted(LogKit::LOW,"\nWARNING: Writing of BLOCKED_LOGS is not implemented yet.\n");
      }

      delete [] warningText;
      delete crava;
    } //end doinversion 
  }
  else // do forward modeling
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBuilding model ...\n");
    crava = new Crava(model, 0);
    LogKit::LogFormatted(LogKit::LOW,"\n               ... model built\n");
  
    // Computing synthetic seismic
    LogKit::LogFormatted(LogKit::LOW,"\nComputing synthetic seismic ..."); 
    crava->computeSyntSeismic(crava->getPostAlpha(),crava->getPostBeta(),crava->getPostRho());
    LogKit::LogFormatted(LogKit::LOW,"                              ... synthetic seismic computed.\n");
	
    delete crava;
  } 
  delete model;

  Timings::setTimeTotal(wall,cpu);
  Timings::reportAll(LogKit::MEDIUM);

  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA closing  ***\n"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA finished ***\n");
  LogKit::EndLog();
  return(0);
}
