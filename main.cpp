#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#include "lib/systemcall.h"
#include "lib/global_def.h"
#include "lib/timekit.hpp"
#include "lib/utils.h"

#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"

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

int main(int argc, char** argv)
{  
  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::SetScreenLog(LogKit::L_LOW);
  LogKit::StartBuffering();

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
  LogKit::LogFormatted(LogKit::LOW,"\n***************************************************************************************************");
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                                                                         *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                    C  R  A  V  A                                        *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                                                                         *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n***************************************************************************************************\n\n");

  std::cout 
    << "Compiled: " << SystemCall::getDate() << "/" << SystemCall::getTime() << "\n"
    << std::endl;

  const char * userName    = SystemCall::getUserName();
  const char * dateAndTime = SystemCall::getCurrentTime();
  const char * hostName    = SystemCall::getHostName();
  LogKit::LogFormatted(LogKit::LOW,"Log written by                             : %s\n",userName);
  LogKit::LogFormatted(LogKit::LOW,"Date and time                              : %s"  ,dateAndTime);
  LogKit::LogFormatted(LogKit::LOW,"Host                                       : %s\n",hostName);
  delete [] userName;
  delete [] dateAndTime;
  delete [] hostName;

  // Parsing modelfile and reading files
  Model * model = new Model(argv[1]);
  if(model->getFailed())
  {
    LogKit::LogFormatted(LogKit::LOW,"\nErrors detected when loading data.\nAborting.\n");
    return(1);
  }

  Crava * crava;

  if(model->getModelSettings()->getGenerateSeismic() == false)
  {
    if (model->getModelSettings()->getDoInversion())
    {
      time_t timestart, timeend;
      time(&timestart);

      crava = new Crava(model);
      
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

      Corr * corr = model->getCorrelations();
      corr->invFFT();
      corr->getPostVariances();
      corr->printPostVariances();
      if((model->getModelSettings()->getOutputFlag() & ModelSettings::CORRELATION) > 0)
      {
        corr->writeFilePostVariances();
        corr->writeFilePostCovGrids(model->getTimeSimbox());
      }       

      FilterWellLogs * filteredlogs = NULL;
      crava->filterLogs(model->getTimeSimboxConstThick(),filteredlogs);
      crava->computeFaciesProb(filteredlogs);
      delete filteredlogs;

      //
      // Temprary placement.  crava.cpp needs a proper restructuring.
      //
      if((model->getModelSettings()->getOutputFlag() & ModelSettings::BLOCKED_WELLS) > 0) {
        WellData ** wells = model->getWells();
        for (int i=0 ; i<model->getModelSettings()->getNumberOfWells() ; i++)
          wells[i]->getBlockedLogsOrigThick()->writeRMSWell(model->getModelSettings());
      }
      if((model->getModelSettings()->getOutputFlag() & ModelSettings::BLOCKED_LOGS) > 0) {
        LogKit::LogFormatted(LogKit::LOW,"\nWARNING: Writing of BLOCKED_LOGS is not implemented yet.\n");
      }

      delete [] warningText;
      delete crava;
    } //end doinversion 
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"\nBuilding model ...\n");
    crava = new Crava(model);
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
