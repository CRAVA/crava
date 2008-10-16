#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#include "lib/systemcall.h"
#include "lib/global_def.h"
#include "nrlib/segy/segy.hpp"
#include "lib/timekit.hpp"

#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/model.h"
#include "src/wavelet.h"
#include "src/crava.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/welldata.h"
#include "src/filterwelllogs.h"

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
  LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
  LogKit::LogFormatted(LogKit::LOW,"\n***                                                                 ***"); 
  LogKit::LogFormatted(LogKit::LOW,"\n***                        C  R  A  V  A                            ***"); 
  LogKit::LogFormatted(LogKit::LOW,"\n***                                                                 ***"); 
  LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");

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
    LogKit::LogFormatted(LogKit::LOW,"\nErrors detected in model file processing.\nAborting.\n");
    return(1);
  }

  Crava * crava;

  if(model->getModelSettings()->getGenerateSeismic() == false)
  {
    if (model->getModelSettings()->getDoInversion())
    {
      time_t timestart, timeend;
      time(&timestart);
      
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
      LogKit::LogFormatted(LogKit::LOW,"\n***                    Building Stochastic Model                     ***"); 
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");

      crava = new Crava(model);
      
      char * warningText = new char[12*MAX_STRING*crava->getNTheta()];
      
      if(crava->getWarning( warningText ) != 0)
       {
         LogKit::LogFormatted(LogKit::LOW,"\nWarning  !!!\n");
         LogKit::LogFormatted(LogKit::LOW,"%s",warningText);
         LogKit::LogFormatted(LogKit::LOW,"\n");
       }
      crava->printEnergyToScreen();
      
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
      LogKit::LogFormatted(LogKit::LOW,"\n***             Posterior model / Performing Inversion              ***"); 
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");
      crava->computePostMeanResidAndFFTCov();
      time(&timeend);
      LogKit::LogFormatted(LogKit::DEBUGLOW,"\nTime elapsed :  %d\n",timeend-timestart);  
      
      if(model->getModelSettings()->getNumberOfSimulations() > 0)
      {
        LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
        LogKit::LogFormatted(LogKit::LOW,"\n***                Simulating from posterior model                  ***"); 
        LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n\n");
        crava->simulate(model->getRandomGen());
      }
      FilterWellLogs * filteredlogs = NULL;
      crava->filterLogs(model->getTimeSimboxConstThick(),
                        filteredlogs);

      // Posterior covariance
      if((model->getModelSettings()->getOutputFlag() & ModelSettings::CORRELATION) > 0)
      {
        LogKit::LogFormatted(LogKit::LOW,"\nPost process ...\n"); 
        crava->computePostCov();
        LogKit::LogFormatted(LogKit::LOW,"\n             ... post prosess ended\n");
        
      }
      
      crava->computeFaciesProb(filteredlogs);

      //
      // Temprary placement.  crava.cpp needs a proper restructuring.
      //
      WellData ** wells = model->getWells();
      for (int i=0 ; i<model->getModelSettings()->getNumberOfWells() ; i++)
        wells[i]->getBlockedLogsPropThick()->writeRMSWell(model->getModelSettings());


      delete [] warningText;
      delete filteredlogs;
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
    crava->computeSyntSeismic(crava->getpostAlpha(),crava->getpostBeta(),crava->getpostRho());
    LogKit::LogFormatted(LogKit::LOW,"                              ... synthetic seismic computed.\n");
	
    delete crava;
  } 
  delete model;

  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA closing  ***\n"); 
  TimeKit::getTime(wall,cpu);
  LogKit::LogFormatted(LogKit::LOW,"\nTotal CPU  time used in CRAVA: %6d seconds", static_cast<int>(cpu));
  LogKit::LogFormatted(LogKit::LOW,"\nTotal Wall time used in CRAVA: %6d seconds\n", static_cast<int>(wall));
  LogKit::LogFormatted(LogKit::LOW,"\n*** CRAVA finished ***\n");

  LogKit::EndLog();
  return(0);
}
