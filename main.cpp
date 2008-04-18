#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#include "lib/systemcall.h"
#include "lib/global_def.h"
#include "lib/log.h"
#include "lib/segy.h"

#include "src/model.h"
#include "src/wavelet.h"
#include "src/crava.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

int main(int argc, char** argv)
{  
  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::initialize(true);

  double wall=0.0, cpu=0.0;
  LogKit::getTime(wall,cpu);
  LogKit::writeLog("\n***********************************************************************");
  LogKit::writeLog("\n***                                                                 ***"); 
  LogKit::writeLog("\n***                        C  R  A  V  A                            ***"); 
  LogKit::writeLog("\n***                                                                 ***"); 
  LogKit::writeLog("\n***********************************************************************\n\n");

  char segyMode[50];
  char bypassCoordScaling[50];

#ifdef SEGY_ISEX
  strcpy(segyMode,"ISEX");
#else
  strcpy(segyMode,"Seisworks/Charisma");
#endif
#ifdef BYPASS_COORDINATE_SCALING
  strcpy(bypassCoordScaling,"yes");
#else
  strcpy(bypassCoordScaling,"no");
#endif

  std::cout 
    << "Compiled: " << SystemCall::getDate() << "/" << SystemCall::getTime() << "\n"
    << std::endl;

  const char * userName    = SystemCall::getUserName();
  const char * dateAndTime = SystemCall::getCurrentTime();
  const char * hostName    = SystemCall::getHostName();
  LogKit::writeLog("Compile-time directives used in this version:\n");
  LogKit::writeLog("  SegY mode: %s\n",segyMode);
  LogKit::writeLog("  Bypass coordinate scaling: %s\n\n",bypassCoordScaling);
  LogKit::writeLog("Log written by                             : %s\n",userName);
  LogKit::writeLog("Date and time                              : %s"  ,dateAndTime);
  LogKit::writeLog("Host                                       : %s\n",hostName);
  delete [] userName;
  delete [] dateAndTime;
  delete [] hostName;

  // Parsing modelfile and reading files
  Model * model = new Model(argv[1]);
  if(model->getFailed())
  {
    LogKit::writeLog("\nErrors detected in model file processing.\nAborting.\n");
    return(1);
  }

  Crava * crava;

  if(model->getModelSettings()->getGenerateSeismic() == false)
  {
    if (model->getModelSettings()->getDoInversion())
    {
      time_t timestart, timeend;
      time(&timestart);
      
      LogKit::writeLog("\n***********************************************************************");
      LogKit::writeLog("\n***                    Building Stochastic Model                     ***"); 
      LogKit::writeLog("\n***********************************************************************\n\n");

      crava = new Crava(model);
      
      char * warningText = new char[12*MAX_STRING*crava->getNTheta()];
      
      if(crava->getWarning( warningText ) != 0)
       {
         LogKit::writeLog("\nWarning  !!!\n");
         LogKit::writeLog("%s",warningText);
         LogKit::writeLog("\n");
       }
      crava->printEnergyToScreen();
      
      time(&timeend);
      LogKit::writeDebugLog("\nTime elapsed :  %d\n",timeend-timestart);  
      LogKit::writeLog("\n***********************************************************************");
      LogKit::writeLog("\n***             Posterior model / Performing Inversion              ***"); 
      LogKit::writeLog("\n***********************************************************************\n\n");
      crava->computePostMeanResidAndFFTCov();
      time(&timeend);
      LogKit::writeDebugLog("\nTime elapsed :  %d\n",timeend-timestart);  
      
      if(model->getModelSettings()->getNumberOfSimulations() > 0)
      {
        LogKit::writeLog("\n***********************************************************************");
        LogKit::writeLog("\n***                Simulating from posterior model                  ***"); 
        LogKit::writeLog("\n***********************************************************************\n\n");
        crava->simulate(model->getRandomGen());
      }
      
      // Posterior covariance
      if((model->getModelSettings()->getOutputFlag() & ModelSettings::CORRELATION) > 0)
      {
        LogKit::writeLog("\nPost process ...\n"); 
        crava->computePostCov();
        LogKit::writeLog("\n             ... post prosess ended\n");
        
      }
      crava->computeFaciesProb();
      delete [] warningText;
      delete crava;
    } //end doinversion 
  }
  else
  {
    LogKit::writeLog("\nBuilding model ...\n");
    crava = new Crava(model);
    LogKit::writeLog("\n               ... model built\n");

    // Computing synthetic seismic
    LogKit::writeLog("\nComputing synthetic seismic ..."); 
    crava->computeSyntSeismic(crava->getpostAlpha(),crava->getpostBeta(),crava->getpostRho());
    LogKit::writeLog("                              ... synthetic seismic computed.\n");
	
    delete crava;
  } 
  delete model;

  LogKit::writeLog("\n*** CRAVA closing  ***\n"); 
  LogKit::getTime(wall,cpu);
  LogKit::writeLog("\nTotal CPU  time used in CRAVA: %6d seconds", static_cast<int>(cpu));
  LogKit::writeLog("\nTotal Wall time used in CRAVA: %6d seconds\n", static_cast<int>(wall));
  LogKit::writeLog("\n*** CRAVA finished ***\n");

  LogKit::terminate();
  return(0);
}
