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
#include "src/timeline.h"

#include "src/seismicparametersholder.h"
#include "src/doinversion.h"

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
    InputFiles * inputFiles;
    XmlModelFile modelFile(argv[1]);
    inputFiles     = modelFile.getInputFiles();
    ModelSettings * modelSettings = modelFile.getModelSettings();

    bool failedModelFile    = false;
    bool failedInputFiles   = false;
    bool failedLoadingModel = false;

    ModelGeneral    * modelGeneral    = NULL;
    ModelAVOStatic  * modelAVOstatic  = NULL;

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
    
    Simbox * timeBGSimbox   = NULL;
    SeismicParametersHolder seismicParameters;
 
    if(!failedModelFile && !failedInputFiles){

      setupStaticModels(modelGeneral, modelAVOstatic, modelSettings, inputFiles, timeBGSimbox);
 
      failedLoadingModel =  modelGeneral    == NULL || modelGeneral->getFailed()   || 
                            modelAVOstatic  == NULL || modelAVOstatic->getFailed();

      if(failedModelFile || failedInputFiles || failedLoadingModel)
        return(1);

      bool failedFirst = false;

      int eventType;
      int eventIndex;
      int oldTime;
      if(modelGeneral->getTimeLine() == NULL) {//Forward modelling.
        bool failed = doFirstAVOInversion(modelSettings, modelGeneral, modelAVOstatic, seismicParameters, inputFiles, 0, timeBGSimbox);
        if(failed == true || errTxt != "")
          return(1);
      }
      else {
        modelGeneral->getTimeLine()->GetNextEvent(eventType, eventIndex, oldTime);
        switch(eventType) {
          case TimeLine::AVO :
            failedFirst = doFirstAVOInversion(modelSettings, modelGeneral, modelAVOstatic, seismicParameters, inputFiles, eventIndex, timeBGSimbox);
            break;
          case TimeLine::TRAVEL_TIME :
          case TimeLine::GRAVITY :
            errTxt += "Error: Asked for inversion type that is not implemented yet.\n";
            break;
          default :
            errTxt += "Error: Unknown inverstion type.\n";
            break;
        }
   
        if(failedFirst == true || errTxt != "")
          return(1);

        if(timeBGSimbox != NULL)
          delete timeBGSimbox;

        int time;
        while(modelGeneral->getTimeLine()->GetNextEvent(eventType, eventIndex, time) == true) {
          //Advance time (time-oldTime);
          bool failed;
          switch(eventType) {
            case TimeLine::AVO :
              failed = doTimeLapseAVOInversion(modelSettings, modelGeneral, modelAVOstatic, inputFiles, seismicParameters, eventIndex);
              break;
            case TimeLine::TRAVEL_TIME :
            case TimeLine::GRAVITY :
              failed = true;
              break;
            default :
              failed = true;
              break;
          }
          if(failed == true)
            return(1);
        }
      }
    }
    else {
      LogKit::SetFileLog(IO::FileLog()+IO::SuffixTextFiles(), modelSettings->getLogLevel());
      LogKit::EndBuffering();
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
    
    delete modelAVOstatic;
    delete modelGeneral;
    delete modelSettings;
    delete inputFiles;  

    Timings::reportTotal();
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Out of memory: " << ba.what() << std::endl;
  }
  LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA closing  ***\n"); 
  LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA finished ***\n");
  LogKit::EndLog();

  return(0);
}
