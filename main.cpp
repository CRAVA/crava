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
    XmlModelFile modelFile(argv[1]);

    InputFiles     * inputFiles     = modelFile.getInputFiles();
    ModelSettings  * modelSettings  = modelFile.getModelSettings();
    ModelGeneral   * modelGeneral   = NULL;
    ModelAVOStatic * modelAVOstatic = NULL;

    if (modelFile.getParsingFailed()) {
      LogKit::SetFileLog(IO::FileLog()+IO::SuffixTextFiles(), modelSettings->getLogLevel());
      LogKit::EndBuffering();
      return(1);
    }
   
    std::string errTxt = inputFiles->addInputPathAndCheckFiles();
    if(errTxt != "") {
      LogKit::WriteHeader("Error opening files");
      LogKit::LogMessage(LogKit::Error, "\n"+errTxt);
      LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
      LogKit::SetFileLog(IO::FileLog()+IO::SuffixTextFiles(), modelSettings->getLogLevel());
      LogKit::EndBuffering();
      return(1);
    }

    Simbox * timeBGSimbox   = NULL;
    SeismicParametersHolder seismicParameters;
 
    setupStaticModels(modelGeneral, 
                      modelAVOstatic, 
                      modelSettings, 
                      inputFiles, 
                      timeBGSimbox);
    
    if(modelGeneral   == NULL || modelGeneral->getFailed()   || 
       modelAVOstatic == NULL || modelAVOstatic->getFailed())
      return(1);
    
    if(modelGeneral->getTimeLine() == NULL) {//Forward modelling.
      bool failed = doFirstAVOInversion(modelSettings, 
                                        modelGeneral, 
                                        modelAVOstatic, 
                                        seismicParameters, 
                                        inputFiles, 
                                        0, 
                                        timeBGSimbox);
      if(failed)
        return(1);
    }
    else {
      int  eventType;
      int  eventIndex;
      int  oldTime;      
      bool failedFirst = false;
      modelGeneral->getTimeLine()->GetNextEvent(eventType, eventIndex, oldTime);
      switch(eventType) {
      case TimeLine::AVO :
        failedFirst = doFirstAVOInversion(modelSettings, 
                                          modelGeneral, 
                                          modelAVOstatic, 
                                          seismicParameters, 
                                          inputFiles, 
                                          eventIndex, 
                                          timeBGSimbox);
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
        if(failed)
          return(1);
      }
    }
    
    if (modelSettings->getDoInversion() && FFTGrid::getMaxAllowedGrids() != FFTGrid::getMaxAllocatedGrids()) {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: A memory requirement inconsistency has been detected:");
      LogKit::LogFormatted(LogKit::Warning,"\n            Maximum number of grids requested  :  %2d",FFTGrid::getMaxAllowedGrids());
      LogKit::LogFormatted(LogKit::Warning,"\n            Maximum number of grids allocated  :  %2d",FFTGrid::getMaxAllocatedGrids());
      TaskList::addTask("The memory usage estimate failed. Please send your XML-model file and\nthe logFile.txt to the CRAVA developers.");
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
