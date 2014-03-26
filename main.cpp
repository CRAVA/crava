/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#if defined(COMPILE_STORM_MODULES_FOR_RMS)
#include <util/precompile.h>

#include <license/feature.h>
#include <license/features.h>
#include <file/appl.h>
#include <config/releaseinfo.h>
#ifdef GCOMPILE_SIMPLELMSYSTEM
#include <license/simplelmsystem.h>
#else
#include <license/flexlm/flexlmsystem.h>
#endif
#endif

#include "lib/timekit.hpp"
#include "lib/utils.h"

#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "rplib/demmodelling.h"

#include "src/program.h"
#include "src/definitions.h"
#include "src/wavelet.h"
#include "src/avoinversion.h"
#include "src/fftgrid.h"
#include "src/gridmapping.h"
#include "src/simbox.h"
//#include "src/welldata.h"
#include "src/timings.h"
#include "src/spatialwellfilter.h"
#include "src/tasklist.h"
#include "src/commondata.h"

#include "src/xmlmodelfile.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/modelgeneral.h"
#include "src/timeline.h"
#include "src/modelgravitystatic.h"

#include "src/seismicparametersholder.h"
#include "src/doinversion.h"

#if defined(COMPILE_STORM_MODULES_FOR_RMS)

class CravaAppl : public IoAppl {
public:
  CravaAppl(const char  *appl_name,
            const char  *appl_full_name,
            Trademark    trademark,
            const char  *version,
            const char  *visible_version,
            VersionType  version_type,
            FileVersion  file_version,
            const char  *resource_name,
            float        license_version);

  //- Methods inherited from Appl
  virtual int  Main(int argc, char **argv);

private:
  float       license_version_;
};


CravaAppl::CravaAppl(const char  *appl_name,
                     const char  *appl_full_name,
                     Trademark    trademark,
                     const char  *version,
                     const char  *visible_version,
                     VersionType  version_type,
                     FileVersion  file_version,
                     const char  *resource_name,
                     float        license_version)
  : IoAppl(appl_name, appl_full_name, trademark, version, visible_version,
           version_type, file_version, resource_name) {

  license_version_ = license_version;

  Appl::Install(this);
}

int
CravaAppl::Main(int /*argc*/, char ** /*argv*/) {
#ifdef GCOMPILE_SIMPLELMSYSTEM
  // Allow any license feature to be checked out (true = allow_checkout):
  SimpleLmLicenseSystem *dummy = new SimpleLmLicenseSystem(true);
#else
  FlexLmLicenseSystem *dummy = new FlexLmLicenseSystem(license_version_);
#endif
  LicenseSystem *license_system = LicenseSystem::Instance();

  Feature& feature = FEATURE_INVERSION_EXE;
  if (feature.GetState() == Feature::NOTEXIST ||
      !license_system->CheckOut(&feature)) {
    fprintf(stderr,"Can't check out licence feature %s\n",feature.GetFullName().toLatin1().constData());
    exit(2);
  }

  return 1;
}

#endif

int main(int argc, char** argv)
{
  bool AddLicenceInformation = true;
#if defined(COMPILE_STORM_MODULES_FOR_RMS)
  AddLicenceInformation = false;
  //* The Appl is static, otherwise it is not deleted upon call to exit().
  //* TODO: checkUsage() must not call exit(), either return a status, or
  //* throw an exception that is caught here.
  //* So that the app can exit gracefully.
  static
  CravaAppl CravaAppl("RMS",
                      "RMS - Seismic inversion module",
                      Appl::NOTRADEMARK,
                      ReleaseInfo::GetRmsMajorMinorVersionNumber(),
                      ReleaseInfo::GetFullVersionNumber(),
                      Appl::RELEASE,
                      FileVersion(0),
                      "SEISMIC_INVERSION",
                      ReleaseInfo::GetLicenseVersion());
  Appl::Instance()->Main(argc, argv);
#endif

  if(0) {
    //test of DEM RPM
    double effective_bulk_modulus2;
    double effective_shear_modulus2;
    double effective_density2;
    DEMTools::DebugTestCalcEffectiveModulus2(effective_bulk_modulus2,
                                             effective_shear_modulus2,
                                             effective_density2);
    double effective_bulk_modulus;
    double effective_shear_modulus;
    double effective_density;
    DEMTools::DebugTestCalcEffectiveModulus4(effective_bulk_modulus,
                                             effective_shear_modulus,
                                             effective_density);
    //float tmp10 = 5.0f;


  }

  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::SetScreenLog(LogKit::L_Low);
  LogKit::StartBuffering();

  Program program( 2,                     // Major version
                   0,                     // Minor version
                   0,                     // Patch number for bug fixes
                   //"",                  // Use empty string "" for release versions
                   " beta",               // Use empty string "" for release versions
                   -1,                    // Validity of licence in days (-1 = infinite)
                   //"NORSAR",            // Who this copy of CRAVA is licensed to
                   "Norsk Regnesentral/Statoil",
                   AddLicenceInformation);

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  try
  {
    XmlModelFile modelFile(argv[1]);
    InputFiles     * inputFiles     = modelFile.getInputFiles();
    ModelSettings  * modelSettings  = modelFile.getModelSettings();
    CommonData     * common_data     = NULL;
    ModelGeneral   * modelGeneral   = NULL;
    ModelAVOStatic * modelAVOstatic = NULL;
    ModelGravityStatic * modelGravityStatic = NULL;

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

    /*------------------------------------------------------------
    READ COMMON DATA AND PERFORM ESTIMATION BASED ON INPUT FILES
    AND MODEL SETTINGS
    -------------------------------------------------------------*/

    common_data = new CommonData(modelSettings, inputFiles);

    //Loop over intervals
    for (int i_interval = 0; i_interval < common_data->GetMultipleIntervalGrid()->GetNIntervals(); i_interval++) {

      modelGeneral       = NULL;
      modelAVOstatic     = NULL;
      modelGravityStatic = NULL;

      std::string interval_text = "";
      if (common_data->GetMultipleIntervalGrid()->GetNIntervals() > 1)
        interval_text = " for interval" + NRLib::ToString(common_data->GetMultipleIntervalGrid()->GetIntervalName(i_interval));
      LogKit::WriteHeader("Setting up models" + interval_text);

      //Priormodell i 3D

      //From NRLib::Grid to FFT-grid, need to fill in padding.
      SeismicParametersHolder seismicParametersInterval;

      const Simbox * simbox = common_data->GetMultipleIntervalGrid()->GetIntervalSimbox(i_interval);

      //Forventningsgrid (fra multisonegrid)
      seismicParametersInterval.setBackgroundParametersInterval(common_data->GetMultipleIntervalGrid()->GetBackgroundParametersInterval(i_interval),
                                                                simbox->GetNXpad(),
                                                                simbox->GetNYpad(),
                                                                simbox->GetNZpad());
      //korrelasjonsgrid (2m)
      float corr_grad_I = 0.0f;
      float corr_grad_J = 0.0f;
      common_data->GetCorrGradIJ(corr_grad_I, corr_grad_J, simbox);

      float dt        = static_cast<float>(simbox->getdz());
      float low_cut   = modelSettings->getLowCut();
      int low_int_cut = int(floor(low_cut*(simbox->GetNZpad()*0.001*dt))); // computes the integer whis corresponds to the low cut frequency.

      seismicParametersInterval.setCorrelationParameters(common_data->GetPriorParamCov(i_interval),
                                                         common_data->GetPriorCorrT(i_interval),
                                                         common_data->GetPriorCorrXY(i_interval),
                                                         low_int_cut,
                                                         corr_grad_I,
                                                         corr_grad_J,
                                                         simbox->getnx(),
                                                         simbox->getny(),
                                                         simbox->getnz(),
                                                         simbox->GetNXpad(),
                                                         simbox->GetNYpad(),
                                                         simbox->GetNZpad());

      //ModelGeneral, modelAVOstatic, modelGravityStatic, (modelTravelTimeStatic?)
      setupStaticModels(modelGeneral,
                        modelAVOstatic,
                        modelGravityStatic,
                        modelSettings,
                        inputFiles,
                        seismicParametersInterval,
                        common_data,
                        i_interval);

      //Loop over dataset
      //i.   ModelAVODynamic
      //ii.  Inversion
      //iii. Move model one time-step ahead

      //Skal ikke v�re en egen doFirstAVOInversion (det som er spesielt i denne skal v�re gjort i CommonData eller flyttes til doAVOInversion)
      //M� skille ut forward-modelling her et sted.

      if(modelGeneral->GetTimeLine() == NULL) {//Forward modelling.


      }
      else {
        int  eventType;
        int  eventIndex;
        double  oldTime;
        bool failedFirst = false;
        modelGeneral->GetTimeLine()->ReSet();
        modelGeneral->GetTimeLine()->GetNextEvent(eventType, eventIndex, oldTime);
        switch(eventType) {

        case TimeLine::AVO :
          failedFirst = doTimeLapseAVOInversion(modelSettings, modelGeneral, modelAVOstatic, common_data, seismicParametersInterval, eventIndex, i_interval);
          break;
        case TimeLine::TRAVEL_TIME :
          failedFirst = doTimeLapseTravelTimeInversion(modelSettings, modelGeneral, inputFiles, eventIndex, seismicParametersInterval);
          break;
        case TimeLine::GRAVITY :
          failedFirst = doTimeLapseGravimetricInversion(modelSettings, modelGeneral, modelGravityStatic, common_data, inputFiles, eventIndex, seismicParametersInterval);
          break;
        default :
          errTxt += "Error: Unknown inverstion type.\n";
          break;
        }

        //H Why must/is the first run be seperate from the loop below?


        if(failedFirst == true || errTxt != "")
          return(1);

        //delete timeBGSimbox;

        double time;
        int time_index = 0;
        while(modelGeneral->GetTimeLine()->GetNextEvent(eventType, eventIndex, time) == true) {
          modelGeneral->AdvanceTime(time_index, seismicParametersInterval, modelSettings);
          time_index++;
          bool failed;
          switch(eventType) {
          case TimeLine::AVO :
            failed = doTimeLapseAVOInversion(modelSettings,
                                             modelGeneral,
                                             modelAVOstatic,
                                             common_data,
                                             seismicParametersInterval,
                                             eventIndex,
                                             i_interval);
            break;
          case TimeLine::TRAVEL_TIME :
            failed = doTimeLapseTravelTimeInversion(modelSettings,
                                                    modelGeneral,
                                                    inputFiles,
                                                    eventIndex,
                                                    seismicParametersInterval);
            break;
          case TimeLine::GRAVITY :
            failed = doTimeLapseGravimetricInversion(modelSettings,
                                                     modelGeneral,
                                                     modelGravityStatic,
                                                     common_data,
                                                     inputFiles,
                                                     eventIndex,
                                                     seismicParametersInterval);
            break;
          default :
            failed = true;
            break;
          }
          if(failed)
            return(1);
        }
      }
    } //interval_loop








    //Simbox * timeBGSimbox = NULL;
    //SeismicParametersHolder seismicParameters;

    //setupStaticModels(modelGeneral,
    //                  modelAVOstatic,
    //                  modelGravityStatic,
    //                  modelSettings,
    //                  inputFiles,
    //                  seismicParameters,
    //                  timeBGSimbox);

    //if(modelGeneral   == NULL || modelGeneral->getFailed()   ||
    //   modelAVOstatic == NULL || modelAVOstatic->getFailed())
    //  return(1);

    //if(modelGeneral->getTimeLine() == NULL) {//Forward modelling.
      //bool failed = doFirstAVOInversion(modelSettings,
                                        //modelGeneral,
                                        //modelAVOstatic,
                                        //seismicParameters,
                                        //inputFiles,
                                        //0,
                                        //timeBGSimbox);
      //if(failed)
      //  return(1);
    //}
    //else {
    //  int  eventType;
    //  int  eventIndex;
    //  double  oldTime;
    //  bool failedFirst = false;
    //  modelGeneral->getTimeLine()->ReSet();
    //  modelGeneral->getTimeLine()->GetNextEvent(eventType, eventIndex, oldTime);
    //  switch(eventType) {

    //  case TimeLine::AVO :
    //    // In case of 4D inversion
    //    if(modelSettings->getDo4DInversion()){
    //      failedFirst           = doTimeLapseAVOInversion(modelSettings, modelGeneral, modelAVOstatic, inputFiles, seismicParameters, eventIndex);

    //    }
        // In case of 3D inversion
        //else{
          //failedFirst = doFirstAVOInversion(modelSettings,
          //                                  modelGeneral,
          //                                  modelAVOstatic,
          //                                  seismicParameters,
          //                                  inputFiles,
          //                                  eventIndex,
          //                                  timeBGSimbox);
      //  }
      //  break;

      //case TimeLine::TRAVEL_TIME :
      //  failedFirst = doTimeLapseTravelTimeInversion(modelSettings,
      //                                               modelGeneral,
      //                                               inputFiles,
      //                                               eventIndex,
      //                                               seismicParameters);
      //  break;

      //case TimeLine::GRAVITY :
      // errTxt += "Error: Asked for 3D gravimetric inversion: Not available.\n";
      //  break;
      //default :
      //  errTxt += "Error: Unknown inverstion type.\n";
      //  break;
      //}

      //if(failedFirst == true || errTxt != "")
      //  return(1);

      //delete timeBGSimbox;

      //double time;
      //int time_index = 0;
      //while(modelGeneral->getTimeLine()->GetNextEvent(eventType, eventIndex, time) == true) {
      //  modelGeneral->advanceTime(time_index, seismicParameters,modelSettings);
      //  time_index++;
      //  bool failed;
      //  switch(eventType) {
      //  case TimeLine::AVO :
      //    failed = doTimeLapseAVOInversion(modelSettings, modelGeneral, modelAVOstatic, inputFiles, seismicParameters, eventIndex);
      //    break;
      //  case TimeLine::TRAVEL_TIME :
      //    failed = doTimeLapseTravelTimeInversion(modelSettings,
      //                                            modelGeneral,
      //                                            inputFiles,
      //                                            eventIndex,
      //                                            seismicParameters);
        //  break;
        //case TimeLine::GRAVITY :
        //  failed = doTimeLapseGravimetricInversion(modelSettings,
        //                                           modelGeneral,
        //                                           modelGravityStatic,
        //                                           inputFiles,
        //                                           eventIndex,
        //                                           seismicParameters);
        //  break;
        //default :
        //  failed = true;
        //  break;
        //}
        //if(failed)
        //  return(1);
    //  }
    //}

    if(modelSettings->getDo4DInversion())
    {
      bool failed;
      if(modelSettings->getDo4DRockPhysicsInversion())
      {

        failed = modelGeneral->Do4DRockPhysicsInversion(modelSettings);

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
    delete common_data;
    delete modelSettings;
    delete inputFiles;

    Timings::reportTotal();
    LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA closing  ***\n");
    LogKit::LogFormatted(LogKit::Low,"\n*** CRAVA finished ***\n");
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Out of memory: " << ba.what() << std::endl;
    std::string error_message = std::string("Out of memory: ") + ba.what() + "\n";
    LogKit::LogMessage(LogKit::Error, error_message);
  }

#if defined(COMPILE_STORM_MODULES_FOR_RMS)
  Feature& feature = FEATURE_INVERSION_EXE;
  LicenseSystem::Instance()->CheckIn(&feature);
#endif

  LogKit::EndLog(); //Debug messages may occur when variables go out of scope above, so this must be the final call.
  return(0);
}
