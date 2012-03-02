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

#include "src/program.h"
#include "src/definitions.h"
#include "src/corr.h"
#include "src/wavelet.h"
#include "src/crava.h"
#include "src/fftgrid.h"
#include "src/gridmapping.h"
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
#if defined(COMPILE_STORM_MODULES_FOR_RMS)
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

  if (argc != 2) {
    printf("Usage: %s modelfile\n",argv[0]);
    exit(1);
  }
  LogKit::SetScreenLog(LogKit::L_Low);
  LogKit::StartBuffering();

  Program program( 1,                     // Major version
                   3,                     // Minor version
                   0,                     // Patch number
                   //"",                  // Use empty string "" for release versions
                   " beta",               // Use empty string "" for release versions
                   -1,                    // Validity of licence in days (-1 = infinite)
                   //"Roxar",
                   //"NORSAR");
                   "Norsk Regnesentral/Statoil");  // Who this copy of CRAVA is licensed to

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  try
  {
    XmlModelFile modelFile(argv[1]);

    InputFiles      * inputFiles      = modelFile.getInputFiles();
    ModelSettings   * modelSettings   = modelFile.getModelSettings();
    ModelGeneral    * modelGeneral    = NULL;
    ModelAVOStatic  * modelAVOstatic  = NULL;
    ModelAVODynamic * modelAVOdynamic = NULL;

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

    // Construct ModelGeneral object first.
    // For each data type, construct the static model class before the dynamic.
    Simbox * timeBGSimbox = NULL;
    modelGeneral    = new ModelGeneral(modelSettings, inputFiles, timeBGSimbox);
    modelAVOstatic  = new ModelAVOStatic(modelSettings,
                                         inputFiles,
                                         modelGeneral->getFailedDetails(),
                                         modelGeneral->getTimeCutMapping(),
                                         modelGeneral->getTimeSimbox(),
                                         timeBGSimbox,
                                         modelGeneral->getTimeSimboxConstThick());
    // Wells are adjusted by ModelAVODynamic constructor.
    modelAVOdynamic = new ModelAVODynamic(modelSettings, inputFiles,
                                          modelGeneral->getFailedDetails(),
                                          modelAVOstatic->getFailedDetails(),
                                          modelGeneral->getTimeSimbox(),
                                          timeBGSimbox,
                                          modelGeneral->getCorrelationDirection(),
                                          modelGeneral->getRandomGen(),
                                          modelGeneral->getTimeDepthMapping(),
                                          modelGeneral->getTimeCutMapping(),
                                          modelAVOstatic->getWaveletEstimInterval(),
                                          modelAVOstatic->getWellMoveInterval(),
                                          modelAVOstatic->getFaciesEstimInterval(),
                                          modelAVOstatic);
    if(timeBGSimbox != NULL)
      delete timeBGSimbox;

    delete inputFiles;

    if (modelGeneral    == NULL || modelGeneral->getFailed()   ||
        modelAVOstatic  == NULL || modelAVOstatic->getFailed() ||
        modelAVOdynamic == NULL || modelAVOdynamic->getFailed())
      return(1);

    Crava * crava = NULL;

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

        if(modelSettings->getRunFromPanel() == false || modelSettings->getEstimateFaciesProb() == true) {
          int activeAngles = 0; //How many dimensions for local noise interpolation? Turn off for now.
          if(modelSettings->getUseLocalNoise()==true)
            activeAngles = modelSettings->getNumberOfAngles();
          spatwellfilter->doFiltering(corr,
                                      modelAVOstatic->getWells(),
                                      modelSettings->getNumberOfWells(),
                                      modelSettings->getNoVsFaciesProb(),
                                      activeAngles,
                                      crava,
                                      modelAVOdynamic->getLocalNoiseScales());
        }

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

    if (modelSettings->getDoInversion() && FFTGrid::getMaxAllowedGrids() != FFTGrid::getMaxAllocatedGrids() && modelSettings->getFileGrid()==false) {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: A memory requirement inconsistency has been detected:");
      LogKit::LogFormatted(LogKit::Warning,"\n            Maximum number of grids requested  :  %2d",FFTGrid::getMaxAllowedGrids());
      LogKit::LogFormatted(LogKit::Warning,"\n            Maximum number of grids allocated  :  %2d",FFTGrid::getMaxAllocatedGrids());
      TaskList::addTask("The memory usage estimate failed. Please send your XML-model file and\nthe logFile.txt to the CRAVA developers.");
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
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Out of memory: " << ba.what() << std::endl;
  }

#if defined(COMPILE_STORM_MODULES_FOR_RMS)
  Feature& feature = FEATURE_INVERSION_EXE;
  LicenseSystem::Instance()->CheckIn(&feature);
#endif

  return(0);
}
