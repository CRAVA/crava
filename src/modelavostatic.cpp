#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelavostatic.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/corr.h"
#include "src/analyzelog.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/welldata.h"
#include "src/blockedlogs.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/waveletfilter.h"
#include "src/tasklist.h"

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "lib/lib_matr.h"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


ModelAVOStatic::ModelAVOStatic(ModelSettings        *& modelSettings,
                               const InputFiles      * inputFiles,
                               std::vector<bool>       failedGeneralDetails,
                               GridMapping           * timeCutMapping,
                               Simbox                * timeSimbox,
                               Simbox               *& timeBGSimbox,
                               Simbox                * timeSimboxConstThick,
                               std::vector<WellData *> wells)
{
  forwardModeling_        = modelSettings->getForwardModeling();
  numberOfWells_          = modelSettings->getNumberOfWells();

  priorFacies_            = NULL;
  priorFaciesProbCubes_   = NULL;
  bool failedSimbox       = failedGeneralDetails[0];

  failed_                 = false;
  bool failedExtraSurf    = false;
  bool failedPriorFacies  = false;
  bool failedRockPhysics  = false;

  bool failedLoadingModel = false;

  std::string errText("");

  if(!failedSimbox)
  {
    if (modelSettings->getForwardModeling() == false)
    {
      //
      // INVERSION/ESTIMATION
      //
      loadExtraSurfaces(waveletEstimInterval_, faciesEstimInterval_, wellMoveInterval_,
                        timeSimbox, inputFiles, errText, failedExtraSurf);

      blockLogs(wells, timeSimbox, timeBGSimbox, timeSimboxConstThick, modelSettings);

      checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
      bool estimationMode = modelSettings->getEstimationMode();
      if (estimationMode == false && !failedExtraSurf)
      {
        Simbox * timeCutSimbox = NULL;
        if (timeCutMapping != NULL)
          timeCutSimbox = timeCutMapping->getSimbox(); // For the got-enough-data test
        else
          timeCutSimbox = timeSimbox;

        processPriorFaciesProb(faciesEstimInterval_,
                               priorFacies_,
                               wells,
                               timeSimbox,
                               timeCutSimbox,
                               modelSettings,
                               failedPriorFacies,
                               errText,
                               inputFiles);
      }
    }
    else // forward modeling
      checkAvailableMemory(timeSimbox, modelSettings, inputFiles);
  }
  failedLoadingModel = failedExtraSurf || failedPriorFacies || failedRockPhysics;

  if (failedLoadingModel) {
    LogKit::WriteHeader("Error(s) while loading data");
    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel;
  failed_details_.push_back(failedExtraSurf);
  failed_details_.push_back(failedPriorFacies);
}


ModelAVOStatic::~ModelAVOStatic(void)
{

  if(waveletEstimInterval_.size() == 2) {
    if (waveletEstimInterval_[0] != NULL)
      delete waveletEstimInterval_[0];
    if (waveletEstimInterval_[1] != NULL)
      delete waveletEstimInterval_[1];
  }

  if(faciesEstimInterval_.size() == 2) {
    if (faciesEstimInterval_[0] != NULL)
      delete faciesEstimInterval_[0];
    if (faciesEstimInterval_[1] != NULL)
      delete faciesEstimInterval_[1];
  }

  if (priorFacies_ != NULL)
    delete [] priorFacies_;

  if(wellMoveInterval_.size() == 2) {
    if (wellMoveInterval_[0] != NULL)
      delete wellMoveInterval_[0];
    if (wellMoveInterval_[1] != NULL)
      delete wellMoveInterval_[1];
  }
}



void
ModelAVOStatic::blockLogs(std::vector<WellData *> & wells,
                          Simbox                  * timeSimbox,
                          Simbox                  * timeBGSimbox,
                          Simbox                  * timeSimboxConstThick,
                          ModelSettings          *& modelSettings)
{
  int     nWells         = modelSettings->getNumberOfWells();

  if(nWells > 0) {
    for (int i=0 ; i<nWells ; i++)
    {
      wells[i]->findMeanVsVp(waveletEstimInterval_);

      wells[i]->setBlockedLogsOrigThick( new BlockedLogs(wells[i], timeSimbox, modelSettings->getRunFromPanel()) );
      wells[i]->setBlockedLogsConstThick( new BlockedLogs(wells[i], timeSimboxConstThick) );
      if (timeBGSimbox==NULL)
        wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeSimbox) ); // Need a copy constructor?
      else
        wells[i]->setBlockedLogsExtendedBG( new BlockedLogs(wells[i], timeBGSimbox) );
    }
  }
}

void ModelAVOStatic::processFaciesInformation(ModelSettings     *& modelSettings,
                                              const InputFiles   * inputFiles,
                                              std::string        & tmpErrText,
                                              int                & error) const
{
  int nFacies = modelSettings->getNumberOfFacies();
  if(nFacies == 0) {
    // Get facies information from rock physics prior model
    std::vector<std::string> rock_name = modelSettings->getRockName();

    nFacies = static_cast<int>(rock_name.size());

    for(int i=0; i<nFacies; i++) {
      modelSettings->addFaciesName(rock_name[i]);
      modelSettings->addFaciesLabel(i);
    }
  }
  else {
    std::vector<std::string> facies_name = modelSettings->getFaciesNames();
    std::vector<std::string> rock_name   = modelSettings->getRockName();

    // Compare names in wells with names given in rock physics prior model
    if(rock_name.size()>0) {
      int nRocks  = static_cast<int>(rock_name.size());
      if(nRocks > nFacies) {
        tmpErrText += "Problem with facies logs. The number of rocks in the rock physics prior model is larger than the number of facies found in the wells.\n";
        error++;
      }
      for(int i=0; i<nFacies; i++) {
        int j=0;
        while(j<nRocks) {
          if(facies_name[i] == rock_name[j])
            break;
          j++;
        }
        if(j == nRocks) {
          tmpErrText += "Problem with facies logs. Facies "+facies_name[i]+" found in a well is not one of the rocks given in rock physics prior model\n";
          error++;
        }
      }
    }
    // Compare names in wells with names given in .xml-file
    if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
    {
      typedef std::map<std::string,float> mapType;
      mapType myMap = modelSettings->getPriorFaciesProb();

      for(int i=0;i<nFacies;i++)
      {
        mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
        if (iter==myMap.end())
        {
          tmpErrText += "Problem with facies logs. Facies "+modelSettings->getFaciesName(i)+" is not one of the facies given in the xml-file.\n";
          error++;
        }
      }
    }
    // Compare names in wells with names given as input in proability cubes
    else if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
    {
      typedef std::map<std::string,std::string> mapType;
      mapType myMap = inputFiles->getPriorFaciesProbFile();

      for(int i=0;i<nFacies;i++)
      {
        mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
        if (iter==myMap.end())
        {
          tmpErrText += "Problem with facies logs. Facies "+modelSettings->getFaciesName(i)+" is not one of the facies given in the xml-file.\n";
          error++;
        }
      }
    }
  }
}

void
ModelAVOStatic::checkAvailableMemory(Simbox        * timeSimbox,
                                   ModelSettings * modelSettings,
                                   const InputFiles    * inputFiles)
{
  LogKit::WriteHeader("Estimating amount of memory needed");
  //
  // Find the size of first seismic volume
  //
  float memOneSeis = 0.0f;
  if (inputFiles->getNumberOfSeismicFiles(0) > 0 && inputFiles->getSeismicFile(0,0) != "") {
    memOneSeis = static_cast<float> (NRLib::FindFileSize(inputFiles->getSeismicFile(0,0)));
  }

  //
  // Find the size of one grid
  //
  FFTGrid * dummyGrid = new FFTGrid(timeSimbox->getnx(),
                                    timeSimbox->getny(),
                                    timeSimbox->getnz(),
                                    modelSettings->getNXpad(),
                                    modelSettings->getNYpad(),
                                    modelSettings->getNZpad());
  long long int gridSizePad = static_cast<long long int>(4)*dummyGrid->getrsize();

  delete dummyGrid;
  dummyGrid = new FFTGrid(timeSimbox->getnx(),
                          timeSimbox->getny(),
                          timeSimbox->getnz(),
                          timeSimbox->getnx(),
                          timeSimbox->getny(),
                          timeSimbox->getnz());
  long long int gridSizeBase = 4*dummyGrid->getrsize();
  delete dummyGrid;
  int nGridParameters  = 3;                                      // Vp + Vs + Rho, padded
  int nGridBackground  = 3;                                      // Vp + Vs + Rho, padded
  int nGridCovariances = 6;                                      // Covariances, padded
  int nGridSeismicData = modelSettings->getNumberOfAngles(0);     // One for each angle stack, padded

  int nGridFacies       = modelSettings->getNumberOfFacies()+1;   // One for each facies, one for undef, unpadded.
  int nGridHistograms   = modelSettings->getNumberOfFacies();     // One for each facies, 2MB.
  int nGridKriging      = 1;                                      // One grid for kriging, unpadded.
  int nGridCompute      = 1;                                      // Computation grid, padded (for convenience)
  int nGridFileMode     = 1;                                      // One grid for intermediate file storage

  int nGrids;
  long long int gridMem;
  if(modelSettings->getForwardModeling() == true) {
    if (modelSettings->getFileGrid())  // Use disk buffering
      nGrids = nGridFileMode;
    else
      nGrids = nGridParameters + 1;

    gridMem = nGrids*gridSizePad;
  }
  else {
    if (modelSettings->getFileGrid()) { // Use disk buffering
      nGrids = nGridFileMode;
      if(modelSettings->getKrigingParameter() > 0) {
        nGrids += nGridKriging;
      }
      if(modelSettings->getNumberOfSimulations() > 0)
        nGrids = nGridParameters;
      if(modelSettings->getUseLocalNoise(0)) {
        nGrids = 2*nGridParameters;
      }

      gridMem = nGrids*gridSizePad;
    }
    else {
      //baseP and baseU are the padded and unpadde grids allocated at each peak.
      int baseP = nGridParameters + nGridCovariances;
      if(modelSettings->getUseLocalNoise(0) == true || (modelSettings->getEstimateFaciesProb() && modelSettings->getFaciesProbRelative()))
        baseP += nGridBackground;
      int baseU = 0;
      if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        baseU += modelSettings->getNumberOfFacies();

      //First peak: At inversion
      int peak1P = baseP + nGridSeismicData; //Need seismic data as well here.
      int peak1U = baseU;

      long long int peakGridMem = peak1P*gridSizePad + peak1U*gridSizeBase; //First peak must be currently largest.
      int peakNGrid   = peak1P;                                             //Also in number of padded grids

      if(modelSettings->getNumberOfSimulations() > 0) { //Second possible peak when simulating.
        int peak2P = baseP + 3; //Three extra parameter grids for simulated parameters.
        if(modelSettings->getUseLocalNoise(0) == true &&
           (modelSettings->getEstimateFaciesProb() == false || modelSettings->getFaciesProbRelative() == false))
          peak2P -= nGridBackground; //Background grids are released before simulation in this case.
        int peak2U = baseU;     //Base level is the same, but may increase.
        bool computeGridUsed = ((modelSettings->getOutputGridsElastic() & (IO::AI + IO::LAMBDARHO + IO::LAMELAMBDA + IO::LAMEMU + IO::MURHO + IO::POISSONRATIO + IO::SI + IO::VPVSRATIO)) > 0);
        if(computeGridUsed == true)
          peak2P += nGridCompute;
        else if(modelSettings->getKrigingParameter() > 0) //Note the else, since this grid will use same memory as computation grid if both are active.
          peak2U += nGridKriging;

        if(peak2P > peakNGrid)
          peakNGrid = peak2P;

        long long int peak2Mem = peak2P*gridSizePad + peak2U*gridSizeBase;
        if(peak2Mem > peakGridMem)
          peakGridMem = peak2Mem;
      }

      if(modelSettings->getEstimateFaciesProb() == true) {//Third possible peak when computing facies prob.
        int peak3P = baseP;                //No extra padded grids, so this one can not peak here.
        int peak3U = baseU + nGridFacies;  //But this one will, and may trigger new memory max.
        if((modelSettings->getOtherOutputFlag() & IO::FACIES_LIKELIHOOD) > 0)
          peak3U += 1; //Also needs to store seismic likelihood.

        long long int peak3Mem = peak3P*gridSizePad + peak3U*gridSizeBase + 2000000*nGridHistograms; //These are 2MB when Vs is used.
        if(peak3Mem > peakGridMem)
          peakGridMem = peak3Mem;
      }
      nGrids  = peakNGrid;
      gridMem = peakGridMem;
    }
  }
  FFTGrid::setMaxAllowedGrids(nGrids);
  if(modelSettings->getDebugFlag()>0)
    FFTGrid::setTerminateOnMaxGrid(true);

  int   workSize    = 2500 + static_cast<int>( 0.65*gridSizePad); //Size of memory used beyond grids.

  float mem0        = 4.0f * workSize;
  float mem1        = static_cast<float>(gridMem);
  float mem2        = static_cast<float>(modelSettings->getNumberOfAngles(0))*gridSizePad + memOneSeis; //Peak memory when reading seismic, overestimated.

  float neededMem   = mem0 + std::max(mem1, mem2);

  float megaBytes   = neededMem/(1024.f*1024.f);
  float gigaBytes   = megaBytes/1024.f;

  LogKit::LogFormatted(LogKit::High,"\nMemory needed for reading seismic data       : %10.2f MB\n",mem2/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding internal grids (%2d): %10.2f MB\n",nGrids, mem1/(1024.f*1024.f));
  LogKit::LogFormatted(LogKit::High,  "Memory needed for holding other entities     : %10.2f MB\n",mem0/(1024.f*1024.f));

  if (megaBytes > 1000.0f)
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f gigaBytes\n",gigaBytes);
  else
    LogKit::LogFormatted(LogKit::Low,"\nMemory needed by CRAVA:  %.1f megaBytes\n",megaBytes);

  if(mem2>mem1)
    LogKit::LogFormatted(LogKit::Low,"\n This estimate is too high because seismic data are cut to fit the internal grid\n");
  if (!modelSettings->getFileGrid()) {
    //
    // Check if we can hold everything in memory.
    //
    modelSettings->setFileGrid(false);
    char ** memchunk  = new char*[nGrids];

    int i = 0;
    try {
      for(i = 0 ; i < nGrids ; i++)
        memchunk[i] = new char[static_cast<size_t>(gridSizePad)];
    }
    catch (std::bad_alloc& ) //Could not allocate memory
    {
      modelSettings->setFileGrid(true);
      LogKit::LogFormatted(LogKit::Low,"Not enough memory to hold all grids. Using file storage.\n");
    }

    for(int j=0 ; j<i ; j++)
      delete [] memchunk[j];
    delete [] memchunk;
  }
}

void ModelAVOStatic::processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
                                            float                       *& priorFacies,
                                            std::vector<WellData *>        wells,
                                            Simbox                       * timeSimbox,
                                            Simbox                       * timeCutSimbox,
                                            ModelSettings                * modelSettings,
                                            bool                         & failed,
                                            std::string                  & errTxt,
                                            const InputFiles             * inputFiles)
{
  if (modelSettings->getEstimateFaciesProb())
  {
    LogKit::WriteHeader("Prior Facies Probabilities");

    int error = 0;
    std::string tmpErrText = "";
    processFaciesInformation(modelSettings,
                             inputFiles,
                             tmpErrText,
                             error);
    if (error>0)
      errTxt += "Prior facies probabilities failed.\n"+tmpErrText;


    int nFacies = modelSettings->getNumberOfFacies();

    if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_WELLS)
    {
      if (nFacies > 0)
      {
        int   nz      = timeSimbox->getnz();
        float dz      = static_cast<float>(timeSimbox->getdz());
        int   nWells  = modelSettings->getNumberOfWells();
        int   nFacies = modelSettings->getNumberOfFacies();
        int   ndata   = nWells*nz;

        int ** faciesCount = new int * [nWells];
        for (int w = 0 ; w < nWells ; w++)
          faciesCount[w] = new int[nFacies];

        for (int w = 0 ; w < nWells ; w++)
          for (int i = 0 ; i < nFacies ; i++)
            faciesCount[w][i] = 0;

        int * faciesLog = new int[ndata];   // NB! *internal* log numbering (0, 1, 2, ...)
        for (int i = 0 ; i < ndata ; i++)
          faciesLog[i] = IMISSING;

        float * vtAlpha   = new float[nz];  // vt = vertical trend
        float * vtBeta    = new float[nz];
        float * vtRho     = new float[nz];
        int   * vtFacies  = new int[nz];

        int nUsedWells = 0;

        for (int w = 0 ; w < nWells ; w++)
        {
          if(wells[w]->getNFacies() > 0) // Well has facies log
          {
            //
            // Note that we use timeSimbox to calculate prior facies probabilities
            // instead of the simbox with parallel top and base surfaces. This
            // will make the prior probabilities slightly different, but that
            // should not be a problem.
            //
            BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
            int nBlocks = bl->getNumberOfBlocks();
            //
            // Set facies data outside facies estimation interval IMISSING
            //
            int * blFaciesLog = new int[nBlocks];
            Utils::copyVector(bl->getFacies(), blFaciesLog, nBlocks);

            if (faciesEstimInterval.size() > 0) {
              const double * xPos  = bl->getXpos();
              const double * yPos  = bl->getYpos();
              const double * zPos  = bl->getZpos();
              for (int i = 0 ; i < nBlocks ; i++) {
                const double zTop  = faciesEstimInterval[0]->GetZ(xPos[i],yPos[i]);
                const double zBase = faciesEstimInterval[1]->GetZ(xPos[i],yPos[i]);
                if ( (zPos[i] - 0.5*dz) < zTop || (zPos[i] + 0.5*dz) > zBase)
                  blFaciesLog[i] = IMISSING;
              }
            }

            bl->getVerticalTrend(bl->getAlpha(),vtAlpha);
            bl->getVerticalTrend(bl->getBeta(),vtBeta);
            bl->getVerticalTrend(bl->getRho(),vtRho);
            bl->getVerticalTrend(blFaciesLog,vtFacies);
            delete [] blFaciesLog;

            for(int i=0 ; i<nz ; i++)
            {
              int facies;
              if(vtAlpha[i] != RMISSING && vtBeta[i] != RMISSING && vtRho[i] != RMISSING)
                facies = vtFacies[i];
              else
                facies = IMISSING;

              faciesLog[w*nz + i] = facies;
              if(facies != IMISSING)
                faciesCount[w][facies]++;
            }
            nUsedWells++;
          }
        }
        delete [] vtAlpha;
        delete [] vtBeta;
        delete [] vtRho;
        delete [] vtFacies;

        if (nUsedWells > 0) {
          //
          // Probabilities
          //
          LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each blocked well: \n");
          LogKit::LogFormatted(LogKit::Low,"\nBlockedWell              ");
          for (int i = 0 ; i < nFacies ; i++)
            LogKit::LogFormatted(LogKit::Low,"%12s ",modelSettings->getFaciesName(i).c_str());
          LogKit::LogFormatted(LogKit::Low,"\n");
          for (int i = 0 ; i < 24+13*nFacies ; i++)
            LogKit::LogFormatted(LogKit::Low,"-");
          LogKit::LogFormatted(LogKit::Low,"\n");
          for (int w = 0 ; w < nWells ; w++)
          {
            if(wells[w]->getNFacies() > 0) // Well has facies log
            {
              float tot = 0.0;
              for (int i = 0 ; i < nFacies ; i++) {
                tot += static_cast<float>(faciesCount[w][i]);
              }

              LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[w]->getWellname().c_str());
              for (int i = 0 ; i < nFacies ; i++) {
                float faciesProb = static_cast<float>(faciesCount[w][i])/tot;
                LogKit::LogFormatted(LogKit::Low," %12.4f",faciesProb);
              }
              LogKit::LogFormatted(LogKit::Low,"\n");
            }
          }
          LogKit::LogFormatted(LogKit::Low,"\n");
          //
          // Counts
          //
          LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each blocked well: \n");

          LogKit::LogFormatted(LogKit::Medium,"\nBlockedWell              ");
          for (int i = 0 ; i < nFacies ; i++)
            LogKit::LogFormatted(LogKit::Medium,"%12s ",modelSettings->getFaciesName(i).c_str());
          LogKit::LogFormatted(LogKit::Medium,"\n");
          for (int i = 0 ; i < 24+13*nFacies ; i++)
            LogKit::LogFormatted(LogKit::Medium,"-");
          LogKit::LogFormatted(LogKit::Medium,"\n");
          for (int w = 0 ; w < nWells ; w++)
          {
            if(wells[w]->getNFacies() > 0)
            {
              float tot = 0.0;
              for (int i = 0 ; i < nFacies ; i++)
                tot += static_cast<float>(faciesCount[w][i]);
              LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[w]->getWellname().c_str());
              for (int i = 0 ; i < nFacies ; i++) {
                LogKit::LogFormatted(LogKit::Medium," %12d",faciesCount[w][i]);
              }
              LogKit::LogFormatted(LogKit::Medium,"\n");
            }
          }
          LogKit::LogFormatted(LogKit::Medium,"\n");

          for (int w = 0 ; w < nWells ; w++)
            delete [] faciesCount[w];
          delete [] faciesCount;

          //
          // Make prior facies probabilities
          //
          float sum = 0.0f;
          int * nData = new int[nFacies];
          for(int i=0 ; i<nFacies ; i++)
            nData[i] = 0;

          for(int i=0 ; i<ndata ; i++) {
            if(faciesLog[i] != IMISSING) {
              nData[faciesLog[i]]++;
            }
          }
          delete [] faciesLog;

          for(int i=0 ; i<nFacies ; i++)
            sum += nData[i];

          if (sum > 0) {
            LogKit::LogFormatted(LogKit::Low,"Facies probabilities based on all blocked wells:\n\n");
            LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
            LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
            priorFacies = new float[nFacies];
            for(int i=0 ; i<nFacies ; i++) {
              priorFacies[i] = float(nData[i])/sum;
              LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",modelSettings->getFaciesName(i).c_str(),priorFacies[i]);
            }
          }
          else {
            LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No valid facies log entries have been found\n");
            modelSettings->setEstimateFaciesProb(false);
            TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");

          }
          delete [] nData;
        }
        else
        {
          LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but there");
          LogKit::LogFormatted(LogKit::Warning,"\n         are no wells with facies available and CRAVA will therefore not");
          LogKit::LogFormatted(LogKit::Warning,"\n         be able to estimate these probabilities...\n");
          modelSettings->setEstimateFaciesProb(false);

          TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
        }
      }
      else
      {
        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but no facies");
        LogKit::LogFormatted(LogKit::Warning,"\n         have been found and CRAVA will therefore not be able to estimate");
        LogKit::LogFormatted(LogKit::Warning,"\n         these probabilities...\n");
        modelSettings->setEstimateFaciesProb(false);
        TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
      }
    }
    else if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
    {
      priorFacies = new float[nFacies];
      typedef std::map<std::string,float> mapType;
      mapType myMap = modelSettings->getPriorFaciesProb();

      for(int i=0;i<nFacies;i++)
      {
        mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));
        if(iter!=myMap.end())
          priorFacies[i] = iter->second;
        else
        {
          LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",modelSettings->getFaciesName(i).c_str());
          modelSettings->setEstimateFaciesProb(false);
          TaskList::addTask("Check that facies " +NRLib::ToString(modelSettings->getFaciesName(i).c_str())+" is given a prior probability in the xml-file");
        }
      }
      LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
      LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
      for(int i=0 ; i<nFacies ; i++) {
        LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",modelSettings->getFaciesName(i).c_str(),priorFacies[i]);
      }

    }
    else if(modelSettings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
    {
      readPriorFaciesProbCubes(inputFiles,
                               modelSettings,
                               priorFaciesProbCubes_,
                               timeSimbox,
                               timeCutSimbox,
                               errTxt,
                               failed);

       typedef std::map<std::string,std::string> mapType;
       mapType myMap = inputFiles->getPriorFaciesProbFile();

       LogKit::LogFormatted(LogKit::Low,"Facies         Probability in file\n");
       LogKit::LogFormatted(LogKit::Low,"----------------------------------\n");
       for(mapType::iterator it=myMap.begin();it!=myMap.end();it++)
         LogKit::LogFormatted(LogKit::Low,"%-15s %10s\n",(it->first).c_str(),(it->second).c_str());

    }
  }
}
void ModelAVOStatic::readPriorFaciesProbCubes(const InputFiles  * inputFiles,
                                              ModelSettings     * modelSettings,
                                              FFTGrid         **& priorFaciesProbCubes,
                                              Simbox            * timeSimbox,
                                              Simbox            * timeCutSimbox,
                                              std::string       & errTxt,
                                              bool              & failed)
{
  int nFacies = modelSettings->getNumberOfFacies();
  priorFaciesProbCubes = new FFTGrid*[nFacies];

  typedef std::map<std::string,std::string> mapType;
  mapType myMap = inputFiles->getPriorFaciesProbFile();
  for(int i=0;i<nFacies;i++)
  {
    mapType::iterator iter = myMap.find(modelSettings->getFaciesName(i));

    if(iter!=myMap.end())
    {
      const std::string & faciesProbFile = iter->second;
      const SegyGeometry      * dummy1 = NULL;
      const TraceHeaderFormat * dummy2 = NULL;
      const float               offset = modelSettings->getSegyOffset(0); //Facies estimation only allowed for one time lapse
      std::string errorText("");
      ModelGeneral::readGridFromFile(faciesProbFile,
                       "priorfaciesprob",
                       offset,
                       priorFaciesProbCubes[i],
                       dummy1,
                       dummy2,
                       FFTGrid::PARAMETER,
                       timeSimbox,
                       timeCutSimbox,
                       modelSettings,
                       errorText,
                       true);
      if(errorText != "")
      {
        errorText += "Reading of file \'"+faciesProbFile+"\' for prior facies probability for facies \'"
                     +modelSettings->getFaciesName(i)+"\' failed\n";
        errTxt += errorText;
        failed = true;
      }
    }
    else
    {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",
                           modelSettings->getFaciesName(i).c_str());
      TaskList::addTask("Check that facies "+NRLib::ToString(modelSettings->getFaciesName(i).c_str())+" is given prior probability in the xml-file");
      modelSettings->setEstimateFaciesProb(false);
      break;
    }
  }
}

void
ModelAVOStatic::loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
                                  std::vector<Surface *> & faciesEstimInterval,
                                  std::vector<Surface *> & wellMoveInterval,
                                  Simbox                 * timeSimbox,
                                  const InputFiles             * inputFiles,
                                  std::string            & errText,
                                  bool                   & failed)
{
  const double x0 = timeSimbox->getx0();
  const double y0 = timeSimbox->gety0();
  const double lx = timeSimbox->getlx();
  const double ly = timeSimbox->getly();
  const int    nx = timeSimbox->getnx();
  const int    ny = timeSimbox->getny();
  //
  // Get wavelet estimation interval
  //
  const std::string & topWEI  = inputFiles->getWaveletEstIntFileTop(0); //Same for all time lapses
  const std::string & baseWEI = inputFiles->getWaveletEstIntFileBase(0);//Same for all time lapses

  if (topWEI != "" && baseWEI != "") {
    waveletEstimInterval.resize(2);
    try {
      if (NRLib::IsNumber(topWEI))
        waveletEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWEI.c_str()));
      else {
        Surface tmpSurf(topWEI);
        waveletEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }

    try {
      if (NRLib::IsNumber(baseWEI))
        waveletEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWEI.c_str()));
      else {
        Surface tmpSurf(baseWEI);
        waveletEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
  //
  // Get facies estimation interval
  //
  const std::string & topFEI  = inputFiles->getFaciesEstIntFile(0);
  const std::string & baseFEI = inputFiles->getFaciesEstIntFile(1);

  if (topFEI != "" && baseFEI != "") {
    faciesEstimInterval.resize(2);
    try {
      if (NRLib::IsNumber(topFEI))
        faciesEstimInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topFEI.c_str()));
      else {
        Surface tmpSurf(topFEI);
        faciesEstimInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }

    try {
      if (NRLib::IsNumber(baseFEI))
        faciesEstimInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseFEI.c_str()));
      else {
        Surface tmpSurf(baseFEI);
        faciesEstimInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
  //
  // Get well move interval
  //
  const std::string & topWMI  = inputFiles->getWellMoveIntFile(0);
  const std::string & baseWMI = inputFiles->getWellMoveIntFile(1);

  if (topWMI != "" && baseWMI != "") {
    wellMoveInterval.resize(2);
    try {
      if (NRLib::IsNumber(topWMI))
        wellMoveInterval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWMI.c_str()));
      else {
        Surface tmpSurf(topWMI);
        wellMoveInterval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }

    try {
      if (NRLib::IsNumber(baseWMI))
        wellMoveInterval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWMI.c_str()));
      else {
        Surface tmpSurf(baseWMI);
        wellMoveInterval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = true;
    }
  }
}

void ModelAVOStatic::writeWells(std::vector<WellData *> wells, ModelSettings * modelSettings) const
{
  int nWells  = modelSettings->getNumberOfWells();
  for(int i=0;i<nWells;i++)
    wells[i]->writeWell(modelSettings->getWellFormatFlag());
}

void ModelAVOStatic::writeBlockedWells(std::vector<WellData *> wells, ModelSettings * modelSettings) const
{
  int nWells  = modelSettings->getNumberOfWells();
  for(int i=0;i<nWells;i++)
    wells[i]->getBlockedLogsOrigThick()->writeWell(modelSettings);
}

void ModelAVOStatic::addSeismicLogs(std::vector<WellData *> wells,
                                    FFTGrid      ** seisCube,
                                    ModelSettings * modelSettings,
                                    int             nAngles)
{
  int nWells  = modelSettings->getNumberOfWells();

    for (int iAngle = 0 ; iAngle < nAngles ; iAngle++)
    {
      seisCube[iAngle]->setAccessMode(FFTGrid::RANDOMACCESS);
      for(int i=0;i<nWells;i++)
        wells[i]->getBlockedLogsOrigThick()->setLogFromGrid(seisCube[iAngle],iAngle,nAngles,"SEISMIC_DATA");
      seisCube[iAngle]->endAccess();
  }
}

void ModelAVOStatic::generateSyntheticSeismic(Wavelet      ** wavelet,
                                              std::vector<WellData *> wells,
                                              float        ** reflectionMatrix,
                                              Simbox        * timeSimbox,
                                              ModelSettings * modelSettings,
                                              int             nAngles)
{
  int nWells  = modelSettings->getNumberOfWells();
  int nzp     = modelSettings->getNZpad();
  int nz      = timeSimbox->getnz();

  int i;

  for( i=0; i<nWells; i++ )
  {
    if( wells[i]->isDeviated() == false )
      wells[i]->getBlockedLogsOrigThick()->generateSyntheticSeismic(reflectionMatrix,nAngles,wavelet,nz,nzp,timeSimbox);
  }
}



void ModelAVOStatic::deleteDynamicWells(std::vector<WellData *> wells,
                                        int                     nWells)
{
  for(int i=0; i<nWells; i++){
    wells[i]->getBlockedLogsConstThick()->deleteDynamicBlockedLogs();
    wells[i]->getBlockedLogsExtendedBG()->deleteDynamicBlockedLogs();
    wells[i]->getBlockedLogsOrigThick()->deleteDynamicBlockedLogs();
  }
}
