#include <time.h>

#include "src/corr.h"
#include "src/crava.h"
#include "src/spatialwellfilter.h"
#include "src/modelsettings.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/inputfiles.h"
#include "src/modelgeneral.h"
#include "src/seismicparametersholder.h"
#include "src/simbox.h"

void setupStaticModels(ModelGeneral   *& modelGeneral,
                       ModelAVOStatic *& modelAVOstatic,
                       ModelSettings   * modelSettings,
                       InputFiles      * inputFiles,
                       Simbox         *& timeBGSimbox)
{
  // Construct ModelGeneral object first.
  // For each data type, construct the static model class before the dynamic.
  modelGeneral    = new ModelGeneral(modelSettings, inputFiles, timeBGSimbox);
  modelAVOstatic  = new ModelAVOStatic(modelSettings,
                                       modelGeneral,
                                       inputFiles,
                                       modelGeneral->getTimeCutMapping(),
                                       modelGeneral->getTimeSimbox(),
                                       timeBGSimbox,
                                       modelGeneral->getTimeSimboxConstThick(),
                                       modelGeneral->getWells());
}


bool doFirstAVOInversion(ModelSettings           * modelSettings,
                         ModelGeneral            * modelGeneral,
                         ModelAVOStatic          * modelAVOstatic,
                         SeismicParametersHolder & seismicParameters,
                         InputFiles              * inputFiles,
                         int                       vintage,
                         Simbox                  * timeBGSimbox)
{

  ModelAVODynamic * modelAVOdynamic = NULL;

  // Wells are adjusted by ModelAVODynamic constructor.
  modelAVOdynamic = new ModelAVODynamic(modelSettings,
                                        inputFiles,
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
                                        modelAVOstatic,
                                        modelGeneral,
                                        vintage);

  bool failedLoadingModel = modelAVOdynamic == NULL || modelAVOdynamic->getFailed();

  if(failedLoadingModel == false){

    Crava             * crava          = NULL;
    SpatialWellFilter * spatwellfilter = NULL;

    if(!modelSettings->getForwardModeling()){
      if (modelSettings->getDoInversion()){
        spatwellfilter = new SpatialWellFilter(modelSettings->getNumberOfWells());
        crava          = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter);
        seismicParameters.setSeismicParameters(crava->getPostAlpha(), crava->getPostBeta(), crava->getPostRho(), crava->getCorrelations());
      }
    }
    else
      crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter);

    delete crava;
    delete spatwellfilter;
  }

  modelAVOstatic->deleteDynamicWells(modelGeneral->getWells(),modelSettings->getNumberOfWells());
  delete modelAVOdynamic;

  return(failedLoadingModel);
}


bool doTimeLapseAVOInversion(ModelSettings           * modelSettings,
                             ModelGeneral            * modelGeneral,
                             ModelAVOStatic          * modelAVOstatic,
                             InputFiles              * inputFiles,
                             SeismicParametersHolder & seismicParameters,
                             int                       vintage)
{
  ModelAVODynamic * modelAVOdynamic = NULL;

  modelAVOdynamic = new ModelAVODynamic(modelSettings,
                                        inputFiles,
                                        modelAVOstatic,
                                        modelGeneral,
                                        seismicParameters,
                                        modelGeneral->getTimeSimbox(),
                                        modelGeneral->getCorrelationDirection(),
                                        modelGeneral->getTimeDepthMapping(),
                                        modelGeneral->getTimeCutMapping(),
                                        vintage);

  bool failedLoadingModel = modelAVOdynamic == NULL || modelAVOdynamic->getFailed();

  if(failedLoadingModel == false){

    Crava             * crava          = NULL;
    SpatialWellFilter * spatwellfilter = NULL;

    crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter);

    delete crava;
  }

  modelAVOstatic->deleteDynamicWells(modelGeneral->getWells(),modelSettings->getNumberOfWells());
  delete modelAVOdynamic;

  return(failedLoadingModel);
}

bool allocate4DGrids(SeismicParametersHolder & seismicParameters, ModelSettings * modelSettings, ModelGeneral * modelGeneral, Simbox * timeBGSimbox){

  // FFTGrids needed for 4D inversion. To be allocated in memory.
  FFTGrid *vp;
  FFTGrid *vs;
  FFTGrid *rho;
  FFTGrid *crCovVpVp;
  FFTGrid *crCovVpVs;
  FFTGrid *crCovVpRho;
  FFTGrid *crCovVsVs;
  FFTGrid *crCovVsRho;
  FFTGrid *crCovRhoRho;

  // Parameters for generating new FFTGrids
  const int nx    = timeBGSimbox->getnx();
  const int ny    = timeBGSimbox->getny();
  const int nz    = timeBGSimbox->getnz();
  const int nxPad = modelSettings->getNXpad();
  const int nyPad = modelSettings->getNYpad();
  const int nzPad = modelSettings->getNZpad();

  // Creating the 4D grids
  vp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  vp->createRealGrid();
  vs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  vs->createRealGrid();
  rho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  rho->createRealGrid();
  crCovVpVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovVpVp->createRealGrid();
  crCovVpVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovVpVs->createRealGrid();
  crCovVpRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovVpRho->createRealGrid();
  crCovVsVs  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovVsVs->createRealGrid();
  crCovVsRho  = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovVsRho->createRealGrid();
  crCovRhoRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  crCovRhoRho->createRealGrid();

  // Merge the allocated 4D grids in the structure seismicParametersHolder
  modelGeneral->get3DPriorFrom4D(seismicParameters, vp, vs, rho, crCovVpVp, crCovVpVs, crCovVpRho, crCovVsVs, crCovVsRho, crCovRhoRho);

  return 1;
}
