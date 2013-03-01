#include <time.h>

#include "src/corr.h"
#include "src/crava.h"
#include "src/spatialwellfilter.h"
#include "src/modelsettings.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/modeltraveltimedynamic.h"
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
                                        vintage,
                                        seismicParameters);

  bool failedLoadingModel = modelAVOdynamic == NULL || modelAVOdynamic->getFailed();

  if(failedLoadingModel == false){

    Crava             * crava          = NULL;
    SpatialWellFilter * spatwellfilter = NULL;

    if(!modelSettings->getForwardModeling()){
      if (modelSettings->getDoInversion()){
        spatwellfilter = new SpatialWellFilter(modelSettings->getNumberOfWells());
        crava          = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter, seismicParameters);
        seismicParameters.setSeismicParameters(crava->getPostAlpha(), crava->getPostBeta(), crava->getPostRho(), crava->getCorrelations());
      }
    }
    else
      crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter, seismicParameters);

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

    crava = new Crava(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, spatwellfilter, seismicParameters);

    delete crava;
  }

  modelAVOstatic->deleteDynamicWells(modelGeneral->getWells(),modelSettings->getNumberOfWells());
  delete modelAVOdynamic;

  return(failedLoadingModel);
}

bool
doTimeLapseTravelTimeInversion(const ModelSettings           * modelSettings,
                               const ModelGeneral            * modelGeneral,
                               const InputFiles              * inputFiles,
                               const int                     & vintage,
                               SeismicParametersHolder       & /*seismicParameters*/)
{
  ModelTravelTimeDynamic * modelTravelTimeDynamic = NULL;

  modelTravelTimeDynamic = new ModelTravelTimeDynamic(modelSettings,
                                                      modelGeneral,
                                                      inputFiles,
                                                      vintage);

  bool failedLoadingModel = modelTravelTimeDynamic == NULL || modelTravelTimeDynamic->getFailed();

  delete modelTravelTimeDynamic;

  return(failedLoadingModel);
}

bool allocate4DGrids(SeismicParametersHolder & seismicParameters, ModelSettings * modelSettings, ModelGeneral * modelGeneral, Simbox * timeSimbox)
{
  // Parameters for generating new FFTGrids
  const int nx    = timeSimbox->getnx();
  const int ny    = timeSimbox->getny();
  const int nz    = timeSimbox->getnz();
  const int nxPad = modelSettings->getNXpad();
  const int nyPad = modelSettings->getNYpad();
  const int nzPad = modelSettings->getNZpad();

  // Allocate grids in seismicParameters. These are the 3 + 6 static grids.
  seismicParameters.allocateGrids(nx, ny, nz, nxPad, nyPad, nzPad);
  
  // Alocate dynamic grids needed for 4D inversion. 3 + 6 + 9 dynamic grids.
  // These grids are held by the class variable state4d_
  modelGeneral->complete4DBackground(nx, ny, nz, nxPad, nyPad, nzPad);

  // Merge the allocated 4D grids in the structure seismicParametersHolder
  modelGeneral->getInitial3DPriorFrom4D(seismicParameters);

  return false;
}
