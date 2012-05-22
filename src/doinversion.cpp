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
                                       modelGeneral->getFailedDetails(),
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
