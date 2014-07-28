#include <time.h>

#include "src/avoinversion.h"
#include "src/spatialwellfilter.h"
#include "src/modelsettings.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/modeltraveltimedynamic.h"
#include "src/modelgravitystatic.h"
#include "src/modelgravitydynamic.h"
#include "src/inputfiles.h"
#include "src/modelgeneral.h"
#include "src/seismicparametersholder.h"
#include "src/simbox.h"
#include "src/gravimetricinversion.h"

void setupStaticModels(ModelGeneral            *& modelGeneral,
                       ModelAVOStatic          *& modelAVOstatic,
                       //ModelGravityStatic      *& modelGravityStatic,
                       ModelSettings            * modelSettings,
                       InputFiles               * inputFiles,
                       SeismicParametersHolder  & seismicParameters,
                       CommonData               * commonData,
                       int                        i_interval)
{
  // Construct ModelGeneral object first.
  // For each data type, construct the static model class before the dynamic.
  modelGeneral    = new ModelGeneral(modelSettings,
                                     inputFiles,
                                     seismicParameters,
                                     commonData,
                                     i_interval);

  modelAVOstatic  = new ModelAVOStatic(modelSettings,
                                       inputFiles,
                                       commonData,
                                       modelGeneral->GetSimbox(),
                                       i_interval);

  // Add some logic to decide if modelGravityStatic should be created. To be done later.
  //H-Debugging
  //modelGravityStatic = new ModelGravityStatic(modelSettings,
  //                                            commonData,
  //                                            modelGeneral->GetSimbox(),
  //                                            i_interval);

  //Add in ModelTravelTimeStatic when ready
}

bool doTimeLapseAVOInversion(ModelSettings           * modelSettings,
                             ModelGeneral            * modelGeneral,
                             ModelAVOStatic          * modelAVOstatic,
                             CommonData              * commonData,
                             SeismicParametersHolder & seismicParameters,
                             int                       vintage,
                             int                       i_interval)
{
  //For intervals: Combination of doFirstAVOInversion and doTimeLapseAVOInversion

  ModelAVODynamic * modelAVOdynamic = NULL;

  // Wells are adjusted by ModelAVODynamic constructor.
  modelAVOdynamic = new ModelAVODynamic(modelSettings,
                                        modelGeneral,
                                        commonData,
                                        seismicParameters,
                                        modelGeneral->GetSimbox(),
                                        vintage,
                                        i_interval);

  bool failedLoadingModel = modelAVOdynamic == NULL || modelAVOdynamic->GetFailed();

  if(failedLoadingModel == false) {

    AVOInversion * avoinversion = new AVOInversion(modelSettings, modelGeneral, modelAVOstatic, modelAVOdynamic, seismicParameters);

    delete avoinversion;
  }

  //delete modelAVOdynamic;

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

bool
doTimeLapseGravimetricInversion(ModelSettings           * modelSettings,
                                ModelGeneral            * modelGeneral,
                                ModelGravityStatic      * modelGravityStatic,
                                CommonData              * commonData,
                                int                     & vintage,
                                SeismicParametersHolder & seismicParameters)
{
  ModelGravityDynamic * modelGravityDynamic = NULL;

  modelGravityDynamic = new ModelGravityDynamic(modelSettings,
                                                modelGravityStatic,
                                                modelGeneral->GetSimbox(),
                                                commonData,
                                                vintage,
                                                seismicParameters);

  bool failedLoadingModel = modelGravityDynamic == NULL || modelGravityDynamic->GetFailed();

   if(failedLoadingModel == false) {

    GravimetricInversion * gravimetricInversion = new GravimetricInversion(modelGeneral,
                                                                           modelGravityStatic,
                                                                           modelGravityDynamic,
                                                                           seismicParameters,
                                                                           modelSettings);

    delete gravimetricInversion;
  }

  delete modelGravityDynamic;

  return(failedLoadingModel);
}
