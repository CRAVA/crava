#ifndef DOINVERSION_H
#define DOINVERSION_H

#include <stdio.h>

class ModelSettings;
class ModelAVODynamic;
class ModelAVOStatic;
class ModelGravityStatic;
class ModelGravityDynamic;
class ModelGeneral;
class InputFiles;
class Simbox;
class SeismicParametersHolder;

void setupStaticModels(ModelGeneral            *& modelGeneral,
                       ModelAVOStatic          *& modelAVOstatic,
                       //ModelGravityStatic      *& modelGravityStatic,
                       ModelSettings            * modelSettings,
                       InputFiles               * inputFiles,
                       SeismicParametersHolder  & seismicParameters,
                       CommonData               * commonData,
                       int                        i_interval);

bool doTimeLapseAVOInversion(ModelSettings           * modelSettings,
                             ModelGeneral            * modelGeneral,
                             ModelAVOStatic          * modelAVOstatic,
                             CommonData              * commonData,
                             SeismicParametersHolder & seismicParameters,
                             //InputFiles              * inputFiles,
                             int                       vintage,
                             int                       i_interval);

bool doTimeLapseTravelTimeInversion(const ModelSettings     * modelSettings,
                                    ModelGeneral            * modelGeneral,
                                    ModelTravelTimeStatic   * modelTravelTimeStatic,
                                    const InputFiles        * inputFiles,
                                    const int               & vintage,
                                    SeismicParametersHolder & seismicParameters);

bool doTimeLapseGravimetricInversion(ModelSettings           * modelSettings,
                                     ModelGeneral            * modelGeneral,
                                     ModelGravityStatic      * modelGravityStatic,
                                     CommonData              * commonData,
                                     //InputFiles              * inputFiles,
                                     int                     & vintage,
                                     SeismicParametersHolder & seismicParameters);


#endif

