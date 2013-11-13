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

//void setupStaticModels(ModelGeneral           *& modelGeneral,
//                       ModelAVOStatic         *& modelAVOstatic,
//                       ModelGravityStatic     *& modelGravityStatic,
//                       ModelSettings           * modelSettings,
//                       InputFiles              * inputFiles,
//                       SeismicParametersHolder & seismicParameters,
//                       Simbox                 *& timeBGSimbox);

void setupStaticModels(ModelGeneral            *& modelGeneral,
                       ModelAVOStatic          *& modelAVOstatic,
                       ModelGravityStatic      *& modelGravityStatic,
                       ModelSettings            * modelSettings,
                       InputFiles               * inputFiles,
                       SeismicParametersHolder  & seismicParameters,
                       CommonData               * commonData,
                       int                        i_interval);

//bool doFirstAVOInversion(ModelSettings           * modelSettings,
//                         ModelGeneral            * modelGeneral,
//                         ModelAVOStatic          * modelAVOStatic,
//                         SeismicParametersHolder & seismicParameters,
//                         InputFiles              * inputFiles,
//                         int                       sortedVintage,
//                         Simbox                  * timeBGSimbox);

bool doFirstAVOInversion(ModelSettings           * modelSettings, //Intervals
                         ModelGeneral            * modelGeneral,
                         ModelAVOStatic          * modelAVOstatic,
                         CommonData              * commonData,
                         SeismicParametersHolder & seismicParameters,
                         InputFiles              * inputFiles,
                         int                       vintage,
                         int                       i_interval);

bool doTimeLapseAVOInversion(ModelSettings           * modelSettings,
                             ModelGeneral            * modelGeneral,
                             ModelAVOStatic          * modelAVOStatic,
                             InputFiles              * inputFiles,
                             SeismicParametersHolder & seismicParameters,
                             int                       sortedVintage);

bool doTimeLapseTravelTimeInversion(const ModelSettings           * modelSettings,
                                    const ModelGeneral            * modelGeneral,
                                    const InputFiles              * inputFiles,
                                    const int                     & vintage,
                                    SeismicParametersHolder       & seismicParameters);

bool doTimeLapseGravimetricInversion(ModelSettings           * modelSettings,
                                     ModelGeneral            * modelGeneral,
                                     ModelGravityStatic      * modelGravityStatic,
                                     InputFiles              * inputFiles,
                                     int                     & vintage,
                                     SeismicParametersHolder & seismicParameters);


#endif

