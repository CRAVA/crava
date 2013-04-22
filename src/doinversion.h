#ifndef DOINVERSION_H
#define DOINVERSION_H

#include <stdio.h>

class ModelSettings;
class ModelAVODynamic;
class ModelAVOStatic;
class ModelGeneral;
class InputFiles;
class Simbox;
class SeismicParametersHolder;

void setupStaticModels(ModelGeneral            *& modelGeneral,
                       ModelAVOStatic          *& modelAVOstatic,
                       ModelSettings            * modelSettings,
                       InputFiles               * inputFiles,
                       SeismicParametersHolder  & seismicParameters,
                       Simbox                  *& timeBGSimbox);

bool doFirstAVOInversion(ModelSettings           * modelSettings,
                         ModelGeneral            * modelGeneral,
                         ModelAVOStatic          * modelAVOStatic,
                         SeismicParametersHolder & seismicParameters,
                         InputFiles              * inputFiles,
                         int                       sortedVintage,
                         Simbox                  * timeBGSimbox);

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


#endif

