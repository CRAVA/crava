#ifndef DOINVERSION_H
#define DOINVERSION_H

#include <stdio.h>

class ModelSettings;
class ModelAVODynamic;
class ModelAVOStatic;
class ModelGeneral;
class InputFiles;
class Simbox;

void setupStaticModels(ModelGeneral   *& modelGeneral,
                       ModelAVOStatic *& modelAVOstatic,
                       ModelSettings   * modelSettings,
                       InputFiles      * inputFiles,
                       Simbox         *& timeBGSimbox);

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


#endif

