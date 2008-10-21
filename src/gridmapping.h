#ifndef GRIDMAPPING_H
#define GRIDMAPPING_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"
#include "src/modelsettings.h"
#include "lib/global_def.h"
#include "src/definitions.h"
#include "nrlib/stormgrid/stormcontgrid.hpp"

class Simbox;
class FFTGrid;
class ModelSettings;
class ModelFile;
class Model;

class GridMapping{
public:
  GridMapping(const Simbox  * simbox, 
              ModelSettings * modelSettings, 
              char          * errText,
              int             format, 
              FFTGrid       * velocity = NULL);
  ~GridMapping();
  StormContGrid * getMapping(void)  const { return mapping_ ;}
  Simbox        * getSimbox(void)   const { return simbox_  ;}
  int             getFormat(void)         { return format_  ;}

  void            setDepthSurfaces(ModelFile     * modelFile, 
                                   ModelSettings * modelSettings, 
                                   bool          & failed, 
                                   char          * errText,
                                   int             nz);
  void            calculateSurfaceFromVelocity(FFTGrid       * velocity, 
                                               const Simbox  * simbox, 
                                               ModelSettings * modelSettings, 
                                               bool          & failed, 
                                               char          * errText);
  void            makeTimeDepthMapping(FFTGrid      * velocity,
                                       const Simbox * timeSimbox);
  void            makeTimeTimeMapping(const Simbox * timeCutSimbox);

private:
  void            makeTimeDepthMapping(StormContGrid *& mapping,
                                       const Simbox * timeSimbox,
                                       const Simbox * depthSimbox,
                                       FFTGrid      * velocity);
  void            setSimbox(ModelSettings *modelSettings, bool &failed, char *errText, int nz);
 
  StormContGrid * mapping_;
  Simbox        * simbox_;
  Surface       * z0Grid_;
  Surface       * z1Grid_;
  int             format_;  // 3 = both, 2 = storm_binary, 1 = storm_ascii, 0 = nothing
  int             surfmissing_;
};
#endif
