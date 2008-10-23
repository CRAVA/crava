#ifndef GRIDMAPPING_H
#define GRIDMAPPING_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"
#include "lib/global_def.h"
#include "src/definitions.h"
#include "nrlib/stormgrid/stormcontgrid.hpp"

class Simbox;
class FFTGrid;

class GridMapping{
public:
  GridMapping(int format);
  ~GridMapping(void);
  StormContGrid * getMapping(void)  const { return mapping_ ;}
  Simbox        * getSimbox(void)   const { return simbox_  ;}
  int             getFormat(void)         { return format_  ;}

  void            setDepthSurfaces(char ** surfFile, 
                                   bool  & failed, 
                                   char  * errText);
  void            calculateSurfaceFromVelocity(FFTGrid      * velocity, 
                                               const Simbox * simbox);
  void            setDepthSimbox(const Simbox * timeSimbox,
                                 int            nz,
                                 bool         & failed,
                                 char         * errText);
  void            makeTimeDepthMapping(FFTGrid      * velocity,
                                       const Simbox * timeSimbox);
  void            makeTimeTimeMapping(const Simbox * timeCutSimbox);

private:
  StormContGrid * mapping_;
  Simbox        * simbox_;
  Surface       * z0Grid_;
  Surface       * z1Grid_;
  int             format_;  // 3 = both, 2 = storm_binary, 1 = storm_ascii, 0 = nothing
};
#endif
