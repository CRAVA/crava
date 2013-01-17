/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef GRIDMAPPING_H
#define GRIDMAPPING_H

#include <stdio.h>

#include "src/definitions.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"

class Simbox;
class FFTGrid;

class GridMapping
{
public:
  GridMapping();
  ~GridMapping(void);
  StormContGrid * getMapping(void)  const { return mapping_ ;}
  Simbox        * getSimbox(void)   const { return simbox_  ;}
  void            setDepthSurfaces(const std::vector<std::string> & surfFile,
                                   bool                           & failed,
                                   std::string                    & errText);
  void            calculateSurfaceFromVelocity(FFTGrid      * velocity,
                                               const Simbox * simbox);
  void            setDepthSimbox(const Simbox * timeSimbox,
                                 int            nz,
                                 int            outputFormat,
                                 bool         & failed,
                                 std::string  & errText);
  void            makeTimeDepthMapping(FFTGrid      * velocity,
                                       const Simbox * timeSimbox);
  void            makeTimeTimeMapping(const Simbox * timeCutSimbox);

  void            setMappingFromVelocity(FFTGrid * velocity, const Simbox * timeSimbox);

  //Please do not renumber the modes below. It is very convenient that TOPGIVEN+BOTTOMGIVEN = BOTHGIVEN.
  enum            surfaceModes{NONEGIVEN = 0, TOPGIVEN = 1, BOTTOMGIVEN = 2, BOTHGIVEN = 3};


private:
  StormContGrid * mapping_;
  Simbox        * simbox_;

  Surface       * z0Grid_;
  Surface       * z1Grid_;

  int             surfaceMode_;
};
#endif
