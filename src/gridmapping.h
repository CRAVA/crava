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
class Model;
class ModelSettings;
class ModelFile;



class GridMapping{
public:
  GridMapping(const Simbox *simbox, ModelFile *modelFile, ModelSettings *modelSettings, bool depthmode, bool &failed, char *errText,bool format, FFTGrid *velocity = NULL);
  ~GridMapping();

 void makeMapping(FFTGrid *velocity = NULL, const Simbox *timeSimbox = NULL);
 void calculateSurfaceFromVelocity(FFTGrid *velocity, const Simbox *simbox, ModelSettings *modelSettings, bool &failed, char *errText);
 StormContGrid *getMapping()const {return mapping_;}
 Simbox *getSimbox()const {return simbox_;}
 bool getFormat() {return format_;}
 
 
private:
  void setSurfaces(ModelFile *modelFile, ModelSettings *modelSettings, bool &failed, char *errText,int surfmissing);
 void setSimbox(ModelSettings *modelSettings, bool &failed, char *errText, int nz);
 
  bool format_;  // 0 = storm_binary, 1 = storm_ascii
  StormContGrid *mapping_;
  Simbox        *simbox_;
  int surfmissing_;
  bool depthmode_;
  Surface *z0Grid_;
  Surface *z1Grid_;


};
#endif
