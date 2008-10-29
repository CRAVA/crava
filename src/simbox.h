#ifndef SIMBOX_H
#define SIMBOX_H

#include <string.h>
#include "src/definitions.h"
#include "nrlib/volume/volume.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/segy/segy.hpp"

class Simbox : public NRLib2::Volume 
{
public:
  Simbox(void);
  Simbox(double x0, double y0, Surface * z0, double lx, 
         double ly, double lz, double rot, double dx, double dy, double dz); //Assumes constant thickness.
  Simbox(const Simbox *simbox);
  ~Simbox();

  int            getIndex(double x, double y, double z)   const;
  int            getClosestZIndex(double x, double y, double z);
  void           getIndexes(double x, double y, int & xInd, int & yInd) const;
  void           getIndexes(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void           getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void           getZInterpolation(double x, double y, double z, 
                                   int & index1, int & index2, double & t) const;
  void           getCoord(int xInd, int yInd, int zInd, double &x, double &y, double &z) const;
  void           getXYCoord(int xInd, int yInd, double &x, double &y) const;
  int            getnx()                         const { return nx_ ;}
  int            getny()                         const { return ny_ ;}
  int            getnz()                         const { return nz_ ;}
  double         getdx()                         const { return dx_ ;}
  double         getdy()                         const { return dy_ ;}
  double         getdz()                         const { return dz_ ;} // Maximum dz (nz is constant).
  double         getlx()                         const { return(GetLX());}
  double         getly()                         const { return(GetLY());}
  double         getlz()                         const { return(GetLZ());} // Maximum thickness
  double         getx0()                         const { return(GetXMin());}
  double         gety0()                         const { return(GetYMin());}
  double         getAngle()                      const { return(GetAngle());}
  int            getIL0()                        const { return inLine0_ ;}
  int            getXL0()                        const { return crossLine0_ ;}
  int            getILStep()                     const { return ilStep_  ;}
  int            getXLStep()                     const { return xlStep_  ;}
  bool           getIsConstantThick()            const { return constThick_;}
  double         getMinRelThick()                const { return minRelThick_;} // Returns minimum relative thickness.
  double         getRelThick(int i, int j)       const;                        // Local relative thickness.
  double         getRelThick(double x, double y) const;                        // Local relative thickness.
  double         getAvgRelThick(void)            const;
  void           getMinMaxZ(double & minZ, double & maxZ) const;
  bool           getILxflag()                    const {return(ILxflag_);}
  int            isInside(double x, double y) const;
  int            insideRectangle(const SegyGeometry *  geometry) const;
  double         getTop(double x, double y) const;
  double         getBot(double x, double y) const;
  char         * getStormHeader(int cubetype, int nx, int ny, int nz, bool flat = false, bool ascii = false) const;
  void           writeTopBotGrids(std::string topname, std::string botname, int outputFormat);
  int            checkError(double lzLimit, char * errText);
  int            setArea(const SegyGeometry * geometry, char * errText);
  void           setSeisLines(int * lineParams);
  void           setDepth(Surface * zref, double zShift, double lz, double dz);
  void           setDepth(Surface * z0, Surface * z1, int nz);
  int            status() const {return(status_);}
  void           externalFailure() {status_ = EXTERNALERROR;}
  void           findIJFromILXL(int IL, int XL, int &i, int &j) const;
  enum           simboxstatus{BOXOK, INTERNALERROR, EXTERNALERROR, EMPTY, NOAREA, NODEPTH};

private:
  double         dx_, dy_, dz_;            // Working resolution.
  int            nx_, ny_, nz_;            // Number of cells in each direction.
  int            status_;                  // Since Simbox may be incomplete or with error
  double         cosrot_, sinrot_;         // Saving time in transformations.
  int            inLine0_, crossLine0_;
  int            ilStep_, xlStep_;
  std::string    topName_;
  std::string    botName_;
  bool           constThick_;
  double         minRelThick_;
  bool           ILxflag_;
};
#endif
