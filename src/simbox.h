/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SIMBOX_H
#define SIMBOX_H

#include <string.h>

#include "nrlib/volume/volume.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "src/definitions.h"
#include "src/io.h"

class Simbox : public NRLib::Volume
{
public:
  Simbox(void);

  Simbox(double                 x0,
         double                 y0,
         const Surface        & z0,
         double                 lx,
         double                 ly,
         double                 lz,
         double                 rot,
         double                 dx,
         double                 dy,
         double                 dz); //Assumes constant thickness.

  Simbox(const Simbox         * simbox,
         const std::string    & interval_name,
         int                    n_layers,
         const Surface        & top_surface,
         const Surface        & bot_surface,
         Surface              * single_corr_surface,
         std::string          & err_text,
         bool                 & failed);

  Simbox(const Simbox         * simbox,
         const std::string    & interval_name,
         int                    n_layers,
         const Surface        & top_surface,
         const Surface        & base_surface,
         std::string          & err_text,
         bool                 & failed,
         const Surface        * top_corr_surface,
         const Surface        * base_corr_surface);

  Simbox(const Simbox         *simbox);

  ~Simbox();

  // GET functions ------------------------------------------------------

  int            getIndex(double x, double y, double z)   const;
  int            getClosestZIndex(double x, double y, double z);
  void           getIndexes(double x, double y, int & xInd, int & yInd) const;
  void           getIndexes(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void           getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void           getZInterpolation(double x, double y, double z,
                                   int & index1, int & index2, double & t) const;
  void           getInterpolationIndexes(double x, double y, double z,
                                         double & xInd, double & yInd, double & zInd) const;
  void           getCoord(int xInd, int yInd, int zInd, double &x, double &y, double &z) const;
  void           getXYCoord(int xInd, int yInd, double &x, double &y) const;
  void           getZCoord(int zInd, double x, double y, double &z) const;
  int            getnx()                         const { return nx_                      ;}
  int            getny()                         const { return ny_                      ;}
  int            getnz()                         const { return nz_                      ;}
  double         getdx()                         const { return dx_                      ;}
  double         getdy()                         const { return dy_                      ;}
  double         getMinDz()                      const { return getdz()*getMinRelThick() ;} // Maximum dz (nz is constant).
  double         getMaxDz()                      const { return getdz()                  ;} // Maximum dz (nz is constant).
  double         getdz()                         const { return dz_                      ;} // Maximum dz (nz is constant).
  double         getdz(int i, int j)             const { return dz_*getRelThick(i,j)     ;}
  double         getlx()                         const { return GetLX()                  ;}
  double         getly()                         const { return GetLY()                  ;}
  double         getlz()                         const { return GetLZ()                  ;} // Maximum thickness
  double         getMaxLz()                      const { return GetLZ()                  ;} // Maximum thickness
  double         getx0()                         const { return GetXMin()                ;}
  double         gety0()                         const { return GetYMin()                ;}
  double         getAngle()                      const { return GetAngle()               ;}
  double         getIL0()                        const { return inLine0_                 ;}
  double         getXL0()                        const { return crossLine0_              ;}
  double         getILStepX()                    const { return ilStepX_                 ;}
  double         getILStepY()                    const { return ilStepY_                 ;}
  double         getXLStepX()                    const { return xlStepX_                 ;}
  double         getXLStepY()                    const { return xlStepY_                 ;}
  bool           getIsConstantThick()            const { return constThick_              ;}
  std::string    GetIntervalName()               const { return interval_name_           ;}
  double         getMinRelThick()                const { return minRelThick_             ;} // Returns minimum relative thickness.
  const Surface* GetTopErodedSurface()           const { return top_eroded_surface_      ;}
  const Surface* GetBaseErodedSurface()          const { return base_eroded_surface_     ;}

  double         getRelThick(int i, int j)       const;                                     // Local relative thickness.
  double         getRelThick(double x, double y) const;                                     // Local relative thickness.
  double         getAvgRelThick(void)            const;
  void           getMinMaxZ(double & minZ, double & maxZ) const;
  double         GetErodedTopZMin()              const { return GetTopErodedSurface()->Min()  ;}
  double         GetErodedTopZMax()              const { return GetTopErodedSurface()->Max()  ;}
  double         GetErodedBotZMin()              const { return GetBaseErodedSurface()->Min() ;}
  double         GetErodedBotZMax()              const { return GetBaseErodedSurface()->Max() ;}
  double         getTopZMin()                    const { return(GetTopZMin(nx_,ny_))          ;}
  double         getTopZMax()                    const { return(GetTopZMax(nx_,ny_))          ;}
  double         getBotZMin()                    const { return(GetBotZMin(nx_,ny_))          ;}
  double         getBotZMax()                    const { return(GetBotZMax(nx_,ny_))          ;}
  int            isInside(double x, double y) const;
  int            insideRectangle(const SegyGeometry *  geometry) const;
  double         getTop(int i, int j) const;
  double         getBot(int i, int j) const;
  double         getTop(double x, double y) const;
  double         getBot(double x, double y) const;
  double         GetTopErodedSurface(int i, int j) const;
  double         GetTopErodedSurface(double x, double y) const;
  double         GetBotErodedSurface(int i, int j) const;
  double         GetBotErodedSurface(double x, double y) const;
  std::string    getStormHeader(int cubetype, int nx, int ny, int nz, bool flat = false, bool ascii = false) const;

  // SET functions ------------------------------------------------------

  void           setTopBotName(const std::string & topname, const std::string & botname, int outputFormat);
  void           SetTopBotErodedSurfaces(const Surface * top_surf, const Surface base_surf);
  bool           setArea(const SegyGeometry * geometry, std::string & errText);
  void           setILXL(const SegyGeometry * geometry);
  void           setDepth(const Surface & zRef, double zShift, double lz, double dz, bool skipCheck = false);
  void           setDepth(const Surface & z0, const Surface & z1, int nz, bool skipCheck = false);
  void           setDepth(const NRLib::Surface<double>& top_surf, const NRLib::Surface<double>& bot_surf, int nz, bool skipCheck);

  // other public functions

  void           writeTopBotGrids(const std::string & topname, const std::string & botname, const std::string & subdir, int outputFormat);
  int            calculateDz(double lzLimit, std::string & errText);
  bool           isAligned(const SegyGeometry * geometry) const; //Checks if IL/XL form geometry maps nicely.
  int            status() const {return(status_);}
  void           externalFailure() {status_ = EXTERNALERROR;}
  void           getMinAndMaxXY(double &xmin, double &xmax, double &ymin, double &ymax) const;
  enum           simboxstatus{BOXOK, INTERNALERROR, EXTERNALERROR, EMPTY, NOAREA, NODEPTH};
  NRLib::Vector  FindPlane(const Surface * surf);
  Surface *      CreatePlaneSurface(const NRLib::Vector & planeParams,
                                    Surface             * templateSurf) const;

private:
  double         dx_, dy_, dz_;            // Working resolution.
  int            nx_, ny_, nz_;            // Number of cells in each direction.
  int            status_;                  // Since Simbox may be incomplete or with error
  double         cosrot_, sinrot_;         // Saving time in transformations.
  std::string    topName_;
  std::string    botName_;
  std::string    interval_name_;

  Surface      * top_eroded_surface_;      // Original top surface before padding (eroded if there are multiple intervals)
  Surface      * base_eroded_surface_;     // Original base surface before padding (eroded if there are multiple intervals)

  //Note: IL/XL information is just carried passively by this class.
  double         inLine0_, crossLine0_;    // XL, IL at origin, not necessarily int.
  double         ilStepX_, ilStepY_;       // Change in XL when moving along x and y
  double         xlStepX_, xlStepY_;       // Change in XL when moving along x and y

  bool           constThick_;
  double         minRelThick_;
};
#endif
