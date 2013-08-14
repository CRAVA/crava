/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef INTERVALSIMBOX_H
#define INTERVALSIMBOX_H

#include "nrlib/volume/volume.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/flens/nrlib_flens.hpp"
#include "src/simbox.h"
#include "src/io.h"

class IntervalSimbox : public NRLib::Volume
{
public:
  IntervalSimbox(void);

  // Constructor with single correlation surface
  IntervalSimbox(const Simbox         * simbox,
                 const std::string    & interval_name,
                 int                    n_layers,
                 const Surface        & top_surface,
                 const Surface        & bot_surface,
                 Surface              * single_corr_surface,
                 std::string          & err_text,
                 bool                 & failed);

  // Constructor with two correlation surfaces
  IntervalSimbox(const Simbox         * simbox,
                 const std::string    & interval_name,
                 int                    nz,
                 const Surface        & top_surface,
                 const Surface        & bot_surface,
                 std::string          & err_text,
                 bool                 & failed,
                 const Surface        * top_corr_surface = NULL,
                 const Surface        * base_corr_surface = NULL);

  IntervalSimbox(double             x0,
                 double             y0,
                 const Surface    & z0,
                 double             lx,
                 double             ly,
                 double             lz,
                 double             rot,
                 double             dx,
                 double             dy,
                 double             dz); //Assumes constant thickness.

  IntervalSimbox(const Simbox   * simbox);

  ~IntervalSimbox();

  // GET FUNCTIONS

  int            GetNx()                         const { return nx_                      ;}
  int            GetNy()                         const { return ny_                      ;}
  int            GetNz()                         const { return nz_                      ;}
  double         GetDx()                         const { return dx_                      ;}
  double         GetDy()                         const { return dy_                      ;}
  double         GetDz()                         const { return dz_                      ;} // Maximum dz (nz is constant).
  double         GetDz(int i, int j)             const { return dz_*GetRelThick(i,j)     ;}
  double         GetMinDz()                      const { return dz_*GetMinRelThick()     ;} // Maximum dz (nz is constant).
  double         GetMaxDz()                      const { return dz_                      ;} // Maximum dz (nz is constant).
  double         GetMaxLz()                      const { return GetLZ()                  ;} // Maximum thickness
  double         GetX0()                         const { return GetXMin()                ;}
  double         GetY0()                         const { return GetYMin()                ;}
  double         GetIL0()                        const { return inline0_                 ;}
  double         GetXL0()                        const { return crossline0_              ;}
  double         GetILStepX()                    const { return il_step_X_               ;}
  double         GetILStepY()                    const { return il_step_Y_               ;}
  double         GetXLStepX()                    const { return xl_step_X_               ;}
  double         GetXLStepY()                    const { return xl_step_Y_               ;}
  bool           GetIsConstantThick()            const { return const_thick_             ;}
  double         GetMinRelThick()                const { return min_rel_thick_           ;} // Returns minimum relative thickness.
  double         GetRelThick(int i, int j)       const;                                     // Local relative thickness.
  double         GetRelThick(double x, double y) const;                                     // Local relative thickness.
  double         GetAvgRelThick(void)            const;
  int            GetStatus()                     const { return(status_)                 ;}
  std::string    GetIntervalName()               const { return interval_name_           ;}
  const Surface & GetBotCorrSurface()            const { return *bot_correlation_surface_;}
  const Surface & GetTopCorrSurface()            const { return *top_correlation_surface_;}

  int            GetIndex(double          x,
                          double          y,
                          double          z)   const;

  int            GetClosestZIndex(double  x,
                                  double  y,
                                  double  z);

  void           GetIndexes(double        x,
                            double        y,
                            int         & x_ind,
                            int         & y_ind) const;

  void           GetIndexes(double        x,
                            double        y,
                            double        z,
                            int         & x_ind,
                            int         & y_ind,
                            int         & z_ind) const;

  void           GetIndexesFull(double    x,
                                double    y,
                                double    z,
                                int     & x_ind,
                                int     & y_ind,
                                int     & z_ind) const;

  void           GetZInterpolation(double       x,
                                   double       y,
                                   double       z,
                                   int        & index_1,
                                   int        & index_2,
                                   double     & t) const;

  void           GetInterpolationIndexes(double         x,
                                         double         y,
                                         double         z,
                                         double       & x_ind,
                                         double       & y_ind,
                                         double       & z_ind) const;

  void           GetCoord(int             x_ind,
                          int             y_ind,
                          int             z_ind,
                          double        & x,
                          double        & y,
                          double        & z) const;

  void           GetXYCoord(int           x_ind,
                            int           y_ind,
                            double      & x,
                            double      & y) const;

  void           GetZCoord(int            z_ind,
                           double         x,
                           double         y,
                           double       & z) const;

  void           GetMinMaxZ(double      & min_z,
                            double      & max_z) const;

  double         GetTop(int               i,
                        int               j) const;

  double         GetBot(int               i,
                        int               j) const;

  double         GetTop(double            x,
                        double            y) const;

  double         GetBot(double            x,
                        double            y) const;

  std::string    GetStormHeader(int       cube_type,
                                int       nx,
                                int       ny,
                                int       nz,
                                bool      flat = false,
                                bool      ascii = false) const;

  void           GetMinAndMaxXY(double  & x_min,
                                double  & x_max,
                                double  & y_min,
                                double  & y_max) const;

  // SET FUNCTIONS

  void           SetCorrelationSurfaces(Surface     * top_corr_surface,
                                        Surface     * bot_corr_surface);

  void           SetTopCorrelationSurface(Surface   * top_corr_surface){ top_correlation_surface_ = top_corr_surface;}

  void           SetBotCorrelationSurface(Surface   * bot_corr_surface){ bot_correlation_surface_ = bot_corr_surface;}

  bool           SetArea(const SegyGeometry   * geometry,
                         std::string          & err_text);

  void           SetILXL(const SegyGeometry   * geometry);

  void           SetDepth(const Surface       & z_ref,
                          double                z_shift,
                          double                lz,
                          double                dz,
                          bool                  skip_check = false);

  void           SetDepth(const Surface       & z0,
                          const Surface       & z1,
                          int                   nz,
                          bool                  skip_check = false);

  void           SetTopBotName(const std::string        & top_name,
                               const std::string        & bot_name,
                               int                        output_format);

  // OTHER FUNCTIONS

  Surface *       CreatePlaneSurface(const NRLib::Vector & planeParams,
                                     Surface             * templateSurf);

  NRLib::Vector  FindPlane(const Surface * surf);

  int            IsInside(double                x,
                          double                y) const;

  int            InsideRectangle(const SegyGeometry     * geometry) const;

  void           WriteTopBotGrids(const std::string     & top_name,
                                  const std::string     & bot_name,
                                  const std::string     & subdir,
                                  int                     output_format);

  int            CalculateDz(double                       lz_limit,
                             std::string                & errText);

  bool           IsAligned(const SegyGeometry           * geometry) const; //Checks if IL/XL form geometry maps nicely.

  void           ExternalFailure() {status_ = EXTERNALERROR;}

  // VARIABLES

  enum           Simboxstatus{BOXOK, INTERNALERROR, EXTERNALERROR, EMPTY, NOAREA, NODEPTH};

private:

  double         dx_, dy_, dz_;            // Working resolution.
  int            nx_, ny_, nz_;            // Number of cells in each direction.
  int            status_;                  // Since Simbox may be incomplete or with error
  double         cosrot_, sinrot_;                  // Saving time in transformations.
  std::string    top_name_;                // Top surface name
  std::string    bot_name_;                // Botton surface name
  std::string    interval_name_;           // Interval name

  //Note: IL/XL information is carried passively by this class.
  double         inline0_, crossline0_;    // XL, IL at origin, not necessarily int.
  double         il_step_X_, il_step_Y_;       // Change in XL when moving along x and y
  double         xl_step_X_, xl_step_Y_;       // Change in XL when moving along x and y

  Surface     * top_correlation_surface_;
  Surface     * bot_correlation_surface_;
  Surface     * single_correlation_surface_;

  bool           const_thick_;
  double         min_rel_thick_;

};

#endif
