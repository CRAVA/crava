/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#include <assert.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "nrlib/volume/volume.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/definitions.h"

//
// Empty constructor
//
Simbox::Simbox(void)
  : Volume(),
  top_eroded_surface_(NULL),
  base_eroded_surface_(NULL)
{
  interval_name_ = "";
  status_      = EMPTY;
  topName_     = "";
  botName_     = "";
  constThick_  = true;
  minRelThick_ = 1.0;
  dz_          = 0;
  inLine0_     = -0.5;
  crossLine0_  = -0.5;
  xlStepX_     = 1;
  xlStepY_     = 0;
  ilStepX_     = 0;
  ilStepY_     = 1;
  grad_x_      = 0;
  grad_y_      = 0;
  lz_eroded_   = 0;
}

//
// Constructor where a simple geometric cube defines the simbox
//
Simbox::Simbox(double x0, double y0, const Surface & z0, double lx,
               double ly, double lz, double rot, double dx, double dy, double dz)
  :  Volume(),
  top_eroded_surface_(NULL),
  base_eroded_surface_(NULL)
{
  interval_name_ = "";
  status_      = BOXOK;
  topName_     = "";
  botName_     = "";
  SetDimensions(x0,y0,lx,ly);
  SetAngle(rot);

  Surface z1(z0);
  z1.Add(lz);
  SetSurfaces(z0,z1); //Automatically sets lz correct in this case.
  lz_eroded_   = lz;

  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  dz_          = dz;
  nx_          = int(0.5+lx/dx_);
  ny_          = int(0.5+ly/dy_);
  nz_          = int(0.5+lz/dz_);
  constThick_  = true;
  minRelThick_ = 1.0;
  inLine0_     = -0.5;
  crossLine0_  = -0.5;
  xlStepX_     =  cosrot_/dx_;
  xlStepY_     =  sinrot_/dx_;
  ilStepX_     = -sinrot_/dy_;
  ilStepY_     =  cosrot_/dy_;
  grad_x_      = 0;
  grad_y_      = 0;
}


//
// Constructor that copies all simbox data from simbox
//
Simbox::Simbox(const Simbox * simbox):
Volume(*simbox),
top_eroded_surface_(NULL),
base_eroded_surface_(NULL)
{

  interval_name_  = "";
  status_         = simbox->status_;
  cosrot_         = cos(GetAngle());
  sinrot_         = sin(GetAngle());
  dx_             = simbox->dx_;
  dy_             = simbox->dy_;
  dz_             = simbox->dz_;
  nx_             = simbox->nx_;
  ny_             = simbox->ny_;
  nz_             = simbox->nz_;
  inLine0_        = simbox->inLine0_;
  crossLine0_     = simbox->crossLine0_;
  ilStepX_        = simbox->ilStepX_;
  ilStepY_        = simbox->ilStepY_;
  xlStepX_        = simbox->xlStepX_;
  xlStepY_        = simbox->xlStepY_;
  lz_eroded_      = 0;
  topName_        = "";
  botName_        = "";
  constThick_     = simbox->constThick_;
  minRelThick_    = simbox->minRelThick_;
  topName_        = simbox->topName_;
  botName_        = simbox->botName_;
  grad_x_         = 0;
  grad_y_         = 0;

  std::string   s = "";
  this->setArea(simbox, static_cast<int>(simbox->getnx()),
                  static_cast<int>(simbox->getny()), s);
  CopyAllPadding(*simbox, simbox->getMaxLz(), s);
  SetErodedSurfaces(simbox->GetTopErodedSurface(), simbox->GetBaseErodedSurface(), true);

}

//
//
//

Simbox::Simbox(const Simbox   & simbox):
Volume(simbox),
top_eroded_surface_(NULL),
base_eroded_surface_(NULL)
{
  std::string s   = "";
  this->setArea(&simbox, static_cast<int>(simbox.getnx()),
                  static_cast<int>(simbox.getny()), s);
  this->CopyAllPadding(simbox, simbox.getMaxLz(), s);
  this->SetErodedSurfaces(simbox.GetTopErodedSurface(), simbox.GetBaseErodedSurface());
  interval_name_  = "";
  status_         = simbox.status_;
  cosrot_         = cos(GetAngle());
  sinrot_         = sin(GetAngle());
  dx_             = simbox.dx_;
  dy_             = simbox.dy_;
  dz_             = simbox.dz_;
  inLine0_        = simbox.inLine0_;
  crossLine0_     = simbox.crossLine0_;
  ilStepX_        = simbox.ilStepX_;
  ilStepY_        = simbox.ilStepY_;
  xlStepX_        = simbox.xlStepX_;
  xlStepY_        = simbox.xlStepY_;
  lz_eroded_      = 0;
  constThick_     = simbox.constThick_;
  minRelThick_    = simbox.minRelThick_;
  topName_        = simbox.topName_;
  botName_        = simbox.botName_;
  grad_x_         = 0;
  grad_y_         = 0;
}


//
// Constructor that copies the area from estimation_simbox and where
// top surface and base_surface both define the inversion interval
// and correlation surfaces
//
Simbox::Simbox(const Simbox             * estimation_simbox,
               const std::string        & interval_name,
               int                        n_layers,
               double                     lz_limit,
               const Surface            & top_surface,
               const Surface            & base_surface,
               std::string              & err_text,
               bool                     & failed)
: Volume(*estimation_simbox),
  top_eroded_surface_(NULL),
  base_eroded_surface_(NULL)
{
  interval_name_  = interval_name;
  status_         = NODEPTH;
  cosrot_         = cos(estimation_simbox->getAngle());
  sinrot_         = sin(estimation_simbox->getAngle());
  dx_             = estimation_simbox->getdx();
  dy_             = estimation_simbox->getdy();
  dz_             = -1;
  nx_             = estimation_simbox->getnx();
  ny_             = estimation_simbox->getny();
  nz_             = n_layers;
  nx_pad_         = nx_;
  ny_pad_         = ny_;
  nz_pad_         = nz_;
  x_pad_fac_      = (nx_pad_ / nx_);
  y_pad_fac_      = (ny_pad_ / ny_);
  z_pad_fac_      = (nz_pad_ / nz_);
  inLine0_        = estimation_simbox->getIL0();
  crossLine0_     = estimation_simbox->getXL0();
  ilStepX_        = estimation_simbox->getILStepX();
  ilStepY_        = estimation_simbox->getILStepY();
  xlStepX_        = estimation_simbox->getXLStepX();
  xlStepY_        = estimation_simbox->getXLStepY();
  constThick_     = estimation_simbox->constThick_;
  minRelThick_    = estimation_simbox->minRelThick_;
  lz_eroded_      = 0;
  topName_        = "";
  botName_        = "";
  grad_x_         = 0;
  grad_y_         = 0;
  setDepth(top_surface, base_surface, n_layers);
  this->calculateDz(lz_limit, err_text);
  SetErodedSurfaces(top_surface, base_surface);

  if (estimation_simbox->CheckSurface(top_surface) == false){
    err_text += "Error: Surface "+ topName_ +" does not cover volume.\n";
    failed = true;
  }
  if(estimation_simbox->CheckSurface(base_surface) == false){
    err_text += "Error: Surface"+ botName_ +"does not cover volume.\n";
    failed  = true;
  }
    // Set Eroded surfaces in the simbox
  if(!failed){
    SetErodedSurfaces(top_surface, base_surface);
    // Set surfaces in the volume
    // This should set the status to BOXOK
    setDepth(top_surface, base_surface, n_layers);
  }
  if(status_!=BOXOK){
    failed = true;
    err_text+="Simbox setup failed: there is a problem with the surfaces.";
  }
}


//
// Constructor for intervals with one correlation direction
//
Simbox::Simbox(const Simbox         * simbox,
               const std::string    & interval_name,
               int                    n_layers,
               double                 lz_limit,
               const Surface        & top_surface,
               const Surface        & bot_surface,
               Surface              * single_corr_surface,
               int                    other_output,
               int                    output_domain,
               int                    output_format,
               std::string          & err_text,
               bool                 & failed)
: Volume(*simbox)
{

  interval_name_  = interval_name;
  status_         = BOXOK;
  cosrot_         = cos(simbox->GetAngle());
  sinrot_         = sin(simbox->GetAngle());
  SetAngle (simbox->GetAngle());
  dx_             = simbox->getdx();
  dy_             = simbox->getdy();
  dz_             = -1;
  nx_             = simbox->getnx();
  ny_             = simbox->getny();
  nz_             = n_layers;
  nx_pad_         = nx_;
  ny_pad_         = ny_;
  nz_pad_         = nz_;
  x_pad_fac_      = (nx_pad_ / nx_);
  y_pad_fac_      = (ny_pad_ / ny_);
  z_pad_fac_      = (nz_pad_ / nz_);
  topName_        = "";
  botName_        = "";
  lz_eroded_      = 0;
  inLine0_        = simbox->getIL0();
  crossLine0_     = simbox->getXL0();
  ilStepX_        = simbox->getILStepX();
  ilStepY_        = simbox->getILStepY();
  xlStepX_        = simbox->getXLStepX();
  xlStepY_        = simbox->getXLStepY();
  setDepth(top_surface, bot_surface, n_layers);
  this->calculateDz(lz_limit, err_text);
  SetErodedSurfaces(top_surface, bot_surface);

  if (simbox->CheckSurface(*single_corr_surface) != true) {
    err_text += "Error: Correlation surface "+single_corr_surface->GetName()+" does not cover volume.\n";
    failed = true;
  }
  // Extend simbox as in ModelGeneral::SetupExtendedTimeSimbox

  NRLib::Vector corr_plane_parameters = FindPlane(single_corr_surface);
  Surface * mean_surface;
  if(single_corr_surface->GetNI() > 2)
  mean_surface = new Surface(*single_corr_surface);
  else {
    mean_surface = new Surface(dynamic_cast<const Surface &>(top_surface));
    if(mean_surface->GetNI() == 2) { //Extend corrSurf to cover other surfaces.
      double minX = mean_surface->GetXMin();
      double maxX = mean_surface->GetXMax();
      double minY = mean_surface->GetYMin();
      double maxY = mean_surface->GetYMax();
      if(minX > single_corr_surface->GetXMin())
        minX = single_corr_surface->GetXMin();
      if(maxX < single_corr_surface->GetXMax())
        maxX = single_corr_surface->GetXMax();
      if(minY > single_corr_surface->GetYMin())
        minY = single_corr_surface->GetYMin();
      if(maxY < single_corr_surface->GetYMax())
        maxY = single_corr_surface->GetYMax();
      single_corr_surface->SetDimensions(minX, minY, maxX-minX, maxY-minY);
    }
  }

  // initialize mean surface to 0
  for(size_t i=0;i<mean_surface->GetN();i++)
    (*mean_surface)(i) = 0;

  // Find mean of top and base surface
  mean_surface->AddNonConform(&GetTopSurface());
  mean_surface->AddNonConform(&GetBotSurface());
  mean_surface->Multiply(0.5);
  NRLib::Vector ref_plane_parameters = FindPlane(mean_surface);

  // tilt the mean plane with the correlation plane
  ref_plane_parameters -= corr_plane_parameters;
  grad_x_ = ref_plane_parameters(1);
  grad_y_ = ref_plane_parameters(2);

  // Create plane from parameters and add the original corr surface
  Surface * ref_plane = CreatePlaneSurface(ref_plane_parameters, mean_surface);
  ref_plane->AddNonConform(single_corr_surface);

  // Create new top surface
  Surface new_top_surface(*ref_plane);
  new_top_surface.SubtractNonConform(&(top_surface));
  double shift_top = new_top_surface.Max();
  shift_top *= -1.0;
  new_top_surface.Add(shift_top);
  new_top_surface.AddNonConform(&(simbox->GetTopSurface()));

  // Create new base surface
  Surface new_base_surface(*ref_plane);
  new_base_surface.SubtractNonConform(&(bot_surface));
  double shift_bot = new_base_surface.Min();
  shift_bot *= -1.0;
  double thick    = shift_bot-shift_top;
  double dz       = getdz();
  int    nz       = int(thick/dz);
  double residual = thick - nz*dz;
  if (residual > 0.0) {
    shift_bot += dz-residual;
    nz++;
  }
  if (nz != n_layers) {
    LogKit::LogFormatted(LogKit::High,"\nNumber of layers in interval "+ interval_name +" increased from %d", n_layers);
    LogKit::LogFormatted(LogKit::High," to %d because of the correlation direction.\n",nz);
  }

  new_base_surface.Add(shift_bot);
  new_base_surface.AddNonConform(&(bot_surface));

  setDepth(new_top_surface, new_base_surface, nz);

  if((other_output & IO::EXTRA_SURFACES) > 0 && (output_domain & IO::TIMEDOMAIN) > 0) {
    std::string top_surf_name  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Extended";
    std::string base_surf_name = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Extended";
    writeTopBotGrids(top_surf_name,
                     base_surf_name,
                     IO::PathToInversionResults(),
                     output_format);
  }

  delete mean_surface;
  delete ref_plane;
}

// Constructor with two correlation surfaces -------------------------
Simbox::Simbox(const Simbox         * simbox,
               const std::string    & interval_name,
               int                    n_layers,
               double                 lz_limit,
               const Surface        & top_surface,
               const Surface        & base_surface,
               const Surface        * top_corr_surface,
               const Surface        * base_corr_surface,
               int                    other_output,
               int                    output_domain,
               int                    output_format,
               std::string          & err_text,
               bool                 & failed)
: Volume(*simbox)
{
  interval_name_  = interval_name;
  status_         = BOXOK;
  cosrot_         = cos(simbox->GetAngle());
  sinrot_         = sin(simbox->GetAngle());
  dx_             = simbox->getdx();
  dy_             = simbox->getdy();
  dz_             = -1;
  nx_             = simbox->getnx();
  ny_             = simbox->getny();
  nz_             = n_layers;
  nx_pad_         = nx_;
  ny_pad_         = ny_;
  nz_pad_         = nz_;
  x_pad_fac_      = (nx_pad_ / nx_);
  y_pad_fac_      = (ny_pad_ / ny_);
  z_pad_fac_      = (nz_pad_ / nz_);
  lz_eroded_      = 0;
  inLine0_        = simbox->getIL0();
  crossLine0_     = simbox->getXL0();
  ilStepX_        = simbox->getILStepX();
  ilStepY_        = simbox->getILStepY();
  xlStepX_        = simbox->getXLStepX();
  xlStepY_        = simbox->getXLStepY();
  setDepth(top_surface, base_surface, n_layers);
  SetErodedSurfaces(top_surface, base_surface);
  this->calculateDz(lz_limit, err_text);

  //
  // Check that the two surfaces do not intersect
  //
  for (size_t i = 0; i < base_surface.GetNI(); i++){
    for (size_t j = 0; j < base_surface.GetNJ(); j++){
      double z_top_corr   = top_corr_surface->GetZ(base_surface.GetX(i), base_surface.GetY(j));
      double z_base_corr  = base_corr_surface->GetZ(base_surface.GetX(i), base_surface.GetY(j));
      if (z_top_corr > z_base_corr){
        err_text += "Error: The top correlation surface crosses the base correlation surface for interval "+ interval_name +".\n";
        failed = true;
      }
    }
  }

  //
  // Check that the correlation surfaces cover the inversion surfaces
  //
  if (simbox->CheckSurface(*top_corr_surface) == false) {
    err_text += "Error: Top correlation surface "+ top_corr_surface->GetName() +"does not cover volume.\n";
    failed = true;
  }
  if (simbox->CheckSurface(*base_corr_surface) == false){
    err_text += "Error: Base correlation surface "+ base_corr_surface->GetName() +"does not cover volume.\n";
    failed = true;
  }

  //
  // Should the corr surfaces have the same resolution?
  //

  if(!failed){
    Surface * mean_corr_surface;
    //NRLib::Vector corr_plane_parameters_top = FindPlane(top_corr_surface);
    //NRLib::Vector corr_plane_parameters_base = FindPlane(base_corr_surface);

    Surface * mean_surface;
    // Use the finest resolution
    double resolution_top   = top_corr_surface->GetDX()*top_corr_surface->GetDY();
    double resolution_base  = base_corr_surface->GetDX()*base_corr_surface->GetDY();
    if(resolution_top != resolution_base){
      if(resolution_top > resolution_base){ // base corr surface has higher resolution
        mean_surface            = new Surface(*base_corr_surface);
        mean_corr_surface       = new Surface(*base_corr_surface);
      }
      else{                                 // top corr surface has the highest resolution
        mean_surface            = new Surface(*top_corr_surface);
        mean_corr_surface       = new Surface(*top_corr_surface);
      }

    }
    else{
      mean_surface              = new Surface(*top_corr_surface);
      mean_corr_surface         = new Surface(*top_corr_surface);
    }

    // Initialize mean surface to 0
    for(size_t i = 0; i < mean_surface->GetN(); i++){
      (*mean_surface)(i)        = 0;
      (*mean_corr_surface)(i)   = 0;
    }

    // Find mean of top and base surface

    mean_surface->AddNonConform(&GetTopSurface());
    mean_surface->AddNonConform(&GetBotSurface());
    mean_surface->Multiply(0.5);
    NRLib::Vector ref_plane_parameters = FindPlane(mean_surface);

    // Find the mean of the top and base correlation surface
    mean_corr_surface->AddNonConform(top_corr_surface);
    mean_corr_surface->AddNonConform(base_corr_surface);
    mean_corr_surface->Multiply(0.5);
    //NRLib::Vector corr_plane_parameters = FindPlane(mean_corr_surface);


    // tilt the mean plane with the correlation plane

    NRLib::Vector corr_plane_parameters_mean = FindPlane(mean_corr_surface);
    ref_plane_parameters -= corr_plane_parameters_mean;
    grad_x_               = ref_plane_parameters(1);
    grad_y_               = ref_plane_parameters(2);

    // Create plane from parameters and add the original corr surfaces

    //Surface * ref_plane_base    = CreatePlaneSurface(ref_plane_parameters, mean_surface);
    Surface * ref_plane     = CreatePlaneSurface(ref_plane_parameters, mean_surface);

    // Create new top surface
    ref_plane->AddNonConform(top_corr_surface);
    ref_plane->AddNonConform(base_corr_surface);
    ref_plane->Multiply(0.5);
    Surface new_top(*ref_plane);
    new_top.SubtractNonConform(&(top_surface));
    double shift_top = new_top.Max();
    shift_top *= -1.0;
    new_top.Add(shift_top);
    new_top.AddNonConform(&(simbox->GetTopSurface()));

    // Create new base surface
    //ref_plane_base->AddNonConform(base_corr_surface);
    Surface new_base(*ref_plane);
    new_base.SubtractNonConform(&(base_surface));
    double shift_bot = new_base.Min();
    shift_bot *= -1.0;

    double thick    = shift_bot-shift_top;
    double dz       = getdz();
    int    nz       = int(thick/dz);
    double residual = thick - nz*dz;
    if (residual > 0.0) {
      shift_bot += dz-residual;
      nz++;
    }

    if (nz != n_layers) {
      LogKit::LogFormatted(LogKit::High,"\nNumber of layers in interval "+ interval_name +" increased from %d", n_layers);
      LogKit::LogFormatted(LogKit::High," to %d in grid created using the correlation direction.\n",nz);
    }

    new_base.Add(shift_bot);
    new_base.AddNonConform(&(base_surface));

    setDepth(new_top, new_base, nz);

    if((other_output & IO::EXTRA_SURFACES) > 0 && (output_domain & IO::TIMEDOMAIN) > 0) {
      std::string top_surf_name  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Extended";
      std::string base_surf_name = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Extended";
      writeTopBotGrids(top_surf_name,
                       base_surf_name,
                       IO::PathToInversionResults(),
                       output_format);
    }

    delete  ref_plane;
    //delete  ref_plane_base;
    delete  mean_surface;
    delete  mean_corr_surface;
  }
}

//
// Destructor
//
Simbox::~Simbox()
{
  delete top_eroded_surface_;
  delete base_eroded_surface_;
}

//
// Copy constructor
//
Simbox & Simbox::operator=(const Simbox   & rhs)
{
  if(this ==  &rhs) return *this;

  // Erik N: Should use copy and swap for exception safety - but not possible
  // for NRLib::Volume ?
  Simbox tmp(rhs);

  std::string   s = "";
  // sets top and base surfaces, nz_, lz_ dz_, nx_pad_, ny_pad_, nz_pad_ and nx_pad_factor_,
  // ny_pad_factor_, nz_pad_factor_
  this->CopyAllPadding(rhs, rhs.getMaxLz(), s);
  // sets   x_min_, y_min_,lx_ and ly_
  this->setArea(&rhs, static_cast<int>(rhs.getnx()),
                  static_cast<int>(rhs.getny()), s);
  // sets lz_eroded_, top and base eroded surfaces
  SetErodedSurfaces(rhs.GetTopErodedSurface(), rhs.GetBaseErodedSurface());

  std::swap(interval_name_, tmp.interval_name_);
  std::swap(status_, tmp.status_);
  this->SetAngle(tmp.GetAngle());
  this->SetTolerance(tmp.GetTolerance());
  std::swap(cosrot_, tmp.cosrot_);
  std::swap(sinrot_, tmp.sinrot_);
  std::swap(dx_, tmp.dx_);
  std::swap(dy_, tmp.dy_);
  std::swap(dz_, tmp.dz_);
  std::swap(inLine0_, tmp.inLine0_);
  std::swap(crossLine0_, tmp.crossLine0_);
  std::swap(ilStepX_, tmp.ilStepX_);
  std::swap(ilStepY_, tmp.ilStepY_);
  std::swap(xlStepX_, tmp.xlStepX_);
  std::swap(xlStepY_, tmp.xlStepY_);
  std::swap(constThick_, tmp.constThick_);
  std::swap(minRelThick_, tmp.minRelThick_);
  std::swap(topName_, tmp.topName_);
  std::swap(botName_, tmp.botName_);
  grad_x_         = 0;
  grad_y_         = 0;

  return *this;
}



//
// --------------------------------------------------------------------------
//

int Simbox::getIndex(double x, double y, double z) const

{
  int index = IMISSING;
  int i, j, k;
  getIndexes(x,y,z,i,j,k);
  if(k != IMISSING && j != IMISSING && i != IMISSING)
    index = int(i+j*nx_+k*nx_*ny_);
  return(index);
}

int
Simbox::getClosestZIndex(double x, double y, double z)
{
  int index = IMISSING;
  int i, j, k;
  getIndexesFull(x,y,z,i,j,k);
  if(i >=0 && i < nx_ && j >=0 && j < ny_)
  {
    if(k < 0)
      k = 0;
    else if(k >= nz_)
      k = nz_-1;
    index = i+j*nx_+k*nx_*ny_;
  }
  return(index);
}

void
Simbox::getIndexes(double x, double y, double z, int & xInd, int & yInd, int & zInd) const
{
  xInd = IMISSING;
  yInd = IMISSING;
  zInd = IMISSING;
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx >= 0 && rx <= GetLX() && ry >= 0 && ry <= GetLY())
  {
    double zBot, zTop = GetTopSurface().GetZ(x,y);
    if(GetTopSurface().IsMissing(zTop) == false)
    {
      zBot = GetBotSurface().GetZ(x,y);
      if(GetBotSurface().IsMissing(zBot) == false &&  z > zTop && z < zBot)
      {
        xInd = int(floor(rx/dx_));
        if(xInd > nx_-1)
          xInd = nx_-1;
        yInd = int(floor(ry/dy_));
        if(yInd > ny_-1)
          yInd = ny_-1;
        zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
        //LogKit::LogFormatted(LogKit::Low,"rx,dx,xInd = %.4f %.4f %d   ry,dy,yInd = %.4f %.4f %d    %d\n",rx,dx_,xInd,ry,dy_,yInd,zInd);
      }
    }
  }
}

void
Simbox::getIndexes(double x, double y, int & xInd, int & yInd) const
{
  xInd = IMISSING;
  yInd = IMISSING;
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx > 0 && rx < GetLX() && ry>0 && ry < GetLY())
  {
    xInd = static_cast<int>(floor(rx/dx_));
    yInd = static_cast<int>(floor(ry/dy_));
  }
}

void
Simbox::getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  xInd = int(floor(rx/dx_));
  yInd = int(floor(ry/dy_));
  zInd = IMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
      zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
  }
}

void
Simbox::getInterpolationIndexes(double x, double y, double z,
                                double & xInd, double & yInd, double & zInd) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  xInd = rx/dx_-0.5;
  yInd = ry/dy_-0.5;
  zInd = RMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
      zInd = static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)-0.5;
  }
}

void
Simbox::getZInterpolation(double x, double y, double z,
                          int & index1, int & index2, double & t) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  int xInd = int(floor(rx/dx_));
  int yInd = int(floor(ry/dy_));
  int zInd2, zInd1;
  index1 = IMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      zInd1 = static_cast<int>(floor((z-zTop)/dz)-0.5); //Find cell center above.
      if(zInd1 >=0 && zInd1 < nz_-1)
      {
        t = (z-zTop)/dz - 0.5 - static_cast<double>(zInd1);
        zInd2 = zInd1+1;
      }
      else
      {
        t = 0;
        if(zInd1 < 0)
          zInd1 = 0;
        else
          zInd1 = nz_-1;
        zInd2 = zInd1;
      }
      index1 = xInd+yInd*nx_+zInd1*nx_*ny_;
      index2 = xInd+yInd*nx_+zInd2*nx_*ny_;
    }
  }
}

bool Simbox::IsPointBetweenOriginalSurfaces(double x, double y, double z) const{
  const NRLib::Surface<double> * top_surf  = &GetTopErodedSurface();
  const NRLib::Surface<double> * base_surf = &GetBaseErodedSurface();
  bool b = false;
  if(isInside(x, y)){
    if (top_surf->GetZ(x,y) <= z && base_surf->GetZ(x,y) > z)
      b = true;
  }
  return b;
}

void
Simbox::getCoord(int xInd, int yInd, int zInd, double &x, double &y, double &z) const
{
  getXYCoord(xInd, yInd, x, y);
  getZCoord(zInd, x, y, z);
}

void
Simbox::getXYCoord(int xInd, int yInd, double &x, double &y) const
{
  double rx = (static_cast<double>(xInd) + 0.5)*dx_;
  double ry = (static_cast<double>(yInd) + 0.5)*dy_;
  x = rx*cosrot_-ry*sinrot_ + GetXMin();
  y = rx*sinrot_+ry*cosrot_ + GetYMin();
}

void
Simbox::getZCoord(int zInd, double x, double y, double &z) const
{
  z = RMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      z = zTop + (static_cast<double>(zInd) + 0.5)*dz;
    }
  }
}

void
Simbox::getMinMaxZ(double &minZ, double &maxZ) const
{
  minZ = GetZMin(nx_,ny_);
  maxZ = GetZMax(nx_,ny_);
}

int
Simbox::isInside(double x, double y) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx < 0 || rx > GetLX() || ry<0 || ry > GetLY())
    return(0);
  else
    return(1);
}

int
Simbox::insideRectangle(const SegyGeometry *  geometry) const
{
  double xr   = geometry->GetX0();
  double yr   = geometry->GetY0();
  double rotr = geometry->GetAngle();
  double lxr  = geometry->Getlx();
  double lyr  = geometry->Getly();
  double dxr  = geometry->GetDx();
  double dyr  = geometry->GetDy();

  // check that incoming rectangle is within simbox +-0.5 grid cells
  int allOk = 1;
  double cosrotr = cos(rotr);
  double sinrotr = sin(rotr);
  double x       = GetXMin();
  double y       = GetYMin();
  double rx      =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  double ry      = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;

  x  = GetXMin()+GetLX()*cosrot_;
  y  = GetYMin()+GetLX()*sinrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;

  x  = GetXMin()-GetLY()*sinrot_;
  y  = GetYMin()+GetLY()*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;

  x  = GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_;
  y  = GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;
  if(rotr<0)
    rotr+=2*NRLib::Pi;

  if (allOk==0) {
    double seisAzimuth = (-1)*rotr*(180/NRLib::Pi);
    double areaAzimuth = (-1)*GetAngle()*(180/NRLib::Pi);
    if (seisAzimuth < 0) seisAzimuth += 360.0;
    if (areaAzimuth < 0) areaAzimuth += 360.0;
    LogKit::LogFormatted(LogKit::Low,"                        x0            y0           lx         ly     azimuth         dx      dy\n");
    LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::Low,"Model area:    %11.2f  %11.2f    %11.2f %11.2f    %8.3f    %7.2f %7.2f\n",
                         GetXMin(), GetYMin(), GetLX(), GetLY(), dx_, dy_, areaAzimuth);
    LogKit::LogFormatted(LogKit::Low,"Seismic area:  %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                         xr, yr, lxr, lyr, dxr, dyr, seisAzimuth);

    LogKit::LogFormatted(LogKit::High,"\nCorner     XY Area                    XY Seismic\n");
    LogKit::LogFormatted(LogKit::High,"-----------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::High,"A %18.2f %11.2f    %11.2f %11.2f\n", GetXMin(),GetYMin(), xr,yr);
    LogKit::LogFormatted(LogKit::High,"B %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()+GetLX()*cosrot_, GetYMin()+GetLX()*sinrot_,
                         xr+lxr*cosrotr, yr+lxr*sinrotr);
    LogKit::LogFormatted(LogKit::High,"C %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()-GetLY()*sinrot_, GetYMin()+GetLY()*cosrot_,
                         xr -lyr*sinrotr, yr +lyr*cosrotr);
    LogKit::LogFormatted(LogKit::High,"D %18.2f %11.2f    %11.2f %11.2f\n",
                         GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_, GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_,
                         xr +lxr*cosrotr-lyr*sinrotr, yr +lxr*sinrotr+lyr*cosrotr);
    //
    // Calculate and write the largest possible AREA based on the (dx, dy, angle) given by user.
    //
    // Not implemented...
  }
  int error = 1 - allOk;
  return error;
}

double
Simbox::getTop(int i, int j) const
{
  double x, y;
  getXYCoord(i,j,x,y);
  double zTop = GetTopSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(zTop))
    zTop = RMISSING;
  return(zTop);
}

double
Simbox::getTop(double x, double y) const
{
  double zTop = GetTopSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(zTop))
    zTop = RMISSING;
  return(zTop);
}

double
Simbox::getBot(int i, int j) const
{
  double x, y;
  getXYCoord(i,j,x,y);
  double zBot = GetBotSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(zBot))
    zBot = RMISSING;
  return(zBot);
}

double
Simbox::getBot(double x, double y) const
{
  double zBot = GetBotSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(zBot))
    zBot = RMISSING;
  return(zBot);
}

double  Simbox::GetTopErodedSurface(int i, int j) const{
  double x, y;
  getXYCoord(i,j,x,y);
  double z_top = GetTopErodedSurface().GetZ(x,y);
  if(GetTopErodedSurface().IsMissing(z_top))
    z_top = RMISSING;
  return z_top;
}

double  Simbox::GetTopErodedSurface(double x, double y) const{
  double z_top = GetTopErodedSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(z_top))
    z_top = RMISSING;
  return(z_top);
}

double  Simbox::GetBotErodedSurface(int i, int j) const{
  double x, y;
  getXYCoord(i,j,x,y);
  double z_base = GetBaseErodedSurface().GetZ(x,y);
  if(GetTopErodedSurface().IsMissing(z_base))
    z_base = RMISSING;
  return z_base;
}

double  Simbox::GetBotErodedSurface(double x, double y) const{
  double z_bot = GetBaseErodedSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(z_bot))
    z_bot = RMISSING;
  return(z_bot);
}

std::string
Simbox::getStormHeader(int cubetype, int nx, int ny, int nz, bool flat, bool ascii) const
{
  if(flat == false)
    assert(topName_ != "");
  std::string header;
  if(ascii == false)
    header = "storm_petro_binary\n";
  else
    header = "storm_petro_ascii\n";

  header += "0 "+NRLib::ToString(cubetype) +" "+ NRLib::ToString(RMISSING,6)+"\n";
  header += "FFTGrid\n";
  if(flat == false)
    header += NRLib::ToString(GetXMin(),6) +" "+ NRLib::ToString(GetLX(),6) +" "+ NRLib::ToString(GetYMin(),6) +" "+ NRLib::ToString(GetLY(),6) +" "+ topName_ +" "+ botName_ +" 0.0 0.0\n";
  else
    header += NRLib::ToString(GetXMin(),6) +" "+ NRLib::ToString(GetLX(),6) +" "+ NRLib::ToString(GetYMin(),6) +" "+ NRLib::ToString(GetLY(),6) +" 0.0 "+ NRLib::ToString(GetLZ(),6)+" 0.0 0.0\n";

  header += NRLib::ToString(GetLZ(),6) +" "+ NRLib::ToString(GetAngle()*180/NRLib::Pi,6)+"\n\n";
  header += NRLib::ToString(nx) +" "+ NRLib::ToString(ny) +" "+ NRLib::ToString(nz)+"\n";
  std::string strHeader(header);

  /*
    ==>
    ==> Code to be used with g++ 4.3.2
    ==>
  std::string strHeader;
  if(ascii == false)
    strHeader += "storm_petro_binary\n";
  else
    strHeader += "storm_petro_ascii\n";

  strHeader += "0 " + NRLib::ToString(cubetype,6) + " " + NRLib::ToString(RMISSING,6) + "\n";
  strHeader += "FFTGrid\n";

  if(flat == false)
    strHeader += NRLib::ToString(GetXMin(),6) + " " + NRLib::ToString(GetLX(),6) + " "
      + NRLib::ToString(GetYMin(),6) + " " + NRLib::ToString(GetLY(),6) + " "
               + topName_ + " "
      + botName_ + " 0.0 0.0\n";
  else
    strHeader += NRLib::ToString(GetXMin(),6) + " " + NRLib::ToString(GetLX(),6) + " "
      + NRLib::ToString(GetYMin(),6) + " " + NRLib::ToString(GetLY(),6) + " "
               + "0.0 "
      + NRLib::ToString(GetLZ(),6)+" 0.0 0.0\n";

  strHeader += NRLib::ToString(GetLZ(),6) + " " + NRLib::ToString(GetAngle()*180/PI,6) + "\n\n";
  strHeader += NRLib::ToString(nx) + " " + NRLib::ToString(ny) + " " + NRLib::ToString(nz) + "\n";
  */

  return(strHeader);
}

void
Simbox::writeTopBotGrids(const std::string & topname,
                         const std::string & botname,
                         const std::string & subdir,
                         int                 outputFormat)
{
  assert(typeid(GetTopSurface()) == typeid(Surface));
  assert(typeid(GetBotSurface()) == typeid(Surface));

  const Surface & wtsurf = dynamic_cast<const Surface &>(GetTopSurface());
  const Surface & wbsurf = dynamic_cast<const Surface &>(GetBotSurface());

  IO::writeSurfaceToFile(wtsurf, topname, subdir, outputFormat);
  IO::writeSurfaceToFile(wbsurf, botname, subdir, outputFormat);
}

void
Simbox::setTopBotName(const std::string & topname,
                      const std::string & botname,
                      int                 outputFormat)
{
  std::string suffix;
  if ((outputFormat & IO::ASCII) > 0 && (outputFormat & IO::STORM) == 0)
    suffix = IO::SuffixAsciiIrapClassic();
  else
    suffix = IO::SuffixStormBinary();

  topName_ = IO::getFilePrefix()+topname+suffix;
  botName_ = IO::getFilePrefix()+botname+suffix;
}

int
Simbox::calculateDz(double lzLimit, std::string & errText)
{
  if(status_ == NODEPTH || status_ == EMPTY)
    status_ = EXTERNALERROR; //At this stage, lack of depth is an error

  if(status_ == EXTERNALERROR || status_ == INTERNALERROR)
    //Earlier internal errors are external for this purpose.
    return(EXTERNALERROR);

  if(status_ == NOAREA)
    return(BOXOK);

  if(dz_ < 0)
  {
    double z0, z1 = 0.0;
    double x, y, rx, ry = 0.5f*dy_;
    double lzCur, lzMin = double(1e+30);
    int i,j;
    for(j=0;j<ny_;j++)
    {
      rx = 0.5f*dx_;
      for(i=0;i<nx_;i++)
      {
        x = rx*cosrot_-ry*sinrot_ + GetXMin();
        y = rx*sinrot_+ry*cosrot_ + GetYMin();
        z0 = GetTopSurface().GetZ(x,y);
        z1 = GetBotSurface().GetZ(x,y);
        if(GetTopSurface().IsMissing(z0) == false && GetBotSurface().IsMissing(z1) == false )
        {
          lzCur = z1 - z0;
          if(lzCur < lzMin)
            lzMin = lzCur;
        }
        rx += dx_;
      }
      ry += dy_;
    }

    if(lzMin < 0.0)
    {
      status_ = INTERNALERROR;
      errText += "-At least parts of the top surface is lower than the base surface. Are surfaces given in wrong order?\n";
    }
    else
    {
      double lzFac = lzMin/GetLZ();
      minRelThick_ = lzFac;
      if(lzFac < lzLimit)
      {
        status_ = INTERNALERROR;
        errText += "-Error with top/bottom grids in interval "+interval_name_+". Minimum thickness should be at least "+NRLib::ToString(lzLimit)+" times maximum, is "+NRLib::ToString(lzFac)+"\n";
      }
      else
      {
        dz_ = GetLZ()/static_cast<double>(nz_);
      }
    }
  }
  return(status_);
}


bool
Simbox::setArea(const SegyGeometry * geometry, std::string & errText)
{
  double x0  = geometry->GetX0();
  double y0  = geometry->GetY0();
  double lx  = geometry->Getlx();
  double ly  = geometry->Getly();
  double rot = geometry->GetAngle();
  double dx  = geometry->GetDx();
  double dy  = geometry->GetDy();

  bool failed = false;

  try
  {
    SetDimensions(x0,y0,lx,ly);
  }
  catch (NRLib::Exception & e)
  {
    errText += "Could not set x0, y0, lx, and ly.\n";
    errText += e.what();
    return true; // Failed
  }
  try
  {
    SetAngle(rot);
  }
  catch (NRLib::Exception & e)
  {
    errText += "Could not set rotation angle.\n";
    errText += e.what();
    failed = true;
    return true; // Failed
  }
  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  nx_          = static_cast<int>(0.5+lx/dx_);
  ny_          = static_cast<int>(0.5+ly/dy_);

  // In case IL/XL information is not available, we fall back
  //  on the following base case values ...
  inLine0_     = -0.5;
  crossLine0_  = -0.5;
  ilStepX_     =  cosrot_/dx_;
  ilStepY_     =  sinrot_/dx_;
  xlStepX_     = -sinrot_/dy_;
  xlStepY_     =  cosrot_/dy_;

  if(status_ == EMPTY)
    status_ = NODEPTH;
  else if(status_ == NOAREA)
    status_ = BOXOK;

  return false; // OK
}

bool
Simbox::setArea(const NRLib::Volume * volume, int nx, int ny, std::string & errText, bool scale)
{
  double scale_value = 1.0;

  if (scale == true) //SGRI
    scale_value = 1000;

  double x0  = volume->GetXMin()*scale_value;
  double y0  = volume->GetYMin()*scale_value;
  double lx  = volume->GetLX()*scale_value;
  double ly  = volume->GetLY()*scale_value;
  double rot = volume->GetAngle();
  double dx  = lx/static_cast<double>(nx);
  double dy  = ly/static_cast<double>(ny);

  bool failed = false;

  try
  {
    SetDimensions(x0,y0,lx,ly);
  }
  catch (NRLib::Exception & e)
  {
    errText += "Could not set x0, y0, lx, and ly.\n";
    errText += e.what();
    return true; // Failed
  }
  try
  {
    SetAngle(rot);
  }
  catch (NRLib::Exception & e)
  {
    errText += "Could not set rotation angle.\n";
    errText += e.what();
    failed = true;
    return true; // Failed
  }
  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  nx_          = static_cast<int>(0.5+lx/dx_);
  ny_          = static_cast<int>(0.5+ly/dy_);

  // In case IL/XL information is not available, we fall back
  //  on the following base case values ...
  inLine0_     = -0.5;
  crossLine0_  = -0.5;
  ilStepX_     =  cosrot_/dx_;
  ilStepY_     =  sinrot_/dx_;
  xlStepX_     = -sinrot_/dy_;
  xlStepY_     =  cosrot_/dy_;

  if(status_ == EMPTY)
    status_ = NODEPTH;
  else if(status_ == NOAREA)
    status_ = BOXOK;

  return false; // OK
}


void
Simbox::setDepth(const Surface & zRef, double zShift, double lz, double dz, bool skipCheck)
{
  Surface zTop(zRef);
  zTop.Add(zShift);
  Surface zBot(zTop);
  zBot.Add(lz);
  SetSurfaces(zTop,zBot,skipCheck);
  dz_ = dz;
  nz_ = int(0.5+lz/dz_);
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;
}

void Simbox::setDepth(const NRLib::Surface<double>& top_surf,
                      const NRLib::Surface<double>& bot_surf, int nz, bool skipCheck)
{
  SetSurfaces(top_surf, bot_surf, skipCheck);
  nz_ = nz;
  dz_ = -1;
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;

  constThick_ = false;
}

void
Simbox::setDepth(const Surface & z0, const Surface & z1, int nz, bool skipCheck)
{
  SetSurfaces(z0, z1, skipCheck);
  nz_ = nz;
  nz_pad_ = nz;
  //dz_ = -1;
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;

  constThick_ = false;
}

void
Simbox::setILXL(const SegyGeometry * geometry)
{
  xlStepX_ = geometry->GetXLStepX();
  xlStepY_ = geometry->GetXLStepY();
  ilStepX_ = geometry->GetILStepX();
  ilStepY_ = geometry->GetILStepY();

  float x0 = static_cast<float>(GetXMin());
  float y0 = static_cast<float>(GetYMin());
  geometry->FindContILXL(x0, y0, inLine0_, crossLine0_); //Sets IL0 ,XL0
}

bool
Simbox::isAligned(const SegyGeometry * geometry) const
{
  double x,y;
  getXYCoord(0, 0, x, y);
  int IL0, XL0;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL0, XL0);
  getXYCoord(1, 0, x, y);
  int ILx1, XLx1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), ILx1, XLx1);
  getXYCoord(0, 1, x, y);
  int ILy1, XLy1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), ILy1, XLy1);
  getXYCoord(nx_-1, 0, x, y);
  int IL1, XL1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL1, XL1);
  getXYCoord(0, ny_-1, x, y);
  int IL2, XL2;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL2, XL2);
  getXYCoord(nx_-1, ny_-1, x, y);
  int IL3, XL3;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL3, XL3);

  int XLdx = XLx1 - XL0;
  int XLdy = XLy1 - XL0;
  int ILdx = ILx1 - IL0;
  int ILdy = ILy1 - IL0;
  if(abs(XLdx*XLdy) > 0 || abs(ILdx*ILdy) > 0)
    return(false); //Moving along one axis lead to change in both il and xl

  int XLndx = XL1 - XL0;
  int XLndy = XL2 - XL0;
  int ILndx = IL1 - IL0;
  int ILndy = IL2 - IL0;
  if(XLndx != (nx_-1)*XLdx || XLndy != (ny_-1)*XLdy ||
     ILndx != (nx_-1)*ILdx || ILndy != (ny_-1)*ILdy)
    return(false); //XL or IL difference at corners not multiple of one-step difference

  if(XL3-XL0 != (nx_-1)*XLdx+(ny_-1)*XLdy || IL3-IL0 != (nx_-1)*ILdx+(ny_-1)*ILdy)
    return(false); //Check final corner, changes failed to match one-step.

  return(true);
}

double
Simbox::getAvgRelThick(void) const
{
  double avgThick = 0.0f;
  for (int i = 0 ; i < nx_ ; i++) {
    for (int j = 0 ; j < ny_ ; j++) {
      avgThick += getRelThick(i, j);
    }
  }
  avgThick /= nx_*ny_;
  return avgThick;
}

double
Simbox::getRelThick(int i, int j) const
{
  double rx = (static_cast<double>(i) + 0.5)*dx_;
  double ry = (static_cast<double>(j) + 0.5)*dy_;
  double x = rx*cosrot_-ry*sinrot_ + GetXMin();
  double y = rx*sinrot_+ry*cosrot_ + GetYMin();
  return(getRelThick(x, y));
}

double Simbox::GetRelThickErodedInterval(double x, double y) const{
  double rel_thick  = 1; //Default value to be used outside grid.
  double z_top      = GetTopErodedSurface().GetZ(x, y);
  double z_base     = GetBaseErodedSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(z_top) == false &&
     GetBotSurface().IsMissing(z_base) == false)
    rel_thick = (z_base-z_top)/GetLZ();
  return(rel_thick);
}

double
Simbox::getRelThick(double x, double y) const
{
  double relThick = 1; //Default value to be used outside grid.
  double zTop = GetTopSurface().GetZ(x,y);
  double zBot = GetBotSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false &&
     GetBotSurface().IsMissing(zBot) == false)
    relThick = (zBot-zTop)/GetLZ();
  return(relThick);
}

void Simbox::getMinAndMaxXY(double &xmin, double &xmax, double &ymin, double &ymax)const
{
  xmin = std::min(GetXMin()+GetLX()*cosrot_, GetXMin());
  xmin = std::min(xmin,GetXMin()-GetLY()*sinrot_);
  xmin = std::min(xmin,GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_);

  xmax = std::max(GetXMin()+GetLX()*cosrot_, GetXMin());
  xmax = std::max(xmax,GetXMin()-GetLY()*sinrot_);
  xmax = std::max(xmax,GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_);

  ymin = std::min(GetYMin()+GetLX()*sinrot_, GetYMin());
  ymin = std::min(ymin,GetYMin()+GetLY()*cosrot_);
  ymin = std::min(ymin,GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_);

  ymax = std::max(GetYMin(),GetYMin()+GetLX()*sinrot_);
  ymax = std::max(ymax,GetYMin()+GetLY()*cosrot_);
  ymax = std::max(ymax,GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_);
}

NRLib::Vector Simbox::FindPlane(const Surface * surf){

  NRLib::SymmetricMatrix A = NRLib::SymmetricZeroMatrix(3);
  NRLib::Vector b(3);
  NRLib::Vector x(3);

  b = 0;

  int nData = 0;

  for(int i=0 ; i<static_cast<int>(surf->GetN()) ; i++) {
    double x, y, z;
    surf->GetXY(i, x, y);
    z = (*surf)(i);
    if(!surf->IsMissing(z)) {
      nData++;
      A(0,1) += x;
      A(0,2) += y;
      A(1,1) += x*x;
      A(1,2) += x*y;
      A(2,2) += y*y;
      b(0)   += z;
      b(1)   += x*z;
      b(2)   += y*z;
    }
  }

  A(0,0) = nData;

  NRLib::CholeskySolve(A, b, x);

  return x;
}

Surface * Simbox::CreatePlaneSurface(const NRLib::Vector & planeParams,
                                     Surface             * templateSurf) const{

  Surface * result = new Surface(*templateSurf);
  for(int i=0;i<static_cast<int>(result->GetN());i++) {
    double x,y;
    result->GetXY(i,x,y);
    (*result)(i) = planeParams(0)+planeParams(1)*x+planeParams(2)*y;
  }
  return(result);
}

void Simbox::SetErodedSurfaces(const NRLib::Surface<double> & top_surf,
                               const NRLib::Surface<double> & bot_surf,
                               bool  skip_check)
{
  if (top_eroded_surface_ != NULL)
    delete top_eroded_surface_;
  if (base_eroded_surface_ != NULL)
    delete base_eroded_surface_;
  top_eroded_surface_  = top_surf.Clone();
  base_eroded_surface_ = bot_surf.Clone();

  if (getlx() > 0 && getly() > 0 && skip_check == false)
    CheckErodedSurfaces();

  lz_eroded_ = RecalculateErodedLZ();
}

bool Simbox::CheckErodedSurfaces() const
{
  if (!CheckSurface(*top_eroded_surface_)) {
    throw NRLib::Exception("The top surface does not cover the volume.");
  }
  if (!CheckSurface(*base_eroded_surface_)) {
    throw NRLib::Exception("The bottom surface does not cover the volume.");
  }
  return true;
}

double Simbox::RecalculateErodedLZ() const
{
  double lz = 0;
  if (this->getlx() > 0.0 || this->getly() > 0.0) { //Only do if area is initialized.
    // An arbitary grid resolution.
    int nx = 100;
    int ny = 100;

    double dx = getlx() / nx;
    double dy = getly() / ny;


    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x, y;
        LocalToGlobalCoord(dx * i, dy * j, x, y);
        lz = std::max(lz_eroded_, base_eroded_surface_->GetZ(x, y) - top_eroded_surface_->GetZ(x, y));
      }
    }
  }
  return lz;
}

void
Simbox::CopyAllPadding(const Simbox & original,
                       double         lz_limit,
                       std::string  & err_txt)
{
  setDepth(original.GetTopSurface(), original.GetBotSurface(), original.getnz(), true);
  calculateDz(lz_limit, err_txt);
  SetNXpad(original.GetNXpad());
  SetNYpad(original.GetNYpad());
  SetNZpad(original.GetNZpad());
  SetXPadFactor(original.GetXPadFactor());
  SetYPadFactor(original.GetYPadFactor());
  SetZPadFactor(original.GetZPadFactor());
}

void
Simbox::SetNoPadding()
{
  SetNXpad(nx_);
  SetNYpad(ny_);
  SetNZpad(nz_);
  SetXPadFactor(0);
  SetYPadFactor(0);
  SetZPadFactor(0);
}

