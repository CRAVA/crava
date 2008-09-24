// $Id$

#include "volume.hpp"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"
#include "../surface/regularsurface.hpp"
#include "../surface/surface.hpp"
#include "../surface/surfaceio.hpp"

using namespace NRLib2;

Volume::Volume()
{
  x_min_       = 0.0;
  y_min_       = 0.0;
  lx_          = 0.0;
  ly_          = 0.0;
  lz_          = 0.0;
  angle_       = 0.0;
  z_top_       = new ConstantSurface(0.0);
  z_bot_       = new ConstantSurface(0.0);
  erosion_top_ = 0;
  erosion_bot_ = 0;
}


Volume::Volume(const Volume & volume)
{
  x_min_ = volume.x_min_;
  y_min_ = volume.y_min_;
  lx_    = volume.lx_;
  ly_    = volume.ly_;
  lz_    = volume.lz_;
  angle_ = volume.angle_;
  if (volume.z_top_ != 0) {
    z_top_ = volume.z_top_->Clone();
  }
  else {
    z_top_ = 0;
  }
  if (volume.z_bot_ != 0) {
    z_bot_ = volume.z_bot_->Clone();
  }
  else {
    z_bot_ = 0;
  }
  if (volume.erosion_top_ != 0) {
    erosion_top_ = volume.erosion_top_->Clone();
  }
  else {
    erosion_top_ = 0;
  }
  if (volume.erosion_bot_ != 0) {
    erosion_bot_ = volume.erosion_bot_->Clone();
  }
  else {
    erosion_bot_ = 0;
  }
}


Volume::~Volume()
{
  delete z_top_;
  delete z_bot_;
  delete erosion_top_;
  delete erosion_bot_;
}


Volume& Volume::operator=(const Volume& rhs)
{
  /// \todo Use "copy and swap" for exception safety.

  if (this == &rhs) return *this;

  x_min_ = rhs.x_min_;
  y_min_ = rhs.y_min_;
  lx_    = rhs.lx_;
  ly_    = rhs.ly_;
  lz_    = rhs.lz_;
  angle_ = rhs.angle_;
  
  delete z_top_;
  if (rhs.z_top_ != 0) {
    z_top_ = rhs.z_top_->Clone();
  }
  else {
    z_top_ = 0;
  }

  delete z_bot_;
  if (rhs.z_bot_ != 0) {
    z_bot_ = rhs.z_bot_->Clone();
  }
  else {
    z_bot_ = 0;
  }

  delete erosion_top_;
  if (rhs.erosion_top_ != 0) {
    erosion_top_ = rhs.erosion_top_->Clone();
  }
  else {
    erosion_top_ = 0;
  }

  delete erosion_bot_;
  if (rhs.erosion_bot_ != 0) {
    erosion_bot_ = rhs.erosion_bot_->Clone();
  }
  else {
    erosion_bot_ = 0;
  }

  return *this;
}


void Volume::SetDimensions(double x_min, double y_min, 
                                double lx, double ly)
{
  x_min_ = x_min;
  y_min_ = y_min;
  lx_ = lx;
  ly_ = ly;
 

  if (!CheckSurface(*z_top_)) {
    throw Exception("The top surface does not cover the volume.");
  }
  if (!CheckSurface(*z_bot_)) {
    throw Exception("The base surface does not cover the volume.");
  }
  if (erosion_top_ != 0 && !CheckSurface(*erosion_top_)) {
    throw Exception("The erosion top surface does not cover the volume.");
  }
  if (erosion_bot_ != 0 && !CheckSurface(*erosion_bot_)) {
    throw Exception("The erosion bottom surface does not cover the volume.");
  }
 lz_ = RecalculateLZ();
}

void Volume::SetAngle(double angle)
{
  angle_ = angle;
  

  if (!CheckSurface(*z_top_)) {
    throw Exception("The top surface does not cover the volume.");
  }
  if (!CheckSurface(*z_bot_)) {
    throw Exception("The bottom surface does not cover the volume.");
  }
  if (erosion_top_ != 0 && !CheckSurface(*erosion_top_)) {
    throw Exception("The erosion top surface does not cover the volume.");
  }
  if (erosion_bot_ != 0 && !CheckSurface(*erosion_bot_)) {
    throw Exception("The erosion bottom surface does not cover the volume.");
  }
  lz_ = RecalculateLZ();
}


void Volume::SetSurfaces(Surface* top_surf,
                         Surface* bot_surf,
                         Surface* erosion_top, 
                         Surface* erosion_bot)
{
  if(lx_ > 0 || ly_ > 0 ) { //Check that area is set.
    if (!CheckSurface(*top_surf)) {
      throw Exception("The top surface does not cover the volume.");
    }
    if (!CheckSurface(*bot_surf)) {
      throw Exception("The bottom surface does not cover the volume.");
    }
    if (erosion_top && !CheckSurface(*erosion_top)) {
      throw Exception("The erosion top surface does not cover the volume.");
    }
    if (erosion_bot && !CheckSurface(*erosion_bot)) {
      throw Exception("The erosion bottom surface does not cover the volume.");
    }
  }

  if (z_top_ != top_surf) {
    delete z_top_;
    z_top_ = top_surf;
  }

  if (z_bot_ != bot_surf) {
    delete z_bot_;
    z_bot_ = bot_surf;
  }

  if (erosion_top_ != erosion_top) {
    delete erosion_top_;
    erosion_top_ = erosion_top;
  }

  if (erosion_bot_ != erosion_bot) {
    delete erosion_bot_;
    erosion_bot_ = erosion_bot;
  }

  lz_ = RecalculateLZ();
}


// Writes surface to file if non-constant. Returns filename or
// surface level if surface is constant.
static std::string WriteSurface(const Surface* surf,
                                const std::string& grid_filename,
                                const std::string& surface_name
) 
{
  if (surf == 0) {
    return "0";
  }

  if (typeid(*surf) == typeid(ConstantSurface)) {
    return (ToString((dynamic_cast<const ConstantSurface*>(surf))->GetZ()));
  }
  else if (typeid(*surf) == typeid(RegularSurface<double>)) {
    std::string filename = grid_filename + surface_name;
    /// \todo Fix this.
    // std::string filename = MainPart(grid_filename) + surface_name + ".s";
    const RegularSurface<double>* rsurf 
      = dynamic_cast<const RegularSurface<double>*>(surf);
    WriteStormBinarySurf(*rsurf, filename);
    return filename;
  }
  else {
    throw Exception("Bug: Trying to write unsupported surface type to file.");
  }
}


void Volume::WriteVolumeToFile(std::ofstream& file, 
                                    const std::string& filename) const
{
  file << x_min_ << " " << lx_ << " " << y_min_ << " " << ly_ << " "
       << WriteSurface(z_top_, filename, "_top") << " "  
       << WriteSurface(z_bot_, filename, "_bot") << " "
       << WriteSurface(erosion_top_, filename, "_erosion_top") << " "  
       << WriteSurface(erosion_bot_, filename, "_erosion_bot") << "\n"
       << GetLZ() << " " << (180.0*angle_)/M_PI << "\n";
}


void Volume::ReadVolumeFromFile(std::ifstream& file, int line, const std::string& path)
{
  std::string token;

  GetNextToken(file, token, line);
  x_min_ = ParseType<double>(token);
  GetNextToken(file, token, line);
  lx_ = ParseType<double>(token);
  GetNextToken(file, token, line);
  y_min_ = ParseType<double>(token);
  GetNextToken(file, token, line);
  ly_ = ParseType<double>(token);

  GetNextToken(file, token, line);
  if (IsType<double>(token)) {
    z_top_ = new ConstantSurface(ParseType<double>(token));  
  } else {
    z_top_ = new RegularSurface<double>(ReadStormBinarySurf(path+token));
    if (!CheckSurface(*z_top_)) {
      throw Exception("The top surface does not fit with the volume.");
    }
  }
  GetNextToken(file, token, line);
  if (IsType<double>(token)) {
    z_bot_ = new ConstantSurface(ParseType<double>(token));  
  } else {
    z_bot_ = new RegularSurface<double>(ReadStormBinarySurf(path+token));
    if (!CheckSurface(*z_bot_)) {
      throw Exception("The bottom surface does not fit with the volume.");
    }
  }
  GetNextToken(file, token, line);
  if (IsType<double>(token)) {
    erosion_top_ = new ConstantSurface(ParseType<double>(token));  
  } else {
    erosion_top_ = new RegularSurface<double>(ReadStormBinarySurf(path+token));
    if (!CheckSurface(*erosion_top_)) {
      throw Exception("The erosion top surface does not fit with the volume.");
    }
  }
  GetNextToken(file, token, line);
  if (IsType<double>(token)) {
    erosion_bot_ = new ConstantSurface(ParseType<double>(token));  
  } else {
    erosion_bot_ = new RegularSurface<double>(ReadStormBinarySurf(path+token));
    if (!CheckSurface(*erosion_bot_)) {
      throw Exception("The erosion bottom surface does not fit with the volume.");
    }
  }

  GetNextToken(file, token, line);
  lz_ = ParseType<double>(token);

  GetNextToken(file, token, line);
  angle_ = (M_PI*ParseType<double>(token))/180.0;
}


void Volume::GlobalToLocalCoord(double global_x, 
                                     double global_y,
                                     double& local_x, 
                                     double& local_y) const
{
  double x_rel = global_x - x_min_;
  double y_rel = global_y - y_min_;

  local_x =   std::cos(angle_)*x_rel + std::sin(angle_)*y_rel;
  local_y = - std::sin(angle_)*x_rel + std::cos(angle_)*y_rel;
}


void Volume::LocalToGlobalCoord(double local_x, 
                                     double local_y,
                                     double& global_x, 
                                     double& global_y) const
{
  global_x = std::cos(angle_)*local_x - std::sin(angle_)*local_y + x_min_;
  global_y = std::sin(angle_)*local_x + std::cos(angle_)*local_y + y_min_;
}


double Volume::RecalculateLZ()
{
  double lz = 0;
  if(lx_ > 0 || ly_ > 0) { //Only do if area is initialized.
    // Just using a arbitary grid resolution.
    int nx = 100;
    int ny = 100;

    double dx = lx_ / nx;
    double dy = ly_ / ny;


    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x, y;
        LocalToGlobalCoord(dx * i, dy * j, x, y);
        lz = std::max(lz, z_bot_->GetZ(x, y) - z_top_->GetZ(x, y));
      }
    }
  }
  return lz;
}


bool Volume::CheckSurface(const Surface& surface) const
{
  std::vector<double> x(4);
  std::vector<double> y(4);
  
  LocalToGlobalCoord(  0,   0, x[0], y[0]);
  LocalToGlobalCoord(lx_,   0, x[1], y[1]);
  LocalToGlobalCoord(  0, ly_, x[2], y[2]);
  LocalToGlobalCoord(lx_, ly_, x[3], y[3]);

  double x_min = *(std::min_element(x.begin(), x.end()));
  double x_max = *(std::min_element(y.begin(), y.end()));
  double y_min = *(std::max_element(x.begin(), x.end()));
  double y_max = *(std::max_element(y.begin(), y.end()));

  return surface.EnclosesRectangle(x_min, y_min, x_max, y_max);
}
int Volume::isInside(double x, double y)
{
  double rx = (x-x_min_)*cos(angle_)+(y-y_min_)*sin(angle_);
  double ry = -(x-x_min_)*sin(angle_) + (y-y_min_)*cos(angle_);
  if(rx < 0 || rx > lx_ || ry<0 || ry > ly_)
    return(0);
  else
    return(1);
}
