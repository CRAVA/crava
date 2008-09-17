// $Id$

#include "stormcontgrid.hpp"

#include <fstream>
#include <math.h>
#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"
#include "../surface/surface.hpp"

using namespace NRLib2;

const float STD_MISSING_CODE = -999.0F;
const std::string format_desc[2] = {"storm_petro_binary", 
                                    "storm_petro_ascii"};


StormContGrid::StormContGrid(int nx, int ny, int nz)
  : Grid<float>(nx, ny, nz, STD_MISSING_CODE)
{
  // Default values
  file_format_ = STORM_BINARY;
  missing_code_ = STD_MISSING_CODE;
  zone_number_ = 0;
  model_file_name_ = "ModelFile";
  variable_name_ = "UNKNOWN";
}

StormContGrid::StormContGrid(const Volume &vol, int nx, int ny, int nz)
:Volume(vol)
{
  file_format_ = STORM_BINARY;
  missing_code_ = STD_MISSING_CODE;
  zone_number_ = 0;
  model_file_name_ = "ModelFile";
  variable_name_ = "UNKNOWN";
  Resize(nx,ny,nz);

}

StormContGrid::StormContGrid(const std::string& filename)
{
  ReadFromFile(filename);
}


void StormContGrid::ReadFromFile(const std::string& filename, bool commonPath)
{
   /// \todo Replace with safe open function.
  std::ifstream file(filename.c_str(), std::ios::out | std::ios::binary);

  std::string path = "";
  if(commonPath == true)
    path = GetPath(filename);

  if (!file) {
    throw new IOError("Error opening " + filename);
  }

  // Current line number
  int line = 0;
  
  // Header
  try {
    std::string token;
    GetNextToken(file, token, line);
    if (token == format_desc[STORM_BINARY]) {
      file_format_ = STORM_BINARY;
    } 
    else if (token == format_desc[STORM_ASCII]) {
      file_format_ = STORM_ASCII;
    }
    else {
      throw FileFormatError("Unknown format: " + token);
    }
  
    GetNextToken(file, token, line);
    zone_number_ = ParseType<int>(token);
    GetNextToken(file, model_file_name_, line);
    GetNextToken(file, token, line);
    missing_code_ = ParseType<double>(token);
    GetNextToken(file, variable_name_, line);

    ReadVolumeFromFile(file, line, path);

    GetNextToken(file, token, line);
    int nx = ParseType<int>(token);
    GetNextToken(file, token, line);
    int ny = ParseType<int>(token);
    GetNextToken(file, token, line);
    int nz = ParseType<int>(token);

    Resize(nx, ny, nz);
    
    switch (file_format_) {
    case STORM_BINARY:
      ReadBinaryFloatArray(file, begin(), GetN());
      break;
    case STORM_ASCII:
      for (int i = 0; i < GetN(); ++i) {
        GetNextToken(file, token, line);
        (*this)(i) = ParseType<float>(token);
      }
      break;
    default:
      throw Exception("Bug in STORM grid parser: unknown fileformat");
    }
 
    try {
      GetNextToken(file, token, line);
      int n_barriers = ParseType<int>(token);
      if (n_barriers != 0) {
        throw FileFormatError("Number of barriers greater than 0 found. "
          "Only grids without barriers are supported.");
      }
    }
    catch (EndOfFile& ) {
      // Number of barriers not present in file.
    }
  }
  catch (EndOfFile& ) {
    throw FileFormatError("Unexcpected end of file found while parsing "
      " \"" + filename + "\"");
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a " 
      "STORM file at line " + ToString(line) + ":" + e.what());
  }
}


void StormContGrid::WriteToFile(const std::string& filename, const std::string& predefinedHeader) const
{
  /// @todo Replace with safe open function.
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
  if (!file) {
    throw new IOError("Error opening " + filename);
  }

  // Header
  if(predefinedHeader == "") {
    file << format_desc[file_format_] << "\n\n"
         << zone_number_ << " " << model_file_name_ << " " 
         << missing_code_ << "\n\n" << variable_name_ << "\n\n" ;
    
    WriteVolumeToFile(file, filename);
    file << "\n";
    file << GetNI() << " " << GetNJ() << " " << GetNK() << "\n";
  }
  else
    file << predefinedHeader;
  // Data
  switch (file_format_) {
  case STORM_BINARY:
    WriteBinaryFloatArray(file, begin(), end());
    break;
  case STORM_ASCII:
    for (const_iterator it = begin(); it != end(); ++it) {
      file << *it << " ";
    }
    break;
  default:
    throw Exception("Unknown fileformat");
  }
  
  // Final 0 (Number of barriers)
  file << 0;
}


/// \todo Common implementation with StormFaciesGrid
int StormContGrid::FindIndex(double x, double y, double z) const 
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);

  int i = static_cast<int>(local_x / GetDX());
  int j = static_cast<int>(local_y / GetDY());
  double z_top = GetTopSurface().GetZ(x, y);
  double z_bot = GetBotSurface().GetZ(x, y);
  double dz = (z_bot - z_top) / GetNK();
  int k;
  if (dz == 0) {
    k = 0;
  }
  else {
    k = static_cast<int>((z - z_top)/dz);
  }

  return GetIndex(i, j, k);
}

double StormContGrid::getValueZInterpolated(double x, double y, double z)const
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);

  int i = static_cast<int>(local_x / GetDX());
  int j = static_cast<int>(local_y / GetDY());
  double z_top = GetTopSurface().GetZ(x, y);
  double z_bot = GetBotSurface().GetZ(x, y);
  double dz = (z_bot - z_top) / GetNK();
  if(i<0 || i>GetNI() || j<0 || j> GetNJ())
    return(missing_code_);

  float value;
  if(z<=z_top+0.5*dz)
  {
    value = (*this)(GetIndex(i,j,0));
   // z = z_top+0.5*dz;
  }
  else if(z>=z_bot+0.5*dz)
  {
    value = (*this)(GetIndex(i,j,GetNK()-1));
   // z = z_bot+0.5*dz;
  }
  else
  {
    int zInd1 = static_cast<int>(floor((z-z_top)/dz)-0.5);
    double t = (z-z_top)/dz - 0.5 - static_cast<double>(zInd1);
    int zInd2 = zInd1+1;
    value = (*this)(GetIndex(i,j,zInd1))*(1-t)+(*this)(GetIndex(i,j,zInd2))*t;
  }
  return(value);

  
}

double StormContGrid::getValueClosestInZ(double x, double y, double z)const
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);

  int i = static_cast<int>(local_x / GetDX());
  int j = static_cast<int>(local_y / GetDY());
  double z_top = GetTopSurface().GetZ(x, y);
  double z_bot = GetBotSurface().GetZ(x, y);
  double dz = (z_bot - z_top) / GetNK();
  if(i<0 || i>GetNI() || j<0 || j> GetNJ())
    return(missing_code_);

  if(z<z_top+0.5*dz)
    z = z_top+0.5*dz;
  else if(z>z_bot+0.5*dz)
    z = z_bot+0.5*dz;
  
  int zInd1 = static_cast<int>(floor((z-z_top)/dz)-0.5);
  float t = (z-z_top)/dz - 0.5 - static_cast<double>(zInd1);
  int zInd2 = zInd1+1;
  float value = (1-t)*(*this)(GetIndex(i,j,zInd1))+t*(*this)(GetIndex(i,j,zInd2));
  return(value);
  

}
double StormContGrid::RecalculateLZ()
{
  double lz = 0;
  for (int i = 0; i < GetNI(); ++i) {
    for (int j = 0; j < GetNJ(); ++j) {
      double x, y;
      LocalToGlobalCoord(i * GetDX(), j * GetDY(), x, y);
      lz = std::max(lz, GetBotSurface().GetZ(x, y) - 
                        GetTopSurface().GetZ(x, y));
    }
  }

  return lz;
}
