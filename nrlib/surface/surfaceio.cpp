// $Id$

#include <fstream>
#include <string>

#include "surfaceio.hpp"

#include "regularsurface.hpp"
#include "surface.hpp"
#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"

using namespace NRLib2;

/// \todo Move to a suitable place
bool Equal(double a, double b)
{
  if ( fabs(a-b) < 0.0001 * a ) {
    return true;
  }
  return false;
}


RegularSurface<double> NRLib2::ReadStormSurf(const std::string& filename)
{
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file) {
    throw new IOError("Error operning " + filename + " for reading.");
  }
  int line = 0;
  std::string token;
  GetNextToken(file, token, line);
  file.close();
  RegularSurface<double> surface;
  if (token == "STORMGRID_BINARY")
    surface = NRLib2::ReadStormBinarySurf(filename);
  else
    surface = NRLib2::ReadIrapClassicAsciiSurf(filename);
  return(surface);
}
  
RegularSurface<double> NRLib2::ReadStormBinarySurf(const std::string& filename)
{
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file) {
    throw new IOError("Error operning " + filename + " for reading.");
  }
  
  int line = 0;
  
  // Header
  try {
    std::string token;
    GetNextToken(file, token, line);
    if (token != "STORMGRID_BINARY") {
      throw FileFormatError("Error reading " + filename + ", file is not "
                            "in STORM binary format.");
    }
    GetNextToken(file, token, line);
    int ni = ParseType<int>(token);
    GetNextToken(file, token, line);
    int nj = ParseType<int>(token);
    GetNextToken(file, token, line);
    double dx = ParseType<double>(token);
    GetNextToken(file, token, line);
    double dy = ParseType<double>(token);
    GetNextToken(file, token, line);
    double x_min = ParseType<double>(token);
    GetNextToken(file, token, line);
    double x_max = ParseType<double>(token);
    GetNextToken(file, token, line);
    double y_min = ParseType<double>(token);
    GetNextToken(file, token, line);
    double y_max = ParseType<double>(token);

    double lx = x_max - x_min;
    double ly = y_max - y_min;

    /// @todo good floating poing comparison.
    if (!Equal(lx/(ni-1), dx)) {
      throw FileFormatError("Inconsistent data in file. dx != lx/nx.");
    }
    if (!Equal(ly/(nj-1), dy)) {
      throw FileFormatError("Inconsistent data in file. dy != ly/ny.");
    }
    
    RegularSurface<double> surface(x_min, y_min, lx, ly, ni, nj);

    ReadBinaryDoubleArray(file, surface.begin(), surface.GetN());

    if (!CheckEndOfFile(file)) {
      throw FileFormatError("File too long.");
    }

    return surface;
  }
  catch (EndOfFile& ) {
    throw FileFormatError("Unexcpected end of file found while parsing "
      " \"" + filename + "\"");
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a " 
      "STORM surface file at line " + ToString(line) + ":" + e.what());
  }
}


RegularSurface<double> NRLib2::ReadIrapClassicAsciiSurf(const std::string& filename)
{
  std::ifstream file(filename.c_str(), std::ios::in);
  if (!file) {
    throw new IOError("Error operning " + filename + " for reading.");
  }
  
  int line = 0;
  std::string token;
  // Header
  try {
    GetNextToken(file, token, line);
    int ni = ParseType<int>(token);
    GetNextToken(file, token, line);
    int nj = ParseType<int>(token);
    GetNextToken(file, token, line);
    double dx = ParseType<double>(token);
    GetNextToken(file, token, line);
    double dy = ParseType<double>(token);
    GetNextToken(file, token, line);
    double x_min = ParseType<double>(token);
    GetNextToken(file, token, line);
    double x_max = ParseType<double>(token);
    GetNextToken(file, token, line);
    double y_min = ParseType<double>(token);
    GetNextToken(file, token, line);
    double y_max = ParseType<double>(token);

    double lx = x_max - x_min;
    double ly = y_max - y_min;

    /// @todo good floating poing comparison.
    if (!Equal(lx/(ni-1), dx)) {
      throw FileFormatError("Inconsistent data in file. dx != lx/nx.");
    }
    if (!Equal(ly/(nj-1), dy)) {
      throw FileFormatError("Inconsistent data in file. dy != ly/ny.");
    }
    
    RegularSurface<double> surface(x_min, y_min, lx, ly, ni, nj);

    int i;
    for(i=0;i<ni*nj;i++)
      file >> surface(i);

    if (!CheckEndOfFile(file)) {
      throw FileFormatError("File too long.");
    }

    return surface;
  }
  catch (EndOfFile& ) {
    throw FileFormatError("Unexcpected end of file found while parsing "
      " \"" + filename + "\"");
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a " 
      "STORM surface file at line " + ToString(line) + ":" + e.what());
  }
}


void NRLib2::WriteIrapClassicAsciiSurf(const RegularSurface<double>& surf, 
                                       const std::string& filename)
{
  /// \todo Replace with safe open function.
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
  if (!file) {
    throw new IOError("Error opening " + filename + " for writing.");
  }

  file << surf.GetNI() << " " << surf.GetNJ() << " " 
       << surf.GetDX() << " " << surf.GetDY() << "\n"
       << surf.GetXMin() << " " << surf.GetXMax() << " "
       << surf.GetYMin() << " " << surf.GetYMax() << "\n";

  int i;
  for(i=0;i<surf.GetN();i++)
    file << surf(i) << " ";
  file.close();
}

void NRLib2::WriteStormBinarySurf(const RegularSurface<double>& surf, 
                                  const std::string& filename)
{
  /// \todo Replace with safe open function.
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
  if (!file) {
    throw new IOError("Error opening " + filename + " for writing.");
  }
  
  file << "STORMGRID_BINARY\n\n"
       << surf.GetNI() << " " << surf.GetNJ() << " " 
       << surf.GetDX() << " " << surf.GetDY() << "\n"
       << surf.GetXMin() << " " << surf.GetXMax() << " "
       << surf.GetYMin() << " " << surf.GetYMax() << "\n";

  WriteBinaryDoubleArray(file, surf.begin(), surf.end());
  file.close();
}

//void WritePointAsciiSurf(const RegularSurface<double>& surf,
//                         const std::string& filename);
