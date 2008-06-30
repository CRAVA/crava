// $Id: surfaceio.hpp 44 2008-04-11 13:48:08Z perroe $

#ifndef NRLIB_SURFACEIO_HPP
#define NRLIB_SURFACEIO_HPP

#include "surface.hpp"
#include "regularsurface.hpp"

namespace NRLib2 {
  //enum FileFormat{
  //  SURF_IRAP_ASCII,
  //  SURF_STORM_BINARY
  //};

  //Surface ReadSurface(const std::string& filename);

  ///This function will determine format, and call either
  ///ReadIrapClassic or ReadStormBinary.
  RegularSurface<double> ReadStormSurf(const std::string& filename);

  RegularSurface<double> ReadIrapClassicAsciiSurf(const std::string& filename);

  RegularSurface<double> ReadStormBinarySurf(const std::string& filename);

  // void WriteIrapClassicAsciiSurf(const RegularSurface<double>& surf, 
  //                               const std::string& filename);

  void WriteStormBinarySurf(const RegularSurface<double>& surf, 
                            const std::string& filename);

  // void WritePointAsciiSurf(const RegularSurface<double>& surf,
  //                         const std::string& filename);
}

#endif // NRLIB_SURFACEIO_HPP
