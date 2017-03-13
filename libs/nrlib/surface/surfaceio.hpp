// $Id: surfaceio.hpp 1232 2014-01-13 12:22:40Z gudmundh $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NRLIB_SURFACEIO_HPP
#define NRLIB_SURFACEIO_HPP

#include <string>
#include <vector>
#include <locale>

namespace NRLib {
  template <class A> class RegularSurface;
  template <class A> class RegularSurfaceRotated;
  template <class A> class Grid2D;

  const double MULT_IRAP_MISSING =   -999.25;
  const double IRAP_MISSING      = 9999900.0;
  const double STORM_MISSING     =    -999.0;

  enum SurfaceFileFormat {
    SURF_UNKNOWN,
    SURF_IRAP_CLASSIC_ASCII,
    SURF_STORM_BINARY,
    SURF_SGRI,
    SURF_RMS_POINTS_ASCII,
    SURF_MULT_ASCII,
    SURF_XYZ_ASCII
    //  SURF_PLAIN_ASCII
    //  SURF_CPS3_ASCII
  };

  /// \brief Find type of file.
  SurfaceFileFormat FindSurfaceFileType(const std::string& filename);

  /// \brief String describing file format
  std::string GetSurfFormatString(SurfaceFileFormat format);

  template <class A>
  void ReadStormBinarySurf(const std::string & filename,
                           RegularSurface<A> & surface);

  template <class A>
  void ReadIrapClassicAsciiSurf(const std::string & filename,
                                RegularSurface<A> & surface,
                                double            & angle);

  template <class A>
  void ReadSgriSurf(const std::string & filename,
                    RegularSurface<A> & surface,
                    double            & angle);

  // If labels is non-empty, the labels of the axes on the file are compared with these. Throws if mismatch.
  std::vector<RegularSurfaceRotated<float> > ReadMultipleSgriSurf(const std::string& filename,
                                                                  const std::vector<std::string> & labels);

  template <class A>
  void ReadXYZAsciiSurf(std::string         filename,
                        RegularSurface<A> & surface,
                        double              x_ref,
                        double              y_ref,
                        double              lx,
                        double              ly,
                        int               * ilxl_area,
                        double              il0_ref,
                        double              xl0_ref,
                        double              il_step_x,
                        double              il_step_y,
                        double              xl_step_x,
                        double              xl_step_y);

  template <class A>
  void ReadMulticolumnAsciiSurf(std::string         filename,
                                RegularSurface<A> & surface,
                                double              x_ref,
                                double              y_ref,
                                double              lx,
                                double              ly,
                                int               * ilxl_area,
                                double              il0_ref,
                                double              xl0_ref);

  template <class A>
  void CreateSurfaceFromILXL(RegularSurface<A>   & surface,
                             std::vector<int>    & il_vec,
                             std::vector<int>    & xl_vec,
                             std::vector<double> & z_vec,
                             int                 * ilxl_area,
                             double              & il0_ref,
                             double              & xl0_ref,
                             double              & x_ref,
                             double              & y_ref,
                             double              & lx,
                             double              & ly);

  template <class A>
  void WriteIrapClassicAsciiSurf(const RegularSurface<A> & surf,
                                 double                    angle,
                                 const std::string       & filename);

  template <class A>
  void WriteStormBinarySurf(const RegularSurface<A> & surf,
                            const std::string       & filename);

  bool FindXYZAsciiLine(const std::string& filename, int & header_start_line);

  bool FindMulticolumnAsciiLine(const std::string& filename, int & header_start_line);

  // void WritePointAsciiSurf(const RegularSurface<double>& surf,
  //                         const std::string& filename);

  namespace NRLibPrivate {
    /// \todo Move to a suitable place
    bool Equal(double a, double b);
  }

} // namespace NRLib

  // ----------- TEMPLATE IMPLEMENTATIONS ----------------------

#include <fstream>
#include <string>

#include "regularsurface.hpp"
#include "surface.hpp"
#include "../exception/exception.hpp"
#include "../math/constants.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"


template<class A>
void NRLib::ReadStormBinarySurf(const std::string & filename,
                                RegularSurface<A> & surface)
{
  std::ifstream file;
  OpenRead(file, filename.c_str(), std::ios::in | std::ios::binary);

  int line = 0;

  // Header
  try {
    std::string token = ReadNext<std::string>(file, line);
    if (token != "STORMGRID_BINARY") {
      throw FileFormatError("Error reading " + filename + ", file is not "
                            "in STORM binary format.");
    }

    int ni       = ReadNext<int>(file, line);
    int nj       = ReadNext<int>(file, line);
    double dx    = ReadNext<double>(file, line);
    double dy    = ReadNext<double>(file, line);
    double x_min = ReadNext<double>(file, line);
    double x_max = ReadNext<double>(file, line);
    double y_min = ReadNext<double>(file, line);
    double y_max = ReadNext<double>(file, line);

    double lx = x_max - x_min;
    double ly = y_max - y_min;

    if (!NRLibPrivate::Equal(lx/(ni-1), dx)) {
      throw FileFormatError("Inconsistent data in file. dx != lx/(nx-1).");
    }
    if (!NRLibPrivate::Equal(ly/(nj-1), dy)) {
      throw FileFormatError("Inconsistent data in file. dy != ly/(ny-1).");
    }

    surface.Resize(ni, nj);
    surface.SetDimensions(x_min, y_min, lx, ly);

    DiscardRestOfLine(file, line, true);
    ReadBinaryDoubleArray(file, surface.begin(), surface.GetN());

    surface.SetMissingValue(static_cast<A>(STORM_MISSING));

    if (!CheckEndOfFile(file)) {
      throw FileFormatError("File too long.");
    }

    surface.SetName(GetStem(filename));
  }
  catch (EndOfFile& ) {
    throw FileFormatError("Unexcpected end of file found while parsing "
                           " \"" + filename + "\"");
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a "
      "STORM surface file at line " + ToString(line) + ":" + e.what() + "\n");
  }
}


template <class A>
void NRLib::ReadIrapClassicAsciiSurf(const std::string & filename,
                                     RegularSurface<A> & surface,
                                     double            & angle)
{
  std::ifstream file;
  OpenRead(file, filename);

  int line = 0;
  // Header
  try {
    ReadNext<int>(file, line);   // -996
    int nj       = ReadNext<int>(file, line);
    double dx    = ReadNext<double>(file, line);
    double dy    = ReadNext<double>(file, line);
    // ----------- line shift --------------
    double x_min = ReadNext<double>(file, line);
    double x_max = ReadNext<double>(file, line);
    double y_min = ReadNext<double>(file, line);
    double y_max = ReadNext<double>(file, line);
    // ----------- line shift --------------
    int ni       = ReadNext<int>(file, line);
    angle        = ReadNext<double>(file, line);
    angle        = NRLib::Degree*angle;
    ReadNext<double>(file, line); // rotation origin - x
    ReadNext<double>(file, line); // rotation origin - y
    // ----------- line shift --------------
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    ReadNext<int>(file, line);
    //double lx = (ni-1)*dx;
    //double ly = (nj-1)*dy;
    double lx = (x_max - x_min)*cos(angle);
    double ly = (y_max - y_min)*cos(angle);


    if (!NRLibPrivate::Equal(lx/(ni-1), dx)) {
      std::string text = "Inconsistent data in file. dx != lx/(nx-1).\n";
      text += "dx        = "+NRLib::ToString(dx,2)+"\n";
      text += "lx        = "+NRLib::ToString(lx,2)+"\n";
      text += "nx        = "+NRLib::ToString(ni,0)+"\n";
      text += "lx/(nx-1) = "+NRLib::ToString(lx/(ni - 1),2);
      throw FileFormatError(text);
    }
    if (!NRLibPrivate::Equal(ly/(nj-1), dy)) {
      std::string text = "Inconsistent data in file. dy != ly/(ny-1).\n";
      text += "dy        = "+NRLib::ToString(dy,2)+"\n";
      text += "ly        = "+NRLib::ToString(ly,2)+"\n";
      text += "ny        = "+NRLib::ToString(nj,0)+"\n";
      text += "ly/(ny-1) = "+NRLib::ToString(ly/(nj - 1),2);
      throw FileFormatError(text);
    }


    surface.Resize(ni, nj);
    surface.SetDimensions(x_min, y_min, lx, ly);

    ReadAsciiArrayFast(file, surface.begin(), surface.GetN());

    surface.SetMissingValue(static_cast<A>(IRAP_MISSING));

    surface.SetName(GetStem(filename));

    if (!CheckEndOfFile(file)) {
      throw FileFormatError("File too long.");
    }
  }
  catch (EndOfFile& ) {
    throw FileFormatError("Unexcpected end of file found while parsing "
                          " \"" + filename + "\"");
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a "
      "IRAP ASCII surface file at line " + ToString(line) + ":" + e.what() + "\n");
  }
}


template<class A>
void  NRLib::ReadSgriSurf(const std::string & filename,
                          RegularSurface<A> & surface,
                          double            & angle)
{
  std::ifstream header_file;
  OpenRead(header_file, filename.c_str(), std::ios::in | std::ios::binary);
  int i;
  std::string tmp_str;
  int dim;
  try {
    //Reading record 1: Version header
    getline(header_file, tmp_str);
    //Reading record 2: Grid dimension
    header_file >> dim;
    if(dim!=2)
      throw Exception("Wrong dimension of Sgri file. We expect a surface, dimension should be 2.\n");

    getline(header_file, tmp_str);
    //Reading record 3 ... 3+dim: Axis labels + grid value label
    std::vector<std::string> axis_labels(dim);
    for (i=0; i<dim; i++)
      getline(header_file, axis_labels[i]);
    if (((axis_labels[0].find("X") == std::string::npos) && (axis_labels[0].find("x") == std::string::npos)) ||
      ((axis_labels[1].find("Y") == std::string::npos) && (axis_labels[1].find("y") == std::string::npos)))
      throw Exception("Wrong axis labels. First axis should be x-axis, second axis should be y-axis.\n");
    // if((axis_labels[0]!="X" && axis_labels[0] !="x") || (axis_labels[1]!="Y" && axis_labels[1]!="y"))
    //  throw Exception("Wrong axis labels. First axis should be x-axis, second axis should be y-axis.\n");
    getline(header_file, tmp_str);
    //int config = IMISSING;

    //Reading record 4+dim: Number of grids
    int n_grid;
    header_file >> n_grid;
    if (n_grid < 1) {
      throw Exception("Error: Number of grids read from sgri file must be >0");
    }
    getline(header_file, tmp_str);
    //Reading record 5+dim ... 5+dim+ngrid-1: Grid labels

    for (i=0; i<n_grid; i++)
      getline(header_file, tmp_str);

    std::vector<float> d_values1(dim);
    std::vector<float> d_values2(dim);
    std::vector<int>   i_values(dim);
    //Reading record 5+dim+ngrid: Scaling factor of grid values
    for (i=0; i<dim; i++)
      header_file >> d_values1[i];
    getline(header_file,tmp_str);
    //Reading record 6+dim+ngrid: Number of samples in each dir.
    for (i=0; i<dim; i++)
      header_file >> i_values[i];
    getline(header_file,tmp_str);
    //Reading record 7+dim+ngrid: Grid sampling in each dir.
    for (i=0; i<dim; i++) {
      header_file >> d_values2[i];
    }
    getline(header_file,tmp_str);
    //Reading record 8+dim+ngrid: First point coord.
    std::vector<float> min_values(dim);
    for (i=0; i<dim; i++)
    {
      header_file >> min_values[i];
    }
    int nx = 1;
    int ny = 1;

    double dx, dy;
    nx      = i_values[0];
    dx     = d_values2[0];
    ny      = i_values[1];
    dy     = d_values2[1];

    if (nx < 1) {
      throw Exception("Error: Number of samples in X-dir must be >= 1.\n");
    }
    if (ny < 1) {
      throw Exception("Error: Number of samples in Y-dir must be >= 1.\n");
    }

    if (dx <= 0.0) {
      throw Exception("Error: Grid sampling in X-dir must be > 0.0.\n");

    }
    if (dy <= 0.0) {
      throw Exception("Error: Grid sampling in Y-dir must be > 0.0.\n");
    }

    double lx = nx*dx;
    double ly = ny*dy;

    double x_min = min_values[0]-0.5*dx; //In regular grid, these are at value;
    double y_min = min_values[1]-0.5*dy; //in sgri, at corner of cell, hence move.

    header_file >> angle;

    surface.Resize(nx, ny, 0.0);
    surface.SetDimensions(x_min, y_min, lx, ly);

    getline(header_file, tmp_str);
    //Reading record 10+dim+ngrid: Undef value
    float missing_code;
    header_file >> missing_code;
    surface.SetMissingValue(missing_code);
    getline(header_file, tmp_str);
    //Reading record 11+dim+ngrid: Filename of binary file
    std::string bin_file_name;
    getline(header_file, tmp_str);
    if (!tmp_str.empty()) {
      std::locale loc;
      int i = 0;
      char c = tmp_str[i];
      while (!std::isspace(c,loc)) {
        i++;
        c = tmp_str[i];
      }
      tmp_str.erase(tmp_str.begin()+i, tmp_str.end());
    }
    if (tmp_str.empty())
      bin_file_name = NRLib::ReplaceExtension(filename, "Sgri");
    else {
      std::string path = GetPath(filename);
      bin_file_name = path + "/" + tmp_str;
    }
    //Reading record 12+dim+ngrid: Complex values
    bool has_complex;
    header_file >> has_complex;
    if (has_complex != 0 ) {
      throw Exception("Error: Can not read Sgri binary file. Complex values?");
    }

    surface.SetName(GetStem(bin_file_name));

    std::ifstream bin_file;
    OpenRead(bin_file, bin_file_name, std::ios::in | std::ios::binary);
    ReadBinaryFloatArray(bin_file, surface.begin(), surface.GetN());
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a "
      "Sgri surface file " + e.what() + "\n");
  }
}

template <class A>
void NRLib::ReadXYZAsciiSurf(std::string         filename,
                             RegularSurface<A> & surface,
                             double              x_ref,
                             double              y_ref,
                             double              lx,
                             double              ly,
                             int               * ilxl_area,
                             double              il0_ref,
                             double              xl0_ref,
                             double              il_step_x,
                             double              il_step_y,
                             double              xl_step_x,
                             double              xl_step_y)
{
  try {
    //Create surface with area corresponding to segy-grid (ilxl_area), but use sampling from surface file
    //il0_ref and xl0_ref are il/xl values at rotation corner, this is used to match IL/XL values from surface with segy

    int header_start_line;
    bool is_xyz_ascii = false;
    is_xyz_ascii = FindXYZAsciiLine(filename, header_start_line);

    if (is_xyz_ascii == false)
      throw Exception("Error: Did not recognize file as an xyz ascii file. We require the ascii file to have three columns X,Y,Z in an arbitrary order.\n");

    std::ifstream file;
    NRLib::OpenRead(file, filename);
    int line = 0;
    int line_number = 1;

    for (int i = 0; i < header_start_line; i++) {
      NRLib::DiscardRestOfLine(file, line, false);
      line_number++;
    }

    int x_index  = -1;
    int y_index  = -1;
    int z_index  = -1;
    if (header_start_line > -1) {
      //Header line, contains X, Y, Z
      std::vector<std::string> variable_names(3);
      variable_names[0] = NRLib::ReadNext<std::string>(file, line);
      variable_names[1] = NRLib::ReadNext<std::string>(file, line);
      variable_names[2] = NRLib::ReadNext<std::string>(file, line);

      for (int i = 0; i < 3; i++) {
        if (NRLib::Uppercase(variable_names[i]) == "X" || NRLib::Uppercase(variable_names[i]) == "UTMX")
          x_index = i;
        if (NRLib::Uppercase(variable_names[i]) == "Y" || NRLib::Uppercase(variable_names[i]) == "UTMY")
          y_index = i;
        if (NRLib::Uppercase(variable_names[i]) == "ATTRIBUTE" || NRLib::Uppercase(variable_names[i]) == "Z" || NRLib::Uppercase(variable_names[i]) == "TWT")
          z_index = i;
      }

      std::string err_txt = "";
      if (x_index == -1)
        err_txt += "Could not find variable name for X in file " + filename +". (Tried X, UTMX).\n";
      if (y_index == -1)
        err_txt += "Could not find variable name for Y in file " + filename +". (Tried Y, UTMY).\n";
      if (z_index == -1)
        err_txt += "Could not find variable name for Attribute in file " + filename +". (Tried Attribute, Z, TWT).\n";

      if (err_txt != "")
        throw Exception("Error when finding header information in " + filename + " :" + err_txt + "We require the xyz ascii file to have three columns (X,Y,Z) in an arbitrary order.");

      line_number++;
    }
    else {
      //No header, assume order is: x, y, z
      x_index  = 0;
      y_index  = 1;
      z_index  = 2;
    }

    std::vector<std::vector<double> > data(3);
    std::string line_string;
    while (NRLib::CheckEndOfFile(file)==false) {

      //Read line
      std::string line_string;
      std::getline(file, line_string);
      std::vector<std::string> tokens = NRLib::GetTokens(line_string);

      if (tokens.size() != 3) {
        throw FileFormatError("Error reading line " + NRLib::ToString(line_number) + " in " + filename + ", expected three numbers.");
      }
      else {
        for (int i = 0; i < 3; i++) {

          if (NRLib::IsNumber(tokens[i].c_str()))
            data[i].push_back(atof(tokens[i].c_str()));
          else
            throw FileFormatError("Error reading line " + NRLib::ToString(line_number) + " in " + filename + ", expected three numbers.");
        }
      }
      line_number++;
    }

    std::vector<double> x_vec = data[x_index];
    std::vector<double> y_vec = data[y_index];
    std::vector<double> z_vec = data[z_index];


    //Create IL, XL vector from x, y
    std::vector<int> il_vec = std::vector<int>(x_vec.size());
    std::vector<int> xl_vec = std::vector<int>(x_vec.size());
    for (int i = 0; i < static_cast<int>(x_vec.size()); i++) {
      il_vec[i] = static_cast<int>(0.5+il0_ref+(x_vec[i]-x_ref)*il_step_x+(y_vec[i]-y_ref)*il_step_y);
      xl_vec[i] = static_cast<int>(0.5+xl0_ref+(x_vec[i]-x_ref)*xl_step_x+(y_vec[i]-y_ref)*xl_step_y);
    }

    CreateSurfaceFromILXL(surface, il_vec, xl_vec, z_vec, ilxl_area, il0_ref, xl0_ref, x_ref, y_ref, lx, ly);

    surface.SetMissingValue(static_cast<A>(MULT_IRAP_MISSING));
    surface.SetName(GetStem(filename));

  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as an "
      "XYZ ASCII file: \n" + e.what() + "\n");
  }

}

template <class A>
void NRLib::ReadMulticolumnAsciiSurf(std::string         filename,
                                     RegularSurface<A> & surface,
                                     double              x_ref,
                                     double              y_ref,
                                     double              lx,
                                     double              ly,
                                     int               * ilxl_area,
                                     double              il0_ref,
                                     double              xl0_ref)
{

  try {
    //Create surface with area corresponding to segy-grid (ilxl_area), but use sampling from surface file
    //il0_ref and xl0_ref are il/xl values at rotation corner, this is used to match IL/XL values from surface with segy

    int header_start_line;
    bool is_multicolumn_ascii = false;
    is_multicolumn_ascii = FindMulticolumnAsciiLine(filename, header_start_line);

    if (is_multicolumn_ascii == false)
      throw Exception("Error: Did not recognize file as a multicolumns ascii file. We require the ascii file to have five columns X,Y,Z,IL and XL in an arbitrary order, each with its own header.\n");

    std::ifstream file;
    NRLib::OpenRead(file, filename);
    int line = 0;
    int line_number = 1;

    for (int i = 0; i < header_start_line; i++) {
      NRLib::DiscardRestOfLine(file, line, false);
      line_number++;
    }

    //Header line, contains X, Y, Z, Inline, Crossline
    std::vector<std::string> variable_names(5);
    variable_names[0] = NRLib::ReadNext<std::string>(file, line);
    variable_names[1] = NRLib::ReadNext<std::string>(file, line);
    variable_names[2] = NRLib::ReadNext<std::string>(file, line);
    variable_names[3] = NRLib::ReadNext<std::string>(file, line);
    variable_names[4] = NRLib::ReadNext<std::string>(file, line);
    line_number++;

    std::vector<std::vector<double> > data(5);
    std::string line_string;
    while (NRLib::CheckEndOfFile(file)==false) {

      //Read line
      std::string line_string;
      std::getline(file, line_string);
      std::vector<std::string> tokens = NRLib::GetTokens(line_string);

      if (tokens.size() != 5) {
        throw FileFormatError("Error reading line " + NRLib::ToString(line_number) + " in " + filename + ", expected five numbers.");
      }
      else {
        for (int i = 0; i < 5; i++) {

          if (NRLib::IsNumber(tokens[i].c_str()))
            data[i].push_back(atof(tokens[i].c_str()));
          else
            throw FileFormatError("Error reading line " + NRLib::ToString(line_number) + " in " + filename + ", expected five numbers.");
        }
      }
      line_number++;
    }


    int il_index = -1;
    int xl_index = -1;
    int x_index  = -1;
    int y_index  = -1;
    int z_index  = -1;
    for (int i = 0; i < 5; i++) {
      if (NRLib::Uppercase(variable_names[i]) == "INLINE" || NRLib::Uppercase(variable_names[i]) == "IL")
        il_index = i;
      if (NRLib::Uppercase(variable_names[i]) == "CROSSLINE" || NRLib::Uppercase(variable_names[i]) == "XL")
        xl_index = i;
      if (NRLib::Uppercase(variable_names[i]) == "X" || NRLib::Uppercase(variable_names[i]) == "UTMX")
        x_index = i;
      if (NRLib::Uppercase(variable_names[i]) == "Y" || NRLib::Uppercase(variable_names[i]) == "UTMY")
        y_index = i;
      if (NRLib::Uppercase(variable_names[i]) == "ATTRIBUTE" || NRLib::Uppercase(variable_names[i]) == "Z" || NRLib::Uppercase(variable_names[i]) == "TWT")
        z_index = i;
    }

    std::string err_txt = "";
    if (il_index == -1)
      err_txt += "Could not find variable name for Inline in file " + filename +". (Tried Inline, IL).\n";
    if (xl_index == -1)
      err_txt += "Could not find variable name for Crossline in file " + filename +". (Tried Crossline, XL).\n";
    if (x_index == -1)
      err_txt += "Could not find variable name for X in file " + filename +". (Tried X, UTMX).\n";
    if (y_index == -1)
      err_txt += "Could not find variable name for Y in file " + filename +". (Tried Y, UTMY).\n";
    if (z_index == -1)
      err_txt += "Could not find variable name for Attribute in file " + filename +". (Tried Attribute, Z, TWT).\n";

    if (err_txt != "")
      throw Exception("Error when finding header information in " + filename + " :" + err_txt + "\n. We require the multicolumn ascii file to have five columns (X,Y,Z,IL,XL) in an arbitrary order.");

    std::vector<double> x_vec = data[x_index];
    std::vector<double> y_vec = data[y_index];
    std::vector<double> z_vec = data[z_index];

    std::vector<int> il_vec = std::vector<int>(x_vec.size());
    std::vector<int> xl_vec = std::vector<int>(x_vec.size());

    for (int i = 0; i < static_cast<int>(x_vec.size()); i++) {
      il_vec[i] = static_cast<int>(data[il_index][i]);
      xl_vec[i] = static_cast<int>(data[xl_index][i]);
    }

    CreateSurfaceFromILXL(surface, il_vec, xl_vec, z_vec, ilxl_area, il0_ref, xl0_ref, x_ref, y_ref, lx, ly);

    surface.SetMissingValue(static_cast<A>(MULT_IRAP_MISSING));
    surface.SetName(GetStem(filename));
  }
  catch (Exception& e) {
    throw FileFormatError("Error parsing \"" + filename + "\" as a "
      "Multicolumn ASCII file: \n" + e.what() + "\n");
  }

}

template <class A>
void NRLib::CreateSurfaceFromILXL(RegularSurface<A>   & surface,
                                  std::vector<int>    & il_vec,
                                  std::vector<int>    & xl_vec,
                                  std::vector<double> & z_vec,
                                  int                 * ilxl_area,
                                  double              & il0_ref,
                                  double              & xl0_ref,
                                  double              & x_ref,
                                  double              & y_ref,
                                  double              & lx,
                                  double              & ly)

{

  //Find min and max IL/XL
  std::vector<int> il_vec_sorted = il_vec;
  std::sort(il_vec_sorted.begin(), il_vec_sorted.end());
  il_vec_sorted.erase(std::unique(il_vec_sorted.begin(), il_vec_sorted.end()), il_vec_sorted.end());

  std::vector<int> xl_vec_sorted = xl_vec;
  std::sort(xl_vec_sorted.begin(), xl_vec_sorted.end());
  xl_vec_sorted.erase(std::unique(xl_vec_sorted.begin(), xl_vec_sorted.end()), xl_vec_sorted.end());

  int ni_file = static_cast<int>(il_vec_sorted.size());
  int nj_file = static_cast<int>(xl_vec_sorted.size());

  int il_min_file = il_vec_sorted[0];
  int il_max_file = il_vec_sorted[ni_file-1];

  int xl_min_file = xl_vec_sorted[0];
  int xl_max_file = xl_vec_sorted[nj_file-1];

  int n = static_cast<int>(il_vec.size());
  A missing = static_cast<A>(NRLib::MULT_IRAP_MISSING);

  NRLib::Grid2D<A> ilxl_grid_file_z(ni_file, nj_file, missing);
  NRLib::Grid2D<A> ilxl_grid_file_il(ni_file, nj_file, missing);
  NRLib::Grid2D<A> ilxl_grid_file_xl(ni_file, nj_file, missing);

  int d_il_file = static_cast<int>((il_max_file - il_min_file) / (ni_file - 1));
  int d_xl_file = static_cast<int>((xl_max_file - xl_min_file) / (nj_file - 1));

  //We require the surfaces to have even sampling
  for (int i = 1; i < static_cast<int>(il_vec_sorted.size()); i++)
    if ((il_vec_sorted[i] - il_vec_sorted[i-1]) != d_il_file)
        throw Exception("Found inconsistencies in the surface sampling. Has the surface even sampling?\n");

  for (int i = 1; i < static_cast<int>(xl_vec_sorted.size()); i++)
    if ((xl_vec_sorted[i] - xl_vec_sorted[i-1]) != d_xl_file)
        throw Exception("Found inconsistencies in the surface sampling. Has the surface even sampling?\n");

  for (int k = 0; k < n; k++) {
    //Local IL/XL
    int il_loc = (static_cast<int>(il_vec[k]) - il_min_file)/d_il_file;
    int xl_loc = (static_cast<int>(xl_vec[k]) - xl_min_file)/d_xl_file;

    ilxl_grid_file_z(il_loc, xl_loc)  = static_cast<A>(z_vec[k]);
    ilxl_grid_file_il(il_loc, xl_loc) = static_cast<A>(il_vec[k]);
    ilxl_grid_file_xl(il_loc, xl_loc) = static_cast<A>(xl_vec[k]);
  }

  //Check consistency between IL/XL sampling in surface and expected ((max-min)/n)
  //We only do this for one path in the center of the surface
  int j_index = nj_file/2;
  for (int i = 1; i < ni_file; i++) {
    if (ilxl_grid_file_il(i-1,j_index) != missing && ilxl_grid_file_il(i,j_index) != missing) {
      int diff_il = static_cast<int>(std::abs(ilxl_grid_file_il(i-1,j_index) - ilxl_grid_file_il(i,j_index)));

      if (diff_il != d_il_file) {
        throw Exception("Found sampling of IL-values in file to be " + NRLib::ToString(diff_il) +
                        ", expected sampling equal to " + NRLib::ToString(d_il_file) + ".\n");
      }

    }
  }

  int i_index = ni_file/2;
  for (int j = 1; j < nj_file; j++) {
    if (ilxl_grid_file_xl(i_index,j-1) != missing && ilxl_grid_file_xl(i_index,j) != missing) {
      int diff_xl = static_cast<int>(std::abs(ilxl_grid_file_xl(i_index,j-1) - ilxl_grid_file_xl(i_index,j)));

      if (diff_xl != d_xl_file) {
        throw Exception("Found sampling of IL-values in file to be " + NRLib::ToString(diff_xl) +
                        ", expected sampling equal to " + NRLib::ToString(d_xl_file) + ".\n");
      }

    }
  }

  int il_min_segy  = ilxl_area[0];
  int il_max_segy  = ilxl_area[1];
  int xl_min_segy  = ilxl_area[2];
  int xl_max_segy  = ilxl_area[3];
  int il_step_segy = ilxl_area[4];
  int xl_step_segy = ilxl_area[5];

  //Create IL/XL surface as large as segy geometry, but with sampling from file
  int n_il = (il_max_segy - il_min_segy)/d_il_file + 1;
  int n_xl = (xl_max_segy - xl_min_segy)/d_xl_file + 1;

  //Find IL/XL of rotation corner
  int il0_segy = static_cast<int>(il0_ref+0.5);
  int xl0_segy = static_cast<int>(xl0_ref+0.5);

  // To ensure that the IL XL we find are existing traces
  if (il0_segy < il_min_segy)
    il0_segy -= (il0_segy - il_min_segy) % il_step_segy;
  else if (il0_segy > il_max_segy)
    il0_segy += (il_max_segy - il0_segy) % il_step_segy;

  if (xl0_segy < xl_min_segy)
    xl0_segy -= (xl0_segy - xl_min_segy) % xl_step_segy;
  else if (xl0_segy > xl_max_segy)
    xl0_segy += (xl_max_segy - xl0_segy) % xl_step_segy;

  NRLib::Grid2D<A> surface_grid(n_il, n_xl, static_cast<A>(NRLib::MULT_IRAP_MISSING));
  int grid_i, grid_j;
  for (int i = 0; i < n_il; i++) {
    for (int j = 0; j < n_xl; j++) {

      //Global IL/XL of surface_grid
      int il_glob = il_min_segy + i*d_il_file;
      int xl_glob = xl_min_segy + j*d_xl_file;

      //Get corresponding IL/XL from file_grid
      int il_loc_file = (il_glob - il_min_file)/d_il_file;
      int xl_loc_file = (xl_glob - xl_min_file)/d_xl_file;

      //If surface is smaller than segy-grid, we set is as missing
      A z;
      if (il_loc_file < 0 || il_loc_file > ni_file-1 || xl_loc_file < 0 || xl_loc_file > nj_file-1)
        z = missing;
      else
        z = ilxl_grid_file_z(il_loc_file, xl_loc_file);

      //Fill in correct corner
      if (il0_segy == il_min_segy && xl0_segy == xl_min_segy) {
        grid_i = i;
        grid_j = j;
      }
      else if (il0_segy == il_max_segy && xl0_segy == xl_max_segy) {
        grid_i = n_il - i - 1;
        grid_j = n_xl - j - 1;
      }
      else if (il0_segy == il_min_segy && xl0_segy == xl_max_segy) {
        grid_i = i;
        grid_j = n_xl - j - 1;
      }
      else { //il0_segy == il_max_segy && xl0_segy == xl_min_segy
        grid_i = n_il - i - 1;
        grid_j = j;
      }

      surface_grid(grid_i, grid_j) = z;
    }
  }

  surface = RegularSurface<A>(x_ref, y_ref, lx, ly, surface_grid);
}


template <class A>
void NRLib::WriteIrapClassicAsciiSurf(const RegularSurface<A> & surf,
                                      double                    angle,
                                      const std::string       & filename)
{
  std::ofstream file;
  OpenWrite(file, filename);

  file << std::fixed
       << std::setprecision(6)
       << -996           << " "
       << surf.GetNJ()   << " "
       << surf.GetDX()   << " "
       << surf.GetDY()   << "\n"
       << std::setprecision(2)
       << surf.GetXMin() << " "
       << surf.GetXMax() << " "
       << surf.GetYMin() << " "
       << surf.GetYMax() << "\n"
       << surf.GetNI()   << " "
       << std::setprecision(6)
       << angle*180/NRLib::Pi << " "
       << std::setprecision(2)
       << surf.GetXMin() << " "
       << surf.GetYMin() << "\n"
       << "   0   0   0   0   0   0   0\n";

  file.precision(6);

  if (surf.GetMissingValue() == IRAP_MISSING) {
    for (size_t i = 0; i < surf.GetN(); i++) {
      file << surf(i) << " ";
      if((i+1) % 6 == 0)
       file << "\n";
    }
  }
  else {
    for (size_t i = 0; i < surf.GetN(); i++) {
      if (surf.IsMissing(surf(i)))
        file << IRAP_MISSING << " ";
      else
        file << surf(i) << " ";
      if((i+1) % 6 == 0)
       file << "\n";
    }
  }
  file.close();
}


template <class A>
void NRLib::WriteStormBinarySurf(const RegularSurface<A> & surf,
                                 const std::string       & filename)
{
  std::ofstream file;
  OpenWrite(file, filename.c_str(), std::ios::out | std::ios::binary);

  file.precision(14);

  file << "STORMGRID_BINARY\n\n"
       << surf.GetNI() << " " << surf.GetNJ() << " "
       << surf.GetDX() << " " << surf.GetDY() << "\n"
       << surf.GetXMin() << " " << surf.GetXMax() << " "
       << surf.GetYMin() << " " << surf.GetYMax() << "\n";

  if (surf.GetMissingValue() == STORM_MISSING) {
    // Purify *sometimes* claims a UMR for the call below. No-one understands why...
    WriteBinaryDoubleArray(file, surf.begin(), surf.end());
  }
  else {
    std::vector<double> data(surf.GetN());
    std::copy(surf.begin(), surf.end(), data.begin());
    std::replace(data.begin(), data.end(), surf.GetMissingValue(), static_cast<A>(STORM_MISSING));
    WriteBinaryDoubleArray(file, data.begin(), data.end());
  }
  file.close();
}

#endif // NRLIB_SURFACEIO_HPP
