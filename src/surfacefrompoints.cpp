/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

//#include "libs/nrlib/surface/regularsurface.hpp"
#include "src/surfacefrompoints.h"

SurfaceFromPoints::SurfaceFromPoints(const std::string & filename)
  : surface_(NULL)
{
  std::vector<double> x;     // x-coordinate
  std::vector<double> y;     // y-coordinate
  std::vector<double> twt;   // two-way-time
  std::vector<int>    il;    // inline
  std::vector<int>    xl;    // cross-line

  ReadFile(filename, x, y, twt, il, xl);

  CreateGriddedSurface(surface_, x, y, twt, il, xl);
}

SurfaceFromPoints::~SurfaceFromPoints(void)
{
  // surface_ is not deleted here and must be deleted by caller
}

void SurfaceFromPoints::ReadFile(const std::string   & filename,
                                 std::vector<double> & x,
                                 std::vector<double> & y,
                                 std::vector<double> & twt,
                                 std::vector<int>    & il,
                                 std::vector<int>    & xl)
{
  int N = 9000000; // Larger than most normal surfaces (3000 x 3000)
  x.reserve(N);
  y.reserve(N);
  twt.reserve(N);
  il.reserve(N);
  xl.reserve(N);

  unsigned int nCols = 5; // Required number of columns in file

  std::ifstream file;
  NRLib::OpenRead(file, filename);

  int                      line = 1;
  std::string              dummy;
  std::string              errTxt;
  std::vector<std::string> tokenLine;

  while (!NRLib::CheckEndOfFile(file))
  {
    std::getline(file, dummy);
    tokenLine = NRLib::GetTokens(dummy);

    if (tokenLine.size() == nCols) {
      double d1 = NRLib::ParseType<double>(tokenLine[0]);
      double d2 = NRLib::ParseType<double>(tokenLine[1]);
      double d3 = NRLib::ParseType<double>(tokenLine[2]);
      double d4 = NRLib::ParseType<double>(tokenLine[3]);
      double d5 = NRLib::ParseType<double>(tokenLine[4]);

      int    i4 = static_cast<int>(d4);
      int    i5 = static_cast<int>(d5);

      x.  push_back( d1 );
      y.  push_back( d2 );
      twt.push_back( d3 );
      xl. push_back( i4 ); // Define XL as coming first
      il. push_back( i5 ); // Define IL as coming second
    }
    else {
      errTxt += "ERROR: There are " + NRLib::ToString(tokenLine.size()) + " columns for line ";
      errTxt +=  NRLib::ToString(line) + ". There should have been 5.\n";
      break;
    }
    line++;
  }
  file.close();

  // No need to resize arrays here as they will be soon tossed away.

  if (errTxt != "") {
    throw NRLib::FileFormatError(errTxt);
  }
}


void SurfaceFromPoints::CreateGriddedSurface(Surface                   *& surface,
                                             const std::vector<double>  & x,
                                             const std::vector<double>  & y,
                                             const std::vector<double>  & twt,
                                             const std::vector<int>     & il,
                                             const std::vector<int>     & xl)
{
  int il_min        = il[0];
  int il_second_min = -1;
  int il_max        = il[0];
  int xl_min        = xl[0];
  int xl_second_min = -1;
  int xl_max        = xl[0];

  //
  // Find minimum and maximum IL- and XL-values
  //
  for (size_t i = 1; i < il.size(); ++i) {
    if (il[i] < il_min) {
      il_second_min = il_min;
      il_min        = il[i];
    }
    else if (il[i] > il_min && (il_second_min == -1 || il[i] < il_second_min)) {
      il_second_min = il[i];
    }
    if (il[i] > il_max) {
      il_max        = il[i];
    }
    if (xl[i] < xl_min) {
      xl_second_min = xl_min;
      xl_min        = xl[i];
    }
    else if (xl[i] > xl_min && (xl_second_min == -1 || xl[i] < xl_second_min)) {
      xl_second_min = xl[i];
    }
    if (xl[i] > xl_max) {
      xl_max        = xl[i];
    }
  }

  int il_inc = 1;
  if (il_second_min != -1) {
    il_inc = il_second_min - il_min;
  }

  int xl_inc = 1;
  if (xl_second_min != -1) {
    xl_inc = xl_second_min - xl_min;
  }

  int n_il = (il_max - il_min) / il_inc + 1;
  int n_xl = (xl_max - xl_min) / xl_inc + 1;

  //
  // Find dimension and orientation of the possibly rotated IL-XL-regular-grid
  //
  double x11 = x[0];     // Assuming IL/XL start at first entry in file;
  double x12 = RMISSING;
  double x21 = RMISSING;
  double x22 = RMISSING;
  double y11 = y[0];     // Assuming IL/XL start at first entry in file;
  double y12 = RMISSING;
  double y21 = RMISSING;
  double y22 = RMISSING;

  for (size_t i = 0; i < il.size(); ++i) {
    if (il[i] == il_max && xl[i] == xl_min) {
      x21 = x[i];
      y21 = y[i];
    }
    else if (il[i] == il_min && xl[i] == xl_max) {
      x12 = x[i];
      y12 = y[i];
    }
    else if (il[i] == il_max && xl[i] == xl_max) {
      x22 = x[i];
      y22 = y[i];
    }
  }

  int d_il = il_max - il_min;
  int d_xl = xl_max - xl_min;

  double x0;                          // Rotation reference point X
  double y0;                          // Rotation reference point Y

  double dx_il;
  double dy_il;
  double dy_xl;

  // IL has been associated with X since this was the case for the
  // first test surface that was received.

  if (x21 > x11 && y12 > y11) {       // Normal case
    x0    = x11;
    y0    = y11;
    dx_il = (x21 - x11)/d_il;
    dy_il = (y21 - y11)/d_il;
    dy_xl = (y12 - y11)/d_xl;
  }
  else if (x21 < x11 && y12 > y11)  { // Swap IL information
    x0    = x21;
    y0    = y21;
    dx_il = (x11 - x21)/d_il;
    dy_il = (y11 - y21)/d_il;
    dy_xl = (y12 - y11)/d_xl;
  }
  else if (x21 > x11 && y12 < y11)  { // Swap XL information
    x0    = x12;
    y0    = y12;
    dx_il = (x21 - x11)/d_il;
    dy_il = (y21 - y11)/d_il;
    dy_xl = (y11 - y12)/d_xl;
    std::cout << "1: Nontested surface ASCII XYZ format configuration\n";
  }
  else {                              // Swap IL and XL information
    x0    = x22;
    y0    = y22;
    dx_il = (x11 - x21)/d_il;
    dy_il = (y11 - y21)/d_il;
    dy_xl = (y11 - y12)/d_xl;
    std::cout << "2: Nontested surface ASCII XYZ format configuration\n";
  }

  double rot  = std::atan(dy_il/dx_il);

  double cosR = std::cos(rot);
  double dx   = dx_il*il_inc/cosR;
  double dy   = dy_xl*xl_inc/cosR;

  double xlen = dx_il*d_il;
  double ylen = dy_xl*d_xl;

  /* NBNB-PAL: Keep for a little while
  std::cout << x11 << " " << y11 << std::endl;
  std::cout << x21 << " " << y21 << std::endl;
  std::cout << x12 << " " << y12 << std::endl;
  std::cout << x22 << " " << y22 << "\n" << std::endl;

  std::cout <<" xlen = " << xlen << std::endl;
  std::cout <<" ylen = " << ylen << "\n" << std::endl;

  std::cout <<" dx_il = " << dx_il << std::endl;
  std::cout <<" dy_il = " << dy_il << std::endl;
  std::cout <<" dy_xl = " << dy_xl << "\n" << std::endl;

  std::cout <<" rot     = " << 360*(rot/6.28)     << "\n" << std::endl;
  */

  //
  // Fill IL-XL surface
  //
  RotatedSurface ilxl_surf(x0, y0, xlen, ylen, n_il, n_xl, rot, RMISSING);

  for (size_t k = 0; k < twt.size(); ++k) {
    int i = (il[k] - il_min) / il_inc;
    int j = (xl[k] - xl_min) / xl_inc;

    if (x21 < x11)
      i = n_il - i - 1;  // X and IL runs the opposite way
    if (y12 < y11)
      j = n_xl - j - 1;  // Y and XL runs the opposite way

    ilxl_surf(i, j) = twt[k];
  }
  //ilxl_surf.WriteToFile("ilxl_surf.irap", NRLib::SURF_IRAP_CLASSIC_ASCII); // For debugging

  double du   = 6.25;                  // Seismic data "unit"
  double dxu  = std::abs(du*floor(dx/du + 0.5));
  double dyu  = std::abs(du*floor(dy/du + 0.5));

  double xMin = *(std::min_element(x.begin(), x.end()));
  double xMax = *(std::max_element(x.begin(), x.end()));
  double yMin = *(std::min_element(y.begin(), y.end()));
  double yMax = *(std::max_element(y.begin(), y.end()));

  double xLen = dxu*ceil((xMax - xMin)/dxu);
  double yLen = dyu*ceil((yMax - yMin)/dyu);

  int    nx   = static_cast<int>(xLen/dxu) + 1;
  int    ny   = static_cast<int>(yLen/dyu) + 1;

  surface = new Surface(xMin, yMin, xLen, yLen, nx, ny, RMISSING);

  for (size_t i = 0 ; i < surface->GetNI() ; ++i) {
    for (size_t j = 0 ; j < surface->GetNJ() ; ++j) {
      double xi = xMin + i*dxu;
      double yj = yMin + j*dyu;
      (*surface)(i, j) = ilxl_surf.GetZ(xi, yj);
    }
  }
  //surface->WriteToFile("surface.irap", NRLib::SURF_IRAP_CLASSIC_ASCII);
}
