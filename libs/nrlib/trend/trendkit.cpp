// $Id: trendkit.cpp 1075 2012-09-19 13:42:16Z georgsen $
#include "trendkit.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/regularsurface.hpp"
#include <fstream>

namespace NRLib {

void ReadTrend1DJason(const std::string   & file_name,
                      std::string         & errText,
                      std::vector<double> & trend1d,
                      double              & s_min,
                      double              & dz)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  bool lineIsComment = true;
  int  line          = 0;
  int  thisLine      = 0;

  while( lineIsComment == true) {

    if(CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }

    ReadNextToken(file,dummyStr,line);
    if (line == thisLine)
      DiscardRestOfLine(file,line,false);
    thisLine = line;
    if((dummyStr[0]!='*') &  (dummyStr[0]!='"')) {
      lineIsComment = false;
    }
  }

  s_min = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file))  {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);
  thisLine = line;

  dz = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file)) {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);
  thisLine = line;

  int nz = ParseType<int>(dummyStr);

  trend1d.resize(nz);

  for(int i=0; i<nz; i++) {
    if (CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }
    ReadNextToken(file,dummyStr,line);

    trend1d[i] = ParseType<double>(dummyStr);
  }
  file.close();
}

//----------------------------------------------------//

void ReadTrend1DPlainAscii(const std::string   & file_name,
                           std::string         & /*errText*/,
                           std::vector<double> & trend1d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  int line = 0;

  while(CheckEndOfFile(file)==false){
    ReadNextToken(file,dummyStr,line);
    trend1d.push_back(ParseType<double>(dummyStr));
  }
  file.close();

}

//----------------------------------------------------//

int
GetTrend1DFileFormat(const std::string & file_name,
                     std::string       & errText)
{
  std::string   dummyStr;
  std::string   targetString;
  std::ifstream file;

  // test for jason file format
  OpenRead(file,file_name);

  int  fileformat    = -1;
  int  line          = 0;
  int  thisLine      = 0;
  bool lineIsComment = true;
  bool commentFound = false;

  while (lineIsComment == true) {
    ReadNextToken(file,dummyStr,line);
    if (CheckEndOfFile(file)) {
      errText += "End of trend file "+file_name+" is premature\n";
      return 0;
    }
    else {
      if (thisLine == line) {
        DiscardRestOfLine(file,line,false);
        thisLine = line;
      }
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
        lineIsComment = false;
      else
        commentFound = true;
    }
  }
  file.close();
  if (IsNumber(dummyStr) && commentFound == true) // not convertable number
    fileformat = 0; //Same file format as WAVELET::JASON
  else
    fileformat = 1; // plain ascii file
  return fileformat;
}

//----------------------------------------------------//

void InterpolateTrend1DValues(const double & xi,
                              const double & xi1,
                              const double & fxi,
                              const double & fxi1,
                              const double & yj,
                              double       & fyj)
{
  double t;

  t   = (yj-xi1)/(xi-xi1);

  fyj = fxi*t + fxi1*(1-t);

}

//----------------------------------------------------//

void ResampleTrend1D(const std::vector<double> & x,
                     const std::vector<double> & fx,
                     const std::vector<double> & y,
                     std::vector<double>       & fy)
{
  // Resample fx with sampling x into fy with sampling y

  int nx = static_cast<int>(x.size());
  int ny = static_cast<int>(y.size());

  int i=0;
  int j=0;

  while(i<nx-1) {

    while(j<ny && y[j]>=x[i] && y[j]<x[i+1]) {

      InterpolateTrend1DValues(x[i],
                               x[i+1],
                               fx[i],
                               fx[i+1],
                               y[j],
                               fy[j]);
      j++;
    }

    i++;
  }
}

//----------------------------------------------------//

RegularSurface<double>
ResampleTrend2D(const RegularSurface<double> & surface,
                const std::vector<double>    & x,
                const std::vector<double>    & y,
                const bool                   & transpose)
{
  int length_x = static_cast<int>(x.size());
  int length_y = static_cast<int>(y.size());

  Grid2D<double> resampled_grid(length_x, length_y);

  for(int i=0; i<length_x; i++) {
    for(int j=0; j<length_y; j++) {
      if(transpose)
        resampled_grid(i,j) = surface.GetZ(y[j],x[i]);
      else
        resampled_grid(i,j) = surface.GetZ(x[i],y[j]);
    }
  }

  double x0 = x[0];
  double y0 = y[0];

  RegularSurface<double> resampled_surface(x0, y0, length_x, length_y, resampled_grid);

  return(resampled_surface);
}

//----------------------------------------------------//

void ReadTrend2DPlainAscii(const std::string     & file_name,
                           std::string           & err_txt,
                           NRLib::Grid2D<double> & trend2d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  int line = 0;

  int ni = ReadNext<int>(file, line);
  int nj = ReadNext<int>(file, line);
  trend2d.Resize(ni, nj);

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      if (!CheckEndOfFile(file))
        trend2d(i,j) = ReadNext<double>(file, line);
      else
        err_txt += "Premature end of file for 2D trend: " + file_name + ".\n";
    }
  }

  file.close();
}


}

