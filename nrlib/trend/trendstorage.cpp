#include "trendstorage.hpp"
#include "trend.hpp"
#include "../iotools/fileio.hpp"
#include <fstream>

TrendStorage::TrendStorage()
{
}

TrendStorage::~TrendStorage()
{
}


TrendConstantStorage::TrendConstantStorage(const double value)
: mean_value_(value)
{
}

TrendConstantStorage::~TrendConstantStorage()
{
}

Trend *
TrendConstantStorage::GenerateTrend(const std::string & /*path*/,
                                    std::string       & /*errTxt*/) const
{
  Trend * trend = new ConstantTrend(mean_value_);
  return trend;
}

Trend1DStorage::Trend1DStorage(const std::string file_name,
                               const std::string reference_parameter)
: file_name_(file_name),
  reference_parameter_(reference_parameter)
{
}

Trend1DStorage::~Trend1DStorage()
{
}

Trend *
Trend1DStorage::GenerateTrend(const std::string & path,
                              std::string       & errTxt) const
{
  std::string file_name = path + file_name_;

  int file_format = getTrend1DFileFormat(file_name, errTxt);

  std::vector<double> trend_values;
  double              s_min;
  double              s_max;

  if(file_format < 0) {
    errTxt += "Invalid 1D trend file\n";
    return(0);
  }
  else {
    readTrend1D(file_name,errTxt,trend_values,s_min,s_max);
    Trend * trend = new Trend1D(trend_values,s_min,s_max);
    return trend;
  }
}

void Trend1DStorage::readTrend1D(const std::string   & file_name,
                                 std::string         & errText,
                                 std::vector<double> & trend1d,
                                 double              & s_min,
                                 double              & s_max)    const
{
  std::ifstream file;
  NRLib::OpenRead(file,file_name);
  std::string dummyStr;
  bool lineIsComment = true;
  int  line          = 0;
  int  thisLine      = 0;

  while( lineIsComment == true) {

    if(NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }

    NRLib::ReadNextToken(file,dummyStr,line);
    if (line == thisLine)
      NRLib::DiscardRestOfLine(file,line,false);
    thisLine = line;
    if((dummyStr[0]!='*') &  (dummyStr[0]!='"')) {
      lineIsComment = false;
    }
  }

  s_min = NRLib::ParseType<double>(dummyStr);

  if (NRLib::CheckEndOfFile(file))  {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  NRLib::ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    NRLib::DiscardRestOfLine(file,line,false);
  thisLine = line;

  double dz = NRLib::ParseType<double>(dummyStr);

  if (NRLib::CheckEndOfFile(file)) {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  NRLib::ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    NRLib::DiscardRestOfLine(file,line,false);
  thisLine = line;

  int nz = NRLib::ParseType<int>(dummyStr);

  s_max = s_min+dz*nz;

  trend1d.resize(nz);

  for(int i=0; i<nz; i++) {
    if (NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }
    NRLib::ReadNextToken(file,dummyStr,line);

    trend1d[i] = NRLib::ParseType<double>(dummyStr);
  }
  file.close();
}

int
Trend1DStorage::getTrend1DFileFormat(const std::string & file_name,
                                     std::string       & errText) const
{
  std::string   dummyStr;
  std::string   targetString;
  std::ifstream file;

  // test for jason file format
  NRLib::OpenRead(file,file_name);

  int  fileformat    = -1;
  int  line          = 0;
  int  thisLine      = 0;
  bool lineIsComment = true;

  while (lineIsComment == true) {
    NRLib::ReadNextToken(file,dummyStr,line);
    if (NRLib::CheckEndOfFile(file)) {
      errText += "End of wavelet file "+file_name+" is premature\n";
      return 0;
    }
    else {
      if (thisLine == line) {
        NRLib::DiscardRestOfLine(file,line,false);
        thisLine = line;
      }
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
        lineIsComment = false;
    }
  }
  file.close();
  if (NRLib::IsNumber(dummyStr)) // not convertable number
    fileformat = 0; //Same file format as WAVELET::JASON

  return fileformat;
}

Trend2DStorage::Trend2DStorage(const std::string file_name,
                               const std::string reference_parameter1,
                               const std::string reference_parameter2)
: file_name_(file_name),
  reference_parameter_one_(reference_parameter1),
  reference_parameter_two_(reference_parameter2)
{
}

Trend2DStorage::~Trend2DStorage()
{
}

Trend *
Trend2DStorage::GenerateTrend(const std::string & /*path*/,
                              std::string       & errTxt) const
{
  errTxt += "2D-trend is not implemented yet\n";

  return(0);
}
