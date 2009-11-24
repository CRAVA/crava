#include <fstream>
#include <math.h>

#include "nrlib/iotools/fileio.hpp"

#include "lib/global_def.h"

#include "src/waveletfilter.h"

using namespace NRLib;

WaveletFilter::WaveletFilter(const std::string & fileName,
                             int               & errCode,
                             char              * errText)
{
  readFile(fileName, errCode, errText);
}

WaveletFilter::~WaveletFilter(void)
{
}


bool WaveletFilter::readFile(const std::string & fileName,
                             int               & errCode,
                             char              * errText)
{
/*  std::ifstream inFile;
  NRLib::OpenRead(inFile, fileName);

  int line = 1;
  resPhi_ = GetNext<unsigned int>(inFile, line);
  DiscardRestOfLine(inFile, line, false);
//  inFile >> resPhi_;
  resPsi_ = GetNext<unsigned int>(inFile, line);
  DiscardRestOfLine(inFile, line, false);
  //  inFile >> resPsi_;
  if (resPhi_ < 1 || resPhi_ > 6) {
    sprintf(errText,"%sResolution for phi in filter file %s must be in [1,6].\n",errText,fileName.c_str()); 
    errCode = 1;
  }
  if (resPsi_ < 1 || resPsi_ > 6 || resPhi_ == 4) {
    sprintf(errText,"%sResolution for psi in filter file %s must be in {1,2,3,5,6}.\n",errText,fileName.c_str()); 
    errCode = 1;
  }

  unsigned int nPhi = 360 / resPhi_;
  unsigned int nPsi = 1 + (90 / resPsi_);

  std::string token;

  ampScale_.Resize(nPhi, nPsi, 1.0);
  try {
    for (unsigned int i=0; i<nPhi; i++)
      for (unsigned int j=0; j<nPsi; j++) {
        GetNextToken(inFile, token, line);
        if (token == "-999") {
          sprintf(errText,"%sToo few values for Amplitude Scaling given in file %s. There should be %f values.\n", 
            errText, fileName.c_str(), nPhi*nPsi);
          errCode = 1;
          return(false);
        }
        else
          ampScale_(i,j) = ParseType<double>(token);
      }
  }
  catch (EndOfFile&) {
    sprintf(errText,"%sUnexcpected end of file found while parsing %s.\n", errText, fileName.c_str());
    errCode = 1;
    return(false);
  }

  GetNextToken(inFile, token, line);
  if (token != "-999") {
    sprintf(errText, "%sParameters must be seperated by -999 in file %s.\n", errText, fileName.c_str());
    errCode = 1;
    return(false);
  }
  DiscardRestOfLine(inFile, line, false);
  
  stretch_.Resize(nPhi, nPsi, 1.0);
  try {
    for (unsigned int i=0; i<nPhi; i++)
      for (unsigned int j=0; j<nPsi; j++) {
        GetNextToken(inFile, token, line);
        if (token == "-999") {
          sprintf(errText,"%sToo few values for Wavelet Stretching given in file %s. There should be %f values.\n", 
            errText, fileName.c_str(), nPhi*nPsi);
          errCode = 1;
          return(false);
        }
        else
          stretch_(i,j) = ParseType<double>(token);
      }
  }
  catch (EndOfFile&) {
    sprintf(errText,"%sUnexcpected end of file found while parsing %s.\n", errText, fileName.c_str());
    errCode = 1;
    return(false);
  }

  GetNextToken(inFile, token, line);
  if (token != "-999") {
    sprintf(errText, "%sParameters must be seperated by -999 in file %s.\n", errText, fileName.c_str());
    errCode = 1;
    return(false);
  }
  DiscardRestOfLine(inFile, line, false);
  
  coverage_.Resize(nPhi, nPsi, 1.0);
  try {
    for (unsigned int i=0; i<nPhi; i++)
      for (unsigned int j=0; j<nPsi; j++) {
        GetNextToken(inFile, token, line);
        if (token == "-999") {
          sprintf(errText,"%sToo few values for Survey Coverage given in file %s. There should be %f values.\n", 
            errText, fileName.c_str(), nPhi*nPsi);
          errCode = 1;
          return(false);
        }
        else
          coverage_(i,j) = ParseType<double>(token);
      }
  }
  catch (EndOfFile&) {
    sprintf(errText,"%sUnexcpected end of file found while parsing %s.\n", errText, fileName.c_str());
    errCode = 1;
    return(false);
  }

  inFile.close();
*/
  return(true);
}

