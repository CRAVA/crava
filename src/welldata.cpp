/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/well/norsarwell.hpp"
#include "nrlib/well/laswell.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/welldata.h"
#include "src/modelsettings.h"
#include "src/io.h"

//----------------------------------------------------------------------------
WellData::WellData(const std::string              & wellFileName,
                   const std::vector<std::string> & logNames,
                   const std::vector<bool>        & inverseVelocity,
                   ModelSettings                  * modelSettings,
                   int                              indicatorFacies,
                   int                              indicatorFiltering,
                   int                              indicatorWavelet,
                   int                              indicatorBGTrend,
                   int                              indicatorRealVs,
                   bool                             porosityLogGiven,
                   bool                             faciesLogGiven)
  : modelSettings_(modelSettings),
    wellname_(""),
    wellfilename_(""),
    xpos_(NULL),
    ypos_(NULL),
    zpos_(NULL),
    md_(NULL),
    alpha_(NULL),
    beta_(NULL),
    rho_(NULL),
    facies_(NULL),
    alpha_background_resolution_(NULL),
    beta_background_resolution_(NULL),
    rho_background_resolution_(NULL),
    alpha_seismic_resolution_(NULL),
    beta_seismic_resolution_(NULL),
    rho_seismic_resolution_(NULL),
    blockedLogsOrigThick_(NULL),
    blockedLogsConstThick_(NULL),
    blockedLogsExtendedBG_(NULL),
    useForFaciesProbabilities_(indicatorFacies),
    useForFiltering_(indicatorFiltering),
    useForWaveletEstimation_(indicatorWavelet),
    useForBackgroundTrend_(indicatorBGTrend),
    realVsLog_(indicatorRealVs)
{
  errTxt_="";
  if(wellFileName.find(".nwh",0) != std::string::npos)
    readNorsarWell(wellFileName, logNames, inverseVelocity, porosityLogGiven, faciesLogGiven);
  else if(wellFileName.find(".las",0) != std::string::npos)
    readLASWell(wellFileName, logNames, inverseVelocity, porosityLogGiven, faciesLogGiven);
  else
    readRMSWell(wellFileName, logNames, inverseVelocity, porosityLogGiven, faciesLogGiven);
}

//----------------------------------------------------------------------------
WellData::~WellData()
{
  delete [] xpos_;
  delete [] ypos_;
  delete [] zpos_;
  delete [] alpha_;
  delete [] beta_;
  delete [] rho_;
  delete [] facies_;
  delete [] md_;
  delete [] porosity_;

  if (nFacies_ > 0)
    delete [] faciesNr_;
  if(alpha_background_resolution_!=NULL)
    delete [] alpha_background_resolution_;
  if(beta_background_resolution_!=NULL)
    delete [] beta_background_resolution_;
  if(rho_background_resolution_!=NULL)
    delete [] rho_background_resolution_;
  if(alpha_seismic_resolution_!=NULL)
    delete [] alpha_seismic_resolution_;
  if(beta_seismic_resolution_!=NULL)
    delete [] beta_seismic_resolution_;
  if(rho_seismic_resolution_!=NULL)
    delete [] rho_seismic_resolution_;
  if (blockedLogsOrigThick_ != NULL)
    delete blockedLogsOrigThick_;
  if (blockedLogsConstThick_ != NULL)
    delete blockedLogsConstThick_;
  if (blockedLogsExtendedBG_ != NULL)
    delete blockedLogsExtendedBG_;
}

//----------------------------------------------------------------------------
void
WellData::readRMSWell(const std::string              & wellFileName,
                      const std::vector<std::string> & logNames,
                      const std::vector<bool>        & inverseVelocity,
                      bool                             porosityLogGiven,
                      bool                             faciesLogGiven)
{
  error_ = 0;
  std::ifstream file;
  NRLib::OpenRead(file, wellFileName);
  int i,j, facies;
  std::string token, dummyStr;
  std::vector<std::string> tokenLine;
  double xpos, ypos, zpos;
  double dummy = RMISSING;
  float alpha, beta, rho, porosity;
  wellfilename_ = wellFileName;
  if(file == 0)
  {
    error_ = 1;
    errTxt_ += "Could not open file "+wellFileName+" for reading.\n";
    //NBNB Incomplete solution, but should never happen.
  }
  int nlog; // number of logs in file
  int line = 0;
  NRLib::DiscardRestOfLine(file,line,false); //First two lines contain info we do not need.
  NRLib::DiscardRestOfLine(file,line,false);
  NRLib::ReadNextToken(file, token, line);
  wellname_ = token;
  xpos0_ = NRLib::ReadNext<double>(file, line);
  ypos0_ = NRLib::ReadNext<double>(file, line);
  NRLib::DiscardRestOfLine(file,line,false);
  nlog   = NRLib::ReadNext<int>(file, line);

  //Start searching for key words.

  int nVar = 4;       // z,alpha,beta,rho

  std::vector<std::string> parameterList(logNames.size());

  bool vpLog = false;
  bool vsLog = false;

  if(logNames[0] != "") // Assume that all lognames are filled present if first is.
  {
    parameterList = logNames;
    if (faciesLogGiven)
      nVar = 5; // add facies
    if (porosityLogGiven)
      nVar = 6; // add porosity
    vpLog = !inverseVelocity[0];
    vsLog = !inverseVelocity[1];
  }
 else
  {
    parameterList[0] = "TWT";
    parameterList[1] = "DT";
    parameterList[2] = "RHOB";
    parameterList[3] = "DTS";
    parameterList[4] = "FACIES";
    nVar++; // add facies
  }
  int * pos = new int[nVar];
  for(i=0;i<nVar;i++)
    pos[i] = IMISSING;

  nFacies_ = 0;
  for(i=0;i<nlog;i++)
  {
    NRLib::ReadNextToken(file,token,line);
    for(j=0;j<nVar;j++)
    {
      if( NRLib::Uppercase(token)==parameterList[j])
      {
        pos[j] = i + 4;
        if(j==4 && faciesLogGiven)
        {
          faciesLogName_ = parameterList[4];
          // facies log - save names
          NRLib::ReadNextToken(file,token,line); // read code word DISC
          if (token != "DISC")
          {
            LogKit::LogFormatted(LogKit::Error,"ERROR: Facies log must be discrete.\n");
            exit(1);
          }
          // Find number of facies
          std::getline(file,dummyStr);
          tokenLine = NRLib::GetTokens(dummyStr);
          nFacies_ = static_cast<int>(tokenLine.size())/2;
        }
      }
    }
    if (token != "DISC")
      NRLib::DiscardRestOfLine(file,line,false);
  }

  std::string missVar = "";
  for(i=0 ; i<nVar ; i++)
  {
    if(pos[i]==IMISSING)
    {
      missVar += parameterList[i];
      error_ = 1;
    }
  }
  if(error_ > 0)
    errTxt_ += "Cannot find log(s) "+missVar+" in well file "+wellfilename_+".\n";


  if(pos[0]==IMISSING)
    timemissing_ = 1;
  else
    timemissing_ = 0;

  // Find nd_, the number of observations in well.
  // Count the number of time observations which is not missing values.

  int nData = 0;
  int legalData = 0;
  while (NRLib::CheckEndOfFile(file)==false)
  {
    nData++;
    std::getline(file,dummyStr);
    tokenLine = NRLib::GetTokens(dummyStr);

    int nrec = tokenLine.size();
    if (nrec != 3 + nlog) {
      std::string text;
      text += std::string("\nERROR: Reading of well \'") + wellFileName + "\' failed for log record ";
      text += NRLib::ToString(nData) + " (not counting header lines).\n";
      text += std::string("Incorrect number of logs found ("+NRLib::ToString(nrec-3)+"). Correct number is "+NRLib::ToString(nlog)+". The record is\n");
      text += dummyStr+"\n";
      LogKit::LogMessage(LogKit::Error,text);
      exit(1);
    }

    double timeValue;
    try {
      timeValue = NRLib::ParseType<double>(tokenLine[3]);
    }
    catch (NRLib::Exception & e) {
      std::string text;
      text += std::string("\nERROR: Reading of well \'") + wellFileName + "\' failed for log record ";
      text += NRLib::ToString(nData) + " (not counting header lines).\n";
      text += std::string("\nERROR message is \'") + e.what() + "\'";
      LogKit::LogMessage(LogKit::Error,text);
      exit(1);
    }
    if(timeValue != WELLMISSING) {
      legalData++;   // Found legal TIME variable
    }
  }
  file.close();
  file.clear();
  nd_ = legalData;
  //
  // Check that the number of logs found for each log entry agrees
  // with the number of logs specified in header.
  //
  // A nicer and faster implementation for this is requested...
  //

  int logEntry = 0;
  NRLib::OpenRead(file, wellFileName);
  line = 0;
  for(i=0;i<4+nlog;i++)
    NRLib::DiscardRestOfLine(file,line,false);
  while (NRLib::CheckEndOfFile(file)==false && error_==0)
  {
    logEntry++;
    int  elements = 0;
    bool lastIsBlank = true;
    std::getline(file,token);

    int l = static_cast<int>(token.length());
    for (i=0; i<l; i++)
    {
      if (token[i] != ' ' && token[i] != '\t' && token[i] != '\r' && token[i] != '\0')
      {
        if (lastIsBlank)
          elements++;
        lastIsBlank = false;
      }
      else
        lastIsBlank = true;
    }
  if(elements != nlog+3)
  {
    error_ = 1;
    errTxt_ += "ERROR for well "+wellname_+": The number of log elements (nlogs="+NRLib::ToString(elements-3)+")\n in line "+NRLib::ToString(logEntry)+" does not match header specifications (nlogs="+NRLib::ToString(nlog)+").\n";
   }
  }
  file.close();
  file.clear();

  //
  // Read logs
  //
  int k;
  if (nFacies_ > 0)
    faciesNr_    = new int[nFacies_];

  NRLib::OpenRead(file, wellFileName);
  line = 0;
  for(i=0;i<4+nlog;i++)
  {
    NRLib::ReadNextToken(file,token,line);
    if (NRLib::Uppercase(token) == parameterList[4])
    {
      NRLib::ReadNextToken(file,token,line); // read code word DISC
      // facies types given here
      for(k=0;k<nFacies_;k++)
      {
        NRLib::ReadNextToken(file,token,line);
        faciesNr_[k] = NRLib::ParseType<int>(token);
        NRLib::ReadNextToken(file,token,line);
        faciesNames_.push_back(token);
      }
    }
    NRLib::DiscardRestOfLine(file,line,false);
  }
  double OPENWORKS_MISSING = -999.25;
  bool wrongMissingValues = false;
  // Read data into variables.
  xpos_     = new double[nd_];
  ypos_     = new double[nd_];
  zpos_     = new double[nd_];
  alpha_    = new float[nd_];
  beta_     = new float[nd_];
  rho_      = new float[nd_];
  facies_   = new int[nd_];   // Always allocate a facies log (for code simplicity)
  porosity_ = new float[nd_];
  alpha     = RMISSING;
  beta      = RMISSING;
  rho       = RMISSING;
  zpos      = RMISSING;
  facies    = IMISSING;
  porosity  = RMISSING;
  k         = -1;
  faciesok_ = 1;
  int legal = 0;
  for(i=0;i<nData;i++)
  {
    xpos  = NRLib::ReadNext<double>(file,line);
    ypos  = NRLib::ReadNext<double>(file,line);
    dummy = NRLib::ReadNext<double>(file,line);
    for(j=4;j<=nlog+3;j++)
    {
      dummy = NRLib::ReadNext<double>(file,line);
      if(j==pos[0])
      {
        //Found TIME variable
        if(dummy != WELLMISSING && dummy != OPENWORKS_MISSING)
        {
          zpos = dummy;
          k++;
          legal = 1;
        }
        else
          legal = 0;
      }
      else if(j==pos[1])
      {
        // Found ALPHA variable
        if(dummy != WELLMISSING)
          if (vpLog)
            alpha = static_cast<float>(dummy);
          else
            alpha = static_cast<float>(304800.0/dummy);
        else
          alpha = RMISSING;
      }
      else if(j==pos[3])
      {
        // Found BETA variable
        if(dummy != WELLMISSING)
          if (vsLog)
            beta = static_cast<float>(dummy);
          else
            beta = static_cast<float>(304800.0/dummy);
        else
          beta = RMISSING;
      }
      else if(j==pos[2])
      {
        //Found RHO variable
        if(dummy != WELLMISSING)
          rho = static_cast<float>(dummy);
        else
          rho = RMISSING;
      }
      else if(nVar > 4 && j==pos[4] && faciesLogGiven)
      {
        //Found facies variable
        if(dummy != WELLMISSING)
          facies = static_cast<int>(dummy);
        else
          facies = IMISSING;
        if(facies!=IMISSING)
        {
          int faciesok = 0;
          for(int kk=0 ; kk<nFacies_ ; kk++)
          {
            if(facies == faciesNr_[kk]) {
              faciesok = 1;
              break;
            }
          }
          if(faciesok == 0)
            faciesok_ = 0;
        }
      }
      else if(nVar>5 && j == pos[5] && porosityLogGiven){
        // Porosity variable
        if(dummy != WELLMISSING)
          porosity = static_cast<float>(dummy);
        else
          porosity  = RMISSING;
      }
    }
    if(legal == 1)
    {
      xpos_[k]   = xpos;
      ypos_[k]   = ypos;
      zpos_[k]   = zpos;
      alpha_[k]  = alpha;
      beta_[k]   = beta;
      rho_[k]    = rho;
      facies_[k] = facies;
      porosity_[k] = porosity;
      if (alpha == OPENWORKS_MISSING)
        wrongMissingValues = true;
      if(beta == OPENWORKS_MISSING)
        wrongMissingValues = true;
      if(rho == OPENWORKS_MISSING)
        wrongMissingValues = true;
    }
  }
  file.close();

  if(wrongMissingValues)
  {
    LogKit::LogFormatted(LogKit::Low,"\nWARNING: There are incorrect missing values \'%f\' in well.\n",OPENWORKS_MISSING,wellfilename_.c_str());
    LogKit::LogFormatted(LogKit::Low,"           CRAVA can handle this, but it is not a legal RMS well file format.\n");
  }
  if (k != nd_-1) {
    LogKit::LogFormatted(LogKit::Low,"ERROR: Well %s is corrupt. nd=%d log entries were expected, but only k=%d were read.\n",
                     wellfilename_.c_str(),nd_,k);
    exit(1);
  }
  delete [] pos;
}


void
WellData::readNorsarWell(const std::string              & wellFileName,
                         const std::vector<std::string> & logNames,
                         const std::vector<bool>        & inverseVelocity,
                         bool                             porosityLogGiven,
                         bool                             faciesLogGiven)
{
  error_ = 0;
  wellfilename_ = wellFileName;
  std::string name = NRLib::RemovePath(wellFileName);
  name = NRLib::ReplaceExtension(name, "");
  wellname_ = name;

  try
  {
    NRLib::NorsarWell well(wellFileName);
    processNRLibWell(well, wellFileName, logNames, inverseVelocity, porosityLogGiven, faciesLogGiven, true);
    xpos0_ = well.GetXPosOrigin()*1000;
    ypos0_ = well.GetYPosOrigin()*1000;
  }
  catch (NRLib::Exception & e) {
    errTxt_ += "Error: " + NRLib::ToString(e.what());
    error_ = 1;
  }
}


void
WellData::readLASWell(const std::string              & wellFileName,
                      const std::vector<std::string> & logNames,
                      const std::vector<bool>        & inverseVelocity,
                      bool                             porosityLogGiven,
                      bool                             faciesLogGiven)
{
  error_ = 0;
  wellfilename_ = wellFileName;
  try
  {
    NRLib::LasWell well(wellFileName);
    processNRLibWell(well, wellFileName, logNames, inverseVelocity, porosityLogGiven, faciesLogGiven, false);
    if(errTxt_ == "") {
      xpos0_ = xpos_[0];
      ypos0_ = ypos_[0];
      wellname_ = well.GetWellName();
    }
  }
  catch (NRLib::Exception & e) {
    errTxt_ += "Error: " + NRLib::ToString(e.what());
    error_ = 1;
  }
}

void
WellData::processNRLibWell(const NRLib::Well              & well,
                           const std::string              & wellFileName,
                           const std::vector<std::string> & logNames,
                           const std::vector<bool>        & inverseVelocity,
                           bool                             porosityLogGiven,
                           bool                             faciesLogGiven,
                           bool                             norsarWell)
{

  int nVar = 4;       // z,alpha,beta,rho

  std::vector<std::string> parameterList;

  bool vpLog = false;
  bool vsLog = false;

  if(logNames[0] != "") // Assume that all lognames are filled present if first is.
  {
    parameterList = logNames;
    if(faciesLogGiven)
      nVar = 5;           // facies
    if(porosityLogGiven)
      nVar = 6;           // porosity (even if the facies log is not defined, nVar is set to 6)
    vpLog = !inverseVelocity[0];
    vsLog = !inverseVelocity[1];
  }
  else
  {
    parameterList[0] = "TWT";
    parameterList[1] = "DT";
    parameterList[2] = "RHOB";
    parameterList[3] = "DTS";
    parameterList[4] = "FACIES";
    nVar = 5;
  }

  int nLogs  = 2+nVar;
  int nExtra = 0;
  if(norsarWell==true)
    nExtra = 1; //MD log, needed for writing.
  std::vector<double> * filler = NULL; //to eliminate warning.
  std::vector<const std::vector<double> *> logs(nLogs+nExtra, filler);

  if(well.HasContLog("UTMX") == true)
    logs[0] = &(well.GetContLog("UTMX"));
  else if(well.HasContLog("EASTING") == true)
    logs[0] = &(well.GetContLog("EASTING"));
  else {
    error_ = 1;
    errTxt_ += "Could not find log 'UTMX' or 'EASTING' in well file "+wellFileName+".\n";
    logs[0] = NULL;
  }

  if(well.HasContLog("UTMY") == true)
    logs[1] = &(well.GetContLog("UTMY"));
  else if(well.HasContLog("NORTHING") == true)
    logs[0] = &(well.GetContLog("NORTHING"));
  else {
    error_ = 1;
    errTxt_ += "Could not find log 'UTMY' or 'NORTHING' in well file "+wellFileName+".\n";
    logs[1] = NULL;
  }

  for(int i=0;i<nVar;i++) {
    if(well.HasContLog(parameterList[i]) == false) {
      logs[2+i] = NULL;
      if(!(i == 4 || i ==5) || logNames[0] != "") { // i == 4, 5 corresponds to facies log and porosity log respectively
        error_ = 1;
        errTxt_ += "Could not find log "+parameterList[i]+" in well file "+wellFileName+".\n";
      }
      else if(i==4)
        nLogs = nLogs-1;
    }
    else
      logs[2+i] = &(well.GetContLog(parameterList[i]));
  }

  //Added MD log.
  int mdLog = nLogs;
  if(norsarWell==true) {
    if(well.HasContLog("MD") == false) {
      error_ = 1;
      errTxt_ += "Could not find log 'MD' in well file "+ wellFileName+ ".\n";
      logs[mdLog] = NULL;
    }
    else
      logs[mdLog] = &(well.GetContLog("MD"));
    nLogs++;
  }

  if(logs[2] == NULL)
    timemissing_ = 1;
  else
    timemissing_ = 0;

  if(error_ == 0) {
    faciesok_ = 1;
    std::vector<int> facCodes;
    nd_ = 0;
    for(size_t i=0;i<logs[2]->size();i++)
      if(well.IsMissing((*logs[2])[i]) == false)
        nd_++;

    xpos_     = new double[nd_];
    ypos_     = new double[nd_];
    zpos_     = new double[nd_];
    alpha_    = new float[nd_];
    beta_     = new float[nd_];
    rho_      = new float[nd_];
    facies_   = new int[nd_];   // Always allocate a facies log (for code simplicity)
    porosity_ = new float[nd_];
    if(norsarWell == true)
      md_       = new double[nd_];
    else
      md_       = NULL;

    double scale = 1.0;
    if(norsarWell == true)
      scale = 1000.0;
    int ind   = 0;
    for(size_t i=0;i<logs[0]->size();i++) {
      // check that the time variable is ok
      if(well.IsMissing((*logs[2])[i]) == false) {
        xpos_[ind]  = (*logs[0])[i]*scale;
        ypos_[ind]  = (*logs[1])[i]*scale;
        zpos_[ind]  = (*logs[2])[i]*scale;
        // vp
        if(!well.IsMissing((*logs[3])[i]))
          alpha_[ind] = static_cast<float>((*logs[3])[i]);
        else
          alpha_[ind] = RMISSING;
        // vs
        if(!well.IsMissing((*logs[5])[i]))
          beta_[ind]  = static_cast<float>((*logs[5])[i]);
        else
          beta_[ind] = RMISSING;
        // density
        if(!well.IsMissing((*logs[4])[i]))
          rho_[ind]   = static_cast<float>((*logs[4])[i]);
        else
          rho_[ind] = RMISSING;
        // facies
        if(mdLog != 6 && nLogs > 6 && logs[6] != NULL &&
          !well.IsMissing((*logs[6])[i])) {
          facies_[ind]  = static_cast<int>((*logs[6])[i]);
          if(find(facCodes.begin(), facCodes.end(), facies_[ind]) == facCodes.end())
            facCodes.push_back(facies_[ind]);
        }
        else
          facies_[ind] = IMISSING;
        // porosity
        if(porosityLogGiven && logs[7] !=NULL && mdLog!=7){
          if(!well.IsMissing((*logs[7])[i]))
            porosity_[ind] = static_cast<float>((*logs[7])[i]);
          else
            porosity_[ind] = RMISSING;
        }
        if(norsarWell == true) {
          if(logs[mdLog] != NULL && !well.IsMissing((*logs[mdLog])[i]))
            md_[ind]   = static_cast<float>((*logs[mdLog])[i]);
          else
            md_[ind] = RMISSING;
        }
        ind++;
      }
    }
    nFacies_ = static_cast<int>(facCodes.size());
    if(nFacies_ > 0) {
      faciesLogName_ = parameterList[4];
      faciesNames_.resize(nFacies_);
      faciesNr_ = new int[nFacies_];
      for(size_t i=0;i<facCodes.size();i++) {
        faciesNames_[i] = NRLib::ToString(facCodes[i]);
        faciesNr_[i] = facCodes[i];
      }
    }
  }
}


//----------------------------------------------------------------------------
void
WellData::writeWell(int wellFormat)
{
  if((wellFormat & IO::RMSWELL) > 0)
    writeRMSWell();
  if((wellFormat & IO::NORSARWELL) > 0)
    writeNorsarWell();
}

void
WellData::writeRMSWell(void)
{
  float maxHz_background = modelSettings_->getMaxHzBackground();
  float maxHz_seismic    = modelSettings_->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib::Substitute(wellname,"/","_");
  NRLib::Substitute(wellname," ","_");
  std::string baseName     = IO::PrefixWells() + wellname + IO::SuffixRmsWells();
  std::string wellFileName = IO::makeFullFileName(IO::PathToWells(), baseName);

  std::ofstream file;
  NRLib::OpenWrite(file, wellFileName);

  unsigned int nLogs = 1+3*3;
  if (nFacies_ > 0)
    nLogs += 1;

  std::vector<std::string> params(3);
  params[0] = "Vp";
  params[1] = "Vs";
  params[2] = "Rho";

  //
  // Write HEADER
  //
  file << "1.0\n"
       << "CRAVA\n"
       << wellname_ << " "
       << std::fixed
       << std::setprecision(2)
       << xpos0_ << " "
       << ypos0_ << "\n"
       << nLogs << "\n"
       << "Time  UNK lin"
       << std::endl;

  for (int i =0 ; i < 3 ; i++) {
    file << params[i] <<                                       "  UNK lin\n"
         << params[i] << static_cast<int>(maxHz_background) << "  UNK lin\n"
         << params[i] << static_cast<int>(maxHz_seismic)    << "  UNK lin\n";
  }

  if (nFacies_ > 0) {
    file << faciesLogName_ << "   DISC ";
    for (int i =0 ; i < nFacies_ ; i++)
      file << " " << faciesNr_[i] << " " << faciesNames_[i];
    file << std::endl;
  }

  for (int i=0 ; i<nd_ ; i++) {
    file << std::right
         << std::fixed
         << std::setprecision(4)
         << std::setw(9) << xpos_[i] << " "
         << std::setw(10)<< ypos_[i] << " "
         << std::setw(7) << zpos_[i] << "  "
         << std::setw(7) << zpos_[i] << "  "
         << std::setprecision(5)
         << std::setw(7) << (alpha_[i]==RMISSING                       ? WELLMISSING : alpha_[i])                       << " "
         << std::setw(7) << (alpha_background_resolution_[i]==RMISSING ? WELLMISSING : alpha_background_resolution_[i]) << " "
         << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING    ? WELLMISSING : alpha_seismic_resolution_[i])    << " "
         << std::setw(7) << (beta_[i]==RMISSING                        ? WELLMISSING : beta_[i])                        << " "
         << std::setw(7) << (beta_background_resolution_[i]==RMISSING  ? WELLMISSING : beta_background_resolution_[i])  << " "
         << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING     ? WELLMISSING : beta_seismic_resolution_[i])     << " "
         << std::setw(7) << (rho_[i]==RMISSING                         ? WELLMISSING : rho_[i])                         << " "
         << std::setw(7) << (rho_background_resolution_[i]==RMISSING   ? WELLMISSING : rho_background_resolution_[i])   << " "
         << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING      ? WELLMISSING : rho_seismic_resolution_[i])      << " ";

      if (nFacies_ > 0)
        file << std::setw(3) << (facies_[i]==IMISSING ? static_cast<int>(WELLMISSING) : facies_[i]);

    file << "\n";
  }
  file.close();
}


void
WellData::writeNorsarWell()
{
  float maxHz_background = modelSettings_->getMaxHzBackground();
  float maxHz_seismic    = modelSettings_->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib::Substitute(wellname,"/","_");
  NRLib::Substitute(wellname," ","_");

  //Handle main file.
  std::string baseName = IO::PrefixWells() + wellname + IO::SuffixNorsarWells();
  std::string fileName = IO::makeFullFileName(IO::PathToWells(), baseName);

  std::ofstream mainFile;
  NRLib::OpenWrite(mainFile, fileName);
  mainFile << std::fixed
           << std::setprecision(2);

  std::vector<double> md(nd_,0);
  if(md_ != NULL)
    md[0] = md_[0];
  double dmax = 0;
  double dmin = 1e+30;
  for(int i=1;i<nd_;i++) {
    double dx = xpos_[i]-xpos_[i-1];
    double dy = ypos_[i]-ypos_[i-1];
    double dz = zpos_[i]-zpos_[i-1];
    double d  = sqrt(dx*dx+dy*dy+dz*dz);
    if(d > dmax)
      dmax = d;
    else if(d<dmin)
      dmin = d;
    if(md_ != NULL)
      md[i] = md_[i];
    else
      md[i] = md[i-1] + d;
  }
  mainFile << "[Version information]\nVERSION 1000\nFORMAT ASCII\n\n";
  mainFile << "[Well information]\n";
  mainFile << "MDMIN      m       " << 0.0f << "\n";
  mainFile << "MDMAX      m       " << md[nd_-1] << "\n";
  mainFile << "MDMINSTEP  m       " << dmin << "\n";
  mainFile << "MDMAXSTEP  m       " << dmax << "\n";
  mainFile << "UTMX       m       " << xpos_[0] << "\n";
  mainFile << "UTMY       m       " << ypos_[0] << "\n";
  mainFile << "EKB        m       " << 0.0f << "\n";
  mainFile << "UNDEFVAL   no_unit " << WELLMISSING << "\n\n";

  mainFile << "[Well track data information]\n";
  mainFile << "NUMMD  " << nd_ << "\n";
  mainFile << "NUMPAR 4\n";
  mainFile << "MD      m\n";
  mainFile << "TWT     s\n";
  mainFile << "UTMX    m\n";
  mainFile << "UTMY    m\n\n";

  std::string logBaseName = IO::PrefixWells() + wellname + IO::SuffixNorsarLog();
  std::string logFileName = IO::makeFullFileName(IO::PathToWells(), logBaseName);
  std::string onlyName    = NRLib::RemovePath(logFileName);

  bool gotFacies = nFacies_ > 0;

  int nLogs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz}
  if (gotFacies)
    nLogs += 1;

  std::vector<std::string> params(3);
  params[0] = "VP";
  params[1] = "VS";
  params[2] = "Rho";

  std::vector<std::string> unit(3);
  unit[0] = "m/s";
  unit[1] = "m/s";
  unit[2] = "g/cm^3";

  mainFile << "[Well log data information]\n";
  mainFile << "LOGNAME log\n";
  mainFile << "IN_FILE " << onlyName << "\n";
  mainFile << "NUMPAR " << nLogs+1 << "\n"; //Also count md.
  mainFile << "NUMLINES " << nd_ << "\n";
  mainFile << "MD      m\n";
  for (int i =0 ; i<3 ; i++) {
    mainFile << params[i] << " " << unit[i] << "\n";
    mainFile << params[i] << static_cast<int>(maxHz_background) << " " << unit[i] << "\n";
    mainFile << params[i] << static_cast<int>(maxHz_seismic)    << " " << unit[i] << "\n";
  }

  if (gotFacies)
    mainFile << "FACIES no_unit\n";
  mainFile.close();


  //Write the two other files.
  std::string baseTrackName = IO::PrefixWells() + wellname + IO::SuffixNorsarTrack();
  std::string trackFileName = IO::makeFullFileName(IO::PathToWells(), baseTrackName);
  std::ofstream trackFile;
  NRLib::OpenWrite(trackFile, trackFileName);
  trackFile << std::right
            << std::fixed
            << std::setprecision(2)
            << "[NORSAR Well Track]\n";

  //Note: logFileName created above, needed in mainFile.
  std::ofstream logFile;
  NRLib::OpenWrite(logFile, logFileName);
  logFile << "[NORSAR Well Log]\n";

  for(int i = 0;i<nd_;i++) {
    trackFile << std::setw(7) << md[i] << " " << std::setw(7) << zpos_[i]
              << " " << std::setw(10)<< xpos_[i] << " " << std::setw(10)<< ypos_[i] << "\n";

    logFile   << std::right << std::fixed << std::setprecision(2)
              << std::setw(7) << md[i] << " "
              << std::setw(7) << (alpha_[i]==RMISSING                       ? WELLMISSING : alpha_[i])                       << " "
              << std::setw(7) << (alpha_background_resolution_[i]==RMISSING ? WELLMISSING : alpha_background_resolution_[i]) << " "
              << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING    ? WELLMISSING : alpha_seismic_resolution_[i])    << " "
              << std::setw(7) << (beta_[i]==RMISSING                        ? WELLMISSING : beta_[i])                        << " "
              << std::setw(7) << (beta_background_resolution_[i]==RMISSING  ? WELLMISSING : beta_background_resolution_[i])  << " "
              << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING     ? WELLMISSING : beta_seismic_resolution_[i])     << " "
              << std::setprecision(5)
              << std::setw(7) << (rho_[i]==RMISSING                         ? WELLMISSING : rho_[i])                         << " "
              << std::setw(7) << (rho_background_resolution_[i]==RMISSING   ? WELLMISSING : rho_background_resolution_[i])   << " "
              << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING      ? WELLMISSING : rho_seismic_resolution_[i])      << " ";
    if (gotFacies)
      logFile << (facies_[i]==IMISSING                                 ? static_cast<int>(WELLMISSING) : facies_[i])      << "  ";
    logFile << "\n";
  }
  trackFile.close();
  logFile.close();
}


//----------------------------------------------------------------------------
int WellData::checkError(std::string & errText)
{
  if(error_ == 0)
    return(0);
  else
  {
    errText = errTxt_;
    return(error_);
  }
}

//----------------------------------------------------------------------------
int WellData::checkSimbox(Simbox * simbox)
{
  bool insideArea = false;
  int error = 1;
  if(timemissing_==0)
  {
    for(int i=0;i<nd_;i++)
    {
      if(simbox->isInside(xpos_[i], ypos_[i]))
      {
        insideArea = true;
        //
        // Correct handling of top and base checking.
        //
        if(zpos_[i] > simbox->getTop(xpos_[i], ypos_[i]) && zpos_[i] < simbox->getBot(xpos_[i], ypos_[i]))
        {
          error = 0;
          break;
        }
      }
    }
  }
  else
    error = 0;

  if (error)
  {
    if (insideArea) {
      LogKit::LogFormatted(LogKit::Low," \nWell "+wellname_+":\n");
      LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is inside inversion area but does not hit the inversion volume)\n");
      LogKit::LogFormatted(LogKit::Low,"           (well-depth: min,max = "+NRLib::ToString(zpos_[0])+","+NRLib::ToString(zpos_[nd_-1])+")\n");
    }
    else {
      LogKit::LogFormatted(LogKit::Low," \nWell "+wellname_+":\n");
      LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is not inside inversion area)\n");
    }
  }
  return(error);
}

//----------------------------------------------------------------------------

int WellData::checkVolume(NRLib::Volume & volume) const
{
  bool insideArea = false;
  int  error      = 1;

  if(timemissing_ == 0) {

    NRLib::Surface<double> & top  = volume.GetTopSurface();
    NRLib::Surface<double> & base = volume.GetBotSurface();

    for(int i=0;i<nd_;i++) {

      if(volume.IsInside(xpos_[i], ypos_[i])) {
        insideArea = true;
        //
        // Correct handling of top and base checking.
        //
        double z_top  = top.GetZ(xpos_[i], ypos_[i]);
        double z_base = base.GetZ(xpos_[i], ypos_[i]);

        if(zpos_[i] > z_top && zpos_[i] < z_base) {
          error = 0;
          break;
        }
      }
    }
  }
  else
    error = 0;

  if (error) {
    if (insideArea) {
      LogKit::LogFormatted(LogKit::Low," \nWell "+wellname_+": ");
      LogKit::LogFormatted(LogKit::Low,"IGNORED. Well does not hit the inversion volume.\n");
    }
    else {
      LogKit::LogFormatted(LogKit::Low," \nWell "+wellname_+": ");
      LogKit::LogFormatted(LogKit::Low,"IGNORED. Well is not inside the inversion area.\n");
    }
  }
  return(error);
}

//----------------------------------------------------------------------------
bool
WellData::removeDuplicateLogEntries(const Simbox * simbox, int & nMerges)
{
  bool debug = false;
  bool monotonous = true; //Check that well does not move unreasonably much upwards.
  float maxUpwardDist = -3.0f*static_cast<float>(simbox->getdz()); //Do not allow upwards movement of more than three cells. Calibrate.

  float minMergeDist = modelSettings_->getMaxMergeDist();

  double * xpos_resampled   = new double[nd_];
  double * ypos_resampled   = new double[nd_];
  double * zpos_resampled   = new double[nd_];   // time step
  float  * alpha_resampled  = new float[nd_];
  float  * beta_resampled   = new float[nd_];
  float  * rho_resampled    = new float[nd_];
  int    * facies_resampled = new int[nd_];      // Always included (for convenience)

  int ii = 0;
  int istart = 0;                                // First element in merge
  int iend;                                      // Last element in merge
  while (istart < nd_)                           // Loop over elements in original log
  {
    if (istart == nd_ - 1)                       // Are we at last element?
      iend = istart;
    else
    {
      iend = istart + 1;                         // Start looking one element ahead
      while (zpos_[iend] - zpos_[istart] < minMergeDist && iend < nd_ - 1) {
        iend++;
        if(monotonous == true && zpos_[iend] - zpos_[istart] < maxUpwardDist)
          monotonous = false;
      }
      iend--;
    }

    bool printToScreen = debug && iend != istart;
    if (printToScreen)
      LogKit::LogFormatted(LogKit::Low,"Merge log entries %d -> %d into new entry %d\n",istart,iend,ii);

    mergeCells("time",zpos_resampled, zpos_ ,ii,istart,iend,printToScreen);
    mergeCells("x   ",xpos_resampled, xpos_ ,ii,istart,iend,printToScreen);
    mergeCells("y   ",ypos_resampled, ypos_ ,ii,istart,iend,printToScreen);
    mergeCells("Vp  ",alpha_resampled,alpha_,ii,istart,iend,printToScreen);
    mergeCells("Vs  ",beta_resampled, beta_ ,ii,istart,iend,printToScreen);
    mergeCells("Rho ",rho_resampled,  rho_  ,ii,istart,iend,printToScreen);
    if(isFaciesLogDefined())
      mergeCellsDiscrete("Facies ",facies_resampled, facies_, ii, istart, iend, printToScreen);

    if (printToScreen)
      LogKit::LogFormatted(LogKit::Low,"\n");

    istart = iend + 1;
    ii++;
  }

  nMerges = 0;
  if (ii != nd_)
  {
    nMerges = nd_-ii;
    LogKit::LogFormatted(LogKit::Low,"   Duplicate log entries merged with neighbour : %d\n",nMerges);
    nd_ = ii;
  }

  delete [] xpos_;            // Delete original logs
  delete [] ypos_;
  delete [] zpos_;
  delete [] alpha_;
  delete [] beta_;
  delete [] rho_;
  delete [] facies_;

  xpos_   = xpos_resampled;   // Make original logs point to resampled logs.
  ypos_   = ypos_resampled;
  zpos_   = zpos_resampled;
  alpha_  = alpha_resampled;
  beta_   = beta_resampled;
  rho_    = rho_resampled;
  facies_ = facies_resampled;
  return(monotonous);
}

//----------------------------------------------------------------------------
void
WellData::mergeCells(const std::string & name, double * pos_resampled, double * pos,
                     int ii, int istart, int iend, bool printToScreen)
{
  int nSample = 0;
  pos_resampled[ii] = 0.0f;

  for (int i = istart ; i < iend + 1 ; i++)
  {
    if (pos[i] != RMISSING)
    {
      pos_resampled[ii] += pos[i];
      nSample++;
      if (printToScreen)
        LogKit::LogFormatted(LogKit::Low,"%s     Old:%d   pos = %.3f\n",name.c_str(),i,pos[i]);
    }
  }
  if (nSample == 0)
    pos_resampled[ii] = RMISSING;
  else
    pos_resampled[ii] /= nSample;
  if (printToScreen)
    LogKit::LogFormatted(LogKit::Low,"%s     New:%d   pos = %.3f\n",name.c_str(),ii,pos_resampled[ii]);
}

//----------------------------------------------------------------------------
void
WellData::mergeCells(const std::string & name, float * log_resampled, float * log, int ii,
                     int istart, int iend, bool printToScreen)
{
  int nSample = 0;
  log_resampled[ii] = 0.0f;

  for (int i = istart ; i < iend + 1 ; i++)
  {
    if (log[i] != RMISSING)
    {
      log_resampled[ii] += log[i];
      nSample++;
      if (printToScreen)
        LogKit::LogFormatted(LogKit::Low,"%s     Old:%d   log = %.3f\n",name.c_str(),i,log[i]);
    }
  }
  if (nSample == 0)
    log_resampled[ii] = RMISSING;
  else
    log_resampled[ii] /= nSample;
  if (printToScreen)
    LogKit::LogFormatted(LogKit::Low,"%s     New:%d   log = %.3f\n",name.c_str(),ii,log_resampled[ii]);
}

//---------------------------------------------------------------------------

void
WellData::mergeCellsDiscrete(const std::string & name, int * log_resampled, int * log, int ii,
                     int istart, int iend, bool printToScreen)
{
  int nSample = 0;
  for (int i = istart ; i < iend + 1 ; i++)
  {
    if (log[i] != RMISSING)
    {
      //log_resampled[ii] += log[i];
      nSample++;
      if (printToScreen)
        LogKit::LogFormatted(LogKit::Low,"%s     Old:%d   log = %d\n",name.c_str(),i,log[i]);
    }
  }
  log_resampled[ii] = log[istart];
 if (printToScreen)
    LogKit::LogFormatted(LogKit::Low,"%s     New:%d   log = %d\n",name.c_str(),ii,log_resampled[ii]);

}


//----------------------------------------------------------------------------
void
WellData::setWrongLogEntriesUndefined(int & count_alpha, int & count_beta, int & count_rho)
{
  //
  // Log values outside the minimum and maximum limits given below
  // are set to RMISSING as such values are very unlikely.
  //
  bool debug = true;

  float alpha_min = modelSettings_->getAlphaMin();
  float alpha_max = modelSettings_->getAlphaMax();
  float beta_min  = modelSettings_->getBetaMin();
  float beta_max  = modelSettings_->getBetaMax();
  float rho_min   = modelSettings_->getRhoMin();
  float rho_max   = modelSettings_->getRhoMax();

  count_alpha = 0;
  count_beta  = 0;
  count_rho   = 0;

  for (int i = 0 ; i < nd_ ; i++)
  {
    if (alpha_[i] != RMISSING && (alpha_[i] < alpha_min  || alpha_[i] > alpha_max))
    {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Vp = %.2f\n",zpos_[i],alpha_[i]);
      alpha_[i] = RMISSING;
      count_alpha++;
    }
    if (beta_[i] != RMISSING && (beta_[i] < beta_min  || beta_[i] > beta_max))
    {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Vs = %.2f\n",zpos_[i],beta_[i]);
      beta_[i] = RMISSING;
      count_beta++;
    }
    if (rho_[i] != RMISSING && (rho_[i] < rho_min || rho_[i] > rho_max))
    {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Rho = %.2f\n",zpos_[i],rho_[i]);
      rho_[i] = RMISSING;
      count_rho++;
    }
  }

  if (count_alpha > 0 || count_beta || count_rho > 0)
    LogKit::LogFormatted(LogKit::Low,"   Log entries have been set undefined (Vp:%d Vs:%d Rho:%d)\n",
                     count_alpha,count_beta,count_rho);
}

//----------------------------------------------------------------------------
void
WellData::filterLogs(void)
{
  float maxHz_background = modelSettings_->getMaxHzBackground();
  float maxHz_seismic    = modelSettings_->getMaxHzSeismic();

  float  * alpha_interpolated = new float[nd_];
  float  * beta_interpolated = new float[nd_];
  float  * rho_interpolated = new float[nd_];

  float  * alpha_resampled = new float[nd_];
  float  * beta_resampled = new float[nd_];
  float  * rho_resampled = new float[nd_];

  float  * alpha_filtered = new float[nd_];
  float  * beta_filtered = new float[nd_];
  float  * rho_filtered = new float[nd_];

  double * time_resampled = new double[nd_];
  double   dt;

  alpha_background_resolution_ = new float[nd_];
  beta_background_resolution_ = new float[nd_];
  rho_background_resolution_ = new float[nd_];

  alpha_seismic_resolution_ = new float[nd_];
  beta_seismic_resolution_ = new float[nd_];
  rho_seismic_resolution_ = new float[nd_];

  //
  // Time
  //
  bool filtered = resampleTime(time_resampled, nd_, dt); //False if well not monotonous in time.

  if(filtered) {
    //
    // Vp
    //
    resampleLog(alpha_resampled, alpha_, zpos_, time_resampled, nd_, dt);         // May generate missing values
    interpolateLog(alpha_interpolated, alpha_resampled, nd_);                     // Interpolate missing values

    applyFilter(alpha_filtered, alpha_interpolated, nd_, dt, maxHz_background);
    resampleLog(alpha_resampled, alpha_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(alpha_background_resolution_, alpha_resampled, nd_);

    applyFilter(alpha_filtered, alpha_interpolated, nd_, dt, maxHz_seismic);
    resampleLog(alpha_resampled, alpha_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(alpha_seismic_resolution_, alpha_resampled, nd_);

    //
    // Vs
    //
    resampleLog(beta_resampled, beta_, zpos_, time_resampled, nd_, dt);
    interpolateLog(beta_interpolated, beta_resampled, nd_);

    applyFilter(beta_filtered, beta_interpolated, nd_, dt, maxHz_background);
    resampleLog(beta_resampled, beta_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(beta_background_resolution_, beta_resampled, nd_);

    applyFilter(beta_filtered, beta_interpolated, nd_, dt, maxHz_seismic);
    resampleLog(beta_resampled, beta_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(beta_seismic_resolution_, beta_resampled, nd_);

    //
    // Rho
    //
    resampleLog(rho_resampled, rho_, zpos_, time_resampled, nd_, dt);
    interpolateLog(rho_interpolated, rho_resampled, nd_);

    applyFilter(rho_filtered, rho_interpolated, nd_, dt, maxHz_background);
    resampleLog(rho_resampled, rho_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(rho_background_resolution_, rho_resampled, nd_);

    applyFilter(rho_filtered, rho_interpolated, nd_, dt, maxHz_seismic);
    resampleLog(rho_resampled, rho_filtered, time_resampled, zpos_, nd_, dt);
    interpolateLog(rho_seismic_resolution_, rho_resampled, nd_);
  }
  else {
    for(int i=0;i<nd_;i++) {
      alpha_background_resolution_[i] = RMISSING;
      alpha_seismic_resolution_[i]    = RMISSING;
      beta_background_resolution_[i]  = RMISSING;
      beta_seismic_resolution_[i]     = RMISSING;
      rho_background_resolution_[i]   = RMISSING;
      rho_seismic_resolution_[i]      = RMISSING;
    }
  }
  delete [] alpha_interpolated;
  delete [] beta_interpolated;
  delete [] rho_interpolated;
  delete [] alpha_resampled;
  delete [] beta_resampled;
  delete [] rho_resampled;
  delete [] alpha_filtered;
  delete [] beta_filtered;
  delete [] rho_filtered;
  delete [] time_resampled;

}

//----------------------------------------------------------------------------
bool
WellData::resampleTime(double * time_resampled,
                       int      nd,
                       double & dt)
{
  //Only resample if monotonous increasing in time.
  double time_begin = zpos_[0];
  double time_end   = zpos_[nd - 1];
  bool   monotonous = true;
  for (int i = 1 ; (i < nd && monotonous == true); i++)
    if(zpos_[i] < zpos_[i-1])
      monotonous = false;

  if(monotonous == false)
    return(false);

  //
  // Make new time scale with constant sampling density
  //

  if (time_begin != RMISSING && time_end != RMISSING)
  {
    dt = (time_end - time_begin)/(nd - 1);            // average sampling density
    for (int i = 0 ; i < nd ; i++)
      time_resampled[i] = time_begin + i*dt;
  }
  else
  {
    LogKit::LogFormatted(LogKit::Warning,"WARNING: First or last time sample is undefined. Cannot estimate average sampling density.\n");
    LogKit::LogFormatted(LogKit::Warning,"         time[first] = %12.2f\n",time_begin);
    LogKit::LogFormatted(LogKit::Warning,"         time[last]  = %12.2f\n",time_end);
    return(false);
  }

  return(true);

  //printf("i time[i] dt  time_resampled[i] dt   %d  %7.3f         %7.3f\n",0,time[0],time_resampled[0]);
  //for (unsigned int i = 1 ; i < nd ; i++)
  //{
  //  printf("i time[i] dt  time_resampled[i] dt   %d  %7.3f %.3f   %7.3f %.3f\n",i,time[i],time[i]-time[i-1],time_resampled[i],dt);
  //}
}

//----------------------------------------------------------------------------
void
WellData::resampleLog(float        * log_resampled,
                      const float  * log,
                      const double * time,
                      const double * time_resampled,
                      int            nd,
                      double         dt)
{
  bool resample_log = true;

  if (!resample_log)
  {
    for (int i = 0 ; i < nd ; i++)
    {
      log_resampled[i] = log[i];
    }
  }

  //
  // Initialise as undefined
  //
  for (int i = 0 ; i < nd ; i++)
  {
    log_resampled[i] = RMISSING;
  }

  //
  // Set end points equal to original log. Set the rest undefined
  //
  log_resampled[0]    = log[0];
  log_resampled[nd-1] = log[nd-1];

  //
  // Paals first-shot interpolation...
  //
  // Aritmethic mean of raw values in interval (i-0.5)*dt < t < (i+0.5)*dt
  //
  int j = 0;
  while (time[j] < time_resampled[1] - 0.5*dt) // Find starting position
    j++;

  for (int i = 1 ; i < nd - 1 ; i++)  // End points are already set
  {
    // Start gathering values
    float        value = 0.0f;
    unsigned int count = 0;

    while ((time[j] < time_resampled[i] + 0.5*dt) && (j<nd-1))
    {
      if (log[j] != RMISSING)
      {
        value += log[j];
        count++;
        //printf("i j value   %d %d  %7.3f\n",i,j,value);

      }
      j++;
    }
    if (count > 0)
      log_resampled[i] = value/count;
  }

  //for (unsigned int i = 0 ; i < nd ; i++)
  //{
  //  printf("i log[i] log_resampled[i]    %d  %7.3f   %7.3f\n",i,log[i],log_resampled[i]);
  //}
}

//----------------------------------------------------------------------------
void
WellData::interpolateLog(float * log_interpolated, const float *log_resampled, int nd)
{
  int i;
  for (i = 0 ; i < nd ; i++) {
    log_interpolated[i] = log_resampled[i];
  }

  // Skip leading RMISSING

  i = 0;
  while (i<nd && log_resampled[i]==RMISSING)
    i++;

  // When the log has intermediate RMISSING... use linear interpolation. Skip trailing RMISSING

  while (i<nd) {
    if (log_resampled[i]==RMISSING) {
      int last_nonmissing = i - 1;                // last defined value (e.g. ..., 2.31, -99999, -99999, ...)
      while (i<nd && log_resampled[i]==RMISSING)  //                               ^^^^
        i++;
      if (i<nd) {
        int   first_nonmissing = i;               // first defined value (e.g. ..., -99999, -99999, 2.29, ...)
        int   j0 = last_nonmissing;               //                                                ^^^^
        int   j1 = first_nonmissing;
        float l0 = log_resampled[j0];
        float l1 = log_resampled[j1];

        float a = (l1 - l0)/float(j1 - j0);

        for (int j = j0 + 1 ; j < j1 ; j++) {
          log_interpolated[j] = a*(j - j0) + l0;
        }
      }
    }
    else
      i++;
  }
  //for (int i=0 ; i<nd ; i++) {
  //  printf("i log_resampled[i] log_interpolated[i]  %d  %.3f  %.3f\n",i,log_resampled[i],log_interpolated[i]);
  //}
}

//----------------------------------------------------------------------------
void
WellData::applyFilter(float * log_filtered, float *log_interpolated, int n_time_samples,
                      double dt_milliseconds, float maxHz)
{
  //
  // Extract nonmissing part of log
  //
  int i=0;
  while (i<n_time_samples && log_interpolated[i]==RMISSING)
    i++;
  int first_nonmissing = i;
  i = n_time_samples - 1;
  while (i>0 && log_interpolated[i]==RMISSING)
    i--;
  int last_nonmissing = i;
  int n_time_samples_defined = last_nonmissing - first_nonmissing + 1;

  for(i=0 ; i < n_time_samples ; i++) {            // Initialise with RMISSING
    log_filtered[i] = RMISSING;
  }

  if (n_time_samples_defined > 0)
  {
    //
    // Setup FFT
    //
    int   nt  = 2*n_time_samples_defined;
    int   cnt = nt/2 + 1;
    int   rnt = 2*cnt;

    fftw_real*    rAmp = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
    fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);

    for (i=0 ; i<n_time_samples_defined ; i++) {          // Array to filter is made symmetric
      rAmp[i]      = log_interpolated[first_nonmissing + i];
      rAmp[nt-i-1] = rAmp[i];
    }

    //for (int i=0 ; i<nt ; i++) {
    //  printf("i=%d, log_interpolated[i]=%7.4f\n",i,rAmp[i]);
    //}

    //
    // Transform to Fourier domain
    //
    rfftwnd_plan p1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_one_real_to_complex(p1, rAmp, cAmp);
    fftwnd_destroy_plan(p1);

    //for (int i=0 ; i<cnt ; i++) {
    //  printf("i=%2d, cAmp.re[i]=%11.4f  cAmp.im[i]=%11.4f\n",i,cAmp[i].re,cAmp[i].im);
    //}

    //
    // Filter using Odd's magic vector...
    //
    float dt  = static_cast<float> (dt_milliseconds/1000.0); // Sampling density in seconds
    float T   = (nt - 1)*dt;                                 // Time sample
    float w   = 1/T;                                         // Lowest frequency that can be extracted from log
    int   N   = int(maxHz/w + 0.5f);                         // Number of elements of Fourier vector to keep

    if(cnt < N+1) {
      LogKit::LogMessage(LogKit::Warning, "Warning: The vertical resolution is too low to allow filtering of well logs to %3.1f Hz.\n");
    }

    float * magic_vector = new float[cnt];
    for(i=0 ; ((i < N+1) && (i < cnt)); i++) {
      magic_vector[i] = 1.0;
    }
    for(;i < cnt ; i++) {
      magic_vector[i] = 0.0;
    }
    for (i=0 ; i<cnt ; i++) {
      cAmp[i].re *= magic_vector[i];
      cAmp[i].im *= magic_vector[i];
    }

    //for (int i=0 ; i<cnt ; i++) {
    //  printf("i=%2d, cAmp.re[i]=%11.4f  cAmp.im[i]=%11.4f\n",i,cAmp[i].re,cAmp[i].im);
    //}

    //
    // Backtransform to time domain
    //
    rfftwnd_plan p2 = rfftwnd_create_plan(1, &nt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_one_complex_to_real(p2, cAmp, rAmp);
    fftwnd_destroy_plan(p2);

    float scale= float(1.0/nt);
    for(i=0 ; i < rnt ; i++) {
      rAmp[i] *= scale;
    }

    //
    // Fill log_filtered[]
    //
    for(i=0 ; i < n_time_samples_defined ; i++) {
      log_filtered[first_nonmissing + i] = rAmp[i];      // Fill with values where defined
    }
    delete [] magic_vector;
    fftw_free(rAmp);

    //for (int i=0 ; i<n_time_samples ; i++) {
    //  printf("i log_interpolated[i] log_filtered[i]  %d  %.3f  %.3f\n",i,log_interpolated[i],log_filtered[i]);
    //}
  }
}

//----------------------------------------------

bool compare(const std::pair<int, float>& i1, const std::pair<int, float>& i2)
{
  return (i1.second < i2.second);
}

//----------------------------------------------------------------------------
void
WellData::findILXLAtStartPosition(void)
{
  const SegyGeometry * full_geometry = modelSettings_->getSeismicDataAreaParameters();

  if (full_geometry != NULL) {
    int IL0, XL0;
    full_geometry->FindILXL(static_cast<float>(xpos0_), static_cast<float>(ypos0_), IL0, XL0);
    LogKit::LogFormatted(LogKit::Low,"   IL/XL at start of well is %d/%d\n",IL0,XL0);
  }
}

//----------------------------------------------------------------------------
void
WellData::lookForSyntheticVsLog(float & rank_correlation)
{
  float corr_threshold = modelSettings_->getMaxRankCorr();

  //
  // Estimate the correlation between Vp and Vs logs. To be able to identify
  // nonlinear relationships between the logs we use rank correlation.
  //
  typedef std::pair<int, float> Item;
  std::vector<Item> sorted_alpha;
  std::vector<Item> sorted_beta;

  //
  // Store Vp and Vs in sortable structs. Note that we can only use nonmissing values.
  //
  for (int i = 0 ; i < nd_ ; i++)
  {
    if (alpha_[i] != RMISSING && beta_[i] != RMISSING)
    {
      sorted_alpha.push_back(Item(i, alpha_[i]));
      sorted_beta.push_back(Item(i,beta_[i]));
    }
  }

  int n = static_cast<int>(sorted_alpha.size());
  if (n > 0)
  {
    //
    // Sort Vp and Vs logs.
    //
    std::sort(sorted_alpha.begin(), sorted_alpha.end(), compare);
    std::sort(sorted_beta.begin(), sorted_beta.end(), compare);

    //
    // Estimate correlation between sorted alpha and beta sorted by alpha
    //
    float mean = float(n)/2.0f; // We start the indexing at 0 rather than 1
    float cov_rank = 0.0;       // Covariance between ranks of alpha and beta (sorted on alpha)
    float var_rank = 0.0;       // Variance in ranks of alpha and beta (which are equal)

    for (int i = 0 ; i < n ; i++)
    {
      var_rank +=(i - mean)*(i - mean);
      for (int j = 0 ; j < n ; j++)
      {

        if (sorted_beta[j].first == sorted_alpha[i].first)
        {
          cov_rank += (j - mean)*(i - mean);
        }
      }
    }
    rank_correlation = cov_rank/var_rank; // Skip division by n-1 in both nominator and denominator

    if (rank_correlation > corr_threshold)
    {
      if(realVsLog_ == ModelSettings::NOTSET) {
        LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. Treating Vs log as synthetic.\n",rank_correlation);
        realVsLog_ = ModelSettings::NO;
      }
      else {
        if(realVsLog_ == ModelSettings::YES)
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f, but well log is defined as real.\n",rank_correlation);
        else
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
      }
    }
    else
    {
      switch(realVsLog_) {
        case ModelSettings::YES :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as real.)\n",rank_correlation);
          break;
        case ModelSettings::NO :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
          break;
        default :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is treated as real.)\n",rank_correlation);
          break;
      }
    }
  }
  else
  {
    LogKit::LogFormatted(LogKit::Low,"   Cannot calculate Vp-Vs rank correlation. One or both logs are empty.\n");
  }

  bool useFilter  = modelSettings_->getUseFilterForFaciesProb();
  bool useVpVsRho = modelSettings_->getNoVsFaciesProb() == false;

  if(useFilter && useVpVsRho && realVsLog_ == ModelSettings::NO && useForFaciesProbabilities_ == ModelSettings::NOTSET)
    useForFaciesProbabilities_ = ModelSettings::NO;
}

//----------------------------------------------------------------------------
void
WellData::calculateDeviation(float  & devAngle,
                             Simbox * timeSimbox)
{
  float maxDevAngle   = modelSettings_->getMaxDevAngle();
  float thr_deviation = float(tan(maxDevAngle*NRLib::Pi/180.0));  // Largest allowed deviation
  float max_deviation =  0.0f;
  float max_dz        = 10.0f;                      // Calculate slope each 10'th millisecond.

  //
  // Find first log entry in simbox
  //
  int iFirst = IMISSING;
  int i;
  for(i=0 ; i < nd_ ; i++)
  {
    if(timeSimbox->isInside(xpos_[i], ypos_[i]))
    {
      if (zpos_[i] > timeSimbox->getTop(xpos_[i], ypos_[i]))
      {
        iFirst = i;
        break;
      }
    }
  }

  if (iFirst != IMISSING) {
    //
    // Find last log entry in simbox
    //
    int iLast = iFirst;
    for(i = iFirst + 1 ; i < nd_ ; i++)
    {
      if(timeSimbox->isInside(xpos_[i], ypos_[i]))
        {
          if (zpos_[i] > timeSimbox->getBot(xpos_[i], ypos_[i]))
            break;
        }
      else
        break;
      iLast = i;
    }

    if (iLast > iFirst) {
      double x0 = xpos_[iFirst];
      double y0 = ypos_[iFirst];
      double z0 = zpos_[iFirst];
      for (int i = iFirst+1 ; i < iLast+1 ; i++) {
        double x1 = xpos_[i];
        double y1 = ypos_[i];
        double z1 = zpos_[i];
        float dz = static_cast<float>(z1 - z0);

        if (dz > max_dz || i == iLast) {
          float deviation = static_cast<float>(sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0))/dz);
          if (deviation > max_deviation) {
            x0 = x1;
            y0 = y1;
            z0 = z1;
            max_deviation = deviation;
          }
        }
      }
    }
    devAngle = static_cast<float>(atan(max_deviation)*180.0/NRLib::Pi);
    LogKit::LogFormatted(LogKit::Low,"   Maximum local deviation is %.1f degrees.",devAngle);

    if (max_deviation > thr_deviation)
    {
      if(useForWaveletEstimation_ == ModelSettings::NOTSET) {
        useForWaveletEstimation_ = ModelSettings::NO;
      }
      isDeviated_ = true;
      LogKit::LogFormatted(LogKit::Low," Well is treated as deviated.\n");
    }
    else
    {
      isDeviated_ = false;
      LogKit::LogFormatted(LogKit::Low,"\n");
    }
  }
  else {
    isDeviated_ = false;
    LogKit::LogFormatted(LogKit::Low,"Well is outside inversion interval. Cannot calculate deviation.\n");
  }
}

//----------------------------------------------------------------------------
void
WellData::countFacies(Simbox *simbox, int * faciesCount)
{
  for (int i=0 ; i < nFacies_ ; i++)
    faciesCount[i] = 0;

  if (nFacies_ > 0) {
    for(int i=0 ; i < nd_ ; i++) {
      if (facies_[i] != IMISSING) {
        //
        // Count facies only when logs
        //
        if (alpha_[i] != RMISSING && beta_[i] != RMISSING && rho_[i] != RMISSING) {
          //
          // Count facies only for part of well inside simbox
          //
          if(simbox->isInside(xpos_[i], ypos_[i])) {
            if(zpos_[i] > simbox->getTop(xpos_[i], ypos_[i]) && zpos_[i] < simbox->getBot(xpos_[i], ypos_[i]))  {
              //
              // Find facies index...
              //
              for (int j=0 ; j<nFacies_ ; j++) {
                if (faciesNr_[j] == facies_[i]) {
                  faciesCount[j]++;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void WellData::getMinMaxFnr(int &min, int&max) const
{
  int i;
  //int premin, premax;
  min = faciesNr_[0];
  max = faciesNr_[0];
  for(i=1;i<nFacies_;i++)
  {
    if(faciesNr_[i]<min)
      min = faciesNr_[i];
    if(faciesNr_[i]>max)
      max = faciesNr_[i];
  }
}

void WellData::moveWell(Simbox * timeSimbox, double deltaX, double deltaY, float kMove)
{

  int    i;
  double deltaZ;
  double topOld, topNew;

  topOld = timeSimbox->getTop(xpos_[0], ypos_[0]);

  for(i=0; i<nd_; i++)
  {
    xpos_[i] = xpos_[i]+deltaX;
    ypos_[i] = ypos_[i]+deltaY;
  }

  topNew = timeSimbox->getTop(xpos_[0], ypos_[0]);

  deltaZ = topNew - topOld + kMove;

  for(i=0; i<nd_; i++)
    zpos_[i] = zpos_[i]+deltaZ;

}

void WellData::findMeanVsVp(void)
{
  meanVsVp_ = 0.0f;
  nVsVp_    = 0;

  for (int i=0 ; i < nd_ ; i++) {
    if (alpha_background_resolution_[i] != RMISSING &&
        beta_background_resolution_[i]  != RMISSING) {
      meanVsVp_ += beta_background_resolution_[i]/alpha_background_resolution_[i];
      nVsVp_    += 1;
    }
  }
  meanVsVp_ /= nVsVp_;

  LogKit::LogFormatted(LogKit::Low,"   Average Vp/Vs ratio is %5.3f.\n",1.0f/meanVsVp_);
}
