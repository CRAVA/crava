#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/lib_misc.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/well/norsarwell.hpp"

#include "src/definitions.h"
#include "src/welldata.h"
#include "src/modelsettings.h"
#include "src/model.h"
#include "src/io.h"

//----------------------------------------------------------------------------
WellData::WellData(const std::string              & wellFileName, 
                   const std::vector<std::string> & logNames,
                   const std::vector<bool>        & inverseVelocity,
                   ModelSettings                  * modelSettings,
                   int                              indicatorFacies,
                   int                              indicatorWavelet,
                   int                              indicatorBGTrend,
                   int                              indicatorRealVs,
                   bool                             faciesLogGiven)
  : modelSettings_(modelSettings),
    wellname_(NULL),
    wellfilename_(NULL),
    xpos_(NULL),
    ypos_(NULL),
    zpos_(NULL),
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
    useForWaveletEstimation_(indicatorWavelet),  
    useForBackgroundTrend_(indicatorBGTrend),
    realVsLog_(indicatorRealVs)
{
  sprintf(errTxt_,"%c",'\0');
  if(wellFileName.find(".nwh",0) != std::string::npos)
    readNorsarWell(wellFileName, logNames, inverseVelocity, faciesLogGiven);
  else
    readRMSWell(wellFileName, logNames, inverseVelocity, faciesLogGiven);
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

  if (nFacies_ > 0)
  {
    delete [] faciesLogName_;
    for (int i = 0 ; i < nFacies_ ; i++)
      if(faciesNames_[i] != NULL)
        delete [] faciesNames_[i];
    delete [] faciesNames_;
    delete [] faciesNr_;
  }
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
  delete [] wellfilename_;
  delete [] wellname_;
}

//----------------------------------------------------------------------------
void
WellData::readRMSWell(const std::string              & wellFileName, 
                      const std::vector<std::string> & logNames, 
                      const std::vector<bool>        & inverseVelocity,
                      bool                             faciesLogGiven)
{
  error_ = 0;
  FILE* file = fopen(wellFileName.c_str(), "r");
  int i,j, facies;
  double xpos, ypos, zpos, dummy;
  float alpha, beta, rho;
  wellfilename_ = new char[MAX_STRING];
  strcpy(wellfilename_,wellFileName.c_str());
  sprintf(errTxt_,"%c",'\0');
  if(file == 0)
  {
    error_ = 1;
    sprintf(errTxt_,"Could not open file %s for reading.\n", wellFileName.c_str());
    //NBNB Incomplete solution, but should never happen.
  }
  int nlog; // number of logs in file

  wellname_ = new char[MAX_STRING];

  //First two lines contain info. we do not need.
  readToEOL(file);
  readToEOL(file);
  fscanf(file, "%s %lf %lf",wellname_,&xpos0_, &ypos0_); // read wellname and positions 
  readToEOL(file);                                       // Line may contain a dummy number 
  fscanf(file,"%d",&nlog);                               // read number of logs 

  char tmpStr[MAX_STRING];

  //Start searching for key words.

  int nVar = 5;       // z,alpha,beta,rho, and facies

  std::vector<std::string> parameterList(5);

  bool vpLog = false;
  bool vsLog = false;

  if(logNames[0] != "") // Assume that all lognames are filled present if first is.
  {
    parameterList = logNames; 
    if (!faciesLogGiven)
      nVar = 4; 
    vpLog = !inverseVelocity[0];
    vsLog = !inverseVelocity[1];
  }
 else
  {
    parameterList[0] = "TWT"; // Names are preliminary, to be changed
    parameterList[1] = "DT";
    parameterList[2] = "RHOB";
    parameterList[3] = "DTS";
    parameterList[4] = "FACIES";
  }
  int * pos = new int[nVar];
  for(i=0;i<nVar;i++)
    pos[i] = IMISSING;

  nFacies_ = 0;
  for(i=0;i<nlog;i++)
  { 
    fscanf(file,"%s",tmpStr);
    for(j=0;j<nVar;j++)
    {
      if(strcmp(uppercase(tmpStr),parameterList[j].c_str())==0)
      {
        pos[j] = i + 4;
        if(j==4)
        {
          faciesLogName_ = new char[MAX_STRING];
          strcpy(faciesLogName_,parameterList[4].c_str());
          // facies log - save names
          fscanf(file,"%s",tmpStr); // read code word DISC
          if(strcmp(tmpStr,"DISC")!=0)
          {
            LogKit::LogFormatted(LogKit::ERROR,"ERROR: Facies log must be discrete.\n");
            exit(1);
          }
          char EOL   = '\n';
          char blank = ' ';
          char rc = static_cast<char>(fgetc(file));
          while(rc == blank)
            rc = static_cast<char>(fgetc(file));
          while(rc != EOL)
          {
            // Read numbers or strings           
            while(rc != blank && rc!=EOL)
              rc = static_cast<char>(fgetc(file));
            // Read blanks
            while(rc == blank && rc!=EOL)
              rc = static_cast<char>(fgetc(file));
            nFacies_++;
          }
          nFacies_ = nFacies_/2;
        }
      }
    }
    if(strcmp(tmpStr,"DISC")!=0) readToEOL(file);
  }

  char missVar[MAX_STRING];
  strcpy(missVar,"");
  for(i=0 ; i<nVar ; i++)
  {
    if(pos[i]==IMISSING)
    {
      sprintf(missVar,"%s %s", missVar, parameterList[i].c_str());
      error_ = 1;
    }
  }
  if(error_ > 0)
    sprintf(errTxt_,"Cannot find log(s) %s in well file %s.\n",missVar,wellfilename_);


  if(pos[0]==IMISSING)
    timemissing_ = 1;
  else
    timemissing_ = 0;

  // Find nd_, the number of observations in well. 
  // Count the number of time observations which is not missing values.

  int nData = 0; 
  int legalData = 0;
  while(fscanf(file,"%s",tmpStr) != EOF)
  {
    nData++;
    fscanf(file," %lf", &dummy); // read y, which we do not need yet.
    fscanf(file,"%lf", &dummy);  // read z which we do not need.
    for(j=4;j<=nlog+3;j++)
    {
      fscanf(file,"%lf",&dummy);
      if(j==pos[0] && dummy != WELLMISSING)
      {
        // Found legal TIME variable
        legalData++;   
      }
    }
  }
  fclose(file);
  nd_ = legalData;
  //
  // Check that the number of logs found for each log entry agrees
  // with the number of logs specified in header.
  //
  // A nicer and faster implementation for this is requested...
  //
  int logEntry = 0;
  file = fopen(wellFileName.c_str(),"r");
  for(i=0;i<4+nlog;i++)
    readToEOL(file);
  while(fscanf(file,"%[^\n]\n",tmpStr) != EOF && error_ == 0)
  {
    logEntry++;
    int  elements = 0;
    bool lastIsBlank = true;
    for (i=0 ; i<=static_cast<int>(strlen(tmpStr)) ; i++)
    {
      if (tmpStr[i] != ' ' && tmpStr[i] != '\t'  && tmpStr[i] != '\r'  && tmpStr[i] != '\0') 
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
    sprintf(errTxt_,"%sERROR for well %s : The number of log elements (nlogs=%d) in line %d does not match header specifications (nlogs=%d).\n",
            errTxt_, wellname_, elements-3, logEntry, nlog);
   }
  }
  fclose(file);

  //
  // Read logs
  //
  int k;
  if (nFacies_ > 0)
  {
    faciesNr_    = new int[nFacies_];
    faciesNames_ = new char*[nFacies_];
    for(i=0 ; i<nFacies_ ; i++) {
      faciesNames_[i] = new char[MAX_STRING];
    }
  }

  file = fopen(wellFileName.c_str(),"r");
  for(i=0;i<4+nlog;i++)
  {
    fscanf(file,"%s",tmpStr);
    if(strcmp(uppercase(tmpStr),parameterList[4].c_str())==0)
    {
     fscanf(file,"%s",tmpStr); // read code word DISC
      // facies types given here
      for(k=0;k<nFacies_;k++)
      {
        fscanf(file,"%d",&(faciesNr_[k]));
        fscanf(file,"%s",faciesNames_[k]);
      }
    }
    readToEOL(file);
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
  alpha     = RMISSING;
  beta      = RMISSING;
  rho       = RMISSING;
  zpos      = RMISSING;
  facies    = IMISSING;
  k         = -1;
  faciesok_ = 1;
  int legal = 0;
  for(i=0;i<nData;i++)
  {
    // fscanf(file, "%lf %lf %lf %f %f %f",&xpos, &ypos, &zpos, &alpha, &beta, &rho);
    fscanf(file, "%lf %lf",&xpos, &ypos);
    fscanf(file,"%lf", &dummy); // read z which we do not need.
    for(j=4;j<=nlog+3;j++)
    {
      fscanf(file,"%lf",&dummy);
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
      else if(j==pos[4])
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
      if (alpha == OPENWORKS_MISSING)
        wrongMissingValues = true;
      if(beta == OPENWORKS_MISSING)
        wrongMissingValues = true;
      if(rho == OPENWORKS_MISSING)
        wrongMissingValues = true;
    }
  }
  fclose(file);

  if(wrongMissingValues)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nWARNING: There are incorrect missing values \'%f\' in well.\n",OPENWORKS_MISSING,wellfilename_);
    LogKit::LogFormatted(LogKit::LOW,"           CRAVA can handle this, but it is not a legal RMS well file format.\n");
  }
  if (k != nd_-1) {
    LogKit::LogFormatted(LogKit::LOW,"ERROR: Well %s is corrupt. nd=%d log entries were expected, but only k=%d were read.\n",
                     wellfilename_,nd_,k);
    exit(1);
  }
  delete [] pos;
}


void
WellData::readNorsarWell(const std::string              & wellFileName, 
                         const std::vector<std::string> & logNames, 
                         const std::vector<bool>        & inverseVelocity,
                         bool                             faciesLogGiven)
{
  error_ = 0;
  wellfilename_ = new char[MAX_STRING];
  strcpy(wellfilename_,wellFileName.c_str());
  wellname_ = new char[MAX_STRING];
  std::string name = NRLib::RemovePath(wellFileName);
  name = NRLib::ReplaceExtension(name, "");
  strcpy(wellname_, name.c_str());
  
  try 
  {
    NRLib::NorsarWell well(wellFileName);

    int nVar = 5;       // z,alpha,beta,rho, and facies

    std::vector<std::string> parameterList;

    bool vpLog = false;
    bool vsLog = false;

    if(logNames[0] != "") // Assume that all lognames are filled present if first is.
    {
      parameterList = logNames; 
      if (!faciesLogGiven)
        nVar = 4; 
      vpLog = !inverseVelocity[0];
      vsLog = !inverseVelocity[1];
    }
   else
    {
      parameterList[0] = "TWT"; // Names are preliminary, to be changed
      parameterList[1] = "DT";
      parameterList[2] = "RHOB";
      parameterList[3] = "DTS";
      parameterList[4] = "FACIES";
    }

    int nLogs = 2+nVar;
    std::vector<std::vector<double> *> logs(nLogs, NULL);
    logs[0] = well.GetContLog("UTMX");
    if(logs[0] == NULL) {
      error_ = 1;
      sprintf(errTxt_,"%sError: Could not find log 'UTMX' in well file '%s'.\n", errTxt_, 
        wellFileName.c_str());
    }
    logs[1] = well.GetContLog("UTMY");
    if(logs[1] == NULL) {
      error_ = 1;
      sprintf(errTxt_,"%sError: Could not find log 'UTMY' in well file '%s'.\n", 
        errTxt_, wellFileName.c_str());
    }
    for(int i=0;i<nVar;i++) {
      logs[2+i] = well.GetContLog(parameterList[i]);
      if(logs[2+i] == NULL) {
        if(i != 4 || logNames[0] != "") {
          error_ = 1;
          sprintf(errTxt_,"%sError: Could not find log '%s' in well file '%s'.\n", 
            errTxt_, parameterList[i].c_str(), wellFileName.c_str());
        }
        else if(i==4)
          nLogs = nLogs-1;
      }
    }

    if(error_ == 0) {
      faciesok_ = 1;
      std::vector<int> facCodes;
      nd_       = static_cast<int>(logs[0]->size());
      xpos_     = new double[nd_];
      ypos_     = new double[nd_];
      zpos_     = new double[nd_];
      alpha_    = new float[nd_];
      beta_     = new float[nd_];
      rho_      = new float[nd_];
      facies_   = new int[nd_];   // Always allocate a facies log (for code simplicity)
      for(int i=0;i<nd_;i++) {
        xpos_[i]  = (*logs[0])[i];
        ypos_[i]  = (*logs[1])[i];
        zpos_[i]  = (*logs[2])[i]*1000;
        if(!well.IsMissing((*logs[3])[i]))
          alpha_[i] = static_cast<float>((*logs[3])[i]);
        else 
          alpha_[i] = RMISSING;
        if(!well.IsMissing((*logs[5])[i]))
          beta_[i]  = static_cast<float>((*logs[5])[i]);
        else 
          beta_[i] = RMISSING;
        if(!well.IsMissing((*logs[4])[i]))
          rho_[i]   = static_cast<float>((*logs[4])[i]);
        else 
          rho_[i] = RMISSING;
        if(nLogs > 6 && logs[6] != NULL && 
          !well.IsMissing((*logs[4])[i])) {
          facies_[i]  = static_cast<int>((*logs[6])[i]);
          if(find(facCodes.begin(), facCodes.end(), facies_[i]) == facCodes.end())
            facCodes.push_back(facies_[i]);
        }
        else
          facies_[i] = IMISSING;
      }
      nFacies_ = facCodes.size();
    }
  }
  catch (NRLib::Exception & e) {
    sprintf(errTxt_,"Error: %s\n",e.what());
    error_ = 1;
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
int WellData::checkError(char *errText)
{
  if(error_ == 0)
    return(0);
  else
  {
    strcpy(errText,errTxt_);
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
    if (insideArea)
      LogKit::LogFormatted(LogKit::LOW,"   IGNORED (well is inside inversion area but does not hit the inversion volume)\n");
    else
      LogKit::LogFormatted(LogKit::LOW,"   IGNORED (well is not inside inversion area)\n");
  }
  return(error);
}

//----------------------------------------------------------------------------
void
WellData::removeDuplicateLogEntries(int & nMerges) 
{
  bool debug = false;
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
      while (zpos_[iend] - zpos_[istart] < minMergeDist && iend < nd_ - 1)
        iend++; 
      iend--;
    }

    bool printToScreen = debug && iend != istart;
    if (printToScreen)
      LogKit::LogFormatted(LogKit::LOW,"Merge log entries %d -> %d into new entry %d\n",istart,iend,ii);

    mergeCells("time",zpos_resampled, zpos_ ,ii,istart,iend,printToScreen);
    mergeCells("x   ",xpos_resampled, xpos_ ,ii,istart,iend,printToScreen);
    mergeCells("y   ",ypos_resampled, ypos_ ,ii,istart,iend,printToScreen);
    mergeCells("Vp  ",alpha_resampled,alpha_,ii,istart,iend,printToScreen);
    mergeCells("Vs  ",beta_resampled, beta_ ,ii,istart,iend,printToScreen);
    mergeCells("Rho ",rho_resampled,  rho_  ,ii,istart,iend,printToScreen);
    if(isFaciesLogDefined())
      mergeCellsDiscrete("Facies ",facies_resampled, facies_, ii, istart, iend, printToScreen);

    if (printToScreen)
      LogKit::LogFormatted(LogKit::LOW,"\n");
    
    istart = iend + 1;
    ii++;
  }    

  nMerges = 0;
  if (ii != nd_)
  {
    nMerges = nd_-ii; 
    LogKit::LogFormatted(LogKit::LOW,"   Duplicate log entries merged with neighbour : %d\n",nMerges);
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
}

//----------------------------------------------------------------------------
void
WellData::mergeCells(const char * name, double * pos_resampled, double * pos, 
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
        LogKit::LogFormatted(LogKit::LOW,"%s     Old:%d   pos = %.3f\n",name,i,pos[i]);
    }
  }
  if (nSample == 0)
    pos_resampled[ii] = RMISSING;
  else
    pos_resampled[ii] /= nSample;
  if (printToScreen)
    LogKit::LogFormatted(LogKit::LOW,"%s     New:%d   pos = %.3f\n",name,ii,pos_resampled[ii]);
}

//----------------------------------------------------------------------------
void
WellData::mergeCells(const char * name, float * log_resampled, float * log, int ii, 
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
        LogKit::LogFormatted(LogKit::LOW,"%s     Old:%d   log = %.3f\n",name,i,log[i]);
    }
  }
  if (nSample == 0)
    log_resampled[ii] = RMISSING;
  else
    log_resampled[ii] /= nSample;
  if (printToScreen)
    LogKit::LogFormatted(LogKit::LOW,"%s     New:%d   log = %.3f\n",name,ii,log_resampled[ii]);
}

//---------------------------------------------------------------------------

void
WellData::mergeCellsDiscrete(const char * name, int * log_resampled, int * log, int ii, 
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
        LogKit::LogFormatted(LogKit::LOW,"%s     Old:%d   log = %d\n",name,i,log[i]);
    }
  }
  log_resampled[ii] = log[istart];
 if (printToScreen)
    LogKit::LogFormatted(LogKit::LOW,"%s     New:%d   log = %d\n",name,ii,log_resampled[ii]);

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
        LogKit::LogFormatted(LogKit::LOW,"   Set undefined:   time = %.2f   Vp = %.2f\n",zpos_[i],alpha_[i]);
      alpha_[i] = RMISSING; 
      count_alpha++;
    }
    if (beta_[i] != RMISSING && (beta_[i] < beta_min  || beta_[i] > beta_max))
    {
      if (debug) 
        LogKit::LogFormatted(LogKit::LOW,"   Set undefined:   time = %.2f   Vs = %.2f\n",zpos_[i],beta_[i]);
      beta_[i] = RMISSING;
      count_beta++;
    }
    if (rho_[i] != RMISSING && (rho_[i] < rho_min || rho_[i] > rho_max))
    {
      if (debug) 
        LogKit::LogFormatted(LogKit::LOW,"   Set undefined:   time = %.2f   Rho = %.2f\n",zpos_[i],rho_[i]);
      rho_[i] = RMISSING;
      count_rho++;
    }
  }

  if (count_alpha > 0 || count_beta || count_rho > 0)
    LogKit::LogFormatted(LogKit::LOW,"   Log entries have been set undefined (Vp:%d Vs:%d Rho:%d)\n",
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
  resampleTime(time_resampled, nd_, dt); 

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

  /*
  for (unsigned int i = 0 ; i < nd_ ; i++)
  {
    printf("i = %d\n",i);
    printf("alpha_[i]                         %.3f\n", alpha_[i] );
    printf("beta_[i]                          %.3f\n", beta_[i] );
    printf("rho_[i]                           %.3f\n", rho_[i] );
    printf("alpha_seismic_resolution_[i]      %.3f\n", alpha_seismic_resolution_[i] );
    printf("beta_seismic_resolution_[i]       %.3f\n", beta_seismic_resolution_[i] );
    printf("rho_seismic_resolution_[i]        %.3f\n", rho_seismic_resolution_[i] );
    printf("alpha_background_resolution_[i]   %.3f\n", alpha_background_resolution_[i] );
    printf("beta_background_resolution_[i]    %.3f\n", beta_background_resolution_[i] );
    printf("rho_background_resolution_[i]     %.3f\n", rho_background_resolution_[i] );
  }
  */

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
void
WellData::resampleTime(double * time_resampled,
                       int      nd, 
                       double & dt) 
{
  //
  // Make new time scale with constant sampling density
  //
  double time_begin = zpos_[0];
  double time_end   = zpos_[nd - 1];
  
  if (time_begin != RMISSING && time_end != RMISSING)
  {
    dt = (time_end - time_begin)/(nd - 1);            // average sampling density
    for (int i = 0 ; i < nd ; i++) 
      time_resampled[i] = time_begin + i*dt;          
  }
  else 
  {
    LogKit::LogFormatted(LogKit::LOW,"ERROR: First or last time sample is undefined. Cannot estimate average sampling density.\n");
    LogKit::LogFormatted(LogKit::LOW,"       time[first] = %12.2f\n",time_begin);
    LogKit::LogFormatted(LogKit::LOW,"       time[last]  = %12.2f\n",time_end);
    exit(1);
  }

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
  
  for(int i=0 ; i < n_time_samples ; i++) {            // Initialise with RMISSING
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
    
    for (int i=0 ; i<n_time_samples_defined ; i++) {          // Array to filter is made symmetric
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

    float * magic_vector = new float[cnt];
    for(int i=0 ; i < N+1 ; i++) {
      magic_vector[i] = 1.0;
    }
    for(int i=N+1 ; i < cnt ; i++) {
      magic_vector[i] = 0.0; 
    }
    for (int i=0 ; i<cnt ; i++) {
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
    for(int i=0 ; i < rnt ; i++) {
      rAmp[i] *= scale;  
    }
    
    //
    // Fill log_filtered[]
    //  
    for(int i=0 ; i < n_time_samples_defined ; i++) { 
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
typedef struct 
{
  float        value;
  unsigned int index;
} Item, *pItem;


//----------------------------------------------
int compare(const void *l, const void *r)
{
  double res = ((pItem)l)->value - ((pItem)r)->value;
  if (res < 0) return -1;
  if (res > 0) return  1;
  return 0;
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
  Item * sorted_alpha = new Item[nd_];
  Item * sorted_beta  = new Item[nd_];

  //
  // Store Vp and Vs in sortable structs. Note that we can only use nonmissing values.
  //
  int n = 0;
  for (int i = 0 ; i < nd_ ; i++) 
  {
    if (alpha_[i] != RMISSING && beta_[i] != RMISSING) 
    {  
      sorted_alpha[n].value = alpha_[i];
      sorted_alpha[n].index = i;
      sorted_beta[n].value  = beta_[i];
      sorted_beta[n].index  = i;
      n++;
    }
  }

  if (n > 0)
  {
    //
    // Sort Vp and Vs logs. 
    //
    qsort(sorted_alpha, n, sizeof(Item), compare);
    qsort(sorted_beta , n, sizeof(Item), compare);
    
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
        if (sorted_beta[j].index == sorted_alpha[i].index) 
        {
          cov_rank += (j - mean)*(i - mean);
        }  
      }
    }
    rank_correlation = cov_rank/var_rank; // Skip division by n-1 in both nominator and denominator
    
    if (rank_correlation > corr_threshold) 
    {
      if(realVsLog_ == ModelSettings::NOTSET) {
        LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f. Treating Vs log as synthetic.\n",rank_correlation);
        realVsLog_ = ModelSettings::NO;
      }
      else {
        if(realVsLog_ == ModelSettings::YES)
          LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f, but well log is defined as real.\n",rank_correlation);
        else
          LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
      }
    }
    else
    {
      switch(realVsLog_) {
        case ModelSettings::YES :
          LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as real.)\n",rank_correlation);
          break;
        case ModelSettings::NO :
          LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
          break;
        default :
          LogKit::LogFormatted(LogKit::LOW,"   Vp-Vs rank correlation is %5.3f. (Well log is treated as real.)\n",rank_correlation);
          break;
      }
    }
  }
  else 
  {
    LogKit::LogFormatted(LogKit::LOW,"   Cannot calculate Vp-Vs rank correlation. One or both logs are empty.\n");
  }
  if(realVsLog_ == ModelSettings::NO && useForFaciesProbabilities_ == ModelSettings::NOTSET)
    useForFaciesProbabilities_ = ModelSettings::NO;

  delete [] sorted_alpha;
  delete [] sorted_beta;
}

//----------------------------------------------------------------------------
void 
WellData::calculateDeviation(float  & devAngle,
                             Simbox * timeSimbox)
{
  float maxDevAngle   = modelSettings_->getMaxDevAngle();
  float thr_deviation = float(tan(maxDevAngle*PI/180.0));  // Largest allowed deviation
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
  double x0 = xpos_[iFirst];
  double y0 = ypos_[iFirst];
  double z0 = zpos_[iFirst];
  for (int i = iFirst+1 ; i < iLast+1 ; i++) 
  {
    double x1 = xpos_[i];
    double y1 = ypos_[i];
    double z1 = zpos_[i];
    float dz = static_cast<float>(z1 - z0);
    if (dz > max_dz)
    {
      float deviation = static_cast<float>(sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0))/(z1 - z0));
      if (deviation > max_deviation) 
      {
        max_deviation = deviation;
      }
    }
  }
  devAngle = static_cast<float>(atan(max_deviation)*180.0/PI);
  LogKit::LogFormatted(LogKit::LOW,"   Maximum local deviation is %.1f degrees.",devAngle);

  if (max_deviation > thr_deviation) 
  {
    if(useForWaveletEstimation_ == ModelSettings::NOTSET)   
      useForWaveletEstimation_ = ModelSettings::NO;    
    isDeviated_ = true;
    LogKit::LogFormatted(LogKit::LOW," Well is treated as deviated.\n");
  }
  else
  {
    isDeviated_ = false;
    LogKit::LogFormatted(LogKit::LOW,"\n");
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

void WellData::getMeanVsVp(float & muA, float & muB)
{
  int i;
  int n = 0;

  muA = 0;
  muB = 0;

  for( i=0; i<nd_; i++ )
  {
    if (alpha_background_resolution_[i] != RMISSING && beta_background_resolution_[i] != RMISSING )
    {
      muA += alpha_background_resolution_[i];
      muB += beta_background_resolution_[i];
      n += 1;
    }
  }
  muA /= n;
  muB /= n;
}