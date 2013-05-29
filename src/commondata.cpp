/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/commondata.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/modelgeneral.h"

#include "nrlib/well/well.hpp"

#include "src/seismicstorage.h"
#include "nrlib/well/norsarwell.hpp"

CommonData::CommonData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles)
{

  //if(readSeismicData(modelSettings,
  //                   inputFiles) == true)
  //  read_seismic_ = true; //True or false if there is no seismic data?

  if(readWellData(modelSettings,
                  inputFiles) == true)
    read_wells_ = true;


}

CommonData::~CommonData() {
}

bool
CommonData::createOuterTemporarySimbox() {

  return true;
}

bool CommonData::readSeismicData(ModelSettings  * modelSettings,
                                 InputFiles     * inputFiles)
{
  //Skip if there is no AVO-seismic.
  int timelaps_seismic_files = 0;
  for(size_t i = 0; i < inputFiles->getTimeLapseSeismicFiles().size(); i++) {
    timelaps_seismic_files += inputFiles->getNumberOfSeismicFiles(i);
  }
  if(timelaps_seismic_files == 0 && inputFiles->getSeismicFiles().size() == 0)
    return(true);
  //Skip if mode is estimation and wavelet/noise is not set for estimation.
  if(modelSettings->getEstimationMode() == true && modelSettings->getEstimateWaveletNoise() == false)
    return(true);

//  const std::vector<std::string> & seismic_files = inputFiles->getSeismicFiles();
  const std::vector<std::vector<std::string> > seismic_timelapse_files = inputFiles->getTimeLapseSeismicFiles();

  int nTimeLapses = modelSettings->getNumberOfTimeLapses(); //Returnerer timeLapseAngle_.size()
//  int nVintages = modelSettings->getNumberOfVintages();
  int numberOfAngles = modelSettings->getNumberOfTimeLapses();
  //bool failed = false;

  for(int thisTimeLapse = 0; thisTimeLapse < nTimeLapses; thisTimeLapse++) {
    if(inputFiles->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) {

      LogKit::WriteHeader("Reading seismic data");

      std::vector<float> angles = modelSettings->getAngle(thisTimeLapse);
      std::vector<float> offset = modelSettings->getLocalSegyOffset(thisTimeLapse);

      for (int i = 0 ; i < numberOfAngles ; i++) {

        std::string filename = inputFiles->getSeismicFile(thisTimeLapse,i);
        int fileType = IO::findGridType(filename);

        if(fileType == IO::SEGY) { //From ModelGeneral::readSegyFile

          SegY * segy = NULL;
          TraceHeaderFormat * format = modelSettings->getTraceHeaderFormat(thisTimeLapse,i);

          if(format == NULL) { //Unknown format
            std::vector<TraceHeaderFormat*> traceHeaderFormats(0);
            if (modelSettings->getTraceHeaderFormat() != NULL)
            {
              traceHeaderFormats.push_back(modelSettings->getTraceHeaderFormat());
            }
            segy = new SegY(filename,
                            offset[i],
                            traceHeaderFormats,
                            true); // Add standard formats to format search
          }
          else //Known format, read directly.
            segy = new SegY(filename, offset[i], *format);

          float guard_zone = modelSettings->getGuardZone();

          //ModelGeneral::checkThatDataCoverGrid ?
          float dz = segy->GetDz();
          float z0 = offset[i] + 0.5f*dz;
          float zn = z0 + (segy->GetNz() - 1)*dz;


          float top_grid = static_cast<float>(estimation_simbox_->getTopZMin());
          float bot_grid = static_cast<float>(estimation_simbox_->getBotZMax());

          float top_guard = top_grid - guard_zone;
          float bot_guard = bot_grid + guard_zone;

          if (top_guard < z0) {
            float z0_new = z0 - ceil((z0 - top_guard)/dz)*dz;
            LogKit::LogFormatted(LogKit::Warning, "\nThere is not enough seismic data above the interval of interest. The seismic data\n"
                                                  " must start at "+NRLib::ToString(z0_new)+"ms (in CRAVA grid) to allow for a "
                                                  + NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n"
                                                  + "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n"
                                                  + "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n"
                                                  + "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n"
                                                  + "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n"
                                                  + "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n"
                                                  + "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n");
          }
          if (bot_guard > zn) {
            float zn_new = zn + ceil((bot_guard - zn)/dz)*dz;
            LogKit::LogFormatted(LogKit::Warning, "\nThere is not enough seismic data below the interval of interest. The seismic data\n"
                                                  "must end at "+NRLib::ToString(zn_new)+"ms (in CRAVA grid) to allow for a "
                                                  + NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n"
                                                  + "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n"
                                                  + "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n"
                                                  + "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n"
                                                  + "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n"
                                                  + "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n"
                                                  + "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n");
          }
          //End checkThatDataCoverGrid

          float padding         = 2*guard_zone;
          bool  relativePadding = false;
          bool  onlyVolume      = true;

          segy->ReadAllTraces(estimation_simbox_, //timeCutSimBox
                              padding,
                              onlyVolume,
                              relativePadding);
          segy->CreateRegularGrid();

          SeismicStorage seismicdata(filename, 0, angles[i], segy);
          seismic_data_.push_back(seismicdata);

        } //SEGY
        else if(fileType == IO::STORM || fileType == IO::SGRI) { //From ModelGeneral::readStormFile
          StormContGrid * stormgrid = NULL;
          bool failed = false;

          try
          {
            stormgrid = new StormContGrid(0,0,0);
            stormgrid->ReadFromFile(filename);
          }
          catch (NRLib::Exception & e)
          {
            LogKit::LogFormatted(LogKit::Warning, "Error when reading storm-file " + filename +": " + NRLib::ToString(e.what()) + "\n");
            failed = true;
          }

          if(failed == false) {
            SeismicStorage seismicdata;

            if(fileType == IO::STORM)
              seismicdata = SeismicStorage(filename, 1, angles[i], stormgrid);
            else
              seismicdata = SeismicStorage(filename, 2, angles[i], stormgrid);

            seismic_data_.push_back(seismicdata);
          }

        } //STORM / SGRI
        else
        {
          LogKit::LogFormatted(LogKit::Warning, "Error when reading file " + filename +". File type not recognized.\n");
        }
      }
    }
  }

  return true;
}

bool CommonData::readWellData(ModelSettings  * modelSettings,
                              InputFiles     * inputFiles)
{

  int nWells = modelSettings->getNumberOfWells();
  if(nWells > 0)
    LogKit::WriteHeader("Reading wells");
  int well_error = 0;
  std::vector<std::string> logNames = modelSettings->getLogNames();
  std::vector<bool> inverseVelocity = modelSettings->getInverseVelocity();
  bool faciesLogGiven = modelSettings->getFaciesLogGiven();

  for(int i=0 ; i<nWells; i++) {
    int error = 0;
    NRLib::Well new_well;
    std::string wellFileName = inputFiles->getWellFile(i);

    std::string well_name;
    std::vector<std::string>  faciesNames;

    std::vector<double> xpos;
    std::vector<double> ypos;
    std::vector<double> zpos;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> rho;
    std::vector<double> md;
    std::vector<int> facies; // Always allocate a facies log (for code simplicity). Keep this ??

    if(wellFileName.find(".nwh",0) != std::string::npos) { //Norsar well, from WellData::readNorsarWell

      std::string name = NRLib::RemovePath(wellFileName);
      name = NRLib::ReplaceExtension(name, "");
      well_name = name;
      int faciesok; // all faciesnumbers read are present in header
      double xpos0;
      double ypos0;

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
          parameterList[0] = "TWT";
          parameterList[1] = "DT";
          parameterList[2] = "RHOB";
          parameterList[3] = "DTS";
          parameterList[4] = "FACIES";
        }

        int nLogs  = 2+nVar;
        int nExtra = 1; //MD log, needed for writing.
        std::vector<double> * filler = NULL; //to eliminate warning.
        std::vector<std::vector<double> *> logs(nLogs+nExtra, filler);
        if(well.HasContLog("UTMX") == false) {
          error = 1;
          LogKit::LogFormatted(LogKit::Error,"Could not find log 'UTMX' in well file "+wellFileName+".\n");
          logs[0] = NULL;
        }
        else
          logs[0] = &(well.GetContLog("UTMX"));
        if(well.HasContLog("UTMY") == false) {
          error = 1;
          LogKit::LogFormatted(LogKit::Error,"Could not find log 'UTMY' in well file "+wellFileName+".\n");
          logs[1] = NULL;
        }
        logs[1] = &(well.GetContLog("UTMY"));
        for(int i=0;i<nVar;i++) {
          if(well.HasContLog(parameterList[i]) == false) {
            logs[2+i] = NULL;
            if(i != 4 || logNames[0] != "") {
              error = 1;
              LogKit::LogFormatted(LogKit::Error,"Could not find log "+parameterList[i]+" in well file "+wellFileName+".\n");
            }
            else if(i==4)
              nLogs = nLogs-1;
          }
          else
            logs[2+i] = &(well.GetContLog(parameterList[i]));
        }

        //Added MD log.
        int mdLog = nLogs;
        if(well.HasContLog("MD") == false) {
          error = 1;
          LogKit::LogFormatted(LogKit::Error,"Could not find log 'MD' in well file "+wellFileName+".\n");
          logs[mdLog] = NULL;
        }
        logs[mdLog] = &(well.GetContLog("MD"));
        nLogs++;

        //if(logs[2] == NULL)
        //  timemissing = 1;
        //else
        //  timemissing =0;
        int nFacies = 0;
        if(error == 0) {
          faciesok = 1;
          std::vector<int> facCodes;
          int nd = 0;
          for(size_t i=0;i<logs[2]->size();i++)
            if(well.IsMissing((*logs[2])[i]) == false)
              nd++;

          for(size_t i=0;i<logs[0]->size();i++) {
            if(well.IsMissing((*logs[2])[i]) == false) {
              xpos.push_back((*logs[0])[i]*1000);
              ypos.push_back((*logs[1])[i]*1000);
              zpos.push_back((*logs[2])[i]*1000);
              if(!well.IsMissing((*logs[3])[i]))
                alpha.push_back(static_cast<double>((*logs[3])[i]));
              else
                alpha.push_back(RMISSING);
              if(!well.IsMissing((*logs[5])[i]))
                beta.push_back(static_cast<float>((*logs[5])[i]));
              else
                beta.push_back(RMISSING);
              if(!well.IsMissing((*logs[4])[i]))
                rho.push_back(static_cast<float>((*logs[4])[i]));
              else
                rho.push_back(RMISSING);
              if(mdLog != 6 && nLogs > 6 && logs[6] != NULL && !well.IsMissing((*logs[6])[i])) {
                facies.push_back(static_cast<int>((*logs[6])[i]));

                if(find(facCodes.begin(), facCodes.end(), facies[facies.size()]) == facCodes.end())
                  facCodes.push_back(facies[facies.size()]);
              }
              else
                facies.push_back(IMISSING);
              if(!well.IsMissing((*logs[mdLog])[i]))
                md.push_back(static_cast<float>((*logs[mdLog])[i]));
              else
                md.push_back(RMISSING);
            }
          }
          nFacies = static_cast<int>(facCodes.size());

        }
        xpos0 = well.GetXPosOrigin()*1000;
        ypos0 = well.GetYPosOrigin()*1000;

      }
      catch (NRLib::Exception & e) {
        LogKit::LogFormatted(LogKit::Error,"Error: " + NRLib::ToString(e.what()));
        error = 1;
      }

    }
    else if(wellFileName.find(".rms",0) != std::string::npos) { //RMS well, from WellData::ReadRMSWell

      double dummy = RMISSING;
      int j,k;
      //int timemissing = 0;
      std::ifstream file;
      NRLib::OpenRead(file, wellFileName);
      std::string token, dummyStr;
      std::vector<std::string> tokenLine;
      //From WellData::ReadRMSWell

      int nlog; // number of logs in file
      int line = 0;
      NRLib::DiscardRestOfLine(file,line,false); //First two lines contain info we do not need.
      NRLib::DiscardRestOfLine(file,line,false);
      NRLib::ReadNextToken(file, token, line);
      well_name = token;
      double xpos0 = NRLib::ReadNext<double>(file, line); //Needed?
      double ypos0 = NRLib::ReadNext<double>(file, line); //Needed
      NRLib::DiscardRestOfLine(file,line,false);
      nlog   = NRLib::ReadNext<int>(file, line);

      int faciesok; // all faciesnumbers read are present in header

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
        parameterList[0] = "TWT";
        parameterList[1] = "DT";
        parameterList[2] = "RHOB"; //DTS before RHOB?
        parameterList[3] = "DTS";
        parameterList[4] = "FACIES";
      }

      int * pos = new int[nVar];
      for(k=0;k<nVar;k++)
        pos[k] = IMISSING;

      int nFacies = 0;
      std::string faciesLogName;
      for(k=0;k<nlog;k++)
      {
        NRLib::ReadNextToken(file,token,line);
        for(j=0;j<nVar;j++)
        {
          if( NRLib::Uppercase(token)==parameterList[j])
          {
            pos[j] = k + 4;
            if(j==4)
            {
              faciesLogName = parameterList[4];
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
              nFacies = static_cast<int>(tokenLine.size())/2;
            }
          }
        }
        if (token != "DISC")
          NRLib::DiscardRestOfLine(file,line,false);
      }

      std::string missVar = "";
      for(k=0 ; k<nVar ; k++)
      {
        if(pos[k]==IMISSING)
        {
          missVar += parameterList[k];
          error = 1;
        }
      }
      if(error > 0)
        LogKit::LogFormatted(LogKit::Error,"Cannot find log(s) "+missVar+" in well file "+wellFileName+".\n");

      //if(pos[0]==IMISSING)
      //  timemissing = 1;
      //else
      //  timemissing = 0;

      // Find nd_, the number of observations in well.
      // Count the number of time observations which is not missing values.

      int nData = 0;
      int legalData = 0;
      int nd;
      while (NRLib::CheckEndOfFile(file)==false)
      {
        nData++;
        try {
          dummy = NRLib::ReadNext<double>(file,line); // Read x which we do not need yet.
          dummy = NRLib::ReadNext<double>(file,line); // Read y which we do not need yet.
          dummy = NRLib::ReadNext<double>(file,line); // Read z which we do not need.

          for(j=4;j<=nlog+3;j++) {
            dummy = NRLib::ReadNext<double>(file,line); // Read z which we do not need.
            if(j==pos[0] && dummy != WELLMISSING) {
              legalData++;   // Found legal TIME variable
            }
          }
        }
        catch (NRLib::IOError e) {
          std::string text;
          text += std::string("\nERROR: Reading of well \'") + wellFileName + "\' failed for log record ";
          text += NRLib::ToString(nData) + " (not counting header lines).\n";
          text += std::string("\nERROR message is \'") + e.what() + "\'";
          LogKit::LogMessage(LogKit::Error,text);
          exit(1);
        }
      }
      file.close();
      file.clear();
      nd = legalData;

      //
      // Check that the number of logs found for each log entry agrees
      // with the number of logs specified in header.
      //
      // A nicer and faster implementation for this is requested...
      //

      int logEntry = 0;
      NRLib::OpenRead(file, wellFileName);
      line = 0;
      for(k=0;k<4+nlog;k++)
        NRLib::DiscardRestOfLine(file,line,false);
      while (NRLib::CheckEndOfFile(file)==false && error==0)
      {
        logEntry++;
        int  elements = 0;
        bool lastIsBlank = true;
        std::getline(file,token);

        int l = static_cast<int>(token.length());
        for (k=0; k<l; k++)
        {
          if (token[k] != ' ' && token[k] != '\t' && token[k] != '\r' && token[k] != '\0')
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
          error = 1;
          LogKit::LogFormatted(LogKit::Error,"ERROR for well "+well_name +": The number of log elements (nlogs="+NRLib::ToString(elements-3)+")\n in line "+NRLib::ToString(logEntry)+" does not match header specifications (nlogs="+NRLib::ToString(nlog)+").\n");
        }
      }
      file.close();
      file.clear();

      //
      // Read logs
      //
      std::vector<int> faciesNr;
      //if (nFacies > 0)
      //  faciesNr.resize(nFacies);

      NRLib::OpenRead(file, wellFileName);
      line = 0;
      for(k=0;k<4+nlog;k++)
      {
        NRLib::ReadNextToken(file,token,line);
        if (NRLib::Uppercase(token) == parameterList[4])
        {
          NRLib::ReadNextToken(file,token,line); // read code word DISC
          // facies types given here
          for(int kk=0;kk<nFacies;kk++)
          {
            NRLib::ReadNextToken(file,token,line);
            //faciesNr[kk] = NRLib::ParseType<int>(token);
            faciesNr.push_back(NRLib::ParseType<int>(token));
            NRLib::ReadNextToken(file,token,line);
            faciesNames.push_back(token);
          }
        }
        NRLib::DiscardRestOfLine(file,line,false);
      }
      double OPENWORKS_MISSING = -999.25;
      bool wrongMissingValues = false;

      double xpos_data = RMISSING;
      double ypos_data = RMISSING;
      double zpos_data = RMISSING;
      double alpha_data = RMISSING;
      double beta_data = RMISSING;
      double rho_data = RMISSING;
      int facies_data = IMISSING;
      k         = -1;
      faciesok = 1;
      int legal = 0;
      for(k=0;k<nData;k++)
      {
        xpos_data  = NRLib::ReadNext<double>(file,line);
        ypos_data  = NRLib::ReadNext<double>(file,line);
        dummy = NRLib::ReadNext<double>(file,line);
        for(j=4;j<=nlog+3;j++)
        {
          dummy = NRLib::ReadNext<double>(file,line);
          if(j==pos[0])
          {
            //Found TIME variable
            if(dummy != WELLMISSING && dummy != OPENWORKS_MISSING)
            {
              zpos_data = dummy;
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
                alpha_data = static_cast<double>(dummy);
              else
                alpha_data = static_cast<double>(304800.0/dummy);
            else
              alpha_data = RMISSING;
          }
          else if(j==pos[3])
          {
            // Found BETA variable
            if(dummy != WELLMISSING)
              if (vsLog)
                beta_data = static_cast<double>(dummy);
              else
                beta_data = static_cast<double>(304800.0/dummy);
            else
              beta_data = RMISSING;
          }
          else if(j==pos[2])
          {
            //Found RHO variable
            if(dummy != WELLMISSING)
              rho_data = static_cast<double>(dummy);
            else
              rho_data = RMISSING;
          }
          else if(nVar > 4 && j==pos[4])
          {
            //Found facies variable
            if(dummy != WELLMISSING)
              facies_data = static_cast<int>(dummy);
            else
              facies_data = IMISSING;
            if(facies_data!=IMISSING)
            {
              int faciesok_check = 0;
              for(int kk=0 ; kk<nFacies ; kk++)
              {
                if(facies_data == faciesNr[kk]) {
                  faciesok_check = 1;
                  break;
                }
              }
              if(faciesok_check == 0)
                faciesok = 0;
            }
          }
        }
        if(legal == 1)
        {
          //Cut against full_inversion_volume
          //double z_top = full_inversion_volume_->GetTopSurface().GetZ(xpos_data, ypos_data);
          //double z_bot = full_inversion_volume_->GetBotSurface().GetZ(xpos_data, ypos_data);

          //if(zpos_data < z_top && zpos_data > z_bot) {
            xpos.push_back(xpos_data);
            ypos.push_back(ypos_data);
            zpos.push_back(zpos_data);
            alpha.push_back(alpha_data);
            beta.push_back(beta_data);
            rho.push_back(rho_data);
            facies.push_back(facies_data);
            if (alpha_data == OPENWORKS_MISSING)
              wrongMissingValues = true;
            if(beta_data == OPENWORKS_MISSING)
              wrongMissingValues = true;
            if(rho_data == OPENWORKS_MISSING)
              wrongMissingValues = true;
          //}
        }
      }
      file.close();
    }

    //Add well logs
    new_well.SetWellName(well_name);
    new_well.AddContLog("xpos", xpos);
    new_well.AddContLog("ypos", ypos);
    new_well.AddContLog("zpos", zpos);
    new_well.AddContLog("alpha", alpha);
    new_well.AddContLog("beta", beta);
    new_well.AddContLog("rho", rho);
    if(md.size() > 0)
      new_well.AddContLog("md", md);

    new_well.AddDiscLog("facies", facies);
    new_well.SetMissing(RMISSING);

    wells_.push_back(new_well);

    well_error += error;
  }

  if(well_error > 0)
    return false;
  else
    return true;
}

bool CommonData::blockWellsForEstimation() {
  return true;

}

bool CommonData::setupReflectionMatrixAndTempWavelet() {
  return true;

}

bool CommonData::optimizeWellLocations() {
  return true;
}

bool CommonData::estimateWaveletShape() {
  return true;

}

bool CommonData::estimatePriorCorrelation() {
  return true;

}

bool CommonData::setupEstimationRockPhysics() {
  return true;

}