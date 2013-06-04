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
#include "src/timeline.h"
#include "src/fftgrid.h"
#include "src/modelgeneral.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "src/seismicstorage.h"
#include "nrlib/well/norsarwell.hpp"

CommonData::CommonData(ModelSettings  * model_settings,
                       InputFiles     * input_files)
:
  outer_temp_simbox_(false),
  read_seismic_(false),
  read_wells_(false),
  block_wells_(false),
  setup_reflection_matrix_(false),
  optimize_well_location_(false),
  wavelet_estimation_shape_(false),
  prior_corr_estimation_(false),
  setup_estimation_rock_physics_(false)
{
  bool failed = false;
  std::string err_text = "";
  //if(readSeismicData(modelSettings,
  //                   inputFiles) == true)
  //  read_seismic_ = true; //True or false if there is no seismic data?

  // 1. set up outer simbox
  outer_temp_simbox_ = createOuterTemporarySimbox(model_settings, input_files, estimation_simbox_, full_inversion_volume_, err_text);
  failed = !outer_temp_simbox_;

  // 3. read well data
  if(!failed){
    read_wells_ = readWellData(model_settings, input_files)
    failed = !read_wells_;
  }

  // 4. block wells for estimation
  // if well position is to be optimized
  if(!failed){
    if (model_settings->getOptimizeWellLocation() || model_settings->getEstimateWaveletNoise() || model_settings->getEstimateCorrelations()){

    }
  }

}

CommonData::~CommonData() {
  //delete estimation_simbox_;
  //delete full_inversion_volume_;
}

bool CommonData::createOuterTemporarySimbox(ModelSettings   * model_settings,
                                            InputFiles      * input_files,
                                            Simbox          & estimation_simbox,
                                            NRLib::Volume   & full_inversion_volume,
                                            std::string     & err_text){

  // parameters
  bool failed = false;
  std::string grid_file("");
  int  area_specification  = model_settings->getAreaSpecification();
  

  // FIND SEGY GEOMETRY FOR INVERSION AREA -------------------------------------------------------------------------

  bool estimation_mode_need_ILXL = model_settings->getEstimationMode() &&
                                (area_specification == ModelSettings::AREA_FROM_GRID_DATA ||
                                 area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM ||
                                 area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE ||
                                (model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0 ||
                                (model_settings->getOutputGridFormat() & IO::SEGY) > 0);

  if(model_settings->getForwardModeling())
    grid_file = input_files->getBackFile(0);    // Get geometry from earth model (Vp)
  else {
    if (model_settings->getEstimationMode() == false || estimation_mode_need_ILXL)
      grid_file = input_files->getSeismicFile(0,0); // Get area from first seismic data volume
  }
  SegyGeometry * ILXL_geometry = NULL; //Geometry with correct XL and IL settings.

  std::string area_type;

  // tilfelle i uten snap to seismic
  if (area_specification == ModelSettings::AREA_FROM_UTM)
  {
    // The geometry is already present in modelSettings (geometry_ or geometry_full_ ? )
    area_type = "Model file";
    LogKit::LogFormatted(LogKit::High,"\nArea information has been taken from model file\n");
  }
  // tilfelle ii uten snap to seismic
  else if(area_specification == ModelSettings::AREA_FROM_SURFACE){
    LogKit::LogFormatted(LogKit::High,"\nFinding area information from surface \'"+input_files->getAreaSurfaceFile()+"\'\n");
    area_type = "Surface";
    RotatedSurface surf(input_files->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    model_settings->setAreaParameters(&geometry);
  }

  else if(area_specification == ModelSettings::AREA_FROM_GRID_DATA         || // tilfelle iii 
          area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM || // tilfelle i med snap to seismic
          area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE){ // tilfelle ii med snap to seismic
    LogKit::LogFormatted(LogKit::High,"\nFinding inversion area from grid data in file \'"+grid_file+"\'\n");
    area_type = "Grid data";
    std::string tmp_err_text;
    SegyGeometry * geometry;
    getGeometryFromGridOnFile(grid_file,
                              model_settings->getTraceHeaderFormat(0,0), //Trace header format is the same for all time lapses
                              geometry,
                              tmp_err_text);
    if(geometry!=NULL){
      //tilfelle iv og v
      if (model_settings->getAreaILXL().size() > 0 || model_settings->getSnapGridToSeismicData()) {
        SegyGeometry * full_geometry = geometry;

        std::vector<int> area_ILXL;
        bool got_area = true;

        //
        // Geometry is given as XY, but we snap it to IL and XL.
        //
        if (model_settings->getSnapGridToSeismicData()) {
          SegyGeometry * template_geometry = NULL;
          if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM) {
            template_geometry = model_settings->getAreaParameters();
          }
          else if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE) {
            RotatedSurface surf(input_files->getAreaSurfaceFile());
            template_geometry = new SegyGeometry(surf);
          }
          else {
            err_text += "CRAVA has been asked to identify a proper ILXL inversion area based\n";
            err_text += "on XY input information, but no UTM coordinates or surface have\n";
            err_text += "been specified in model file.\n";
            got_area = false;
          }
          if (got_area) {
            area_ILXL = full_geometry->findAreaILXL(template_geometry);
          }
        }
        else {
          area_ILXL = model_settings->getAreaILXL();
        }

        if (got_area) {
          try {
            bool interpolated, aligned;
            geometry = full_geometry->GetILXLSubGeometry(area_ILXL, interpolated, aligned);

            std::string text;
            if(interpolated == true) {
              if(aligned == true) {
                text  = "Check IL/XL specification: Specified IL- or XL-step is not an integer multiple\n";
                text += "   of those found in the seismic data. Furthermore, the distance between first\n";
                text += "   and last XL and/or IL does not match the step size.\n";
                TaskList::addTask(text);
              }
              else {
                text  = "Check IL/XL specification: Specified IL- or XL-step is not an integer multiple\n";
                text += "   of those found in the seismic data.\n";
                TaskList::addTask(text);
              }
            }
            else if(aligned == true) {
              text  = "Check IL/XL specification: Either start or end of IL and/or XL interval does not\n";
              text += "   align with IL/XL in seismic data, or end IL and/or XL is not an integer multiple\n";
              text += "   of steps away from the start.\n";
              TaskList::addTask(text);
            }
          }
          catch (NRLib::Exception & e) {
            err_text += "Error: "+std::string(e.what());
            geometry->WriteILXL(true);
            geometry = NULL;
            failed = true;
          }
        }
        delete full_geometry;
      }
      else {
        geometry->WriteILXL();
      }
      if(!failed) {
        model_settings->setAreaParameters(geometry);
        ILXL_geometry = geometry;
      }
    }
    else {
      err_text += tmp_err_text;
      failed = true;
    }
  }

  // SET THE AREA FOR THE ESTIMATION SIMBOX

  if(!failed){
    const SegyGeometry * area_params = model_settings->getAreaParameters();
    failed = estimation_simbox.setArea(area_params, err_text);
    if(failed)
    {
      writeAreas(area_params,&estimation_simbox_,area_type);
      err_text += "The specified AREA extends outside the surface(s).\n";
    }
  }

  // SET TOP AND BASE SURFACES FOR THE FORWARD/INVERSION INTERVALS -----------------------------------------------

  // if multiple intervals
  if(model_settings->getIntervalNames().size() > 0){
    setSurfacesMultipleIntervals(estimation_simbox, full_inversion_volume, input_files, model_settings, err_text, failed);
  }
  // single interval described by either one or two surfaces
  else{
    setSurfacesSingleInterval(estimation_simbox, full_inversion_volume, input_files->getTimeSurfFiles(), model_settings, err_text, failed);
  }

  return (!failed);
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


          float top_grid = static_cast<float>(estimation_simbox_.getTopZMin());
          float bot_grid = static_cast<float>(estimation_simbox_.getBotZMax());

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


          //DOES NOT WORK WITH "Simbox estimation_simbox_"
          //segy->ReadAllTraces(estimation_simbox_, //timeCutSimBox
          //                    padding,
          //                    onlyVolume,
          //                    relativePadding);
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
  
bool CommonData::setupEstimationRockPhysics(){
  return true;
}

int
CommonData::computeTime(int year, int month, int day) const
{
  if(year == IMISSING)
    return(0);

  int delta_year = year-1900; //Ok baseyear.
  int time = 365*delta_year+delta_year/4; //Leap years.
  if(month == IMISSING)
    time += 182;
  else {
    std::vector<int> acc_days(12,0);
    acc_days[1]  = acc_days[0]  + 31;
    acc_days[2]  = acc_days[1]  + 28;
    acc_days[3]  = acc_days[2]  + 31;
    acc_days[4]  = acc_days[3]  + 30;
    acc_days[5]  = acc_days[4]  + 31;
    acc_days[6]  = acc_days[5]  + 30;
    acc_days[7]  = acc_days[6]  + 31;
    acc_days[8]  = acc_days[7]  + 31;
    acc_days[9]  = acc_days[8]  + 30;
    acc_days[10] = acc_days[9]  + 31;
    acc_days[11] = acc_days[10] + 30;

    time += acc_days[month];

    if(day == IMISSING)
      time += 15;
    else
      time += day;
  }
  return(time);
}

void
CommonData::getGeometryFromGridOnFile(const std::string           grid_file,
                                      const TraceHeaderFormat   * thf,
                                      SegyGeometry             *& geometry,
                                      std::string               & err_text)
{
  geometry = NULL;

  if(grid_file != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int file_type = IO::findGridType(grid_file);
    if(file_type == IO::CRAVA) {
      geometry = geometryFromCravaFile(grid_file);
    }
    else if(file_type == IO::SEGY) {
      try
      {
        geometry = SegY::FindGridGeometry(grid_file, thf);
      }
      catch (NRLib::Exception & e)
      {
        err_text = e.what();
      }
    }
    else if(file_type == IO::STORM)
      geometry = geometryFromStormFile(grid_file, err_text);
    else if(file_type==IO::SGRI) {
      bool scale = true;
      geometry = geometryFromStormFile(grid_file, err_text, scale);
    }
    else {
      err_text = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    err_text = "Cannot get geometry from file. The file name is empty.\n";
  }
}

SegyGeometry * CommonData::geometryFromCravaFile(const std::string & file_name)
{
  std::ifstream bin_file;
  NRLib::OpenRead(bin_file, file_name, std::ios::in | std::ios::binary);

  std::string file_type;
  getline(bin_file,file_type);

  double x0      = NRLib::ReadBinaryDouble(bin_file);
  double y0      = NRLib::ReadBinaryDouble(bin_file);
  double dx      = NRLib::ReadBinaryDouble(bin_file);
  double dy      = NRLib::ReadBinaryDouble(bin_file);
  int    nx      = NRLib::ReadBinaryInt(bin_file);
  int    ny      = NRLib::ReadBinaryInt(bin_file);
  double IL0     = NRLib::ReadBinaryDouble(bin_file);
  double XL0     = NRLib::ReadBinaryDouble(bin_file);
  double il_step_X = NRLib::ReadBinaryDouble(bin_file);
  double il_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_X = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double rot     = NRLib::ReadBinaryDouble(bin_file);

  bin_file.close();

  SegyGeometry * geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                             IL0, XL0, il_step_X, il_step_Y,
                                             xl_step_X, xl_step_Y, rot);
  return(geometry);
}

SegyGeometry * CommonData::geometryFromStormFile(const std::string & file_name,
                                                 std::string       & err_text,
                                                 bool                scale)
{
  SegyGeometry  * geometry  = NULL;
  StormContGrid * storm_grid = NULL;
  std::string     tmp_err_text;
  float scale_hor;
  if(scale==false)
  {
    scale_hor = 1.0;
  }
  else //from sgri file
  {
    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scale_hor  = 1000.0;
  }
  try
  {
    storm_grid = new StormContGrid(0,0,0);
    storm_grid->ReadFromFile(file_name);
    storm_grid->SetMissingCode(RMISSING);
  }
  catch (NRLib::Exception & e)
  {
    tmp_err_text = e.what();
  }

  if (tmp_err_text == "") {
    double x0      = storm_grid->GetXMin()*scale_hor;
    double y0      = storm_grid->GetYMin()*scale_hor;
    double dx      = storm_grid->GetDX()*scale_hor;
    double dy      = storm_grid->GetDY()*scale_hor;
    int    nx      = static_cast<int>(storm_grid->GetNI());
    int    ny      = static_cast<int>(storm_grid->GetNJ());
    double rot     = storm_grid->GetAngle();
    double IL0     = 0.0;  ///< Dummy value since information is not contained in format
    double XL0     = 0.0;  ///< Dummy value since information is not contained in format
    double IL_step_X =   1;  ///< Dummy value since information is not contained in format
    double IL_step_Y =   1;  ///< Dummy value since information is not contained in format
    double XL_step_X =   1;  ///< Dummy value since information is not contained in format
    double XL_step_Y =   1;  ///< Dummy value since information is not contained in format
    geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                IL0, XL0, IL_step_X, IL_step_Y,
                                XL_step_X, XL_step_Y, rot);
  }
  else {
    err_text += tmp_err_text;
  }

  if (storm_grid != NULL)
    delete storm_grid;

  return(geometry);
}

void CommonData::writeAreas(const SegyGeometry * area_params,
                            Simbox             * time_simbox,
                            std::string        & text)
{
  double area_x0   = area_params->GetX0();
  double area_y0   = area_params->GetY0();
  double area_lx   = area_params->Getlx();
  double area_ly   = area_params->Getly();
  double area_dx   = area_params->GetDx();
  double area_dy   = area_params->GetDy();
  double area_rot  = area_params->GetAngle();
  double area_xmin = RMISSING;
  double area_xmax = RMISSING;
  double area_ymin = RMISSING;
  double area_ymax = RMISSING;

  findSmallestSurfaceGeometry(area_x0, area_y0, area_lx, area_ly, area_rot,
                              area_xmin, area_ymin, area_xmax, area_ymax);

  LogKit::LogFormatted(LogKit::Low,"\nThe top and/or base time surfaces do not cover the area specified by the "+text);
  LogKit::LogFormatted(LogKit::Low,"\nPlease extrapolate surfaces or specify a smaller AREA in the model file.\n");
  LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx        ly     azimuth          dx      dy\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
  double azimuth = (-1)*area_rot*(180.0/M_PI);
  if (azimuth < 0)
    azimuth += 360.0;
  LogKit::LogFormatted(LogKit::Low,"Model area       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n\n",
                       area_x0, area_y0, area_lx, area_ly, azimuth, area_dx, area_dy);

  LogKit::LogFormatted(LogKit::Low,"Area                    xmin         xmax           ymin        ymax\n");
  LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %11.2f %11.2f\n",
                       text.c_str(),area_xmin, area_xmax, area_ymin, area_ymax);
  const NRLib::Surface<double> & top  = time_simbox->GetTopSurface();
  const NRLib::Surface<double> & base = time_simbox->GetBotSurface();
  LogKit::LogFormatted(LogKit::Low,"Top surface      %11.2f  %11.2f    %11.2f %11.2f\n",
                       top.GetXMin(), top.GetXMax(), top.GetYMin(), top.GetYMax());
  LogKit::LogFormatted(LogKit::Low,"Base surface     %11.2f  %11.2f    %11.2f %11.2f\n",
                       base.GetXMin(), base.GetXMax(), base.GetYMin(), base.GetYMax());
}


void CommonData::findSmallestSurfaceGeometry(const double   x0,
                                             const double   y0,
                                             const double   lx,
                                             const double   ly,
                                             const double   rot,
                                             double       & x_min,
                                             double       & y_min,
                                             double       & x_max,
                                             double       & y_max)
{
  x_min = x0 - ly*sin(rot);
  x_max = x0 + lx*cos(rot);
  y_min = y0;
  y_max = y0 + lx*sin(rot) + ly*cos(rot);
  if (rot < 0) {
    x_min = x0;
    x_max = x0 + lx*cos(rot) - ly*sin(rot);
    y_min = y0 + lx*sin(rot);
    y_max = y0 + ly*cos(rot);
  }
}

void CommonData::setSurfacesSingleInterval(Simbox                           & estimation_simbox,
                                           NRLib::Volume                    & full_inversion_volume,
                                           const std::vector<std::string>   & surf_file,
                                           ModelSettings                    * model_settings,
                                           std::string                      & err_text,
                                           bool                             & failed)
{
  const std::string & top_surface_file_name = surf_file[0];

  //bool   generate_seismic    = model_settings->getForwardModeling();
  //bool   estimation_mode     = model_settings->getEstimationMode();
  //bool   generate_background = model_settings->getGenerateBackground();
  bool   parallel_surfaces   = model_settings->getParallelTimeSurfaces();
  //int    nz                 = model_settings->getTimeNz();
  //int    output_format       = model_settings->getOutputGridFormat();
  //int    output_domain       = model_settings->getOutputGridDomain();
  //int    output_grids_elastic = model_settings->getOutputGridsElastic();
  //int    output_grids_other   = model_settings->getOutputGridsOther();
  //int    output_grids_seismic = model_settings->getOutputGridsSeismic();
  double d_top               = model_settings->getTimeDTop();
  double lz                 = model_settings->getTimeLz();
  //double dz                 = model_settings->getTimeDz();

  Surface * top_surface = NULL;
  Surface * base_surface = NULL;

  Surface * top_surface_flat = NULL;
  Surface * base_surface_flat = NULL;

  try {
    double x_min, x_max;
    double y_min, y_max;
    findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                  estimation_simbox.getlx(), estimation_simbox.getly(),
                                  estimation_simbox.getAngle(),
                                  x_min,y_min,x_max,y_max);
    if (NRLib::IsNumber(top_surface_file_name)) {
      // Find the smallest surface that covers the simbox. For simplicity
      // we use only four nodes (nx=ny=2).
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(top_surface_file_name.c_str()));
    }
    else {
      top_surface = new Surface(top_surface_file_name);
    }
  }
  catch (NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if(!failed) {
    if(parallel_surfaces) { //Only one reference surface
      //simbox->setDepth(*top_surface_flat, d_top, lz, dz, model_settings->getRunFromPanel());
      top_surface->Add(d_top);
      base_surface = new Surface(*top_surface);
      base_surface->Add(lz);
      //full_inversion_volume->SetSurfaces(*top_surface, *base_surface, model_settings->getRunFromPanel());
    }
    else { //Two reference surfaces
      const std::string & base_surface_file_name = surf_file[1];
      try {
        if (NRLib::IsNumber(base_surface_file_name)) {
          // Find the smallest surface that covers the simbox. For simplicity
          // we use only four nodes (nx=ny=2).
          double x_min, x_max;
          double y_min, y_max;
          findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                      estimation_simbox.getlx(), estimation_simbox.getly(),
                                      estimation_simbox.getAngle(),
                                      x_min, y_min, x_max, y_max);
          base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(base_surface_file_name.c_str()));
        }
        else {
          base_surface = new Surface(base_surface_file_name);
        }
      }
      catch (NRLib::Exception & e) {
        err_text += e.what();
        failed = true;
      }
      
    }
    /*
    if (!failed) {
      if((output_domain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        simbox->setTopBotName(topSurf,baseSurf,output_format);
        if (generate_seismic) {
          simbox->writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   output_format);
        }
        else if (!estimation_mode){
          if (output_grids_elastic > 0 || output_grids_other > 0 || output_grids_seismic > 0)
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToInversionResults(),
                                     output_format);
        }
        if((output_format & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((output_grids_elastic & IO::BACKGROUND) > 0 ||
              (output_grids_elastic & IO::BACKGROUND_TREND) > 0 ||
              (estimation_mode && generate_background)) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToBackground(),
                                     output_format);
          }
          if ((output_grids_other & IO::CORRELATION) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     output_format);
          }
          if ((output_grids_seismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     output_format);
          }
          if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     output_format);
          }
        }
      }
    } */
  }
  if(!failed){
    try{ //initialize full_inversion_volume and set the flat top and base surfaces of the simbox

      full_inversion_volume = NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(2,2);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(2,2);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(), 
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        2, 2, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        2, 2, z_max);
      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, 2, model_settings->getRunFromPanel());
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }
  delete top_surface;
  delete base_surface;
  delete top_surface_flat;
  delete base_surface_flat;
}

void CommonData::setSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
                                              NRLib::Volume                  & full_inversion_volume,
                                              const InputFiles               * input_files,
                                              const ModelSettings            * model_settings,
                                              std::string                    & err_text,
                                              bool                           & failed){

  // Get interval surface data ------------------------------------------------------------------------------

  const std::vector<std::string>            interval_names                         =  model_settings->getIntervalNames();
  unsigned int                              n_intervals                            =  interval_names.size();
  const std::map<std::string, std::string>  interval_base_time_surfaces            =  input_files->getIntervalBaseTimeSurfaces();
  // implicit assumption that the top surface is the first one given and the base surface is the last one
  const std::string                         top_surface_file_name                  =  input_files->getTimeSurfFile(0);
  const std::string                         base_surface_file_name                 =  interval_base_time_surfaces.find( interval_names[n_intervals-1])->second;

  Surface * top_surface = NULL; 
  Surface * base_surface = NULL;

  Surface * top_surface_flat = NULL;
  Surface * base_surface_flat = NULL;

  try{
    if(NRLib::IsNumber(top_surface_file_name)){
      double x_min, x_max;
      double y_min, y_max;
      findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                  estimation_simbox.getlx(), estimation_simbox.getly(),
                                  estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(top_surface_file_name.c_str()));
    }
    else{
      top_surface = new Surface(top_surface_file_name);
    }

  }
  catch(NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if(!failed){
    try{
      if(NRLib::IsNumber(base_surface_file_name)){
        double x_min, x_max;
        double y_min, y_max;
        findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                    estimation_simbox.getlx(), estimation_simbox.getly(),
                                    estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
        base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(base_surface_file_name.c_str()));
      }
      else{
        base_surface = new Surface(base_surface_file_name);
      }
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }
  if(!failed){
    try{ // initialize full_inversion_volume and set the flat top and base surfaces of the simbox
      full_inversion_volume =  NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(2,2);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(2,2);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(), 
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        2, 2, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        2, 2, z_max);
      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, 2, model_settings->getRunFromPanel());
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }

  delete top_surface;
  delete base_surface;
  delete top_surface_flat;
  delete base_surface_flat;
  
}
