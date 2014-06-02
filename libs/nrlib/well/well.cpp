// $Id: well.cpp 882 2011-09-23 13:10:16Z perroe $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include <cassert>
#include <fstream>
#include <string>

#include "well.hpp"
#include "norsarwell.hpp"
#include "rmswell.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/surface/surface.hpp"
//#include "src/definitions.h"

using namespace NRLib;

Well::Well()
{
  well_rmissing_ = -999.0;
  well_imissing_ = -999;
}


Well::Well(const std::string & name,
           double              rmissing,
           int                 imissing)
{
  well_name_ = name;
  well_rmissing_ = rmissing;
  well_imissing_ = imissing;
}


Well::Well(const std::string & file_name,
           bool              & read_ok,
           const std::string & facies_log)
{
  ReadWell(file_name, read_ok, facies_log);
}


Well::Well(const std::map<std::string,std::vector<double> > & cont_log,
           const std::map<std::string,std::vector<int> >    & disc_log,
           const std::string                                & well_name)
{
  cont_log_       = cont_log;
  disc_log_       = disc_log;
  well_name_      = well_name;
  well_rmissing_  = -999.0;
  well_imissing_  = -999;
}

Well::~Well(){

}


void
Well::ReadWell(const std::string              & file_name,
               bool                           & read_ok,
               const std::string              & facies_log)
{
  if(file_name.find(".nwh",0) != std::string::npos) {
    NRLib::NorsarWell well(file_name);
    well_name_ = well.GetWellName();
    cont_log_  = well.GetContLog();
    n_data_    = well.GetNData();
    disc_log_  = well.GetDiscLog();
    read_ok    = true;

    //Norsar wells has its own missing value in the well-log
    SetMissing(well.GetContMissing());

    // assume that facies logs from norsar wells are not used
    has_facies_log_ = false;
  }
  else if(file_name.find(".rms",0) != std::string::npos || file_name.find(".txt",0) != std::string::npos) {
    NRLib::RMSWell well(file_name);
    well_name_ = well.GetWellName();
    n_data_    = well.GetNData();
    cont_log_  = well.GetContLog();
    disc_log_  = well.GetDiscLog();
    x_pos0_    = well.GetXPos0();
    y_pos0_    = well.GetYPos0();

    std::map<int, std::string> facies_map;
    if(well.HasDiscLog(facies_log))
      facies_map = well.GetDiscNames(facies_log);

    if(facies_map.size() > 0){
      has_facies_log_ = true;
      facies_map_ = facies_map;
    }
    else{
      has_facies_log_ = false;
    }
    read_ok = true;
  }
  else
    read_ok = false;
}


void
Well::AddContLog(const std::string& name, const std::vector<double>& log)
{
  cont_log_[name] = log;
}

void
Well::AddContLogSeismicResolution(const std::string& name, const std::vector<double>& log)
{
  cont_log_seismic_resolution_[name] = log;
}

void
Well::AddContLogBackgroundResolution(const std::string& name, const std::vector<double>& log)
{
  cont_log_background_resolution_[name] = log;
}


bool Well::HasDiscLog(const std::string& name) const{
  std::map<std::string, std::vector<int> >::const_iterator it = disc_log_.find(name);
  if(it != disc_log_.end()){
    return true;
  }
  else{
    return false;
  }
}

bool
Well::HasContLog(const std::string& name) const
{
  //std::map<std::string, std::vector<double> >::const_iterator item = cont_log_.find(name);
  if(cont_log_.find(name) != cont_log_.end())
    return true;
  else
    return false;
}


std::vector<double>&
Well::GetContLog(const std::string& name)
{
  std::map<std::string, std::vector<double> >::iterator item = cont_log_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}

std::vector<double>&
Well::GetContLogSeismicResolution(const std::string& name)
{
  std::map<std::string, std::vector<double> >::iterator item = cont_log_seismic_resolution_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}

std::vector<double>&
Well::GetContLogBackgroundResolution(const std::string& name)
{
  std::map<std::string, std::vector<double> >::iterator item = cont_log_background_resolution_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}


const std::vector<double>&
Well::GetContLog(const std::string& name) const
{
  std::map<std::string, std::vector<double> >::const_iterator item = cont_log_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}

const std::vector<double>&
Well::GetContLogSeismicResolution(const std::string& name) const
{
  std::map<std::string, std::vector<double> >::const_iterator item = cont_log_seismic_resolution_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}

const std::vector<double>&
Well::GetContLogBackgroundResolution(const std::string& name) const
{
  std::map<std::string, std::vector<double> >::const_iterator item = cont_log_background_resolution_.find(name);
  assert(item != cont_log_.end());
  return item->second;
}


void
Well::RemoveContLog(const std::string& name)
{
  cont_log_.erase(name);
}


void
Well::AddDiscLog(const std::string& name, const std::vector<int>& log)
{
  disc_log_[name] = log;
}


std::vector<int> &
Well::GetDiscLog(const std::string& name)
{
  std::map<std::string, std::vector<int> >::iterator item = disc_log_.find(name);
  assert(item != disc_log_.end());
  return item->second;
}


const std::vector<int> &
Well::GetDiscLog(const std::string& name) const
{
  std::map<std::string, std::vector<int> >::const_iterator item = disc_log_.find(name);
  assert(item != disc_log_.end());
  return item->second;
}


void
Well::RemoveDiscLog(const std::string& name)
{
  disc_log_.erase(name);
}


void Well::SetWellName(const std::string& wellname)
{
  well_name_ = wellname;
}


bool Well::IsMissing(double x)const
{
  if (x == well_rmissing_)
    return true;
  else
    return false;
}

bool Well::IsMissing(int n)const
{
  if (n == well_imissing_)
    return true;
  else
    return false;
}


int Well::GetDiscValue(size_t index, const std::string& logname) const
{
  std::map<std::string, std::vector<int> >::const_iterator item = disc_log_.find(logname);
  assert(item != disc_log_.end());
  const std::vector<int>& log = item->second;
  assert(index < log.size());
  return log[index];
}


double Well::GetContValue(size_t index, const std::string& logname) const
{
  std::map<std::string,std::vector<double> >::const_iterator item = cont_log_.find(logname);
  assert(item != cont_log_.end());
  const std::vector<double>& log = item->second;
  assert(index < log.size());
  return log[index];
}


void Well::SetDiscValue(int value, size_t index, const std::string& logname)
{
  std::map<std::string,std::vector<int> >::iterator item = disc_log_.find(logname);
  assert(item != disc_log_.end());
  assert(index < item->second.size());
  (item->second)[index] = value;
}


void Well::SetContValue(double value, size_t index, const std::string& logname)
{
  std::map<std::string,std::vector<double> >::iterator item = cont_log_.find(logname);
  assert(item != cont_log_.end());
  assert(index < item->second.size());
  (item->second)[index] = value;
}


size_t Well::GetNlog() const
{
  return cont_log_.size() + disc_log_.size();
}


size_t Well::GetNContLog() const
{
  return cont_log_.size();
}


size_t Well::GetContLogLength(const std::string& logname) const
{
  std::map<std::string,std::vector<double> >::const_iterator item = cont_log_.find(logname);
  assert(item != cont_log_.end());

  return (item->second).size();
}

int Well::CheckStormgrid(NRLib::StormContGrid & stormgrid) const
{
  bool inside_area = false;
  int  error       = 1;

  //if(timemissing_ == 0) { //Timemissing not stored when reading logs in commondata ReadWell?

    NRLib::Surface<double> & top  = stormgrid.GetTopSurface();
    NRLib::Surface<double> & base = stormgrid.GetBotSurface();

    const std::vector<double> & x_pos = GetContLog().find("X_pos")->second;
    const std::vector<double> & y_pos = GetContLog().find("Y_pos")->second;
    const std::vector<double> & z_pos = GetContLog().find("Z_pos")->second;

    for(unsigned int i=0; i < n_data_nonmissing_; i++) { //nd_

      if(stormgrid.IsInside(x_pos[i], y_pos[i])) {
        inside_area = true;
        //
        // Correct handling of top and base checking.
        //
        double z_top  = top.GetZ(x_pos[i], y_pos[i]);
        double z_base = base.GetZ(x_pos[i], y_pos[i]);

        if(z_pos[i] > z_top && z_pos[i] < z_base) {
          error = 0;
          break;
        }
      }
    }
  //}
  //else
  //  error = 0;

  if (error) {
    if (inside_area) {
      LogKit::LogFormatted(LogKit::Low," \nWell "+well_name_+": ");
      LogKit::LogFormatted(LogKit::Low,"IGNORED. Well does not hit the inversion volume.\n");
    }
    else {
      LogKit::LogFormatted(LogKit::Low," \nWell "+well_name_+": ");
      LogKit::LogFormatted(LogKit::Low,"IGNORED. Well is not inside the inversion area.\n");
    }
  }
  return(error);
}


int Well::CheckSimbox(Simbox * simbox) const
{
  bool inside_area = false;
  int error        = 1;

  const std::vector<double> & x_pos = GetContLog().find("X_pos")->second;
  const std::vector<double> & y_pos = GetContLog().find("Y_pos")->second;
  const std::vector<double> & z_pos = GetContLog().find("Z_pos")->second;

  for(size_t i = 0; i < n_data_nonmissing_; i++) {
    if(simbox->isInside(x_pos[i], y_pos[i])) {
      inside_area = true;
      //
      // Correct handling of top and base checking.
      //
      if(z_pos[i] > simbox->getTop(x_pos[i], y_pos[i]) && z_pos[i] < simbox->getBot(x_pos[i], y_pos[i])) {
        error = 0;
        break;
      }
    }
  }

  if (error) {
    if (inside_area) {
      LogKit::LogFormatted(LogKit::Low," \nWell "+well_name_+":\n");
      LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is inside inversion area but does not hit the inversion volume)\n");
      LogKit::LogFormatted(LogKit::Low,"           (well-depth: min,max = "+NRLib::ToString(z_pos[0])+","+NRLib::ToString(z_pos[n_data_nonmissing_-1])+")\n");
    }
    else {
      LogKit::LogFormatted(LogKit::Low," \nWell "+well_name_+":\n");
      LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is not inside inversion area)\n");
    }
  }
  return(error);
}

//----------------------------------------------------------------------------
void Well::WriteWell(int         well_format,
                     const float max_hz_background,
                     const float max_hz_seismic) const
{
  if ((well_format & IO::RMSWELL) > 0)
    WriteRMSWell(max_hz_background,
                 max_hz_seismic);
  if ((well_format & IO::NORSARWELL) > 0)
    WriteNorsarWell(max_hz_background,
                    max_hz_seismic);
}

void Well::WriteRMSWell(const float max_hz_background,
                        const float max_hz_seismic) const
{

  std::string well_name(well_name_);
  NRLib::Substitute(well_name,"/","_");
  NRLib::Substitute(well_name," ","_");
  std::string base_name      = IO::PrefixWells() + well_name + IO::SuffixRmsWells();
  std::string well_file_name = IO::makeFullFileName(IO::PathToWells(), base_name);

  std::ofstream file;
  NRLib::OpenWrite(file, well_file_name);

  unsigned int n_logs = 1+3*3;
  //if (nFacies_ > 0)
  if (has_facies_log_)
    n_logs += 1;

  std::vector<std::string> params(3);
  params[0] = "Vp";
  params[1] = "Vs";
  params[2] = "Rho";

  //
  // Write HEADER
  //
  file << "1.0\n"
       << "CRAVA\n"
       << well_name_ << " "
       << std::fixed
       << std::setprecision(2)
       << x_pos0_ << " "
       << y_pos0_ << "\n"
       << n_logs << "\n"
       << "Time  UNK lin"
       << std::endl;

  for (int i =0 ; i < 3 ; i++) {
    file << params[i] <<                                       "  UNK lin\n"
         << params[i] << static_cast<int>(max_hz_background) << "  UNK lin\n"
         << params[i] << static_cast<int>(max_hz_seismic)    << "  UNK lin\n";
  }

  if (has_facies_log_) {
    file << "FaciesLog  DISC "; //faciesLogName_

    for (size_t i = 0; i < facies_map_.size(); i++)
      file << " " << i << " " << facies_map_.find(i)->second;
    file << "\n";
  }

  const std::vector<double> & x_pos = cont_log_.find("X_pos")->second;
  const std::vector<double> & y_pos = cont_log_.find("Y_pos")->second;
  const std::vector<double> & z_pos = cont_log_.find("Z_pos")->second;

  const std::vector<double> & vp  = cont_log_.find("Vp")->second;
  const std::vector<double> & vs  = cont_log_.find("Vs")->second;
  const std::vector<double> & rho = cont_log_.find("Rho")->second;

  const std::vector<double> & vp_seismic_resolution  = cont_log_seismic_resolution_.find("Vp")->second;
  const std::vector<double> & vs_seismic_resolution  = cont_log_seismic_resolution_.find("Vs")->second;
  const std::vector<double> & rho_seismic_resolution = cont_log_seismic_resolution_.find("Rho")->second;

  const std::vector<double> & vp_background_resolution  = cont_log_background_resolution_.find("Vp")->second;
  const std::vector<double> & vs_background_resolution  = cont_log_background_resolution_.find("Vs")->second;
  const std::vector<double> & rho_background_resolution = cont_log_background_resolution_.find("Rho")->second;

  const std::vector<int> & facies = disc_log_.find("Facies")->second;

  for (size_t i = 0; i < n_data_; i++) {
    file << std::right
         << std::fixed
         << std::setprecision(4)
         << std::setw(9) << x_pos[i] << " "
         << std::setw(10)<< y_pos[i] << " "
         << std::setw(7) << z_pos[i] << "  "
         << std::setw(7) << z_pos[i] << "  "
         << std::setprecision(5)
         << std::setw(7) << (vp[i]==RMISSING                        ? WELLMISSING : vp[i])                        << " "
         << std::setw(7) << (vp_background_resolution[i]==RMISSING  ? WELLMISSING : vp_background_resolution[i])  << " "
         << std::setw(7) << (vp_seismic_resolution[i]==RMISSING     ? WELLMISSING : vp_seismic_resolution[i])     << " "
         << std::setw(7) << (vs[i]==RMISSING                        ? WELLMISSING : vs[i])                        << " "
         << std::setw(7) << (vs_background_resolution[i]==RMISSING  ? WELLMISSING : vs_background_resolution[i])  << " "
         << std::setw(7) << (vs_seismic_resolution[i]==RMISSING     ? WELLMISSING : vs_seismic_resolution[i])     << " "
         << std::setw(7) << (rho[i]==RMISSING                       ? WELLMISSING : rho[i])                       << " "
         << std::setw(7) << (rho_background_resolution[i]==RMISSING ? WELLMISSING : rho_background_resolution[i]) << " "
         << std::setw(7) << (rho_seismic_resolution[i]==RMISSING    ? WELLMISSING : rho_seismic_resolution[i])    << " ";

      //if (nFacies_ > 0)
      if (has_facies_log_)
        file << std::setw(3) << (facies[i]==IMISSING ? static_cast<int>(WELLMISSING) : facies[i]);

    file << "\n";
  }
  file.close();
}

void Well::WriteNorsarWell(const float max_hz_background,
                           const float max_hz_seismic) const
{
  std::string well_name(well_name_);
  NRLib::Substitute(well_name,"/","_");
  NRLib::Substitute(well_name," ","_");

  //Handle main file.
  std::string base_name = IO::PrefixWells() + well_name + IO::SuffixNorsarWells();
  std::string file_name = IO::makeFullFileName(IO::PathToWells(), base_name);

  const std::vector<double> & x_pos = cont_log_.find("X_pos")->second;
  const std::vector<double> & y_pos = cont_log_.find("Y_pos")->second;
  const std::vector<double> & z_pos = cont_log_.find("Z_pos")->second;

  std::ofstream main_file;
  NRLib::OpenWrite(main_file, file_name);
  main_file << std::fixed
            << std::setprecision(2);

  std::vector<double> md_log;
  bool has_md_log = false;
  if (cont_log_.find("MD") != cont_log_.end()) {
    has_md_log = true;
    md_log = cont_log_.find("MD")->second;
  }

  std::vector<double> md(n_data_, 0); //nd_
  if (has_md_log)
    md[0] = md_log[0];
  double dmax = 0;
  double dmin = 1e+30;
  for (size_t i = 1; i < n_data_; i++) {
    double dx = x_pos[i]-x_pos[i-1];
    double dy = y_pos[i]-y_pos[i-1];
    double dz = z_pos[i]-z_pos[i-1];
    double d  = sqrt(dx*dx+dy*dy+dz*dz);
    if (d > dmax)
      dmax = d;
    else if (d<dmin)
      dmin = d;
    if (has_md_log)
      md[i] = md_log[i];
    else
      md[i] = md[i-1] + d;
  }
  main_file << "[Version information]\nVERSION 1000\nFORMAT ASCII\n\n";
  main_file << "[Well information]\n";
  main_file << "MDMIN      m       " << 0.0f << "\n";
  main_file << "MDMAX      m       " << md[n_data_-1] << "\n";
  main_file << "MDMINSTEP  m       " << dmin << "\n";
  main_file << "MDMAXSTEP  m       " << dmax << "\n";
  main_file << "UTMX       m       " << x_pos[0] << "\n";
  main_file << "UTMY       m       " << y_pos[0] << "\n";
  main_file << "EKB        m       " << 0.0f << "\n";
  main_file << "UNDEFVAL   no_unit " << WELLMISSING << "\n\n";

  main_file << "[Well track data information]\n";
  main_file << "NUMMD  " << n_data_ << "\n";
  main_file << "NUMPAR 4\n";
  main_file << "MD      m\n";
  main_file << "TWT     s\n";
  main_file << "UTMX    m\n";
  main_file << "UTMY    m\n\n";

  std::string log_base_name = IO::PrefixWells() + well_name + IO::SuffixNorsarLog();
  std::string log_file_name = IO::makeFullFileName(IO::PathToWells(), log_base_name);
  std::string only_name     = NRLib::RemovePath(log_file_name);

  //bool got_facies = nFacies_ > 0;

  int n_logs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz}
  if (has_facies_log_)
    n_logs += 1;

  std::vector<std::string> params(3);
  params[0] = "VP";
  params[1] = "VS";
  params[2] = "Rho";

  std::vector<std::string> unit(3);
  unit[0] = "m/s";
  unit[1] = "m/s";
  unit[2] = "g/cm^3";

  main_file << "[Well log data information]\n";
  main_file << "LOGNAME log\n";
  main_file << "IN_FILE " << only_name << "\n";
  main_file << "NUMPAR " << n_logs+1 << "\n"; //Also count md.
  main_file << "NUMLINES " << n_data_ << "\n";
  main_file << "MD      m\n";
  for (int i =0 ; i<3 ; i++) {
    main_file << params[i] << " " << unit[i] << "\n";
    main_file << params[i] << static_cast<int>(max_hz_background) << " " << unit[i] << "\n";
    main_file << params[i] << static_cast<int>(max_hz_seismic)    << " " << unit[i] << "\n";
  }

  if (has_facies_log_)
    main_file << "FACIES no_unit\n";
  main_file.close();


  //Write the two other files.
  std::string base_track_name = IO::PrefixWells() + well_name + IO::SuffixNorsarTrack();
  std::string track_file_name = IO::makeFullFileName(IO::PathToWells(), base_track_name);
  std::ofstream track_file;
  NRLib::OpenWrite(track_file, track_file_name);
  track_file << std::right
             << std::fixed
             << std::setprecision(2)
             << "[NORSAR Well Track]\n";

  //Note: logFileName created above, needed in mainFile.
  std::ofstream log_file;
  NRLib::OpenWrite(log_file, log_file_name);
  log_file << "[NORSAR Well Log]\n";

  const std::vector<double> & vp  = cont_log_.find("Vp")->second;
  const std::vector<double> & vs  = cont_log_.find("Vs")->second;
  const std::vector<double> & rho = cont_log_.find("Rho")->second;

  const std::vector<double> & vp_seismic_resolution  = cont_log_seismic_resolution_.find("Vp")->second;
  const std::vector<double> & vs_seismic_resolution  = cont_log_seismic_resolution_.find("Vs")->second;
  const std::vector<double> & rho_seismic_resolution = cont_log_seismic_resolution_.find("Rho")->second;

  const std::vector<double> & vp_background_resolution  = cont_log_background_resolution_.find("Vp")->second;
  const std::vector<double> & vs_background_resolution  = cont_log_background_resolution_.find("Vs")->second;
  const std::vector<double> & rho_background_resolution = cont_log_background_resolution_.find("Rho")->second;

  const std::vector<int> & facies = disc_log_.find("Facies")->second;

  for (size_t i = 0; i < n_data_; i++) {
    track_file << std::setw(7) << md[i] << " " << std::setw(7) << z_pos[i]
               << " " << std::setw(10)<< x_pos[i] << " " << std::setw(10)<< y_pos[i] << "\n";

    log_file   << std::right << std::fixed << std::setprecision(2)
               << std::setw(7) << md[i] << " "
               << std::setw(7) << (vp[i]==RMISSING                        ? WELLMISSING : vp[i])                        << " "
               << std::setw(7) << (vp_background_resolution[i]==RMISSING  ? WELLMISSING : vp_background_resolution[i])  << " "
               << std::setw(7) << (vp_seismic_resolution[i]==RMISSING     ? WELLMISSING : vp_seismic_resolution[i])     << " "
               << std::setw(7) << (vs[i]==RMISSING                        ? WELLMISSING : vs[i])                        << " "
               << std::setw(7) << (vs_background_resolution[i]==RMISSING  ? WELLMISSING : vs_background_resolution[i])  << " "
               << std::setw(7) << (vs_seismic_resolution[i]==RMISSING     ? WELLMISSING : vs_seismic_resolution[i])     << " "
               << std::setprecision(5)
               << std::setw(7) << (rho[i]==RMISSING                       ? WELLMISSING : rho[i])                       << " "
               << std::setw(7) << (rho_background_resolution[i]==RMISSING ? WELLMISSING : rho_background_resolution[i]) << " "
               << std::setw(7) << (rho_seismic_resolution[i]==RMISSING    ? WELLMISSING : rho_seismic_resolution[i])    << " ";
    if (has_facies_log_)
      log_file << (facies[i]==IMISSING                                 ? static_cast<int>(WELLMISSING) : facies[i])     << "  ";
    log_file << "\n";
  }
  track_file.close();
  log_file.close();
}

