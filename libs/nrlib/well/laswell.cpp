// $Id: laswell.cpp 1071 2012-09-18 11:42:51Z perroe $

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

#include "laswell.hpp"

#include <cassert>
#include <fstream>
#include <string>
#include <vector>

#include "../iotools/fileio.hpp"

using namespace NRLib;


LasWell::LasWell(const std::string & name,
                 size_t              length,
                 double              start_depth,
                 double              depth_increment,
                 DepthType           depth_type,
                 const std::string & depth_unit,
                 const std::string & company_name,
                 const std::string & field_name,
                 const std::string & location,
                 const std::string & date,
                 const std::string & unique_well_identifier)
  : Well(name, -999.25),
    version_("2.0"),
    wrap_(false),
    depth_unit_(depth_unit),
    start_depth_(start_depth),
    stop_depth_(start_depth + (length-1)*depth_increment),
    depth_increment_(depth_increment),
    company_name_(company_name),
    field_name_(field_name),
    location_(location),
    county_("County"),
    state_("State"),
    country_("Country"),
    service_company_("Service company"),
    date_(date),
    unique_well_identifier_(unique_well_identifier)
{
  std::string depth_mnemonic;
  if (depth_type == LasWell::Depth)
    depth_mnemonic = "DEPT";
  else
    depth_mnemonic = "TIME";

  std::vector<double> depth_log(length);
  for (size_t i = 0; i < length; ++i) {
    depth_log[i] = start_depth + i*depth_increment;
  }

  AddLog(depth_mnemonic, depth_unit, "", depth_log);
}


void LasWell::AddLog(const std::string         & name,
                     const std::string         & units,
                     const std::string         & comment,
                     const std::vector<double> & log)
{
  Well::AddContLog(name, log);
  log_name_.push_back(name);
  log_units_.push_back(units);
  log_comment_.push_back(comment);
}


void LasWell::WriteToFile(const std::string              & filename,
                          const std::vector<std::string> & comment_header)
{
  std::ofstream file;
  OpenWrite(file, filename);

  // Comment header
  for (size_t i = 0; i < comment_header.size(); ++i) {
    file << "# " << comment_header[i] << "\n";
  }
  file << "\n";

  // Version information
  file << "~Version information\n";
  WriteLasLine(file, "VERS", "", version_, "");
  WriteLasLine(file, "WRAP", "", (wrap_ ? "YES" : "NO"), "");
  file << "\n";

  // Well information
  file << "~Well information\n";
  WriteLasLine(file, "STRT", depth_unit_, start_depth_, "");
  WriteLasLine(file, "STOP", depth_unit_, stop_depth_, "");
  WriteLasLine(file, "STEP", depth_unit_, depth_increment_, "");
  WriteLasLine(file, "NULL", "", GetContMissing(), "");
  WriteLasLine(file, "COMP", "", company_name_, "");
  WriteLasLine(file, "WELL", "", GetWellName(), "");
  WriteLasLine(file, "FLD" , "", field_name_, "");
  WriteLasLine(file, "LOC" , "", location_, "");
  WriteLasLine(file, "CNTY", "", county_, "");
  WriteLasLine(file, "STAT", "", state_, "");
  WriteLasLine(file, "CTRY", "", country_, "");
  WriteLasLine(file, "SRCV", "", service_company_, "");
  WriteLasLine(file, "DATE", "", date_, "");
  WriteLasLine(file, "UWI" , "", unique_well_identifier_, "");
  file << "\n";

  if (GetNContLog() == 0) {
    // No log data in file.
    return;
  }

  // Curve information
  // LAS only supports continuous logs...

  file << "~Curve information\n";
  for (size_t i = 0; i < GetNContLog(); ++i) {
    WriteLasLine(file, log_name_[i], log_units_[i], "", log_comment_[i]);
  }
  file << "\n";

  // Data section
  file << "~Ascii\n";

  std::vector<const std::vector<double> *> logs(GetNContLog());
  for (size_t i = 0; i < GetNContLog(); ++i) {
    logs[i] = &(GetContLog(log_name_[i]));
  }

  // We don't support wrapped output yet.
  assert(wrap_ == false);

  file.precision(3);
  file.setf(std::ios_base::fixed);

  for (size_t i = 0; i < logs[0]->size(); ++i) {
    for (size_t j = 0; j < logs.size(); ++j) {
      file << (*logs[j])[i] << " ";
    }
    file << "\n";
  }
}


void LasWell::WriteLasLine(std::ofstream     & file,
                           const std::string & mnemonic,
                           const std::string & units,
                           const std::string & data,
                           const std::string & description)
{
  file << mnemonic << " ." << units << " " << data << " : " << description << "\n";
}


void LasWell::WriteLasLine(std::ofstream     & file,
                           const std::string & mnemonic,
                           const std::string & units,
                           double              data,
                           const std::string & description)
{
  file.precision(3);
  file.setf(std::ios_base::fixed);

  file << mnemonic << " ." << units << " " << data << " : " << description << "\n";
}

