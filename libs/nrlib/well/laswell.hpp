// $Id: laswell.hpp 1071 2012-09-18 11:42:51Z perroe $

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

#ifndef NRLIB_WELL_LASWELL_HPP
#define NRLIB_WELL_LASWELL_HPP

#include <string>
#include <vector>

#include "well.hpp"

namespace NRLib {

/// Support for the LAS well format.
/// Currently only LAS 2.0 is supported.
/// \sa{http://www.cwls.org/las_info.php}
class LasWell : public Well {
public:
  enum DepthType {
    Depth,
    Time
  };

  /// Constructor with relevant parameters.
  LasWell(const std::string & name,
          size_t              length,
          double              start_depth,
          double              depth_increment,
          DepthType           depth_type,
          const std::string & depth_unit,
          const std::string & company_name,
          const std::string & field_name,
          const std::string & location,
          const std::string & date,
          const std::string & unique_well_identifier);

  void AddLog(const std::string         & name,
              const std::string         & units,
              const std::string         & comment,
              const std::vector<double> & log);

  void WriteToFile(const std::string              & filename,
                   const std::vector<std::string> & comment_header);
private:
  void WriteLasLine(std::ofstream     & file,
                    const std::string & mnemonic,
                    const std::string & units,
                    const std::string & data,
                    const std::string & description);

  void WriteLasLine(std::ofstream     & file,
                    const std::string & mnemonic,
                    const std::string & units,
                    double              data,
                    const std::string & description);

  /// Version ID. Possible version ID's include 1.2, 2.0 and 3.0.
  std::string version_;

  /// True if wrap around mode is used. If wrap mode is used the depth value will be
  /// on its own line and all lines of data will be no longer than 80 characters (including CR+LF).
  bool wrap_;

  /// Possible depth units:
  /// For logs in depth: M (meter), F (feet) or FT (feet)
  /// For logs in time:  S (seconds), MS (milliseconds), etc.
  std::string depth_unit_;

  /// Depth of first depth sample.
  double start_depth_;
  /// Depth of last depth sample.
  double stop_depth_;
  /// Depth increment. 0 if the increment is not constant.
  double depth_increment_;

  std::string company_name_;
  std::string field_name_;
  std::string location_;
  std::string county_;
  std::string state_;
  std::string country_;
  std::string service_company_;
  std::string date_;
  std::string unique_well_identifier_;

  std::vector<std::string> log_name_;
  std::vector<std::string> log_units_;
  std::vector<std::string> log_comment_;
};

}

#endif // NRLIB_WELL_LASWELL_HPP
