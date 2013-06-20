// $Id: well.hpp 883 2011-09-26 09:17:05Z perroe $

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

#ifndef NRLIB_WELL_HPP
#define NRLIB_WELL_HPP

#include <vector>
#include <sstream>
#include <map>


namespace NRLib {
  class Well{
  public:
    /// Default constructor
    Well();

    /// Construct well with given name and no logs.
    Well(const std::string & name,
         double              rmissing = -999.0,
         int                 imissing = -999);

    /// Construct well from file
    Well(const std::string & file_name,
         bool & read_ok);

    /// Constructor
    /// \param[in] cont_log Continuous logs
    /// \param[in] disc_log Discrete logs
    /// \param[in] well_name
    Well(const std::map<std::string, std::vector<double> > & cont_log,
         const std::map<std::string, std::vector<int> >    & disc_log,
         const std::string                                 & well_name);

    /// Destructor
    ~Well();

    /// Check existence of continuous log
    bool HasContLog(const std::string& name) const;

    /// Return continuous logs
    std::vector<double> & GetContLog(const std::string& name);
    /// Return continuous logs
    const std::vector<double> & GetContLog(const std::string& name) const;

    /// Return discrete logs
    std::vector<int> & GetDiscLog(const std::string& name);
    /// Return discrete logs
    const std::vector<int> & GetDiscLog(const std::string& name) const;

    /// Read a well from a file-name
    void ReadWell(const std::string&file_name,
                  bool read_ok);

    /// Add a continuous log
    /// Replaces the log if there is already a log with the given name.
    void AddContLog(const std::string& name, const std::vector<double>& log);
    /// Remove continuous log
    /// Does nothing if there is no log with the given name.
    void RemoveContLog(const std::string& name);
    /// Add discrete log
    /// Replaces the log if there is already a log with the given name.
    void AddDiscLog(const std::string& name, const std::vector<int>& log);
    /// Remove discrete log
    /// Does nothing if there is no log with the given name.
    void RemoveDiscLog(const std::string& name);
    ///
    void SetWellName(const std::string& wellname);
    ///
    const std::string& GetWellName() const { return well_name_; };
    /// Return true if x is missing
    bool IsMissing(double x) const;
    /// Return true if n is missing
    bool IsMissing(int n) const;
    /// Return cont. missing value
    double GetContMissing() const { return(well_rmissing_); }
    /// Return disc. missing value
    int GetIntMissing() const { return(well_imissing_); }
    /// Set missing values
    void SetMissing(double value) {well_rmissing_ = value; well_imissing_ = static_cast<int>(value);}
    /// Return discrete value at position index in log with name logname
    /// Returns missing if there is no log with the given name, or index is out of range.
    int GetDiscValue(size_t index, const std::string& logname) const;
    /// Return continuous value at position index in log with name logname
    /// Returns missing if there is no log with the given name, or index is out of range.
    double GetContValue(size_t index, const std::string& logname) const;
    /// Set value at position index in log with name logname
    void SetDiscValue(int value, size_t index, const std::string& logname);
    /// Set value at position index in log with name logname
    void SetContValue(double value, size_t index, const std::string& logname);
    /// Return total number of logs
    size_t GetNlog() const;
    /// Return number of discrete logs
    size_t GetNContLog() const;
    /// Return length of log with name logname
    size_t GetContLogLength(const std::string& logname) const;
    /// Return all continuous logs
    const std::map<std::string,std::vector<double> > & GetContLog() const { return cont_log_; };
    /// Return all discrete logs
    const std::map<std::string,std::vector<int> > & GetDiscLog() const { return disc_log_; };

  private:
    /// Continuous logs
    std::map<std::string,std::vector<double> > cont_log_;
    /// Discrete logs
    std::map<std::string,std::vector<int> >    disc_log_;
    /// Name of well
    std::string well_name_;

    /// Missing value for continous logs.
    double well_rmissing_;
    /// Missing value for discrete logs.
    int    well_imissing_;
  };

}

#endif
