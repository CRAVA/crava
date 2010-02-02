#include <time.h>

#include "src/definitions.h"
#include "src/program.h"
#include <time.h>
#include "lib/systemcall.h"

#include "nrlib/iotools/logkit.hpp"

Program::Program(const unsigned int  major,
                 const unsigned int  minor,
                 const unsigned int  patch,
                 const std::string & extra_text,
                 const int           licence_days,
                 const std::string & licensed_to)
  : major_(major),
    minor_(minor),
    patch_(patch),
    extra_text_(extra_text),
    licence_days_(licence_days),
    licensed_to_(licensed_to)
{
  std::string version = NRLib::ToString(major)+"."+NRLib::ToString(minor)+"."+NRLib::ToString(patch)+extra_text;
  std::string blanks(39 - version.size(), ' '); 

  LogKit::LogFormatted(LogKit::LOW,"\n***************************************************************************************************");
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                                                                         *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                   C R A V A  -  "+version+blanks+                      "*****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                                                                         *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                     Copyright (c) 2003 by Norsk Regnesentral/Statoil                    *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n*****                                                                                         *****"); 
  LogKit::LogFormatted(LogKit::LOW,"\n***************************************************************************************************\n\n");

  CheckForLicenceExpiration(licence_days_,
                            licensed_to_);

  LogKit::LogFormatted(LogKit::LOW,"Log written by                             : "+SystemCall::getUserName()+"\n");
  LogKit::LogFormatted(LogKit::LOW,"Date and time                              : "+SystemCall::getCurrentTime());
  LogKit::LogFormatted(LogKit::LOW,"Host                                       : "+SystemCall::getHostName()+"\n");
}

Program::~Program(void) 
{
}

void
Program::CheckForLicenceExpiration(const int           licence_days,
                                   const std::string & licensed_to) const
{
  if (licence_days >= 0) {
    time_t now = time(0);
    int days_since_compilation = static_cast<int>((now - TimeOfCompilation()) / 86400);

    if (days_since_compilation > licence_days) {
      LogKit::LogFormatted(LogKit::ERROR,"The CRAVA licence have expired. Please contact the Norwegian Computing Center to obtain a new one.\n");
      std::exit(EXIT_FAILURE);
    }
    else if (days_since_compilation < 0) {
      LogKit::LogFormatted(LogKit::ERROR,"The CRAVA licence is not valid yet. Are you tampering with the system clock?\n");
      LogKit::LogFormatted(LogKit::ERROR,"Please contact the Norwegian Computing Center to obtain a new one.\n");
      std::exit(EXIT_FAILURE);
    }
    else {
      int days_left_of_licence = licence_days - days_since_compilation;
      LogKit::LogFormatted(LogKit::ERROR,"License expires in                         : %d days\n",days_left_of_licence);
      LogKit::LogFormatted(LogKit::ERROR,"Licensed to                                : %s\n\n",licensed_to.c_str());
    }
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"Compiled: %s/%s\n\n",SystemCall::getDate().c_str(),SystemCall::getTime().c_str());
  }
}

time_t 
Program::TimeOfCompilation(void) const
{
  /*
    This function uses the gcc macros __TIME__ and __DATE__
    and is not necessarily portable to less developed
    compilers.
  */
  std::string ctm = __TIME__;
  std::string cdt = __DATE__;

  std::string mon = cdt.substr(0,3);
  std::string d   = cdt.substr(4,2);
  std::string y   = cdt.substr(7,4);

  struct tm ct;

  if (mon == "Jan")
    ct.tm_mon = 0;
  else if (mon == "Feb")
    ct.tm_mon = 1;
  else if (mon == "Mar")
    ct.tm_mon = 2;
  else if (mon == "Apr")
    ct.tm_mon = 3;
  else if (mon == "May")
    ct.tm_mon = 4;
  else if (mon == "Jun")
    ct.tm_mon = 5;
  else if (mon == "Jul")
    ct.tm_mon = 6;
  else if (mon == "Aug")
    ct.tm_mon = 7;
  else if (mon == "Sep")
    ct.tm_mon = 8;
  else if (mon == "Oct")
    ct.tm_mon = 9;
  else if (mon == "Nov")
    ct.tm_mon = 10;
  else if (mon == "Dec")
    ct.tm_mon = 11;

  ct.tm_mday  = atoi(d.c_str());
  ct.tm_year  = atoi(y.c_str()) - 1900;

  ct.tm_hour  = atoi(ctm.substr(0,2).c_str());
  ct.tm_min   = atoi(ctm.substr(3,2).c_str());
  ct.tm_sec   = atoi(ctm.substr(6,2).c_str());
  ct.tm_wday  = 0 ; // initialized but not used
  ct.tm_yday  = 0;  // initialized but not used
  ct.tm_isdst = 0;  // initialized but not used

  time_t num_secs = mktime(&ct);  // no umr!
  /*
    Transform to number of seconds since 1970.
    There is no umr error here! This is a known purify bug.
  */

  return num_secs;
}
