/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include<ostream>
#include<vector>
#include<string>

#include "src/tasklist.h"
#include "src/io.h"
#include "nrlib/iotools/logkit.hpp"

std::vector<std::string> TaskList::task_(0);

void TaskList::viewAllTasks(bool useFile)
{
  size_t i;
  size_t size = task_.size();

  NRLib::LogKit::WriteHeader("Suggested tasks");

  if (size > 0) {
    if(useFile) {
      std::string fName = IO::makeFullFileName("",IO::FileTasks()+IO::SuffixTextFiles());
      LogKit::SetFileLog(fName,NRLib::LogKit::Low, 1, true);
    }
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,1,"\n");
    for (i=0; i<size; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,1,"%d: %s \n", i+1, task_[i].c_str());
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,1,"\nNo tasks have been suggested.\n");
  }
}
