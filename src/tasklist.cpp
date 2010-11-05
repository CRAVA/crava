#include<ostream>
#include<vector>
#include<string>

#include "src/tasklist.h"
#include "nrlib/iotools/logkit.hpp"

std::vector<std::string> TaskList::task_(0);

void TaskList::viewAllTasks(void)
{
  size_t i;
  size_t size = task_.size();

  if (size > 0)
  {
    NRLib::LogKit::WriteHeader("Suggested tasks");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n");
    for (i=0; i<size; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"%d: %s \n", i+1, task_[i].c_str());
  }
}
