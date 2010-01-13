#include<ostream>
#include<vector>
#include<string>

#include "src/tasklist.h"
#include "nrlib/iotools/logkit.hpp"
#include "lib/utils.h"

std::vector<std::string> TaskList::task_(0);

void TaskList::viewAllTasks(void)
{
  int i;
  int size = task_.size();

  if (size > 0)
  {
    Utils::writeHeader("Suggested tasks");
    LogKit::LogFormatted(LogKit::LOW,"\n");
    for (i=0; i<size; i++)
      LogKit::LogFormatted(LogKit::LOW,"%d: %s \n", i+1, task_[i].c_str());
  }
}
