#include "src/timeline.h"

TimeLine::TimeLine() 
{
  current_time_  = time_.begin();
  current_event_ = event_type_.begin();
  current_index_ = event_data_index_.begin();
}

TimeLine::~TimeLine() {}

void
TimeLine::AddEvent(int time, int event_type, int event_data_index)
{
  std::list<int>::iterator time_it  = time_.begin();
  std::list<int>::iterator event_it = event_type_.begin();
  std::list<int>::iterator index_it = event_data_index_.begin();

  while(time_it != time_.end() && *time_it <= time) {
    ++time_it;
    ++event_it;
    ++index_it;
  }

  time_.insert(time_it, time);
  event_type_.insert(event_it, event_type);
  event_data_index_.insert(index_it, event_data_index);

  current_time_  = time_.begin();
  current_event_ = event_type_.begin();
  current_index_ = event_data_index_.begin();
}

bool 
TimeLine::GetNextEvent(int & event_type, 
                       int & event_data_index, 
                       int & delta_time)
{
  if(current_time_ != time_.end()) {
    event_type       = *current_event_;
    event_data_index = *current_index_;
    if(current_time_ == time_.begin())
      delta_time = 0;
    else {
      std::list<int>::const_iterator prev_time = current_time_;
      --prev_time;
      delta_time = *current_time_ - *prev_time;
    }
    ++current_time_;
    ++current_event_;
    ++current_index_;
    return(true);
  }
  else
    return(false);
}
