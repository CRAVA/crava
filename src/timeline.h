#ifndef TIMELINE_H
#define TIMELINE_H

#include<list>

class TimeLine {
public:
  TimeLine();

  ~TimeLine();

  void AddEvent(int time, int event_type, int event_data_index);                 //Resets GetNext.
  bool GetNextEvent(int & event_type, int & event_data_index, int & delta_time); //Returns false if no  more events.
                                                                                 //Not const, since it advances iterators.
  void ReSet();
  void GetAllTimes(std::list<int> & time) const {time = time_;}

  enum event_types{AVO, TRAVEL_TIME, GRAVITY};

private:

  std::list<int> time_;             //Measured in days. Parallel indexing of the three tables.
  std::list<int> event_type_;
  std::list<int> event_data_index_; //Contains index number for data for this event type.

  std::list<int>::const_iterator current_time_;
  std::list<int>::const_iterator current_event_;
  std::list<int>::const_iterator current_index_;
};

#endif
