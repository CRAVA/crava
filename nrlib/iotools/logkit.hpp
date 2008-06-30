// $Id:  $

#ifndef NRLIB_LOGKIT_H
#define NRLIB_LOGKIT_H

namespace NRLib2 { 

/// Kit for logging of messages from program.
class LogKit {
public:
  ///Philosophy:
  ///A message has a flag, determining the type of message, and a phase,
  ///determining the stage in the program. Each stream has a flag (which may
  ///be a combination of basic flags) for each phase. If phase and flag 
  ///matches, the message is sent to the stream. A message without phase is
  ///sent to all streams that would send it in at least one phase.
  ///The system can be used without bothering with phases. All NRLib logging
  ///is phase 0 BASIC.
  
  enum MessageFlags {ERROR = 1, WARNING = 2, BASIC = 4, DETAILED = 8, 
                     PEDANTIC = 16, DEBUG = 32, DETAILEDDEBUG = 64};

  ///Set a stream that logs independent of phase.
  void setStream(ostream * logstream, int flag, 
                 bool includeNRLibLogging = true);

  ///Set a full phase dependent stream
  void setStream(ostream * logstream, const std::vector<int> & flags);

  ///Set single-phase stream
  void setStream(ostream * logstream, int flag, int phase);


  ///Send message independent of phase
  void logMessage(const string & message, int flag);

  ///Send message in given phase
  void logMessage(const string & message, int flag, int phase);


  ///Close streams
  void endLog();


private:
  static std::vector<LogStream *> logstreams_;
};

class LogStream {
public:
  LogStream(ostream * logstream, int flag);
  LogStream(ostream * logstream, const std::vector<int> & flags);

  void logMessage(std::string & message, int flag);
  void logMessage(std::string & message, int flag, int phase);

private:
  ostream * logstream_;
  std::vector<int> flags_;
  int maxFlag_;
}

#endif

