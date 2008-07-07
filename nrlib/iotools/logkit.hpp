// $Id:  $

#ifndef NRLIB_LOGKIT_H
#define NRLIB_LOGKIT_H

#include<ostream>
#include<vector>
#include<string>

namespace NRLib2 { 

class LogStream;
struct BufferMessage;

/// Kit for logging of messages from program.
class LogKit {
public:
  ///Philosophy:
  ///A message has a level, determining the type of message, and a phase,
  ///determining the stage in the program. Each stream has a level (which may
  ///be a combination of basic levels) for each phase. If phase and flag 
  ///matches, the message is sent to the stream. A message without phase is
  ///sent to all streams that would send it in at least one phase.
  ///
  ///The system can be used without bothering with phases. All NRLib logging
  ///is phase 0 LOW.
  ///
  
  ///Symbols for use when sending message level and parsing exact levels.
  enum MessageLevels {ERROR = 1, WARNING = 2, LOW = 4, MEDIUM = 8, HIGH = 16, DEBUGLOW = 32, DEBUGHIGH = 64};

  ///Symbols for use when parsing given level and lower.
  enum LimitLevels {L_ERROR = 1, L_WARNING = 3, L_LOW = 7, L_MEDIUM = 15,
                    L_HIGH = 31, L_DEBUGLOW = 63, L_DEBUGHIGH = 127};

  ///Set a file that logs independent of phase.
  static void SetFileLog(const std::string & fileName, int levels, 
                         bool includeNRLibLogging = true);

  ///Set a full phase dependent file log
  static void SetFileLog(const std::string & fileName, 
                         const std::vector<int> & levels);

  ///Set single-phase file log, useful for debugging given phase.
  static void SetFileLog(const std::string & fileName, int levels, int phase);


  ///Set a screen log independent of phase.
  static void SetScreenLog(int levels, bool includeNRLibLogging = true);

  ///Set a full phase dependent screen log
  static void SetScreenLog(const std::vector<int> & levels);


  ///Send message independent of phase
  static void LogMessage(int level, const std::string & message);

  ///Send message in given phase
  static void LogMessage(int level, int phase, const std::string & message);

  ///Send message as c-style format string and arguments.
  static void LogFormatted(int level, std::string format, ...);

  ///Send message as c-style format string and arguments.
  static void LogFormatted(int level, int phase, std::string format, ...);


  ///Close streams
  static void EndLog();

  ///Buffering allows temporary storage of messages for sending to files
  ///opened later. When a file log is opened, the buffer is dumped to it.
  ///EndBuffering should be called once all files are opened.
  static void StartBuffering();
  static void EndBuffering();


private:
  static std::vector<LogStream *> logstreams_;
  static int screenLog_; //Remembers which log is screen.
  static std::vector<BufferMessage *> * buffer_;

  static void SendToBuffer(int level, int phase, const std::string & message);
  static void DumpBuffer(LogStream * logstream);
};

///Class LogStream is for internal use in LogKit only.
class LogStream {
public:
  ///Convention: logstream = NULL means cout.
  LogStream(std::ostream * logstream, int level);
  LogStream(std::ostream * logstream, const std::vector<int> & levels);
  ~LogStream();

  void LogMessage(int level, const std::string & message);
  void LogMessage(int level, int phase, const std::string & message);

private:
  std::ostream * logstream_;
  std::vector<int> levels_;
  int fullLevel_;
  bool deleteStream;
};

struct BufferMessage {
  std::string text_;
  int         phase_;
  int         level_;
};

}
#endif

