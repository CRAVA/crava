#ifndef LOGKIT_H
#define LOGKIT_H

#include <stdio.h>
#include <time.h>

class LogKit 
{
public:
  static void     initialize(bool screen, bool noFile = false);
  static void     setDebug(int level, char * fName);
  static void     setLogFileName(char * fName);
  static void     terminate();
  static void     writeLog(char * format, ...);
  static void     writeDebugLog(char * format, ...);
  static void     writeDebugLog(int level, char * format, ...);
  static void     writeLog(const char * format, ...);                 // NBNB-PAL: I have added 'const char *'
  static void     writeDebugLog(const char * format, ...);            // versions of the logging function to
  static void     writeDebugLog(int level, const char * format, ...); // avoid warnings with new g++ compiler
  static int      getDebugLevel() {return(debugLevel_);}
  static void     getTime(double& wall, double& cpu);
  static void     markTime();
  static double   getPassedTime();
  static void     setFilePrefix(char * filePrefix);
  static char   * makeFullFileName(const char * name, const char * postfix = NULL);


private:
  static char   * filePrefix_;  // Prefix (including path) for all output files
  static bool     screendump_;
  static bool     noLogFile_;
  static FILE   * file_;
  static FILE   * debugFile_;
  static char   * message_;
  static int      debugLevel_;
  static char   * buffer_;      //Holds information before output file is set.
  static clock_t  timeMark_;    //Reference time.
};

#endif
