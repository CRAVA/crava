#include <time.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "lib/log.h"
#include "lib/global_def.h"

void
LogKit::initialize(bool screen, bool noFile)
{
  screendump_ = screen;
  noLogFile_  = noFile;
  if(message_ == NULL)
    message_  = new char[MAX_STRING];
  if(buffer_ == NULL && noFile == false) {
    buffer_ = new char[1];
    buffer_[0] = 0;
  }
  debugLevel_ = 0;
}

void
LogKit::terminate()
{
  if(file_ != NULL)
  {
    fclose(file_);
    file_ = NULL;
  }
  if(debugFile_ != NULL)
  {
    fclose(debugFile_);
    debugFile_ = NULL;
  }
  screendump_ = false;

  delete [] message_;
  message_ = NULL;

  delete [] filePrefix_;
  filePrefix_ = NULL;

  debugLevel_ = 0;
}

void
LogKit::setFilePrefix(char * filePrefix)               
{
  sprintf(filePrefix_,filePrefix);
}

char *
LogKit::makeFullFileName(const char * name, const char * postfix)
{
  char * result;
  if(postfix != NULL)
  {
    result = new char[strlen(name)+strlen(filePrefix_)+strlen(postfix)+1];
    sprintf(result,"%s%s%s",filePrefix_, name, postfix);
  }
  else
  {
    result = new char[strlen(name)+strlen(filePrefix_)+1];
    sprintf(result,"%s%s",filePrefix_, name);
  }
  return(result);
}

void
LogKit::setLogFileName(char * fName)
{
  if(file_ == NULL) {
    file_ = fopen(fName, "w");
    if(file_ == NULL) {
      noLogFile_ = true;
      printf("Failed to open file %s. No log file created.\n", fName);
    }
    else {
      noLogFile_ = false;
      if(buffer_ != NULL)
        fprintf(file_,"%s",buffer_);
      delete [] buffer_;
      buffer_ = NULL;
    }
  }
  else
    writeLog("Log file already given, rejecting new file name %s.\n\n",fName);
}

void
LogKit::setDebug(int level, char * fName)
{
  debugLevel_ = level;
  if(debugFile_ != NULL)
  {
    fclose(debugFile_);
    debugFile_ = NULL;
  }
  if(fName != NULL)
  {
    debugFile_ = fopen(fName, "w");
    if(debugFile_ == NULL)
      printf("Failed to open file %s. No debug log file created.\n", fName);
  }
}

void
LogKit::writeLog(const char * format, ...)
{
  if(message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(file_ != NULL)
    {
      fprintf(file_, message_);
      fflush(file_);
    }
    else if(noLogFile_ == false) { //Buffer messages until logfilename is ok.
      char * tmp = new char[strlen(message_)+strlen(buffer_)+1];
      sprintf(tmp,"%s%s",buffer_,message_);
      delete [] buffer_;
      buffer_ = tmp;
    }
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::writeLog(char * format, ...)
{
  if(message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(file_ != NULL)
    {
      fprintf(file_, message_);
      fflush(file_);
    }
    else if(noLogFile_ == false) { //Buffer messages until logfilename is ok.
      char * tmp = new char[strlen(message_)+strlen(buffer_)+1];
      sprintf(tmp,"%s%s",buffer_,message_);
      delete [] buffer_;
      buffer_ = tmp;
    }
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::writeDebugLog(char * format, ...)
{
  if(debugLevel_ > 0 && message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::writeDebugLog(const char * format, ...)
{
  if(debugLevel_ > 0 && message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::writeDebugLog(int level, char * format, ...)
{
  if(debugLevel_ >= level && message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::writeDebugLog(int level, const char * format, ...)
{
  if(debugLevel_ >= level && message_ != NULL)
  {
    va_list ap;

    va_start(ap, format);
    vsprintf(message_, format, ap);
    if(debugFile_ != NULL)
    {
      fprintf(debugFile_, message_);
      fflush(debugFile_);
    }
    if(screendump_ == true)
      printf(message_);
    va_end(ap);
  }
}

void
LogKit::getTime(double& wall, double& cpu)
{
  time_t  tmpwall = time(NULL);
  clock_t tmpcpu  = clock();

  if (tmpwall == (time_t)(-1))  fprintf (stderr, "Unable to get time()\n");
  if (tmpcpu  == (clock_t)(-1)) fprintf (stderr, "Unable to get clock()\n");

  if (wall==0 && cpu==0) 
  {
    wall = (double) tmpwall; 
    cpu  = (double) tmpcpu; 
  }
  else 
  {
    wall = (double) tmpwall - wall;
    cpu  = ((double) tmpcpu  - cpu)/CLOCKS_PER_SEC;
  }
}


void
LogKit::markTime()
{
  timeMark_  = clock();
}


double
LogKit::getPassedTime()
{
  return(double(clock() - timeMark_)/CLOCKS_PER_SEC);
}


bool    LogKit::screendump_ = true;
bool    LogKit::noLogFile_  = false;
FILE *  LogKit::file_       = NULL;
FILE *  LogKit::debugFile_  = NULL;
char *  LogKit::filePrefix_ = new char[MAX_STRING];
char *  LogKit::message_    = new char[MAX_STRING];
int     LogKit::debugLevel_ = 0;
char *  LogKit::buffer_     = NULL;
clock_t LogKit::timeMark_   = 0;


