// $Id$

#include <fstream>
#include <iostream>
#include <stdarg.h>

#include "logkit.hpp"
#include "../exception/exception.hpp"

using namespace NRLib2;

std::vector<LogStream*> LogKit::logstreams_(0);
int LogKit::screenLog_ = -1;
std::vector<BufferMessage *> * LogKit::buffer_ = NULL;

void
LogKit::SetFileLog(const std::string & fileName, int levels, 
                   bool includeNRLibLogging) 
{
  std::ofstream * file = new std::ofstream(fileName.c_str());
  if (!(*file)) {
    throw new IOError("Error opening " + fileName);
  }
  LogStream * curStream;
  if(includeNRLibLogging == true)
    curStream = new LogStream(file, levels);
  else {
    std::vector<int> phaseLevels;
    phaseLevels.push_back(0); //No logging in NRLib, phase 0.
    phaseLevels.push_back(levels); //This will be used for all other phases.
    curStream = new LogStream(file, phaseLevels);
  }
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}

void
LogKit::SetFileLog(const std::string & fileName, const std::vector<int> & levels) {
  std::ofstream * file = new std::ofstream(fileName.c_str());
  if (!(*file)) {
    throw new IOError("Error opening " + fileName);
  }
  LogStream * curStream = new LogStream(file, levels);
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}

void 
LogKit::SetFileLog(const std::string & fileName, int levels, int phase) {
  std::ofstream * file = new std::ofstream(fileName.c_str());
  if (!(*file)) {
    throw new IOError("Error opening " + fileName);
  }
  std::vector<int> phaseLevels(1000,0);
  phaseLevels[phase] = levels;
  LogStream * curStream = new LogStream(file, phaseLevels);
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}

void
LogKit::SetScreenLog(int levels, bool includeNRLibLogging) 
{
  LogStream * curStream;
  if(includeNRLibLogging == true)
    curStream = new LogStream(NULL, levels);
  else {
    std::vector<int> phaseLevels;
    phaseLevels.push_back(0); //No logging in NRLib, phase 0.
    phaseLevels.push_back(levels); //This will be used for all other phases.
    curStream = new LogStream(NULL, phaseLevels);
  }
  if(screenLog_ < 0) {
    screenLog_ = int(logstreams_.size());
    logstreams_.push_back(curStream);
  }
  else {
    delete logstreams_[screenLog_];
    logstreams_[screenLog_] = curStream;
  }
}

void
LogKit::SetScreenLog(const std::vector<int> & levels) {
  LogStream * curStream = new LogStream(NULL, levels);
  if(screenLog_ < 0) {
    screenLog_ = int(logstreams_.size());
    logstreams_.push_back(curStream);
  }
  else {
    delete logstreams_[screenLog_];
    logstreams_[screenLog_] = curStream;
  }
}


void
LogKit::LogMessage(int level, const std::string & message) {
  unsigned int i;
  for(i=0;i<logstreams_.size();i++)
    logstreams_[i]->LogMessage(level, message);
  SendToBuffer(level,-1,message);
}

void
LogKit::LogMessage(int level, int phase, const std::string & message) {
  unsigned int i;
  for(i=0;i<logstreams_.size();i++)
    logstreams_[i]->LogMessage(level, phase, message);
  SendToBuffer(level,phase,message);
}

void
LogKit::LogFormatted(int level, std::string format, ...) {
  va_list ap;
  char message[1000];
  va_start(ap, format);
  vsprintf(message, format.c_str(), ap);
  va_end(ap);
  LogMessage(level, std::string(message));
}

void
LogKit::LogFormatted(int level, int phase, std::string format, ...) {
  va_list ap;
  char message[1000];
  va_start(ap, format);
  vsprintf(message, format.c_str(), ap);
  va_end(ap);
  LogMessage(level, phase, std::string(message));
}


void
LogKit::EndLog() {
  unsigned int i;
  for(i=0;i<logstreams_.size();i++)
    delete logstreams_[i];

  if(buffer_ != NULL)
    EndBuffering(); //Also deletes buffer.
}

void
LogKit::StartBuffering() {
  buffer_ = new std::vector<BufferMessage *>;
}

void
LogKit::EndBuffering() {
  if(buffer_ != NULL) {
    for(unsigned int i=0;i<buffer_->size();i++) {
      delete (*buffer_)[i];
      (*buffer_)[i] = NULL;
    }
    delete buffer_;
    buffer_ = NULL;
  }
}

void
LogKit::SendToBuffer(int level, int phase, const std::string & message) {
  if(buffer_ != NULL) {
    BufferMessage * bm = new BufferMessage;
    bm->level_ = level;
    bm->phase_ = phase;
    bm->text_  = message;
    buffer_->push_back(bm);
  }
}

void
LogKit::DumpBuffer(LogStream *logstream) {
  if(buffer_ != NULL) {
    for(unsigned int i=0;i<buffer_->size();i++) {
      if((*buffer_)[i]->phase_ < 0)
        logstream->LogMessage((*buffer_)[i]->level_, (*buffer_)[i]->text_);
      else
        logstream->LogMessage((*buffer_)[i]->level_, (*buffer_)[i]->phase_, (*buffer_)[i]->text_);
    }
  }
}



LogStream::LogStream(std::ostream * logstream, int level) {
  fullLevel_ = level;
  if(logstream != NULL) {
    logstream_ = logstream;
    deleteStream = true;
  }
  else {
    logstream_ = &(std::cout);
    deleteStream = false;
  }
}

LogStream::LogStream(std::ostream * logstream, const std::vector<int> & levels) {
  unsigned int i;
  fullLevel_ = 0;
  for(i=0;i<levels.size();i++) {
    fullLevel_ = (fullLevel_ | levels[i]);
    levels_.push_back(levels[i]);
  }
  if(logstream != NULL) {
    logstream_ = logstream;
    deleteStream = true;
  }
  else {
    logstream_ = &(std::cout);
    deleteStream = false;
  }
}

LogStream::~LogStream() {
  if(deleteStream == true) {
    delete logstream_;
  }
}

void
LogStream::LogMessage(int level, const std::string & message) {
  if((level & fullLevel_) > 0) {
    *logstream_ << message;
    logstream_->flush();
  }
}

void
LogStream::LogMessage(int level, int phase, const std::string & message) {
  if(phase < int(levels_.size())) {
    if((level & levels_[phase]) > 0) {
      *logstream_ << message;
      logstream_->flush();
    }
  }
  else if((level & fullLevel_) > 0) {
    *logstream_ << message;
    logstream_->flush();
  }
}


