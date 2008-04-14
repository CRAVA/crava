// $Id: exception.hpp 30 2008-04-04 12:39:42Z perroe $

#ifndef NRLIB_EXCEPTION_HPP
#define NRLIB_EXCEPTION_HPP

#include <string>

namespace NRLib2 {

class Exception 
{ 
public:
  explicit Exception(const std::string& msg = "") : msg_(msg) { }
  virtual ~Exception() {}
  virtual std::string what() {return msg_;}
private:
  std::string msg_;
};

class IndexOutOfRange : public Exception 
{
public:
  explicit IndexOutOfRange(const std::string& msg = "") 
    : Exception(msg) {}
};

class IOError : public Exception
{
public:
  explicit IOError(const std::string& msg = "") 
    : Exception(msg) {}
};

class FileFormatError : public IOError
{
public:
  explicit FileFormatError(const std::string& msg = "") 
    : IOError(msg) {}
};

class EndOfFile : public IOError
{
public:
  explicit EndOfFile(const std::string& msg = "")
    : IOError(msg) {}
};

}

#endif
