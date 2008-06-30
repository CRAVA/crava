// $Id: exception.hpp 43 2008-04-11 13:47:10Z perroe $

#ifndef NRLIB_EXCEPTION_HPP
#define NRLIB_EXCEPTION_HPP

#include <exception>
#include <string>

namespace NRLib2 {

class Exception : public std::exception
{ 
public:
  explicit Exception(const std::string& msg = "") : msg_(msg) { }
  virtual ~Exception() throw() {}
  virtual const char * what() const throw() {return msg_.c_str();}
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
