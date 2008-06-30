// $Id: stringtools.hpp 83 2008-06-16 12:00:18Z perroe $

#ifndef NRLIB_STRINGTOOLS_HPP
#define NRLIB_STRINGTOOLS_HPP

#include <string>
#include <vector>
#include <sstream>

#include "../exception/exception.hpp"

namespace NRLib2 {
  // std::vector<std::string> GetTokens(const std::string& s);

  /// OBS: Negative values are recognized as unsigned integers on Windows.
  /// Works OK on Linux. (gcc 4.*)
  template <typename T>
  bool IsType(const std::string& s);


  /// OBS: Negative values are recognized as unsigned integers on Windows.
  /// Works OK on Linux. (gcc 4.*)
  template <typename T>
  T ParseType(const std::string& s);

  template <typename T>
  std::string ToString(const T obj);
}

/// @todo Use correct exceptions.

template <typename T>
bool NRLib2::IsType(const std::string& s)
{
  std::istringstream i(s);
  T x;
  char c;
  if (!(i >> x) || i.get(c))
    return false;
  return true;
}


template <typename T>
T NRLib2::ParseType(const std::string& s)
{
  std::istringstream i(s);
  T x;
  char c;
  if (!(i >> x))
    throw Exception("Failed to convert \"" + s + "\" to " + typeid(T).name());
  if (i.get(c)) 
    throw Exception("Could not convert whole \"" + s + "\" to " + typeid(T).name());
  return x;
}


template <typename T>
std::string NRLib2::ToString(const T obj)
{
  std::ostringstream o;
  if (!(o << obj)) {
    throw Exception("Bad conversion.");
  }
  return o.str();
}


#endif // NRLIB_STRINGTOOLS_HPP
