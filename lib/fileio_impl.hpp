// $Id: fileio_impl.hpp 30 2008-04-04 12:39:42Z perroe $

#ifndef NRLIB_FILEIO_IMPL_HPP
#define NRLIB_FILEIO_IMPL_HPP

#include <sstream>

#include "exception.hpp"

#include "fileio_internal.hpp"

/// @file Do not include directly, include fileio.hpp instead.

/// @todo Use correct exceptions.

template <typename iterator>
void NRLib2::WriteBinaryFloatArray(std::ostream& stream, 
                                   iterator begin,
                                   iterator end,
                                   NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename iterator::difference_type n_char = 4*std::distance(begin, end);
  char* buffer = new char[n_char];
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      writeIEEEFloatBE(&buffer[4*i], *begin);
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      writeIEEEFloatLE(&buffer[4*i], *begin);
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(buffer, static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}

template <typename iterator>
iterator NRLib2::ReadBinaryFloatArray(std::istream& stream, 
                                      iterator begin,
                                      int n,
                                      NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  char* buffer = new char[4*n];
  float f;

  if (!stream.read(buffer, 4*n)) {
    throw Exception("Error reading from stream.");
  }  
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      parseIEEEFloatBE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      parseIEEEFloatLE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  return begin;
}

template <typename iterator>
void NRLib2::WriteBinaryDoubleArray(std::ostream& stream, 
                                    iterator begin,
                                    iterator end,
                                    NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename iterator::difference_type n_char = 8*std::distance(begin, end);
  char* buffer = new char[n_char];
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      writeIEEEDoubleBE(&buffer[8*i], *begin);
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      writeIEEEDoubleLE(&buffer[8*i], *begin);
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(buffer, static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}

template <typename iterator>
iterator NRLib2::ReadBinaryDoubleArray(std::istream& stream, 
                                       iterator begin,
                                       int n,
                                       NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  char* buffer = new char[8*n];
  double d;

  if (!stream.read(buffer, 8*n)) {
    throw Exception("Error reading from stream.");
  }  
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      parseIEEEDoubleBE(&buffer[8*i], d);
      *begin = d;
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      parseIEEEDoubleLE(&buffer[8*i], d);
      *begin = d;
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  return begin;
}

#endif // NRLIB_FILEIO_IMPL_HPP
