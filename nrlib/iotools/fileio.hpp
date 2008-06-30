// $Id: fileio.hpp 86 2008-06-18 11:41:15Z perroe $

#ifndef NRLIB_FILEIO_HPP
#define NRLIB_FILEIO_HPP

#include <string>
#include <vector>
#include <iostream>

#include "../exception/exception.hpp"

namespace NRLib2 {
  enum Endianess {END_LITTLE_ENDIAN, END_BIG_ENDIAN};

  /// \brief Check if end of file is reached.
  /// \return True if end of file is reached, else false.
  bool CheckEndOfFile(std::istream& stream);

  /// \brief Read next white-space seperated token from file.
  ///        Updates line number when new line is read in.
  /// \note  Since the non-used part of the line is stored in 
  ///        a static variable it is not thread proof.
  /// \todo  Replace with a class with the same functionality.
  void GetNextToken(std::istream& stream, std::string& s, int& line);

  // ---------------------------------
  // 2-byte integer
  // ---------------------------------

  /// \brief Write a 2-byte integer to a binary file.
  void WriteBinaryShort(std::ostream& stream, 
                        short s,
                        Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read a 2-byte integer from a binary file.
  short ReadBinaryShort(std::istream& stream, 
                        Endianess file_format = END_BIG_ENDIAN);

  /// \brief Write an array of binary 2-byte integers.
  template<typename I>
  void WriteBinaryShortArray(std::ostream& stream,
                             I begin,
                             I end,
                             Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read an array of binary 2-byte integers.
  template<typename I>
  I ReadBinaryShortArray(std::istream& stream,
                         I begin,
                         int n,
                         Endianess file_format = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte integer
  // ---------------------------------

  /// \brief Write a 4-byte integer to a binary file.
  void WriteBinaryInt(std::ostream& stream, 
                      int i,
                      Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read a 4-byte integer from a binary file.
  int ReadBinaryInt(std::istream& stream, 
                    Endianess file_format = END_BIG_ENDIAN);

  /// \brief Write an array of binary 4-byte integers.
  template<typename I>
  void WriteBinaryIntArray(std::ostream& stream,
                           I begin,
                           I end,
                           Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read an array of binary 4-byte integers.
  template<typename I>
  I ReadBinaryIntArray(std::istream& stream,
                       I begin,
                       int n,
                       Endianess file_format = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte IEEE floating point number
  // ---------------------------------

  /// \brief Write a 4-byte float on standard IEEE format.
  void WriteBinaryFloat(std::ostream& stream, 
                        float f,
                        Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read a 4-byte float on standard IEEE format.
  float ReadBinaryFloat(std::istream& stream, 
                        Endianess file_format = END_BIG_ENDIAN);

  /// \brief Write an array of 4-byte floats on standard IEEE format.
  template<typename I>
  void WriteBinaryFloatArray(std::ostream& stream,
                             I begin,
                             I end,
                             Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read an array of 4-byte floats on standard IEEE format.
  template<typename I>
  I ReadBinaryFloatArray(std::istream& stream,
                         I begin,
                         int n,
                         Endianess file_format = END_BIG_ENDIAN);

  // ---------------------------------
  // 8-byte IEEE floating point number
  // ---------------------------------

  /// \brief Write a 8-byte float on standard IEEE format.
  void WriteBinaryDouble(std::ostream& stream, 
                         double d,
                         Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read a 8-byte float on standard IEEE format.
  double ReadBinaryDouble(std::istream& stream, 
                          Endianess file_format = END_BIG_ENDIAN);

  /// \brief Write an array of 8-byte floats on standard IEEE format.
  template<typename I>
  void WriteBinaryDoubleArray(std::ostream& stream,
                              I begin,
                              I end,
                              Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read an array of 8-byte floats on standard IEEE format.
  template<typename I>
  I ReadBinaryDoubleArray(std::istream& stream,
                          I begin,
                          int n,
                          Endianess file_format = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte IBM floating point number
  // ---------------------------------

  /// \brief Write a 4-byte float on standard IEEE format.
  void WriteBinaryIbmFloat(std::ostream& stream, 
                           float f,
                           Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read a 4-byte float on standard IEEE format.
  float ReadBinaryIbmFloat(std::istream& stream, 
                           Endianess file_format = END_BIG_ENDIAN);

  /// \brief Write an array of 4-byte floats on standard IEEE format.
  template<typename I>
  void WriteBinaryIbmFloatArray(std::ostream& stream,
                                I begin,
                                I end,
                                Endianess file_format = END_BIG_ENDIAN);

  /// \brief Read an array of 4-byte floats on standard IEEE format.
  template<typename I>
  I ReadBinaryIbmFloatArray(std::istream& stream,
                            I begin,
                            int n,
                            Endianess file_format = END_BIG_ENDIAN);

} // namespace NRLib2


namespace NRLib2_internal {
  /// \todo Use stdint.h if available.
  // typedef unsigned int uint32_t;
  // typedef unsigned long long uint64_t;
  
  union FloatAsInt {
    /*uint32_t*/ unsigned int ui;
    float f;
  };

  union DoubleAsLongLong {
    /*uint64_t*/ unsigned long long ull;
    double d;
  };

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseUInt16BE(const char* buffer, /*uint16_t*/ unsigned short& ui);

  /// Parse unsigned 32-bit integer from little-endian buffer.
  inline void ParseUInt16LE(const char* buffer, /*uint16_t*/ unsigned short& ui);

  /// Write unsigned 32-bit integer to big-endian buffer.
  inline void WriteUInt16BE(char* buffer, /*uint16_t*/ unsigned short us);

  /// Write unsigned 32-bit integer to little-endian buffer.
  inline void WriteUInt16LE(char* buffer, /*uint16_t*/ unsigned short us);

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseUInt32BE(const char* buffer, /*uint32_t*/ unsigned int& ui);

  /// Parse unsigned 32-bit integer from little-endian buffer.
  inline void ParseUInt32LE(const char* buffer, /*uint32_t*/ unsigned int& ui);

  /// Write unsigned 32-bit integer to big-endian buffer.
  inline void WriteUInt32BE(char* buffer, /*uint32_t*/ unsigned int ui);

  /// Write unsigned 32-bit integer to little-endian buffer.
  inline void WriteUInt32LE(char* buffer, /*uint32_t*/ unsigned int ui);

  /// Parse IEEE single-precision float from big-endian buffer. 
  inline void ParseIEEEFloatBE(const char* buffer, float& f);

  /// Parse IEEE single-precision float from little-endian buffer. 
  inline void ParseIEEEFloatLE(const char* buffer, float& f);

  /// Write IEEE single-precision float to big-endian buffer. 
  inline void WriteIEEEFloatBE(char* buffer, float f);

  /// Write IEEE single-precision float to little-endian buffer. 
  inline void WriteIEEEFloatLE(char* buffer, float f);

  /// Parse IEEE double-precision float from big-endian buffer. 
  inline void ParseIEEEDoubleBE(const char* buffer, double& d);

  /// Parse IEEE double-precision float from little-endian buffer. 
  inline void ParseIEEEDoubleLE(const char* buffer, double& d);

  /// Write IEEE double-precision float to big-endian buffer. 
  inline void WriteIEEEDoubleBE(char* buffer, double d);

  /// Write IEEE double-precision float to little-endian buffer. 
  inline void WriteIEEEDoubleLE(char* buffer, double d);

  /// Parse IEEE double-precision float from big-endian buffer. 
  inline void ParseIBMFloatBE(const char* buffer, float& f);

  /// Parse IEEE double-precision float from little-endian buffer. 
  inline void ParseIBMFloatLE(const char* buffer, float& f);

  /// Write IEEE double-precision float to big-endian buffer. 
  inline void WriteIBMFloatBE(char* buffer, float f);

  /// Write IEEE double-precision float to little-endian buffer. 
  inline void WriteIBMFloatLE(char* buffer, float f);
} // namespace NRLib2_internal


template <typename I>
void NRLib2::WriteBinaryShortArray(std::ostream& stream, 
                                   I begin,
                                   I end,
                                   NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename I::difference_type n_char = 2*std::distance(begin, end);
  std::vector<char> buffer(n_char);
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteUInt16BE(&buffer[2*i], static_cast<unsigned short>(*begin));
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteUInt16LE(&buffer[2*i], static_cast<unsigned short>(*begin));
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(&buffer[0], static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename I>
I NRLib2::ReadBinaryShortArray(std::istream& stream, 
                               I begin,
                               int n,
                               NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  std::vector<char> buffer(2*n);
  unsigned short us;

  if (!stream.read(&buffer[0], 2*n)) {
    throw Exception("Error reading from stream.");
  }   
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseUInt16BE(&buffer[2*i], us);
      *begin = static_cast<typename std::iterator_traits<I>::value_type>(us);
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseUInt16LE(&buffer[2*i], us);
      *begin = static_cast<typename std::iterator_traits<I>::value_type>(us);
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return begin;
}


template <typename I>
void NRLib2::WriteBinaryIntArray(std::ostream& stream, 
                                 I begin,
                                 I end,
                                 NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename I::difference_type n_char = 4*std::distance(begin, end);
  std::vector<char> buffer(n_char);
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteUInt32BE(&buffer[4*i], static_cast<unsigned int>(*begin));
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteUInt32LE(&buffer[4*i], static_cast<unsigned int>(*begin));
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(&buffer[0], static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename I>
I NRLib2::ReadBinaryIntArray(std::istream& stream, 
                             I begin,
                             int n,
                             NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  std::vector<char> buffer(4*n);
  unsigned int ui;

  if (!stream.read(&buffer[0], 4*n)) {
    throw Exception("Error reading from stream.");
  }   
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseUInt32BE(&buffer[4*i], ui);
      *begin = static_cast<typename std::iterator_traits<I>::value_type>(ui);
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseUInt32LE(&buffer[4*i], ui);
      *begin = static_cast<typename std::iterator_traits<I>::value_type>(ui);
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return begin;
}


template <typename I>
void NRLib2::WriteBinaryFloatArray(std::ostream& stream, 
                                   I begin,
                                   I end,
                                   NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename I::difference_type n_char = 4*std::distance(begin, end);
  std::vector<char> buffer(n_char);
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIEEEFloatBE(&buffer[4*i], *begin);
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIEEEFloatLE(&buffer[4*i], *begin);
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(&buffer[0], static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename I>
I NRLib2::ReadBinaryFloatArray(std::istream& stream, 
                               I begin,
                               int n,
                               NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  std::vector<char> buffer(4*n);
  float f;

  if (!stream.read(&buffer[0], 4*n)) {
    throw Exception("Error reading from stream.");
  }   
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIEEEFloatBE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIEEEFloatLE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return begin;
}


template <typename I>
void NRLib2::WriteBinaryDoubleArray(std::ostream& stream, 
                                    I begin,
                                    I end,
                                    NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename I::difference_type n_char = 8*std::distance(begin, end);
  std::vector<char> buffer(n_char);
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIEEEDoubleBE(&buffer[8*i], *begin);
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIEEEDoubleLE(&buffer[8*i], *begin);
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(&buffer[0], static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename I>
I NRLib2::ReadBinaryDoubleArray(std::istream& stream, 
                                I begin,
                                int n,
                                NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  std::vector<char> buffer(8*n);
  double d;

  if (!stream.read(&buffer[0], 8*n)) {
    throw Exception("Error reading from stream.");
  }  
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIEEEDoubleBE(&buffer[8*i], d);
      *begin = d;
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIEEEDoubleLE(&buffer[8*i], d);
      *begin = d;
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  return begin;
}


template <typename I>
void NRLib2::WriteBinaryIbmFloatArray(std::ostream& stream, 
                                      I begin,
                                      I end,
                                      NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  typename I::difference_type n_char = 4*std::distance(begin, end);
  std::vector<char> buffer(n_char);
 
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIBMFloatBE(&buffer[4*i], *begin);
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; begin != end; ++begin, ++i) {
      WriteIBMFloatLE(&buffer[4*i], *begin);
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }
  
  if (!stream.write(&buffer[0], static_cast<std::streamsize>(n_char))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename I>
I NRLib2::ReadBinaryIbmFloatArray(std::istream& stream, 
                                  I begin,
                                  int n,
                                  NRLib2::Endianess file_format)
{
  using namespace NRLib2_internal;

  std::vector<char> buffer(4*n);
  float f;

  if (!stream.read(&buffer[0], 4*n)) {
    throw Exception("Error reading from stream.");
  }   
  
  switch (file_format) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIBMFloatBE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n; ++i) {
      ParseIBMFloatLE(&buffer[4*i], f);
      *begin = f;
      ++begin;
    }
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return begin;
}


static const unsigned short smasks[2] = {0x00ff, 0xff00};

static const unsigned int masks[4] = 
{0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000};

static const unsigned long long lmasks[8] = 
{0x00000000000000ffULL, 0x000000000000ff00ULL, 
 0x0000000000ff0000ULL, 0x00000000ff000000ULL,
 0x000000ff00000000ULL, 0x0000ff0000000000ULL,
 0x00ff000000000000ULL, 0xff00000000000000ULL};


// -----------------  16 bit integer ---------------------

void NRLib2_internal::ParseUInt16BE(const char* buffer,
                                    /*uint16_t*/ unsigned short& us)
{
  us = static_cast<unsigned char>(buffer[0]);
  us <<= 8;
  us |= static_cast<unsigned char>(buffer[1]);
}


void NRLib2_internal::ParseUInt16LE(const char* buffer,
                                    /*uint32_t*/ unsigned short& us)
{
  us = static_cast<unsigned char>(buffer[1]);
  us <<= 8;
  us |= static_cast<unsigned char>(buffer[0]);
}


void NRLib2_internal::WriteUInt16BE(char* buffer, 
                                    /*uint32_t*/ unsigned short us)
{
  buffer[0] = static_cast<char>( (us & smasks[1]) >> 8 );
  buffer[1] = static_cast<char>( us & smasks[0] );
}


void NRLib2_internal::WriteUInt16LE(char* buffer, 
                                    /*uint32_t*/ unsigned short us)
{
  buffer[0] = static_cast<char>(us & smasks[0]);
  buffer[1] = static_cast<char>( (us & smasks[1]) >> 8 );
}


// -----------------  32 bit integer ---------------------


void NRLib2_internal::ParseUInt32BE(const char* buffer,
                                    /*uint32_t*/ unsigned int& ui)
{
  ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    ui <<= 8;
    ui |= static_cast<unsigned char>(buffer[i]);
  }
}


void NRLib2_internal::ParseUInt32LE(const char* buffer,
                                    /*uint32_t*/ unsigned int& ui)
{
  ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    ui <<= 8;
    ui |= static_cast<unsigned char>(buffer[3-i]);
  }
}


void NRLib2_internal::WriteUInt32BE(char* buffer, 
                                    /*uint32_t*/ unsigned int ui)
{
  buffer[0] = static_cast<char>( (ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( ui & masks[0] );
}


void NRLib2_internal::WriteUInt32LE(char* buffer, 
                                    /*uint32_t*/ unsigned int ui)
{
  buffer[0] = static_cast<char>(ui & masks[0]);
  buffer[1] = static_cast<char>( (ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (ui & masks[3]) >> 24 );
}


// -----------------  32 bit IEEE floating point -------------------


// Big endian file format.
void NRLib2_internal::ParseIEEEFloatBE(const char* buffer, float& f)
{
  FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[i]);
  }
  f = tmp.f;
}

// Little endian file format.
void NRLib2_internal::ParseIEEEFloatLE(const char* buffer, float& f)
{
  FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[3-i]);
  }
  f = tmp.f;
}


// Big endian file format.
void NRLib2_internal::WriteIEEEFloatBE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  buffer[0] = static_cast<char>( (tmp.ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( tmp.ui & masks[0] );
}

// Little endian file format.
void NRLib2_internal::WriteIEEEFloatLE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  buffer[0] = static_cast<char>(tmp.ui & masks[0]);
  buffer[1] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (tmp.ui & masks[3]) >> 24 );
}


// -----------------  64 bit IEEE floating point --------------------


// Big endian file format.
void NRLib2_internal::ParseIEEEDoubleBE(const char* buffer, double& d)
{
  DoubleAsLongLong tmp;
  tmp.ull = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 8; ++i) {
    tmp.ull <<= 8;
    tmp.ull |= static_cast<unsigned char>(buffer[i]);
  }
  d = tmp.d;
}

// Little endian file format.
void NRLib2_internal::ParseIEEEDoubleLE(const char* buffer, double& d)
{
  DoubleAsLongLong tmp;
  tmp.ull = static_cast<unsigned char>(buffer[7]);
  for (int i = 1; i < 8; ++i) {
    tmp.ull <<= 8;
    tmp.ull |= static_cast<unsigned char>(buffer[7-i]);
  }
  d = tmp.d;
}


// Big endian file format.
void NRLib2_internal::WriteIEEEDoubleBE(char* buffer, double d)
{
  DoubleAsLongLong tmp;  
  tmp.d = d;
  buffer[0] = static_cast<char>( (tmp.ull & lmasks[7]) >> 56);
  buffer[1] = static_cast<char>( (tmp.ull & lmasks[6]) >> 48);
  buffer[2] = static_cast<char>( (tmp.ull & lmasks[5]) >> 40);
  buffer[3] = static_cast<char>( (tmp.ull & lmasks[4]) >> 32);
  buffer[4] = static_cast<char>( (tmp.ull & lmasks[3]) >> 24);
  buffer[5] = static_cast<char>( (tmp.ull & lmasks[2]) >> 16);
  buffer[6] = static_cast<char>( (tmp.ull & lmasks[1]) >> 8 );
  buffer[7] = static_cast<char>(  tmp.ull & lmasks[0] );
}

// Little endian file format.
void NRLib2_internal::WriteIEEEDoubleLE(char* buffer, double d)
{
  DoubleAsLongLong tmp;  
  tmp.d = d;
  buffer[0] = static_cast<char>(  tmp.ull & lmasks[0]);
  buffer[1] = static_cast<char>( (tmp.ull & lmasks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ull & lmasks[2]) >> 16);
  buffer[3] = static_cast<char>( (tmp.ull & lmasks[3]) >> 24);
  buffer[4] = static_cast<char>( (tmp.ull & lmasks[4]) >> 32);
  buffer[5] = static_cast<char>( (tmp.ull & lmasks[5]) >> 40);
  buffer[6] = static_cast<char>( (tmp.ull & lmasks[6]) >> 48);
  buffer[7] = static_cast<char>( (tmp.ull & lmasks[7]) >> 56);
}


// -----------------  32 bit IBM floating point -------------------

static unsigned int IEEEMAX  = 0x7FFFFFFF;
static unsigned int IEMAXIB  = 0x611FFFFF;
static unsigned int IEMINIB  = 0x21200000;


/// Converts Ieee float (represented as 32-bit int) to IBM float.
static inline void Ibm2Ieee(/*uint32_t*/ unsigned int& in)
{
	static int it[8] = { 0x21800000, 0x21400000, 0x21000000, 0x21000000,
			    0x20c00000, 0x20c00000, 0x20c00000, 0x20c00000 };
	static int mt[8] = { 8, 4, 2, 2, 1, 1, 1, 1 };
	/*uint32_t*/ unsigned int manthi, iexp, inabs;
	int ix;

  manthi = in & 0x00ffffff;
  ix     = manthi >> 21;
	iexp   = ( ( in & 0x7f000000 ) - it[ix] ) << 1;
	manthi = manthi * mt[ix] + iexp;
	inabs  = in & 0x7fffffff;
	if ( inabs > IEMAXIB ) manthi = IEEEMAX;
	manthi = manthi | ( in & 0x80000000 );
	in = ( inabs < IEMINIB ) ? 0 : manthi;
}


/// Converts Ieee float (represented as 32-bit int) to IBM float.
static inline void Ieee2Ibm ( /*uint32_t*/ unsigned int& in)
{
	static int it[4] = { 0x21200000, 0x21400000, 0x21800000, 0x22100000 };
	static int mt[4] = { 2, 4, 8, 1 };
	/*uint32_t*/ unsigned int manthi, iexp, ix;

  ix     = ( in & 0x01800000 ) >> 23;
	iexp   = ( ( in & 0x7e000000 ) >> 1 ) + it[ix];
	manthi = ( mt[ix] * ( in & 0x007fffff) ) >> 3;
	manthi = (manthi + iexp) | ( in & 0x80000000 );
	in     = ( in & 0x7fffffff ) ? manthi : 0;
}


// Big endian file format.
void NRLib2_internal::ParseIBMFloatBE(const char* buffer, float& f)
{
  FloatAsInt tmp; 
  tmp.ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[i]);
  }

  Ibm2Ieee(tmp.ui);
  f = tmp.f;
}


// Little endian file format.
void NRLib2_internal::ParseIBMFloatLE(const char* buffer, float& f)
{
  FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[3-i]);
  }

  Ibm2Ieee(tmp.ui);
  f = tmp.f;
}


// Big endian file format.
void NRLib2_internal::WriteIBMFloatBE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  Ieee2Ibm(tmp.ui);
  buffer[0] = static_cast<char>( (tmp.ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( tmp.ui & masks[0] );
}


// Little endian file format.
void NRLib2_internal::WriteIBMFloatLE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  Ieee2Ibm(tmp.ui);
  buffer[0] = static_cast<char>(tmp.ui & masks[0]);
  buffer[1] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (tmp.ui & masks[3]) >> 24 );
}


#endif // NRLIB_FILEIO_HPP
