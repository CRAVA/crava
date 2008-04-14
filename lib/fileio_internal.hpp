// $Id: fileio_internal.hpp 35 2008-04-07 12:01:31Z perroe $

#ifndef NRLIB_FILEIO_INTERNAL_HPP
#define NRLIB_FILEIO_INTERNAL_HPP

#include "fileio.hpp"

#include <cstdio>

// #include <ctypes>

/// @file Internal implementation for file IO.
/// Should not be included directly, include fileio.hpp instead.
/// Not intended to be used outside NRLib2.

static const unsigned int masks[4] = 
{0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000};

static const unsigned long long lmasks[8] = 
{0x00000000000000ff, 0x000000000000ff00, 
 0x0000000000ff0000, 0x00000000ff000000,
 0x000000ff00000000, 0x0000ff0000000000,
 0x00ff000000000000, 0xff00000000000000};

namespace NRLib2_internal {
  union FloatAsInt {
    /*uint32_t*/ unsigned int ui;
    float f;
  };

  union DoubleAsLongLong {
    /* uint64_t */ unsigned long long ull;
    double d;
  };

  /// Parse IEEE single-precision float from big-endian buffer. 
  inline void parseIEEEFloatBE(const char* buffer, float& f);

  /// Parse IEEE single-precision float from little-endian buffer. 
  inline void parseIEEEFloatLE(const char* buffer, float& f);

  /// Write IEEE single-precision float to big-endian buffer. 
  inline void writeIEEEFloatBE(char* buffer, float f);

  /// Write IEEE single-precision float to little-endian buffer. 
  inline void writeIEEEFloatLE(char* buffer, float f);

  /// Parse IEEE double-precision float from big-endian buffer. 
  inline void parseIEEEDoubleBE(const char* buffer, double& d);

  /// Parse IEEE double-precision float from little-endian buffer. 
  inline void parseIEEEDoubleLE(const char* buffer, double& d);

  /// Write IEEE double-precision float to big-endian buffer. 
  inline void writeIEEEDoubleBE(char* buffer, double d);

  /// Write IEEE double-precision float to little-endian buffer. 
  inline void writeIEEEDoubleLE(char* buffer, double d);

} // namespace NRLib2_internal


// Big endian file format.
void NRLib2_internal::parseIEEEFloatBE(const char* buffer, float& f)
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
void NRLib2_internal::parseIEEEFloatLE(const char* buffer, float& f)
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
void NRLib2_internal::writeIEEEFloatBE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  buffer[0] = static_cast<char>( (tmp.ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( tmp.ui & masks[0] );
}

// Little endian file format.
void NRLib2_internal::writeIEEEFloatLE(char* buffer, float f)
{
  FloatAsInt tmp;  
  tmp.f = f;
  buffer[0] = static_cast<char>(tmp.ui & masks[0]);
  buffer[1] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (tmp.ui & masks[3]) >> 24 );
}

// Big endian file format.
void NRLib2_internal::parseIEEEDoubleBE(const char* buffer, double& d)
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
void NRLib2_internal::parseIEEEDoubleLE(const char* buffer, double& d)
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
void NRLib2_internal::writeIEEEDoubleBE(char* buffer, double d)
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
void NRLib2_internal::writeIEEEDoubleLE(char* buffer, double d)
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


#endif // NRLIB_FILEIO_INTERNAL_HPP 
