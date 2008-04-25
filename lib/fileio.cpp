// $Id: fileio.cpp 29 2008-04-04 12:30:38Z perroe $

#include "fileio.hpp"
#include "fileio_internal.hpp"

//#include "../exception/exception.hpp"

using namespace NRLib2_internal;
using namespace NRLib2;

bool NRLib2::CheckEndOfFile(std::istream& stream)
{
  stream.peek();
  if (stream.eof()) {
    return true;
  }
  return false;
}


void NRLib2::WriteBinaryFloat(std::ostream& stream, 
                              float f, 
                              Endianess file_format)
{
  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    writeIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    writeIEEEFloatLE(buffer, f);
    break;
 // default:
 //   throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
 //   throw Exception("Error writing to stream.");
  }
}

float NRLib2::ReadBinaryFloat(std::istream& stream, 
                              Endianess file_format)
{
  float f = -99999.0f;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
 //   throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    parseIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    parseIEEEFloatLE(buffer, f);
    break;
//  default:
//    throw Exception("Invalid file format.");
  }

  return f;
}

void NRLib2::WriteBinaryDouble(std::ostream& stream, 
                               double d, 
                               Endianess file_format)
{
  char buffer[8];

  switch (file_format) {
  case END_BIG_ENDIAN:
    writeIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    writeIEEEDoubleLE(buffer, d);
    break;
//  default:
//    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
//    throw Exception("Error writing to stream.");
  }
}

double NRLib2::ReadBinaryDouble(std::istream& stream, 
                                Endianess file_format)
{
  double d = -99999.0;
  char buffer[8];

  if (!stream.read(buffer, 8)) {
//    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    parseIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    parseIEEEDoubleLE(buffer, d);
    break;
//  default:
//    throw Exception("Invalid file format.");
  }

  return d;
}
