//$Id$

#include "fileio.hpp"

#include "../exception/exception.hpp"

#include <sstream>
#include <locale>

using namespace NRLib2_internal;
using namespace NRLib2;

void NRLib2::GetNextToken(std::istream& stream, 
                          std::string&  s, 
                          int&          line_num)
{
  std::locale loc = std::locale();
  char c(' ');
  while(std::isspace(c,loc) == true && stream.eof() == false) {
    stream.get(c);
    if(c == '\n')
      line_num++;
  }
  s = "";
  while(isspace(c) == false  && stream.eof() == false) {
    s = s+c;
    stream.get(c);
  }
  while(std::isspace(c,loc) == true && c != '\n' && stream.eof() == false) {
    stream.get(c);
    if(c == '\n')
      line_num++;
  }
  if(c != '\n' && stream.eof() == false) //Do not insert anything at eof.
    stream.putback(c);
}


bool NRLib2::CheckEndOfFile(std::istream& stream)
{
  std::string token;
  int line = 0;
  GetNextToken(stream,token,line); //Read possible whitespace at end of file.
  stream.peek();
  if (token == "" && stream.eof()) {
    return true;
  }
  return false;
}


void NRLib2::WriteBinaryShort(std::ostream& stream, 
                              short s,
                              Endianess file_format)
{
  char buffer[2];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteUInt16BE(buffer, static_cast<unsigned short>(s));
    break;
  case END_LITTLE_ENDIAN:
    WriteUInt16LE(buffer, static_cast<unsigned short>(s));
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 2)) {
    throw Exception("Error writing to stream.");
  }
}


short NRLib2::ReadBinaryShort(std::istream& stream, 
                              Endianess file_format)
{
  unsigned short us;
  char buffer[2];

  if (!stream.read(buffer, 2)) {
    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseUInt16BE(buffer, us);
    break;
  case END_LITTLE_ENDIAN:
    ParseUInt16LE(buffer, us);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return static_cast<short>(us);
}


void NRLib2::WriteBinaryInt(std::ostream& stream, 
                            int i,
                            Endianess file_format)
{
  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteUInt32BE(buffer, static_cast<unsigned int>(i));
    break;
  case END_LITTLE_ENDIAN:
    WriteUInt32LE(buffer, static_cast<unsigned int>(i));
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


int NRLib2::ReadBinaryInt(std::istream& stream, 
                          Endianess file_format)
{
  unsigned int ui;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseUInt32BE(buffer, ui);
    break;
  case END_LITTLE_ENDIAN:
    ParseUInt32LE(buffer, ui);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return static_cast<int>(ui);
}


void NRLib2::WriteBinaryFloat(std::ostream& stream, 
                              float f, 
                              Endianess file_format)
{
  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    WriteIEEEFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


float NRLib2::ReadBinaryFloat(std::istream& stream, 
                              Endianess file_format)
{
  float f;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    ParseIEEEFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
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
    WriteIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    WriteIEEEDoubleLE(buffer, d);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 8)) {
    throw Exception("Error writing to stream.");
  }
}


double NRLib2::ReadBinaryDouble(std::istream& stream, 
                                Endianess file_format)
{
  double d;
  char buffer[8];

  if (!stream.read(buffer, 8)) {
    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    ParseIEEEDoubleLE(buffer, d);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return d;
}


void NRLib2::WriteBinaryIbmFloat(std::ostream& stream, 
                                 float f, 
                                 Endianess file_format)
{
  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteIBMFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    WriteIBMFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


float NRLib2::ReadBinaryIbmFloat(std::istream& stream, 
                                 Endianess file_format)
{
  float f;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    throw Exception("Error reading from stream.");
  }  

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIBMFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    ParseIBMFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return f;
}
