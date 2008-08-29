// $Id$
#ifndef COMMONHEADER_HPP
#define COMMONHEADER_HPP

#include <cstdio>
#include <iostream>
#include <fstream>
#include "segy.hpp"

class SegyGeometry;

/// The textual header of a SEG Y file.
/// Also known as the EBCDIC header.
class TextualHeader
{
public:
  /// Default constructor.
  TextualHeader();

  /// The old standard header.
  static TextualHeader standardHeader();

  /// Set a line in the textual header.
  /// @param[in] lineNo Line number.
  /// @param[in] text   Text for the given line. ASCII encoded. Max 75 char.
  /// @return 0 on success, -1 if error in input, 1 if the text was truncated.
  int setLine(int lineNo, const char* text);

  /// Write header to file.
  /// @param[in] file  Output file.
  void write(std::ostream& file) const;

private:
  /// output buffer in EBCDIC encoding.
  unsigned char buffer_[3200];
  //const char * buffer_;
};


class BinaryHeader
{
public:
  /// Constructor
  BinaryHeader(std::istream& file);
  // constructor for writing
  BinaryHeader();
/// Update variables
  void Update(std::istream& file);
  void write(std::ostream& file, SegyGeometry *geometry, double dz, int nz) const;
  short getFormat() {return(format_);}
  int getLino() {return(lino_);}
  short getHns() {return(hns_);}
  short getHdt() {return(hdt_);}
private:
  short format_;
  int lino_;
  short hns_;
  short hdt_;
};

#endif
