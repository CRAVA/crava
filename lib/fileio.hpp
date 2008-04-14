// $Id: fileio.hpp 30 2008-04-04 12:39:42Z perroe $

#ifndef NRLIB_FILEIO_HPP
#define NRLIB_FILEIO_HPP

#include <string>

namespace NRLib2 {
  enum Endianess {END_LITTLE_ENDIAN, END_BIG_ENDIAN};

  /// \brief Check if end of file is reached.
  /// \return True if end of file is reached, else false.
  bool CheckEndOfFile(std::istream& stream);

  // ---------------------------------
  // 4-byte IEEE floating point number
  // ---------------------------------

  /// \brief Write a 4-byte float on standard IEEE format.
  void WriteBinaryFloat(std::ostream& stream, 
                        float f,
                        Endianess file_format = END_BIG_ENDIAN);

  /// Read a 4-byte float on standard IEEE format.
  float ReadBinaryFloat(std::istream& stream, 
                        Endianess file_format = END_BIG_ENDIAN);

  /// Write an array of 4-byte floats on standard IEEE format.
  template<typename iterator>
  void WriteBinaryFloatArray(std::ostream& stream,
                             iterator begin,
                             iterator end,
                             Endianess file_format = END_BIG_ENDIAN);

  /// Read an array of 4-byte floats on standard IEEE format.
  template<typename iterator>
  iterator ReadBinaryFloatArray(std::istream& stream,
                                iterator begin,
                                int n,
                                Endianess file_format = END_BIG_ENDIAN);

  // ---------------------------------
  // 8-byte IEEE floating point number
  // ---------------------------------

  /// Write a 8-byte float on standard IEEE format.
  void WriteBinaryDouble(std::ostream& stream, 
                         double d,
                         Endianess file_format = END_BIG_ENDIAN);

  /// Read a 8-byte float on standard IEEE format.
  double ReadBinaryDouble(std::istream& stream, 
                          Endianess file_format = END_BIG_ENDIAN);

  /// Write an array of 8-byte floats on standard IEEE format.
  template<typename iterator>
  void WriteBinaryDoubleArray(std::ostream& stream,
                              iterator begin,
                              iterator end,
                              Endianess file_format = END_BIG_ENDIAN);

  /// Read an array of 8-byte floats on standard IEEE format.
  template<typename iterator>
  iterator ReadBinaryDoubleArray(std::istream& stream,
                                 iterator begin,
                                 int n,
                                 Endianess file_format = END_BIG_ENDIAN);

//   void WriteBinaryInt32(std::ostream& stream, 
//                         float f,
//                         Endianess file_format = BIG_ENDIAN);

//   float ReadBinaryInt32(std::istream& stream, 
//                         Endianess file_format = BIG_ENDIAN);

//   template<typename iterator>
//   WriteBinaryInt32Array(std::ostream& stream,
//                        iterator begin,
//                        iterator end,
//                        Endianess file_format = BIG_ENDIAN);

//   template<typename iterator>
//   ReadBinaryInt32Array(std::ostream& stream,
//                        iterator begin,
//                        int n,
//                        Endianess file_format = BIG_ENDIAN);

}

#include "fileio_impl.hpp"

#endif // NRLIB_FILEIO_HPP
