#ifndef LIB_MISC_H
#define LIB_MISC_H

#include <stdio.h>
#ifdef __cplusplus
extern "C"
{
#endif


  int    equalRealDoubles( double value1, double value2, int precision);
  char * uppercase(char * text);
  int    findEnd(char * seek, int start, char * find);
  int    findClosestFactorableNumber(int leastint);
  void   readUntilStop(int pos, char * in, char  * out, const char stop);
  void   readToEOL(FILE * file);
  int    isNumber(char * input);
  void   swapDoubles(double * data, int nData);
  void   swap4Bytes(char * data, int nData);
  void   swap2Bytes(char * data, int nData);


#ifdef __cplusplus
}
#endif

#endif
