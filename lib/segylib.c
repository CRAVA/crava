#include "lib/segylib.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  (c) Copyright IBM European Petroleum Application Center., unpublished work.
**      Sample code, no guarantees whatsoever.
*F
*F   gen/Ieee2Ibm.C
*F
*F   PURPOSE
*F
*F   <Describe purpose of the function or functions, and describe the logical
*F   the relationship of the functions in this file.>
*F
*F   SPECIAL REQUIREMENTS AND NOTES
*F
**   
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**/

#define IEEEMAX  0x7FFFFFFF
#define IEMAXIB  0x611FFFFF
#define IEMINIB  0x21200000


typedef unsigned int UInt4;


void
pb_Ibm2Ieee (
             void	*src,
             float	*dsT,
             int	nel)
{
  static int it[8] = { 0x21800000, 0x21400000, 0x21000000, 0x21000000,
    0x20c00000, 0x20c00000, 0x20c00000, 0x20c00000 };
  static int mt[8] = { 8, 4, 2, 2, 1, 1, 1, 1 };
  UInt4 manthi, iexp, inabs, *in, *out;
  int i, ix;

  in = (UInt4 *)src; out = (UInt4 *)dsT;
  for ( i=0; i<nel; i++ )
  {
    manthi = in[i] & 0x00ffffff;
    ix     = manthi >> 21;
    iexp   = ( ( in[i] & 0x7f000000 ) - it[ix] ) << 1;
    manthi = manthi * mt[ix] + iexp;
    inabs  = in[i] & 0x7fffffff;
    if ( inabs > IEMAXIB ) manthi = IEEEMAX;
    manthi = manthi | ( in[i] & 0x80000000 );
    out[i] = ( inabs < IEMINIB ) ? 0 : manthi;
  }
}


void
pb_Ieee2Ibm (
             void	*src,
             float	*dsT,
             int	nel
             )
{
  static int it[4] = { 0x21200000, 0x21400000, 0x21800000, 0x22100000 };
  static int mt[4] = { 2, 4, 8, 1 };
  UInt4 manthi, iexp, ix, *in, *out;
  int i;

  in = (UInt4 *)src; out = (UInt4 *)dsT;
  for ( i=0; i<nel; i++ ) 
  {
    ix     = ( in[i] & 0x01800000 ) >> 23;
    iexp   = ( ( in[i] & 0x7e000000 ) >> 1 ) + it[ix];
    manthi = ( mt[ix] * ( in[i] & 0x007fffff) ) >> 3;
    manthi = (manthi + iexp) | ( in[i] & 0x80000000 );
    out[i] = ( in[i] & 0x7fffffff ) ? manthi : 0;
  }
}


