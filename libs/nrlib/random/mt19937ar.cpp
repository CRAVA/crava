// $Id: mt19937ar.cpp 910 2011-10-26 09:48:12Z perroe $


#include "mt19937ar.hpp"

/**
   @file
   Where the Mersenne Twister is set up.
*/

#ifdef M64 // 64 bit machine
#include "mt19937_64.hpp"
#else
#include "mt19937_32.hpp"
#endif

#
