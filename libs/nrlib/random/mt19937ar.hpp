// $Id: mt19937ar.hpp 910 2011-10-26 09:48:12Z perroe $

#ifndef math_mt19937ar_h
#define math_mt19937ar_h

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Copyright (C) 2005, Mutsuo Saito
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/





/**
   @file
   Declares random number generator.
*/


/**
   Namespace for the Mersenne Twister.
   A very fast random number generator
   of period 2^(19937) - 1.
   @author Makoto Matsumoto
   @author Takuji Nishimura
   @author Mutsuo Saito
   @note Since 2001/4/6 the authors let the Mersenne Twister
   be used freely for any purpose, including commercial use, see
   <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/elicense.html">http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/elicense.html</a>.
   @attention These functions are not meant to be used directly.
   Use the Random interface class in the NRLib namespace
*/


namespace MT19937 {

#ifndef M64
  /// 32 bit seed.
  typedef long unsigned int seed_type;
  typedef int key_type;
#else
  /// 64 bit seed.
  typedef long long unsigned int seed_type;
  typedef long long unsigned int key_type;
#endif


  /**
     Initializer.
     Initializes the Mersenne Twister with a positive seed.
     @param[in] s seed
  */
  void init_genrand(seed_type s);



  /** Initializer.
      Initializes the Mersenne Twister with an array.
      @param[in] init_key is the array for initializing keys.
      @param[in] key_length is its length.
  */
  void init_by_array(seed_type* init_key,
             key_type key_length);

#ifndef M64
  /**
     Random unsigned long integer.
  */
  unsigned long genrand_int32(void);

  /**
     Random long integer.
  */
  long genrand_int31(void);
#else
 /**
     Random unsigned long long integer.
  */
  unsigned long long genrand_int64(void);

  /**
     Random long long integer.
  */

  long long genrand_int63(void);
#endif

  /**
     Random real. Generates a random number on the [0,1] interval.
     @author Isaku Wada
  */
  double genrand_real1(void);

  /**
     Random real. Generates a random number on the [0,1) interval.
     @author Isaku Wada
  */
  double genrand_real2(void);

  /**
     Random real. Generates a random number on the (0,1) interval.
     @author Isaku Wada
  */
  double genrand_real3(void);

  /*
     Random real. Generates generates a random number on [0,1)
     with 53-bit resolution.
     @author Isaku Wada
  */
  //  double genrand_res53(void);


}


//#ifndef DOXYGEN_SKIP
//#define math_mt19937ar_end
//#elif !defined math_mt19937ar_end
//#error math_mt19937ar.h is part of a cyclic dependency structure
//#endif
#endif

