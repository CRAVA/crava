// $Id: randomgenerator.cpp 1175 2013-05-27 07:41:03Z perroe $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "randomgenerator.hpp"
#include "../exception/exception.hpp"

#include <ctime>
#include <cmath>
#include <fstream>

using namespace NRLib;

RandomGenerator::RandomGenerator()
{
  is_initialized_ = false;
  start_seed_     = 0;
}

RandomGenerator::~RandomGenerator()
{
}


void
RandomGenerator::Initialize()
{
  unsigned long seed = static_cast<unsigned long>(time(0));
  InitializeMT(seed);
  start_seed_ = DrawUint32();
  is_initialized_ = true;
  InitializeMT(start_seed_);
}

void
RandomGenerator::Initialize(unsigned long int seed)
{
  start_seed_ = seed;
  is_initialized_ = true;
  InitializeMT(start_seed_);
}


double
RandomGenerator::Norm01()
{
  double u, u1, u2, u3;
  double c, x, v1, v2, w, s, t;

  u = Unif01();
  if (u<0.8638)

  {
    u1 = Unif01();
    u2 = Unif01();
    u3 = Unif01();
    x = 2.0*(u1+u2+u3)-3.0;
  }
  else if (u<0.9745)
  {
    u1 = Unif01();
    u2 = Unif01();
    x = 1.5*(u1+u2-1.0);
  }
  else if (u<0.9973002039)
  {
    do
    {
      u1 = Unif01();
      u2 = Unif01();
      x = 6.0*u1-3.0;
    }
    while(0.358*u2>g(x));
  }
  else
  {
    do
    {
      do
      {
        u1 = Unif01();
        u2 = Unif01();
        v1 = 2.0*u1-1.0;
        v2 = 2.0*u2-1.0;
        w = v1*v1+v2*v2;
      }
      while(w>=1.0);
      c = sqrt((9.0-2.0*log(w))/w);
      s = c*v1;
      t = c*v2;
    }
    while((fabs(s)<=3.0) && (fabs(t)<=3.0));
    if (fabs(s)>3.0)
      x = s;
    else
      x = t;
  }
  return x;
}


unsigned long
RandomGenerator::GetStartSeed()
{
  if (!is_initialized_) {
    throw Exception("Random number generator is not initalized.");
  }
  return start_seed_;
}



double
RandomGenerator::g(double x)
{
  double a = 17.49731196;
  double b = 2.36785163;
  double c = 2.15787544;
  double absx = fabs(x);
  double result = a*exp(-x*x/2.0);
  if (absx<1.0)
    result -= 2.0*b*(3.0-x*x)+c*(1.5-absx);
  else if (absx<1.5)
    result -= b*(3.0-absx)*(3.0-absx)+c*(1.5-absx);
  else
    result -=b*(3.0-absx)*(3.0-absx);
  return result;
}


void
RandomGenerator::InitializeMT(unsigned long seed)
{
  dsfmt_init_gen_rand(&dsfmt, seed);
}

