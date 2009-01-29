#include <math.h>
#include <stdio.h>

#include "lib/random.h"

RandomGen::RandomGen(unsigned int seed)
{
  seed_     = seed;
  seedfile_ = "";
}


RandomGen::RandomGen(const std::string & filename)
{
  seedfile_ = filename;
  FILE * file = fopen(seedfile_.c_str(),"r");
  fscanf(file,"%u",&seed_);
  fclose(file);
}

RandomGen::~RandomGen()
{
  if(seedfile_ != "")
    writeSeedFile(seedfile_);
}

//
// Marsaglia-Bray's method, see Ripley, p. 84.
//
double RandomGen::rnorm01()
{
  double u, u1, u2, u3;
  double c, x, v1, v2, w, s, t;

  u = unif01();
  if(u<0.8638)

  {
    u1 = unif01();
    u2 = unif01();
    u3 = unif01();
    x = 2.0*(u1+u2+u3)-3.0;
  }
  else if(u<0.9745)
  {
    u1 = unif01();
    u2 = unif01();
    x = 1.5*(u1+u2-1.0);
  }
  else if(u<0.9973002039)
  {
    do
    {
      u1 = unif01();
      u2 = unif01();
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
        u1 = unif01();
        u2 = unif01();
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
    if(fabs(s)>3.0)
      x = s;
    else
      x = t;
  }
  return x;


}

double RandomGen::g(double x)
{
  double a = 17.49731196;
  double b = 2.36785163;
  double c = 2.15787544;
  double absx = fabs(x);
  double result = a*exp(-x*x/2.0);
  if(absx<1.0)
    result -= 2.0*b*(3.0-x*x)+c*(1.5-absx);
  else if(absx<1.5)
    result -= b*(3.0-absx)*(3.0-absx)+c*(1.5-absx);
  else
    result -=b*(3.0-absx)*(3.0-absx);
  return result;

}

double RandomGen::unif01()
{
  double x;
  seed_ = MULTIPLIER * seed_ +SHIFT;
  x = static_cast<double>(seed_) *INVMOD;
  return x;
}

int RandomGen::writeSeedFile(const std::string & filename) const
{
  FILE *file;
  int error = 0;
  file = fopen(filename.c_str(),"w");
  if(file == NULL)
    error = 1;
  else
  {
    fprintf(file,"%d",seed_);
    fclose(file);
  }
  return(error);
}

unsigned int RandomGen::seed_ = 23665; //Default seed
