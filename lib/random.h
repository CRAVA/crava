#ifndef RANDOM_H
#define RANDOM_H

#include <string.h>

static const int MULTIPLIER = 69069;
static const int SHIFT = 1;
static const int MODULUS = 256*256*256*64;
static const double INVMOD = ((double) 1 / ((double) MODULUS))/((double) 4); 

class RandomGen{
public:
  RandomGen(unsigned int seed);
  RandomGen(char *filename); //NB: Validity of filename must be externally checked
  ~RandomGen();
  static double rnorm01();
  int writeSeedFile(char *filename) const;
static double unif01();

private:
  static double g(double x);
  
  static unsigned int seed_;
  char * seedfile_;

};
#endif
