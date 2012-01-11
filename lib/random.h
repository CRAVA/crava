#ifndef RANDOM_H
#define RANDOM_H

#include <string.h>
#include <string>

static const int MULTIPLIER = 69069;
static const int SHIFT = 1;
static const int MODULUS = 256*256*256*64;
static const double INVMOD = (static_cast<double>(1) / (static_cast<double>(MODULUS)))/(static_cast<double>(4));
//WAS: static const double INVMOD = ((double) 1 / ((double) MODULUS))/((double) 4);

class RandomGen{
public:
  RandomGen(unsigned int seed);
  RandomGen(const std::string & filename); //NB: Validity of filename must be externally checked
  ~RandomGen();

  int writeSeedFile(const std::string & filename) const;

  static double rnorm01();
  static double unif01();

private:
  static double g(double x);

  static unsigned int seed_;
  std::string         seedfile_;

};
#endif
