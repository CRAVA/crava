#ifndef RPLIB_ROCK_H
#define RPLIB_ROCK_H

#include <assert.h>
#include <vector>


// Abstract rock class.
class Rock {
public:
  Rock();
  virtual ~Rock();

  // Assignment operator, not yet implemented.
  /*Rock& operator=(const Rock& rhs);*/

  virtual Rock * Clone() const = 0;

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const = 0;

  // Rock is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Input parameters:
  //      delta_time : the set of previous and present incremental time steps
  //      rock : the set of previous rock samples
  // Recommended in implementation: assert(delta_time.size() == rock.size() + 1);
  virtual Rock * Evolve(const std::vector<int>         & delta_time,
                        const std::vector< Rock * >    & rock) const = 0;


protected:

};

#endif
