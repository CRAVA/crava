#ifndef CRAVA_SRC_PROGRAM_H
#define CRAVA_SRC_PROGRAM_H

#include <string>

/**
   \class Program program.h crava/src/program.h
   Class to hold global program data such 
   as version numbers, purpose, and the like.
   Version numbers follow the gnu system 
   (major.minor.patch).
*/

class Program
{

public:
  Program(unsigned int major,
          unsigned int minor,
          unsigned int patch,
          int          licence_days);
  
  ~Program(void); 

  const unsigned int   GetMajor(void)        const { return major_        ;}  
  const unsigned int   GetMinor(void)        const { return minor_        ;} 
  const unsigned int   GetPatch(void)        const { return patch_        ;}
  const int            GetLicenceDays(void)  const { return licence_days_ ;}
  
private:      
  void                 CheckForLicenceExpiration(int licence_days) const;
  time_t               TimeOfCompilation(void) const;

  const unsigned int   major_;           ///< Major version number  
  const unsigned int   minor_;           ///< Minor version number
  const unsigned int   patch_;           ///< Patch number
  const int            licence_days_;    ///< Program validity in days (-1 = infinity)
};

#ifndef DOXYGEN_SKIP 
#define  ACTION_END
#elif !defined ACTION_END
#error program.h is part of a cyclic dependency structure
#endif
#endif
