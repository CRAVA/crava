
#ifndef UTILS_H
#define UTILS_H

class Utils
{
public:
  static void    writeTitler(const char * text);

  static void    copyVector(const float * from,
                            float       * to,
                            int           ndim);
  static void    copyMatrix(const float ** from,
                            float       ** to,
                            int            ndim1,
                            int            ndim2);

  static void    writeVector(float * vector,
                             int     ndim);
  static void    writeVector(double * vector,
                             int      ndim);

  static void    writeMatrix(float ** matrix,
                             int      ndim1,
                             int      ndim2);  
  static void    writeMatrix(double ** matrix,
                             int       ndim1,
                             int       ndim2);  
};

#endif
