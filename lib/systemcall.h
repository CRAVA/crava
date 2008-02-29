#ifndef SYSTEM_H
#define SYSTEM_H

#include<string>

class SystemCall
{
public:
  static const char        * getHostName(void);
  static const char        * getUserName(void);
  static const char        * getCurrentTime(void); 

  static const std::string   getTime(void) { return __TIME__ ; }
  static const std::string   getDate(void) { return __DATE__  ;} 
};
#endif

