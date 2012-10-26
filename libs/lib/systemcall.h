/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SYSTEM_H
#define SYSTEM_H

#include<string>

class SystemCall
{
public:
  static const std::string  getHostName(void);
  static const std::string  getUserName(void);
  static const std::string  getCurrentTime(void);

  static const std::string  getTime(void) { return __TIME__ ; }
  static const std::string  getDate(void) { return __DATE__  ;}
};
#endif

