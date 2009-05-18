#include "lib/systemcall.h"
#include "lib/global_def.h"
#include "nrlib/iotools/logkit.hpp"

#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
#include <winsock2.h>
#include <windows.h>
#else
#include <unistd.h>
#include <pwd.h>
#endif

#include <ctime>
#include <string.h>

const std::string
SystemCall::getUserName()
{
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
  std::string strUserName; 
  char * userName = new char[MAX_STRING]; 
  DWORD nUserName = sizeof(userName);
  if (GetUserName(userName, &nUserName)) 
  {
    if (userName == NULL) 
      strUserName("Empty user name found");
  }
  else 
  {
    strUserName("User name not found");
  }
  delete [] userName;
#else
  struct passwd on_the_stack; 
  struct passwd* pw = 0;  
  uid_t  uid    = geteuid();
  size_t size   = sysconf(_SC_GETPW_R_SIZE_MAX);
  char * buffer = new char[size];
#ifdef NO_GETPWUID_R_PTR  
  pw = getpwuid_r(uid, &on_the_stack, buffer, size);  
#else
  getpwuid_r(uid, &on_the_stack, buffer, size, &pw);  
#endif
  std::string strUserName(pw->pw_name); 
  delete [] buffer;
#endif
  return strUserName;
}

const std::string
SystemCall::getHostName()
{
  char * hostname = new char[MAX_STRING];
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
  // WORD wVersionRequested;
  // WSADATA wsaData;
  int err=1;
  // wVersionRequested = MAKEWORD( 2, 2 );
  // err = WSAStartup( wVersionRequested, &wsaData );
  // if ( err == 0 )
  //   gethostname(hostname, MAX_STRING);
  // WSACleanup( );
#else 
  int err = gethostname(hostname, MAX_STRING);
#endif

  std::string strHostname;
  
  if ( err == 0 )
    strHostname = std::string(hostname); 
  else
    strHostname = std::string("*Not set*");

  delete [] hostname;

  return strHostname;
}

const std::string
SystemCall::getCurrentTime(void) 
{ 
  time_t raw; 
  time(&raw); 
  
#if _MSC_VER >= 1400
  const size_t size = 200;
  char * buffer = new char[MAX_STRING];
  ctime_s(buffer, size, &raw);
  std::string strBuffer = std::string(buffer);
  delete [] buffer;
#else
  std::string strBuffer = std::string(ctime(&raw));
#endif
  
  return strBuffer;
} 
