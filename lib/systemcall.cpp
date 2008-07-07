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


const char *
SystemCall::getUserName()
{
  char * userName = new char[MAX_STRING]; 
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
    DWORD nUserName = sizeof(userName);
    if (GetUserName(userName, &nUserName)) 
    {
      if (userName == NULL) 
        sprintf(userName,"Empty user name found\n");
    }
    else 
      sprintf(userName,"User name not found.\n");
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
    sprintf(userName,pw->pw_name); 
    delete [] buffer;
#endif
    return userName;
}

const char *
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
  if ( err != 0 ) {
    sprintf(hostname,"*Not set*");
  }
  return hostname;
}

const char *
SystemCall::getCurrentTime(void) 
{ 
  time_t raw; 
  time(&raw); 

  char * buffer = new char[MAX_STRING];

#if _MSC_VER >= 1400
  const size_t size = 200;
  ctime_s(buffer, size, &raw);
#else
  sprintf(buffer,ctime(&raw));
#endif
  return buffer;
} 
