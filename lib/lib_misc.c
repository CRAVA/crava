#include <string.h>
#include <math.h>
#include <ctype.h>
#include "lib/lib_misc.h"
#include "lib/global_def.h"


int equalRealDoubles( double value1, double value2, int precision)
{
  double difference;
  difference = fabs(value1 - value2);

  if ( difference < pow(0.1,precision)*0.5 )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}



char *uppercase(char * text)
{
  int i=0;
  while(text[i] != '\0')
  {
    text[i] = (char)toupper(text[i]);
    i++;
  }
  return text;
}

int
findEnd(char * seek, int start, char * find)
{
  int i = start;
  int end = (int)(strlen(seek));
  int ok = (int)(strlen(find));
  int flag = 0;
  while(flag < ok && i < end)
  {
    if(seek[i] == find[flag])
      flag++;
    else
      flag = 0;
    i++;
  }
  if(flag == ok)
    return(i-1);
  else
    return(-1);
}

int 
findClosestFactorableNumber(int leastint)
{
  int i,j,k,l,m,n;
  int factor   =       1;

  int maxant2    = (int) (ceil((double) log((float)(leastint)) / log(2.0f) ));
  int maxant3    = (int) (ceil((double) log((float)(leastint)) / log(3.0f) ));
  int maxant5    = (int) (ceil((double) log((float)(leastint)) / log(5.0f) ));
  int maxant7    = (int) (ceil((double) log((float)(leastint)) / log(7.0f) ));
  int maxant11   = (int) 0;
  int maxant13   = (int) 0;
  int closestprod= (int) pow(2.0f,maxant2);

  // kan forbedres ved aa trekke fra i endepunktene.i for lokkene
  for(i=0;i<maxant2+1;i++)
    for(j=0;j<maxant3+1;j++)
      for(k=0;k<maxant5+1;k++)
        for(l=0;l<maxant7+1;l++)
          for(m=0;m<maxant11+1;m++)
            for(n=maxant11;n<maxant13+1;n++)
            {
              factor = (int) (pow(2.0f,i)*pow(3.0f,j)*pow(5.0f,k)*
                pow(7.0f,l)*pow(11.0f,m)*pow(13.0f,n));
              if ((factor >=  leastint) &&  (factor <  closestprod))
              {
                closestprod=factor;
              }
            }
            return closestprod;
}


void
readUntilStop(int pos, char * in, char * out  ,const char stop)
{ 
  //  reads a string from position pos+1 
  //  until the terminating char given in stop 
  //  returns the string   excluding the stop

  int   i = 0; 

  while( ( (out[i]=in[pos+1+i]) != stop ) && in[pos+1+i] != '\0' )
  {
    i++;
  }
}

void
readToEOL(FILE * file)
{
  char rc = (char)fgetc(file);
  while(rc != '\n' && rc != EOF)
    rc = (char)fgetc(file);
}

int
isNumber(char * input)
{
  if((input[0] >= '0' && input[0] <= '9') ||
    (input[0] == '-' && input[1] >= '0' && input[1] <= '9'))
    return(1);
  else
    return(0);
}


void
swapDoubles(double * data, int nData)
{
  char tmp;
  int i, j;
  char * cdata = (char *) data;
  int ds = sizeof(double);
  int hds = ds/2;
  for(i=0;i<nData*ds;i+=ds)
  {
    for(j=0;j<hds;j++)
    {
      tmp = cdata[i+j];
      cdata[i+j] = cdata[i+ds-1-j];
      cdata[i+ds-1-j] = tmp;
    }
  }
}


void
swap4Bytes(char * data, int nData)
{
  char tmp;
  int i;
  for(i=0;i<nData*4;i+=4)
  {
    tmp = data[i];
    data[i] = data[i+3];
    data[i+3] = tmp;
    tmp = data[i+1];
    data[i+1] = data[i+2];
    data[i+2] = tmp;
  }
}

void
swap2Bytes(char * data, int nData)
{
  char tmp;
  int i;
  for(i=0;i<nData*2;i+=2)
  {
    tmp = data[i];
    data[i] = data[i+1];
    data[i+1] = tmp;
  }
}


