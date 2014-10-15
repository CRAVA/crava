/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <string.h>

#include "lib/utils.h"
#include "src/definitions.h"

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "nrlib/iotools/logkit.hpp"

//------------------------------------------------------------
void
Utils::writeTitler(const std::string & text)
{
  const int length = static_cast<int>(text.length());
  LogKit::LogFormatted(LogKit::Low,"\n%s\n",text.c_str());
  for (int i = 0 ; i < length ; i++)
    LogKit::LogFormatted(LogKit::Low,"-");
  LogKit::LogFormatted(LogKit::Low,"\n");
}

//------------------------------------------------------------
void
Utils::copyVector(const int * from,
                  int       * to,
                  int         ndim)
{
  for (int i = 0 ; i < ndim ; i++)
    to[i] = from[i];
}

//------------------------------------------------------------
void
Utils::copyVector(const float * from,
                  float       * to,
                  int           ndim)
{
  for (int i = 0 ; i < ndim ; i++)
    to[i] = from[i];
}

//------------------------------------------------------------
void
Utils::copyMatrix(const float ** from,
                  float       ** to,
                  int            ndim1,
                  int            ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++)
    for (int j = 0 ; j < ndim2 ; j++)
      to[i][j] = from[i][j];
}

//------------------------------------------------------------
void
Utils::writeVectorToFile(const std::string & fileName,
                         float             * vector,
                         int                  ndim)
{
  std::ofstream fout;
  NRLib::OpenWrite(fout,fileName);
  fout << std::setprecision(6);
  for (int i = 0 ; i < ndim ; i++) {
    fout << std::setw(6)  << i + 1
         << std::setw(12) << vector[i] << "\n";
  }
  fout << std::endl;
  fout.close();
}

//------------------------------------------------------------
void
Utils::writeVectorToFile(const std::string        & fileName,
                         const std::vector<float> & vector)
{
  std::ofstream fout;
  NRLib::OpenWrite(fout,fileName);
  fout << std::setprecision(6);
  for (size_t i = 0 ; i < vector.size() ; i++) {
    fout << std::setw(6)  << i + 1
         << std::setw(12) << vector[i] << "\n";
  }
  fout << std::endl;
  fout.close();
}

//------------------------------------------------------------
void
Utils::writeVector(float * vector,
                   int     ndim)
{
  for (int i = 0 ; i < ndim ; i++) {
    LogKit::LogFormatted(LogKit::Low,"%10.6f ",vector[i]);
  }
  LogKit::LogFormatted(LogKit::Low,"\n");
}

//------------------------------------------------------------
void
Utils::writeVector(double * vector,
                   int      ndim)
{
  for (int i = 0 ; i < ndim ; i++) {
    LogKit::LogFormatted(LogKit::Low,"%10.6f ",vector[i]);
  }
  LogKit::LogFormatted(LogKit::Low,"\n");
}

//------------------------------------------------------------
void
Utils::writeMatrix(float ** matrix,
                   int      ndim1,
                   int      ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++) {
    for (int j = 0 ; j < ndim2 ; j++) {
      LogKit::LogFormatted(LogKit::Low,"%10.6f ",matrix[i][j]);
    }
    LogKit::LogFormatted(LogKit::Low,"\n");
  }
}

//------------------------------------------------------------
void
Utils::writeMatrix(double ** matrix,
                   int       ndim1,
                   int       ndim2)
{
  for (int i = 0 ; i < ndim1 ; i++) {
    for (int j = 0 ; j < ndim2 ; j++) {
      LogKit::LogFormatted(LogKit::Low,"%10.6f ",matrix[i][j]);
    }
    LogKit::LogFormatted(LogKit::Low,"\n");
  }
}

//------------------------------------------------------------
void
Utils::fft(fftw_real* rAmp,fftw_complex* cAmp,int nt)
{
  rfftwnd_plan p1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_real_to_complex(p1, rAmp, cAmp);
  fftwnd_destroy_plan(p1);
}

//------------------------------------------------------------
void
Utils::fftInv(fftw_complex* cAmp,fftw_real* rAmp,int nt)
{
  rfftwnd_plan p2 = rfftwnd_create_plan(1, &nt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_complex_to_real(p2, cAmp, rAmp);
  fftwnd_destroy_plan(p2);
  double sf = 1.0/double(nt);
  for(int i=0;i<nt;i++)
    rAmp[i]*=fftw_real(sf);
}
//-----------------------------------------------------------
int
Utils::findEnd(std::string & seek, int start, std::string & find)
{
  size_t i    = start;
  size_t end  = seek.length();
  size_t ok   = find.length();
  size_t flag = 0;

  while(flag < ok && i < end)
  {
    if(seek[i] == find[flag])
      flag++;
    else
      flag = 0;
    i++;
  }
  if(flag == ok)
    return(static_cast<int>(i)-1);
  else
    return(-1);
}
//-------------------------------------------------------------
void
Utils::readUntilStop(int pos, std::string & in, std::string & out, std::string read)
{
  /*
    reads a string from position pos+1
    until the terminating char given in stop
    returns the string excluding the stop
  */

  size_t stop = in.find(read,pos);
  if (stop != std::string::npos)
    stop = in.length();
  out = "";
  size_t i=0;
  while (i < stop-pos-1)
  {
    out.push_back(in[pos+1+i]);
    i++;
  }
}

void
Utils::ShiftTrace(std::vector<double> & trace,
                  bool        shift_up)
{
  std::vector<fftw_real> tmp_trace(trace.size());
  for(size_t i=0;i<trace.size();i++)
    tmp_trace[i] = static_cast<fftw_real>(trace[i]);
  ShiftTrace(&(tmp_trace[0]), trace.size(), shift_up);
  for(size_t i=0;i<trace.size();i++)
    trace[i] = static_cast<double>(tmp_trace[i]);
}



void
Utils::ShiftTrace(std::vector<fftw_real> & trace,
                  bool        shift_up)
{
  ShiftTrace(&(trace[0]), trace.size(), shift_up);
}


void
Utils::ShiftTrace(fftw_real * trace,
                  size_t      n_data,
                  bool        shift_up)
{
  size_t n_small   = 2*n_data;
  size_t fft_small = ((n_small / 2) + 1)*2;
  size_t n_large   = 2*n_small;
  size_t fft_large = ((n_large / 2) + 1)*2;

  std::vector<fftw_real> s_trace(fft_small,0);
  fftw_complex * s_trace_c = reinterpret_cast<fftw_complex*>(&(s_trace[0]));

  size_t i;
  for(i=0;i<n_data;i++)
    s_trace[i] = trace[i];
  size_t change = n_data+n_data/2;

  for(;i<change;i++)
    s_trace[i] = s_trace[n_data-1]*(1+cos(NRLib::Pi*static_cast<double>(i-n_data)/static_cast<double>(change-1-n_data)));
  for(;i<n_small;i++)
    s_trace[i] = s_trace[0]*(1+cos(NRLib::Pi*static_cast<double>(n_small-1-i)/static_cast<double>(n_small-1-change)));

  Utils::fft(&(s_trace[0]), s_trace_c, n_small);

  std::vector<fftw_real> l_trace(fft_large,0);
  fftw_complex * l_trace_c = reinterpret_cast<fftw_complex*>(&(l_trace[0]));
  for(size_t i=0;i<fft_small;i++)
    l_trace[i] = s_trace[i];

  Utils::fftInv(l_trace_c,&(l_trace[0]), n_large);

  if(shift_up == true){
    for(size_t i=0;i<n_data;i++)
      trace[i] = 2*l_trace[2*i+1];
  }
  else {
    for(size_t i=1;i<n_data;i++)
      trace[i] = 2*l_trace[2*i-1];
    trace[0] = l_trace[n_large-1];
  }
}
