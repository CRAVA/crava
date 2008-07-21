#include <stdio.h>

#include "src/corr.h"
#include "src/definitions.h"
#include "src/model.h"
#include "nrlib/iotools/logkit.hpp"

Corr::Corr(float **pointVar0, float **Var0, float *CorrT, int n, float dt, 
           Surface * CorrXY)
{
  Var0_ = Var0;
  pointVar0_ = pointVar0;
  CorrT_ = CorrT;
  n_ = n;
  dt_ = dt;
  CorrXY_ = CorrXY;
}

Corr::~Corr(void)
{
  int i;
  for(i=0;i<3;i++)
    delete [] Var0_[i];
  for(i=0;i<3;i++)
    delete [] pointVar0_[i];
  delete [] pointVar0_;
  delete [] Var0_;
  delete [] CorrT_;
  delete CorrXY_;
}

const float* Corr::getCorrT(int &n, float &dt) const
{
  n = n_;
  dt = dt_;
  return CorrT_;

}

const float** Corr::getVar0() const
{
  return (const float **) Var0_;	
}

const NRLib2::RegularSurface<double>* Corr::getCorrXY() const
{
  return (const NRLib2::RegularSurface<double>*) CorrXY_;
}

void Corr::dumpResult() const
{
  char * filename1= ModelSettings::makeFullFileName("PriorVar0",".dat");
  char * filename2= ModelSettings::makeFullFileName("PriorCorrTUnFiltered",".dat");
  char * filename3= ModelSettings::makeFullFileName("PriorCorrXY",".irap");
  int i,j;
  FILE *file1 = fopen(filename1, "w");
  FILE *file2 = fopen(filename2, "w");
  FILE *file3 = fopen(filename3, "w");

  //NBNB Snakk med Ragnar om du saknar dette veldig.
  //irapgridWritept(file3,CorrXY_); 
  fclose(file3);
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      fprintf(file1,"%f ",Var0_[i][j]);
    }
    fprintf(file1,"\n");
  }
  fclose(file1);
  fprintf(file2,"%f\n",dt_);
  for(i=0;i<n_;i++)
  {
    fprintf(file2,"%f\n",CorrT_[i]);
  }
  fclose(file2);
  delete [] filename1;
  delete [] filename2;
  delete [] filename3;
}

void
Corr::setVar0(float ** var0)
{
  if(Var0_ != NULL) {
    int i;
    for(i=0;i<3;i++)
      delete [] Var0_[i];
    delete [] Var0_;
  }
  Var0_ = var0;
}

void Corr::printVariancesToScreen()
{
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"                         ln Vp     ln Vs    ln Rho         \n");
  LogKit::LogFormatted(LogKit::LOW,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"Parameter variances:   %.1e   %.1e   %.1e (used by program)\n",Var0_[0][0],Var0_[1][1],Var0_[2][2]);
  LogKit::LogFormatted(LogKit::LOW,"Well-log  variances:   %.1e   %.1e   %.1e                  \n",pointVar0_[0][0],pointVar0_[1][1],pointVar0_[2][2]);

  float corr01 = Var0_[0][1]/(sqrt(Var0_[0][0]*Var0_[1][1]));
  float corr02 = Var0_[0][2]/(sqrt(Var0_[0][0]*Var0_[2][2]));
  float corr12 = Var0_[1][2]/(sqrt(Var0_[1][1]*Var0_[2][2]));
  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"Parameter correlations:\n");
  LogKit::LogFormatted(LogKit::LOW,"           ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::LOW,"-------------------------------------\n");
  LogKit::LogFormatted(LogKit::LOW,"ln Vp      %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::LOW,"ln Vs                %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::LOW,"ln Rho                         %5.2f \n",1.0f);
}
