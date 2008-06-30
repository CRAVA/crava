#include "corr.h"
#include "model.h"
#include <stdio.h>

using namespace NRLib2;

#include "lib/log.h"
Corr::Corr(float **pointVar0, float **Var0, float *CorrT, int n, float dt, 
           NRLib2::RegularSurface<double> * CorrXY)
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

const RegularSurface<double>* Corr::getCorrXY() const
{
  return (const RegularSurface<double>*) CorrXY_;
}

void Corr::dumpResult() const
{
  char * filename1= LogKit::makeFullFileName("PriorVar0",".dat");
  char * filename2= LogKit::makeFullFileName("PriorCorrTUnFiltered",".dat");
  char * filename3= LogKit::makeFullFileName("PriorCorrXY",".irap");
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
  LogKit::writeLog("\n");
  LogKit::writeLog("                         ln Vp     ln Vs    ln Rho         \n");
  LogKit::writeLog("--------------------------------------------------------------------\n");
  LogKit::writeLog("Parameter variances:   %.1e   %.1e   %.1e (used by program)\n",Var0_[0][0],Var0_[1][1],Var0_[2][2]);
  LogKit::writeLog("Well-log  variances:   %.1e   %.1e   %.1e                  \n",pointVar0_[0][0],pointVar0_[1][1],pointVar0_[2][2]);

  float corr01 = Var0_[0][1]/(sqrt(Var0_[0][0]*Var0_[1][1]));
  float corr02 = Var0_[0][2]/(sqrt(Var0_[0][0]*Var0_[2][2]));
  float corr12 = Var0_[1][2]/(sqrt(Var0_[1][1]*Var0_[2][2]));
  LogKit::writeLog("\n");
  LogKit::writeLog("Parameter correlations:\n");
  LogKit::writeLog("           ln Vp     ln Vs    ln Rho \n");
  LogKit::writeLog("-------------------------------------\n");
  LogKit::writeLog("ln Vp      %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::writeLog("ln Vs                %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::writeLog("ln Rho                         %5.2f \n",1.0f);
}
