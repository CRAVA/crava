#include <stdlib.h>
#include <math.h>

#include "lib/global_def.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/krigingdata3d.h"
#include "src/welldata.h"
#include "src/bwellpt.h"
#include "src/model.h"
#include "src/definitions.h"

//---------------------------------------------------------------------
KrigingData3D::KrigingData3D(int ntot) 
  : data_(NULL),
    nd_(0)
{
  data_ = new CBWellPt * [ntot];
}

//---------------------------------------------------------------------
KrigingData3D::KrigingData3D(WellData ** wells,
                             int         nWells,
                             int         type) 
  : data_(NULL),
    nd_(0)
{
  //
  // Find total number of data
  //
  int ntot = 0;
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsPropThick()->getNumberOfBlocks();
    ntot += nBlocks;
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }

  data_ = new CBWellPt * [ntot];

  //
  // Fill with data
  //
  for (int w = 0 ; w < nWells ; w++) {

    BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
    const int nBlocks = bl->getNumberOfBlocks();
   
    const float * alpha;
    const float * beta;
    const float * rho;
    if (type == 1) {
      alpha = bl->getAlpha();
      beta  = bl->getBeta();
      rho   = bl->getRho();
    }
    else if (type == 2) {
      alpha = bl->getAlphaHighCutBackground();
      beta  = bl->getBetaHighCutBackground();
      rho   = bl->getRhoHighCutBackground();
    }
    else if (type == 3) {
      alpha = bl->getAlphaHighCutSeismic();
      beta  = bl->getBetaHighCutSeismic();
      rho   = bl->getRhoHighCutSeismic();
    }
    else if (type == 4) {
      alpha = bl->getAlphaSeismicResolution();
      beta  = bl->getBetaSeismicResolution();
      rho   = bl->getRhoSeismicResolution();
    }
    else {
      LogKit::LogFormatted(LogKit::LOW,"ERROR: Undefined log type %d\n",type);
      exit(1);
    }
    
    const int * ipos = bl->getIpos();
    const int * jpos = bl->getJpos();
    const int * kpos = bl->getKpos();
    
    addData(alpha, beta, rho, ipos, jpos, kpos, nBlocks);
  } 
  divide();
}

//---------------------------------------------------------------------
KrigingData3D::~KrigingData3D(void) 
{
  for (int i = 0 ; i < nd_ ; i++)
    delete data_[i];
  delete [] data_;
  data_ = NULL;
}

//---------------------------------------------------------------------
void
KrigingData3D::addData(const float * alpha,
                       const float * beta,
                       const float * rho,
                       const int   * ipos,
                       const int   * jpos,
                       const int   * kpos,
                       const int     nd)
{
  for (int m = 0 ; m < nd ; m++) 
  {
    const int i = ipos[m];
    const int j = jpos[m];
    const int k = kpos[m];

    if (alpha[m] != RMISSING || beta[m] != RMISSING || rho[m] != RMISSING) 
    {
      float a,b,r;
      
      if (alpha[m] != RMISSING)
        a = exp(alpha[m]);
      else 
        a = RMISSING;
      
      if (beta[m] != RMISSING)
        b = exp(beta[m]);
      else 
        b = RMISSING;
          
      if (rho[m] != RMISSING)
        r = exp(rho[m]);
      else 
        r = RMISSING;

      int index = gotBlock(i,j,k);
      if (index == -1) 
      {
        data_[nd_] = new CBWellPt(i, j, k);
        data_[nd_]->AddLog(a, b, r);
        nd_++;
      }
      else
      {
        data_[index]->AddLog(a, b, r);
        //
        // This is not a problem, but normally two wells should not share the same
        // set of blocks (with the possible exception of side-tracks).
        //
        LogKit::LogFormatted(LogKit::LOW,"\nNOTE: Blocked log with indices i,j,k = %d,%d,%d has been referred to twice. This is not a problem\n",i,j,k);
        LogKit::LogFormatted(LogKit::LOW,"\n      but may indicate that the grid is too coarse, or that you have wells with side-tracks?\n");
      }
    }
  }
}

//---------------------------------------------------------------------
int
KrigingData3D::gotBlock(int i, int j, int k) const
{
  for (int b = 0 ; b < nd_ ; b++)
  {
    if (data_[b]->GetI() == i && 
        data_[b]->GetJ() == j && 
        data_[b]->GetK() == k) 
      return b;
  }
  return -1;
}

//---------------------------------------------------------------------
void
KrigingData3D::divide(void)
{ 
  //
  // This method is only required if two wells share the same block,
  // which is not very likely in our case. Possibly therefore, the NOTE 
  // given above is better changed to a WARNING...
  //
  for (int b = 0 ; b < nd_ ; b++)
    data_[b]->Divide();
}

//---------------------------------------------------------------------
void 
KrigingData3D::writeToFile(const std::string type)
{ 
  std::string fileName = ModelSettings::makeFullFileName(type+"_KrigingData"+".dat");

  FILE * file = fopen(fileName.c_str(), "w");

  fprintf(file,"  i   j   k      alpha      beta       rho\n");
  fprintf(file,"------------------------------------------\n");

  for (int m = 0 ; m < nd_ ; m++) 
  {
    float  alpha, beta, rho;
    bool   hasAlpha, hasBeta, hasRho;
    int    i, j, k;
    data_[m]->GetIJK(i, j, k);
    data_[m]->IsValidObs(hasAlpha, hasBeta, hasRho);
    data_[m]->GetAlphaBetaRho(alpha, beta, rho);  

    if (hasAlpha)
      alpha = exp(alpha);
    else 
      alpha = WELLMISSING;
    if (hasBeta)  
      beta = exp(beta);
    else
      beta = WELLMISSING;
    if (hasRho)
      rho = exp(rho);
    else
      rho = WELLMISSING;

    fprintf(file,"%3d %3d %3d   %8.2f  %8.2f  %8.5f\n",i,j,k,alpha,beta,rho);
  }
  fclose(file);
}

