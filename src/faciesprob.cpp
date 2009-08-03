#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "lib/lib_matr.h"
#include "lib/random.h"
#include "lib/kriging1d.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/welldata.h"
#include "src/faciesprob.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/corr.h"
#include "src/crava.h"
#include "src/simbox.h"
#include "src/blockedlogs.h"
#include "src/model.h"
#include "src/definitions.h"
#include "src/filterwelllogs.h"
#include "src/spatialwellfilter.h"

FaciesProb::FaciesProb(Corr          * correlations, 
                       ModelSettings * modelSettings,
                       int             filegrid, 
                       FFTGrid       * bgAlpha, 
                       FFTGrid       * bgBeta, 
                       FFTGrid       * bgRho, 
                       float           p_undef, 
                       float         * priorFacies,
                       FFTGrid       **priorFaciesCubes)
{
  modelSettings_ = modelSettings;
  fileGrid_      = filegrid;
  nzp_           = bgAlpha->getNzp();
  nz_            = bgAlpha->getNz();
  p_undefined_   = p_undef;
  priorFacies_   = priorFacies;
  priorFaciesCubes_ = priorFaciesCubes;
  

  bgAlpha_       = new FFTGrid(bgAlpha);
  bgBeta_        = new FFTGrid(bgBeta);
  bgRho_         = new FFTGrid(bgRho);

  sigma0_        = correlations->getPriorVar0();
  sigmaPost_     = correlations->getPostVar0();
}

FaciesProb::~FaciesProb()
{
  delete [] priorFacies_;

  for(int i=0;i<nFacies_;i++)
  {
    delete density_[i]; 
    delete faciesProb_[i];
  }
  delete [] density_;
  delete [] faciesProb_;
  delete [] nData_;
  
  delete bgAlpha_;
  delete bgBeta_;
  delete bgRho_;

  delete [] alphafiltered_;
  delete [] betafiltered_;
  delete [] rhofiltered_;
  delete [] alphablock_;
  delete [] betablock_;
  delete [] rhoblock_;
  delete [] facieslog_;

  for(int i=0;i<3;i++)
    delete [] sigmae_[i];
  delete [] sigmae_;
}

float**      
FaciesProb::makeFaciesHistAndSetPriorProb(float* alpha, float* beta, float* rho,int* faciesL)
{
  float **hist=new float*[nFacies_];
  int i,j,k,l;
  for(i=0;i<nFacies_;i++)
  {
    hist[i] = new float[nbins_*nbins_*nbinsr_];
    for(l=0;l<nbinsr_;l++)
    {
      for(k=0;k<nbins_;k++)
      {
        for(j=0;j<nbins_;j++)
          hist[i][j+k*nbins_+l*nbins_*nbins_] = 0.0;
      }
    }
  }
  
  for(i=0;i<nFacies_;i++)
    nData_[i] = 0;

  int facies;
  for(i=0;i<ndata_;i++)
  {
    if(faciesL[i]!=IMISSING)
    {
      j=k=l=0;
      j = int (floor((alpha[i]-alphamin_)/dalpha_));
      if(j<0)
        j = 0;
      if(j>nbins_-1)
        j = nbins_-1;
      k = int (floor((beta[i]-betamin_)/dbeta_));
      if(k<0)
        k = 0;
      if(k>nbins_-1)
        k = nbins_-1;
      l = int (floor((rho[i]-rhomin_)/drho_));
      if(l<0)
        l = 0;
      if(l>nbinsr_-1)
        l = nbinsr_-1;
      
      facies = faciesL[i];
      nData_[facies]++;
      hist[facies][j+k*nbins_+l*nbins_*nbins_]+=1.0;
    }
  }
 
 // float sum = 0.0f;
 // for(i=0;i<nFacies_;i++)
 //   sum+=nData_[i];

//  for(i=0;i<nFacies_;i++)
 //    priorFacies_[i] = float(nData_[i])/sum;

  for(i=0;i<nFacies_;i++)
  {
    double nf = 1.0/double(nData_[i]);
    for(j=0;j<nbins_*nbins_*nbinsr_;j++)
       hist[i][j] = float(double(hist[i][j])*double(nf));
  }
  return hist;
}

void       
FaciesProb::makeFaciesDens(int nfac)
{ 
  int i,j,k,l;
  nFacies_ = nfac;

  //priorFacies_ = new float[nFacies_];
  float kstda, kstdb, kstdr, hopt;
  nbins_  = 100;
  nbinsr_ = 50;

  nobs_   = 0;
  nData_  = new int[nFacies_];
  float *smooth = new float[nbins_*nbins_*nbinsr_];
  for(i=0;i<ndata_;i++)
  {
    if(facieslog_[i]!=IMISSING)
      nobs_++;
  }
 
  // Make bins.
  calculateVariances(alphafiltered_,betafiltered_,rhofiltered_,facieslog_);//sets varAlpha_ etc....
  hopt  = static_cast<float>(pow(4.0/7,1.0/7)*pow(static_cast<double>(nobs_),-1.0/7));
  kstda = hopt*sqrt(varAlpha_);
  kstdb = hopt*sqrt(varBeta_);
  kstdr = hopt*sqrt(varRho_);

 
  if(sigmae_[0][0]<kstda*kstda)
    sigmae_[0][0] = kstda*kstda;
  if(sigmae_[1][1]<kstdb*kstdb)
    sigmae_[1][1] = kstdb*kstdb;
  if(sigmae_[2][2]<kstdr*kstdr)
    sigmae_[2][2] = kstdr*kstdr;

  // invert sigmae_
  double **sigmaeinv = new double *[3];
  for(i=0;i<3;i++)
    sigmaeinv[i] = new double [3];

  for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        if(i==j)
          sigmaeinv[i][j] = 1.0;
        else
          sigmaeinv[i][j] = 0.0;
  lib_matrCholR(3, sigmae_);
  lib_matrAXeqBMatR(3, sigmae_, sigmaeinv, 3);


  getMinMax(alphafiltered_,betafiltered_,rhofiltered_,facieslog_);// sets alpamin_etc....
  alphamin_ -= 5.0f*kstda;
  alphamax_ += 5.0f*kstda;
  betamin_  -= 5.0f*kstdb;
  betamax_  += 5.0f*kstdb;
  rhomin_   -= 5.0f*kstdr;
  rhomax_   += 5.0f*kstdr;

  dalpha_ = (alphamax_-alphamin_)/nbins_;
  dbeta_  = (betamax_-betamin_)/nbins_;
  drho_   = (rhomax_-rhomin_)/nbinsr_;

  float** hist=makeFaciesHistAndSetPriorProb(alphafiltered_,betafiltered_,rhofiltered_,facieslog_);

  int jj,jjj,kk,kkk,ll,lll;
  lll=2;
  
  float sum = 0.0f;
  for(l=0;l<nbinsr_;l++)
  {
    kkk=2;
    if(l<=nbinsr_/2)
      ll = l;
    else
    {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0;k<nbins_;k++)
    {
      jjj=2;
      if(k<=nbins_/2)
        kk=k;
      else
      {
        kk = -(k-kkk);
        kkk+=2;
      }
      for(j=0;j<nbins_;j++)
      {
        if(j<=nbins_/2)
          jj=j;
        else
        {
          jj = -(j-jjj);
          jjj+=2;
        }
       // smooth[j+k*nbins_+l*nbins_*nbins_] = exp(-0.5f*(jj*dalpha_*jj*dalpha_/(kstda*kstda)+kk*dbeta_*kk*dbeta_/(kstdb*kstdb)+ll*drho_*ll*drho_/(kstdr*kstdr)));
        smooth[j+k*nbins_+l*nbins_*nbins_] = float(exp(-0.5f*(jj*dalpha_*jj*dalpha_*sigmaeinv[0][0]+kk*dbeta_*kk*dbeta_*sigmaeinv[1][1]+ll*drho_*ll*drho_*sigmaeinv[2][2]+
                                                        2*jj*dalpha_*kk*dbeta_*sigmaeinv[1][0]+2*jj*dalpha_*ll*drho_*sigmaeinv[2][0]+2*kk*dbeta_*ll*drho_*sigmaeinv[2][1])));
        sum = sum+smooth[j+k*nbins_+l*nbins_*nbins_];
      }
    }
  }
  // normalize smoother
  for(l=0;l<nbinsr_;l++)
    for(k=0;k<nbins_;k++)
      for(j=0;j<nbins_;j++)
        smooth[j+k*nbins_+l*nbins_*nbins_]/=sum;


  density_ = new FFTGrid*[nFacies_];
  FFTGrid *smoother;
  smoother = new FFTGrid(nbins_,nbins_,nbinsr_,nbins_,nbins_,nbinsr_);
  smoother->fillInFromArray(smooth);
 // smoother->writeAsciiRaw("smooth");
  smoother->fftInPlace();

  for(i=0;i<nFacies_;i++)
  {
    density_[i] = new FFTGrid(nbins_,nbins_,nbinsr_,nbins_,nbins_,nbinsr_); //No padding because alphamin/alphamax are far from observed values
    density_[i]->fillInFromArray(hist[i]);
    //char* fN=new char[MAX_STRING];
    //sprintf(fN,"Hist%d",i);
    //density_[i]->writeAsciiRaw(fN);
    //delete [] fN;
    density_[i]->fftInPlace();
    density_[i]->multiply(smoother);
    density_[i]->invFFTInPlace();
    density_[i]->multiplyByScalar(float(sqrt(double(nbins_*nbins_*nbinsr_))));
    if(ModelSettings::getDebugLevel()>=1)
    {
      char* fN1=new char[MAX_STRING];
      sprintf(fN1,"Dens_%d",i);
      density_[i]->writeAsciiRaw(fN1);
      delete [] fN1;
    }
  }

  delete smoother;
  delete [] smooth;
  for (int i = 0 ; i < nFacies_ ; i++)
    delete [] hist[i];
  delete [] hist;
  for(i=0;i<3;i++)
    delete [] sigmaeinv[i];
  delete [] sigmaeinv;
}

float**      
FaciesProb::makeFaciesHistAndSetPriorProb2(float* alpha, float* beta, float* rho,int* faciesL)
{
  double** sigma0Inv=new double*[3];
  double** sigma0=new double*[3];
  int i,j,k,l;

  for(i=0;i<3;i++)
  {
    sigma0[i]=new double[3];
    sigma0Inv[i]=new double[3];
    for(j=0;j<3;j++)
    {
      sigma0[i][j]=double(sigma0_[i][j]);
      sigma0Inv[i][j] = 0.0;
    }
     sigma0Inv[i][i] = 1.0;
  }
  
  sigma0[0][0] = MAXIM( sigma0[0][0],varAlpha_);
  sigma0[1][1] = MAXIM( sigma0[0][0],varBeta_);
  sigma0[2][2] = MAXIM( sigma0[0][0],varRho_);

  
  float **hist=new float*[nFacies_];

  for(i=0;i<nFacies_;i++)
  {
    hist[i] = new float[nbins_*nbins_*nbinsr_];
    for(l=0;l<nbinsr_;l++)
    {
      for(k=0;k<nbins_;k++)
      {
        for(j=0;j<nbins_;j++)
          hist[i][j+k*nbins_+l*nbins_*nbins_] = 0.0;
      }
    }
  }
  
  for(i=0;i<nFacies_;i++)
    nData_[i] = 0;

  lib_matrCholR(3, sigma0);
  lib_matrAXeqBMatR(3, sigma0, sigma0Inv, 3);

  int facies;
  for(i=0;i<ndata_;i++)
  {
    if(faciesL[i]!=IMISSING)
    {
      j=k=l=0;
      j = int (floor((alpha[i]-alphamin_)/dalpha_));
      if(j<0)
        j = 0;
      if(j>nbins_-1)
        j = nbins_-1;
      k = int (floor((beta[i]-betamin_)/dbeta_));
      if(k<0)
        k = 0;
      if(k>nbins_-1)
        k = nbins_-1;
      l = int (floor((rho[i]-rhomin_)/drho_));
      if(l<0)
        l = 0;
      if(l>nbinsr_-1)
        l = nbinsr_-1;
      
      facies = faciesL[i];
      nData_[facies]++;
      
      double * vec = new double[3];
      vec[0]=alphablock_[i]-meanA_;vec[1]=betablock_[i]-meanB_;vec[2]=rhoblock_[i]-meanR_;
      int c1,c2;
      double mahDist =0.0;

      for(c1=0;c1<3;c1++)
        for(c2=0;c2<3;c2++)
          mahDist+=vec[c1]*sigma0Inv[c1][c2]*vec[c2];
      delete [] vec;

      hist[facies][j+k*nbins_+l*nbins_*nbins_]+= float(exp(0.5*mahDist));// add divide by prior......;
    }
  }
  float sum = 0.0f;
  for(i=0;i<nFacies_;i++)
    sum+=nData_[i];

  for(i=0;i<nFacies_;i++)
     priorFacies_[i] = float(nData_[i])/sum;

  // Normalizes hist to be ratios of denseities. not multiplied with ...
  for(i=0;i<nFacies_;i++)
  {
    double nf = 1.0/double(nData_[i]);
    for(j=0;j<nbins_*nbins_*nbinsr_;j++)
       hist[i][j] =  float(double(hist[i][j])*double(nf));
  }
 for(i=0;i<3;i++)
  {
    delete [] sigma0[i];
     delete [] sigma0Inv[i];
  }  
   delete []  sigma0;
   delete []  sigma0Inv;

  return hist;
}

void       
FaciesProb::makeFaciesDens2(int nfac)
{
  int i,j,k,l;
  nFacies_ = nfac;
   priorFacies_ = new float[nFacies_];
  float kstda, kstdb, kstdr, hopt;
  nbins_  = 150;
  nbinsr_ = 100;
  nobs_   = 0;
  nData_  = new int[nFacies_];
  float *smooth = new float[nbins_*nbins_*nbinsr_];
  for(i=0;i<ndata_;i++)
  {
    if(facieslog_[i]!=IMISSING)
      nobs_++;
  }

  // Make bins.
  calculateVariances(alphablock_,betablock_,rhoblock_,facieslog_);

  hopt  = static_cast<float>(pow(4.0/7,1.0/7)*pow(static_cast<double>(nobs_),-1.0/7));
  kstda = hopt*sqrt(sigma0_[0][0]);
  kstdb = hopt*sqrt(sigma0_[1][1]);
  kstdr = hopt*sqrt(sigma0_[2][2]);

  double** sigmaSmooth=new double*[3];
  double** sigmaSmoothInv=new double*[3];
  
  for(i=0;i<3;i++)
  {
    sigmaSmooth[i]=new double[3];
    sigmaSmoothInv[i]=new double[3];
    for(j=0;j<3;j++)
    {
      sigmaSmooth[i][j]= double(sigmaPost_[i][j]);
      sigmaSmoothInv[i][j]=0.0;
    }
     sigmaSmoothInv[i][i] = 1.0;
  }

  if((kstda*kstda)> sigmaSmooth[0][0])
    sigmaSmooth[0][0]=double(kstda*kstda);

  if((kstdb*kstdb)> sigmaSmooth[1][1])
    sigmaSmooth[1][1]=double(kstdb*kstdb);

  if((kstdr*kstdr)> sigmaSmooth[2][2])
    sigmaSmooth[2][2]=double(kstdr*kstdr);
  
  getMinMax(alphablock_,betablock_,rhoblock_,facieslog_);
  alphamin_ -= 4.0f*float(sqrt(sigmaSmooth[0][0]));
  alphamax_ += 4.0f*float(sqrt(sigmaSmooth[0][0]));
  betamin_  -= 4.0f*float(sqrt(sigmaSmooth[1][1]));
  betamax_  += 4.0f*float(sqrt(sigmaSmooth[1][1]));
  rhomin_   -= 4.0f*float(sqrt(sigmaSmooth[2][2]));
  rhomax_   += 4.0f*float(sqrt(sigmaSmooth[2][2]));

  dalpha_ = (alphamax_-alphamin_)/nbins_;
  dbeta_  = (betamax_-betamin_)/nbins_;
  drho_   = (rhomax_-rhomin_)/nbinsr_;

  float**hist= makeFaciesHistAndSetPriorProb2(alphablock_,betablock_,rhoblock_,facieslog_);

  lib_matrCholR(3, sigmaSmooth); // NB "destroys" sigmaSmooth
  lib_matrAXeqBMatR(3, sigmaSmooth, sigmaSmoothInv, 3);

  int jj,jjj,kk,kkk,ll,lll;
  lll=2;
  float sum = 0;
  for(l=0;l<nbinsr_;l++)
  {
    kkk=2;
    if(l<=nbinsr_/2)
      ll = l;
    else
    {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0;k<nbins_;k++)
    {
      jjj=2;
      if(k<=nbins_/2)
        kk=k;
      else
      {
        kk = -(k-kkk);
        kkk+=2;
      }
      for(j=0;j<nbins_;j++)
      {
        if(j<=nbins_/2)
          jj=j;
        else
        {
          jj = -(j-jjj);
          jjj+=2;
        }
        double* vec =new double[3];
        vec[0]=jj*dalpha_;vec[1]=kk*dbeta_;vec[2]=ll*drho_;
        int c1,c2;
        double mahDist =0.0;
        for(c1=0;c1<3;c1++)
          for(c2=0;c2<3;c2++)
            mahDist+=vec[c1]*sigmaSmoothInv[c1][c2]*vec[c2];
        delete [] vec;
        smooth[j+k*nbins_+l*nbins_*nbins_] = float(exp(-0.5*mahDist));
        sum = sum+smooth[j+k*nbins_+l*nbins_*nbins_];
      }
    }
  }
  // normalize smoother
  for(l=0;l<nbinsr_;l++)
    for(k=0;k<nbins_;k++)
      for(j=0;j<nbins_;j++)
        smooth[j+k*nbins_+l*nbins_*nbins_]/=sum;

  density_ = new FFTGrid*[nFacies_];
  FFTGrid *smoother;
  smoother = new FFTGrid(nbins_,nbins_,nbinsr_,nbins_,nbins_,nbinsr_);
  smoother->fillInFromArray(smooth);
  //smoother->writeAsciiRaw("smooth_PC");
  smoother->fftInPlace();

  for(i=0;i<nFacies_;i++)
  {
    density_[i] = new FFTGrid(nbins_,nbins_,nbinsr_,nbins_,nbins_,nbinsr_); //No padding because alphamin/alphamax are far from observed values
    density_[i]->fillInFromArray(hist[i]);
    //char* fN=new char[MAX_STRING];
    //sprintf(fN,"Hist_PC%d",i);
    //density_[i]->writeAsciiRaw(fN);
    density_[i]->fftInPlace();
    density_[i]->multiply(smoother);
    density_[i]->invFFTInPlace();
    density_[i]->multiplyByScalar(float(sqrt(double(nbins_*nbins_*nbinsr_))));
   //sprintf(fN,"Dens_PC%d",i);
    //density_[i]->writeAsciiRaw(fN);
    //delete [] fN;
  }

  delete smoother;
  delete [] smooth;
  delete [] hist;

  for(i=0;i<3;i++)
  {
    delete [] sigmaSmooth[i];
     delete [] sigmaSmoothInv[i];
  }  
   delete [] sigmaSmooth;
   delete [] sigmaSmoothInv;
}


void FaciesProb::makeFaciesProb(int nfac, FFTGrid *postAlpha, FFTGrid *postBeta, FFTGrid *postRho)
{
  makeFaciesDens(nfac);
  if(priorFaciesCubes_!=NULL)
    normalizeCubes(priorFaciesCubes_);
  calculateFaciesProb(postAlpha, postBeta, postRho);
}

float FaciesProb::findDensity(float alpha, float beta, float rho, int facies)
{
  int j1,k1,l1;
  int j2,k2,l2;
  float wj,wk,wl;
  j1=k1=l1=j2=k2=l2=0;
  j1 =  int (floor((alpha-alphamin_+dalpha_*0.5)/dalpha_));
  if(j1<0)
    j1 = 0;
  if(j1>nbins_-1)
    j1 = nbins_-1;
  j2=j1+1;
  if(j2>nbins_-1)
    j2 = nbins_-1; 
  wj= MINIM(1,MAXIM(0,(alpha-alphamin_+dalpha_*0.5f-j1*dalpha_)/dalpha_));

  k1  = int (floor((beta-betamin_+dbeta_*0.5f)/dbeta_));
  if(k1<0)
    k1 = 0;
  if(k1>nbins_-1)
    k1 = nbins_-1;
  k2=k1+1;
  if(k2>nbins_-1)
    k2 = nbins_-1; 
  wk= MINIM(1,MAXIM(0,(beta-betamin_+dbeta_*0.5f-k1*dbeta_)/dbeta_));

  l1 = int (floor((rho-rhomin_+drho_*0.5f)/drho_));
  if(l1<0)
    l1 = 0;
  if(l1>nbinsr_-1)
    l1 = nbinsr_-1;
  l2=l1+1;
  if(l2>nbinsr_-1)
    l2 = nbinsr_-1; 
  wl= MINIM(1,MAXIM(0,(rho-rhomin_+drho_*0.5f-l1*drho_)/drho_));

  float value1 = MAXIM(0,density_[facies]->getRealValue(j1,k1,l1));
  float value2 = MAXIM(0,density_[facies]->getRealValue(j1,k1,l2));
  float value3 = MAXIM(0,density_[facies]->getRealValue(j1,k2,l1));
  float value4 = MAXIM(0,density_[facies]->getRealValue(j1,k2,l2));
  float value5 = MAXIM(0,density_[facies]->getRealValue(j2,k1,l1));
  float value6 = MAXIM(0,density_[facies]->getRealValue(j2,k1,l2));
  float value7 = MAXIM(0,density_[facies]->getRealValue(j2,k2,l1));
  float value8 = MAXIM(0,density_[facies]->getRealValue(j2,k2,l2));
  
  float value=0;
  value += (1.0f-wj)*(1.0f-wk)*(1.0f-wl)*value1;
  value += (1.0f-wj)*(1.0f-wk)*(     wl)*value2;
  value += (1.0f-wj)*(     wk)*(1.0f-wl)*value3;
  value += (1.0f-wj)*(     wk)*(     wl)*value4;
  value += (     wj)*(1.0f-wk)*(1.0f-wl)*value5;
  value += (     wj)*(1.0f-wk)*(     wl)*value6;
  value += (     wj)*(     wk)*(1.0f-wl)*value7;
  value += (     wj)*(     wk)*(     wl)*value8;

  if (value > 0.0)
    return(value);
  else
    return 0.0;
}

void FaciesProb::calculateConditionalFaciesProb(WellData **wells, int nWells)
{
  //
  // Get the needed blocked logs
  //
  BlockedLogs ** bw = new BlockedLogs * [nWells];  
  int totBlocks = 0;
  int count = 0;

  for (int i = 0 ; i < nWells ; i++)
  {
    if(wells[i]->getUseForFaciesProbabilities())
    { 
      bw[count] = wells[i]->getBlockedLogsConstThick();
      totBlocks += bw[count]->getNumberOfBlocks();
      count++;
    }
  }
  int nonDeviatedWells = count;

  //
  // Put all blocked facies logs in one vector
  //
  int ** BWfacies = new int * [nonDeviatedWells]; 
  for (int i = 0 ; i < nonDeviatedWells ; i++)
  {
    const int   nBlocks    = bw[i]->getNumberOfBlocks();
    const int * BWfacies_i = bw[i]->getFacies();
    BWfacies[i] = new int[nBlocks]; 
    for (int b = 0 ; b < nBlocks ; b++)
    {
      BWfacies[i][b] = BWfacies_i[b];
      //LogKit::LogFormatted(LogKit::LOW,"ib, BWfacies[i][b] = %d %d\n",ib,BWfacies[ib]);
    }
  }
   
  //
  // Block facies probabilities and put them in one vector
  //
  float *** BWfaciesProb = new float ** [nFacies_];
  for(int f = 0 ; f < nFacies_ ; f++)
  {
    BWfaciesProb[f] = new float * [nonDeviatedWells];
    for (int i = 0 ; i < nonDeviatedWells ; i++)
    {
      const int   nBlocks = bw[i]->getNumberOfBlocks();
      const int * ipos    = bw[i]->getIpos();
      const int * jpos    = bw[i]->getJpos();
      const int * kpos    = bw[i]->getKpos();
      BWfaciesProb[f][i] = new float[nBlocks];
      for(int b = 0 ; b < nBlocks ; b++)
      {
        BWfaciesProb[f][i][b] = faciesProb_[f]->getRealValue(ipos[b],jpos[b],kpos[b]);
        //LogKit::LogFormatted(LogKit::LOW,"f,ib, BWfaciesProb[f][i][b] = %d %d %.5f\n",f,ib,BWfaciesProb[f][ib]);
      }
    }
  }

  //
  // Sum probabilities for a given facies
  //
  float *** sumProb = new float ** [nFacies_];
  int   *** numProb = new int ** [nFacies_];
  for(int f1 = 0 ; f1 < nFacies_ ; f1++)
  {
    sumProb[f1] = new float * [nFacies_];
    numProb[f1] = new int * [nFacies_];
    for(int f2 = 0 ; f2 < nFacies_ ; f2++)
    {
      sumProb[f1][f2] = new float[nonDeviatedWells];
      numProb[f1][f2] = new int[nonDeviatedWells];
      for (int i = 0 ; i < nonDeviatedWells ; i++)
      {
        sumProb[f1][f2][i] = 0.0f;
        numProb[f1][f2][i] = 0;
        for(int b = 0 ; b < bw[i]->getNumberOfBlocks() ; b++)
        {
          if (BWfacies[i][b] == f1) {
            //LogKit::LogFormatted(LogKit::LOW,"f1,f2 = %d,%d   ib = %d    BWfacies[i][b] = %d   sumProb = %.5f\n",f1,f2,ib,BWfacies[i][b],sumProb);
            sumProb[f1][f2][i] += BWfaciesProb[f2][i][b];
            numProb[f1][f2][i] += 1;
          }
        }
      }
    }
  }
  for (int i = 0 ; i < nonDeviatedWells ; i++)
    delete [] BWfacies[i];
  delete [] BWfacies;
   
  for(int f = 0 ; f < nFacies_ ; f++) {
    for (int i = 0 ; i < nonDeviatedWells ; i++)
      delete [] BWfaciesProb[f][i];
    delete [] BWfaciesProb[f];
  }
  delete [] BWfaciesProb;


  float ** condFaciesProb  = new float * [nFacies_];
  for(int f = 0 ; f < nFacies_ ; f++)
    condFaciesProb[f] = new float[nFacies_];

  float * totProb = new float[nFacies_];

  //
  // Estimate P( facies2 | facies1 ) for each well
  //
  for (int i = 0 ; i < nonDeviatedWells ; i++)
  {
    for(int f1 = 0 ; f1 < nFacies_ ; f1++)
    {
      totProb[f1] = 0.0f;
      for(int f2 = 0 ; f2 < nFacies_ ; f2++)
      {
        if (numProb[f1][f2][i] > 0) 
          condFaciesProb[f1][f2] = sumProb[f1][f2][i]/numProb[f1][f2][i];
        else
          condFaciesProb[f1][f2] = 0.0f;
        totProb[f1] += condFaciesProb[f1][f2];
      }
    }
    LogKit::LogFormatted(LogKit::LOW,"\nWell: %s\n",bw[i]->getWellname());
    LogKit::LogFormatted(LogKit::LOW,"\nFacies      |");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings_->getFaciesName(f).c_str());
    LogKit::LogFormatted(LogKit::LOW,"\n------------+");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::LOW,"------------",f);
    LogKit::LogFormatted(LogKit::LOW,"\n");
    for(int f1=0 ; f1 < nFacies_ ; f1++)
    {
      LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings_->getFaciesName(f1).c_str());
      for(int f2=0 ; f2 < nFacies_ ; f2++)
      {
        LogKit::LogFormatted(LogKit::LOW," %11.3f",condFaciesProb[f2][f1]);
      }
      LogKit::LogFormatted(LogKit::LOW,"\n");
    }
    LogKit::LogFormatted(LogKit::LOW,"Total Prob  |");
    for(int f=0 ; f < nFacies_ ; f++)
    {
      LogKit::LogFormatted(LogKit::LOW," %11.3f",totProb[f]);
    }
    LogKit::LogFormatted(LogKit::LOW,"\n");
  }

  //
  // Estimate P( facies2 | facies1 )
  //
  LogKit::LogFormatted(LogKit::HIGH,"\nThe table below gives the mean conditional probability of finding one of");
  LogKit::LogFormatted(LogKit::HIGH,"\nthe facies specified in the left column when one of the facies specified");
  LogKit::LogFormatted(LogKit::HIGH,"\nin the top row are observed in well logs ==> P(A|B)\n");
  for(int f1 = 0 ; f1 < nFacies_ ; f1++)
  {
    totProb[f1] = 0.0f;
    for(int f2 = 0 ; f2 < nFacies_ ; f2++)
    {
      float sumFaciesProb = 0.0f;
      int   numFaciesProb =   0;
      for (int i = 0 ; i < nonDeviatedWells ; i++)
      {
        if (numProb[f1][f2][i] > 0)
        {
          sumFaciesProb += sumProb[f1][f2][i];
          numFaciesProb += numProb[f1][f2][i];
        }
      }
      if (numFaciesProb > 0)
        condFaciesProb[f1][f2] = sumFaciesProb/numFaciesProb;
      else
        condFaciesProb[f1][f2] = 0.0f;
      totProb[f1] += condFaciesProb[f1][f2];
    }
  }

  for(int f1 = 0 ; f1 < nFacies_ ; f1++) {
    for(int f2 = 0 ; f2 < nFacies_ ; f2++) {
      delete [] sumProb[f1][f2];
      delete [] numProb[f1][f2];
    }
    delete [] sumProb[f1];
    delete [] numProb[f1];
  }
  delete [] sumProb;
  delete [] numProb;

  LogKit::LogFormatted(LogKit::LOW,"\nFor all wells:\n");
  LogKit::LogFormatted(LogKit::LOW,"\nFacies      |");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings_->getFaciesName(f).c_str());
  LogKit::LogFormatted(LogKit::LOW,"\n------------+");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::LOW,"------------",f);
  LogKit::LogFormatted(LogKit::LOW,"\n");
  for(int f1=0 ; f1 < nFacies_ ; f1++)
  {
    LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings_->getFaciesName(f1).c_str());
    for(int f2=0 ; f2 < nFacies_ ; f2++)
      LogKit::LogFormatted(LogKit::LOW," %11.3f",condFaciesProb[f2][f1]);
    LogKit::LogFormatted(LogKit::LOW,"\n");
  }
  LogKit::LogFormatted(LogKit::LOW,"Total Prob  |");
  for(int f=0 ; f < nFacies_ ; f++)
  {
    LogKit::LogFormatted(LogKit::LOW," %11.3f",totProb[f]);
  }
  LogKit::LogFormatted(LogKit::LOW,"\n");

  delete [] totProb;

  for(int i=0 ; i < nFacies_ ; i++)
    delete [] condFaciesProb[i];
  delete [] condFaciesProb;
  // We cannot delete blocked logs here, only the array of blocked logs.
  delete [] bw;
}

void FaciesProb::calculateFaciesProb(FFTGrid *alphagrid, FFTGrid *betagrid, FFTGrid *rhogrid)
{
  float * value = new float[nFacies_];
  int i,j,k,l;
  int nx, ny, nz, rnxp, nyp, nzp, smallrnxp;
  float alpha, beta, rho, sum;
  rnxp = alphagrid->getRNxp();
  nyp  = alphagrid->getNyp();
  nzp  = alphagrid->getNzp();
  nx   = alphagrid->getNx();
  ny   = alphagrid->getNy();
  nz   = alphagrid->getNz();
  faciesProb_ = new FFTGrid*[nFacies_]; 
  alphagrid->setAccessMode(FFTGrid::READ);
  betagrid->setAccessMode(FFTGrid::READ);
  rhogrid->setAccessMode(FFTGrid::READ);
  for(i=0;i<nFacies_;i++)
  {
    if(fileGrid_==1)
      faciesProb_[i] = new FFTFileGrid(nx, ny, nz, nx, ny, nz);
    else
      faciesProb_[i] = new FFTGrid(nx, ny, nz, nx, ny, nz);
    faciesProb_[i]->setAccessMode(FFTGrid::WRITE);
    faciesProb_[i]->createRealGrid();
  }
  smallrnxp = faciesProb_[0]->getRNxp();
  float help;
  for(i=0;i<nzp;i++)
  {
    for(j=0;j<nyp;j++)
      for(k=0;k<rnxp;k++)
      {
        // printf("%d %d %d\n",i,j,k);
        alpha = alphagrid->getNextReal();
        beta = betagrid->getNextReal();
        rho = rhogrid->getNextReal();
        // alpha = alphagrid->getRealValue(k,j,i);
        // beta = betagrid->getRealValue(k,j,i);
        // rho = rhogrid->getRealValue(k,j,i);
        if(k<smallrnxp && j<ny && i<nz)
        {
          
          sum = p_undefined_/(nbins_*nbins_*nbinsr_);
          for(l=0;l<nFacies_;l++)
          {
            if(priorFaciesCubes_!=NULL)
              value[l] = priorFaciesCubes_[l]->getRealValue(k,j,i)*findDensity(alpha, beta, rho, l);
            else
              value[l] = priorFacies_[l]*findDensity(alpha, beta, rho, l); 
            sum = sum+value[l];
          }
          for(l=0;l<nFacies_;l++)
          {
            help = value[l]/sum;
            //faciesProb_[l]->setRealValue(k,j,i,help);
            if(k<nx)
              faciesProb_[l]->setNextReal(help);
            else
              faciesProb_[l]->setNextReal(RMISSING);
          }
        }
      }
  }
  delete [] value;
}


void FaciesProb::getMinMax(float* alpha,float* beta,float* rho,int* facies)
{
  int i;
  alphamin_ = 1000.0;
  alphamax_ = -1000.0;
  betamin_ = 1000.0;
  betamax_ = -1000.0;
  rhomin_ = 1000.0;
  rhomax_ = -1000.0;
  for(i=0;i<ndata_;i++)
  {
    if(alpha[i]<alphamin_ && facies[i]!=IMISSING)
      alphamin_ = alpha[i];
    if(alpha[i]>alphamax_ && facies[i] != IMISSING)
      alphamax_ = alpha[i];
    if(beta[i]<betamin_&& facies[i] != IMISSING)
      betamin_ = beta[i];
    if(beta[i]>betamax_&& facies[i] != IMISSING)
      betamax_ = beta[i];
    if(rhofiltered_[i]<rhomin_&& facies[i] != IMISSING)
      rhomin_ = rho[i];
    if(rho[i]>rhomax_&& facies[i] != IMISSING)
      rhomax_ = rho[i];
  }
}

void FaciesProb::calculateVariances(float* alpha,float* beta,float* rho,int* facies)
{
  int i;
  //bool validA, validB, validR;
  int nA, nB, nR;
  nA = nB = nR = 0;
  meanA_ = meanB_ = meanR_ = 0.0;
  for(i=0;i<ndata_;i++)
  {
    if(facies[i]!=IMISSING)
    {
      meanA_+=alpha[i];
      nA++;
      meanB_+=beta[i];
      nB++;
      meanR_+=rho[i];
    nR++;
    }
  }
  meanA_/=nA;
  meanB_/=nB;
  meanR_/=nR;
  varAlpha_ = varBeta_ = varRho_ = 0.0;

  for(i=0;i<ndata_;i++)
  {
    if(facies[i]!=IMISSING)
    {
      varAlpha_+=pow(alpha[i]-meanA_,2);
      varBeta_+=pow(beta[i]-meanB_,2);
      varRho_+=pow(rho[i]-meanR_,2);
    }
  }
  varAlpha_ /= nA-1;
  varBeta_ /= nB-1;
  varRho_ /= nR-1;
}


void FaciesProb::setNeededLogs(FilterWellLogs * filteredLogs,
                               WellData      ** wells,
                               int              nWells,
                               int              nz,
                               RandomGen      * random)
{
  float ** vtAlphaFiltered = filteredLogs->getVtAlphaFiltered(); 
  float ** vtBetaFiltered  = filteredLogs->getVtBetaFiltered(); 
  float ** vtRhoFiltered   = filteredLogs->getVtRhoFiltered(); 

  float ** vtAlpha         = filteredLogs->getVtAlpha(); 
  float ** vtBeta          = filteredLogs->getVtBeta(); 
  float ** vtRho           = filteredLogs->getVtRho(); 

  int    * vtFacies        = new int[nz];

  ndata_ = nWells*nz;

  alphafiltered_ = new float[ndata_];
  betafiltered_  = new float[ndata_];
  rhofiltered_   = new float[ndata_];
  alphablock_    = new float[ndata_];
  betablock_     = new float[ndata_];
  rhoblock_      = new float[ndata_];
  facieslog_     = new int[ndata_];

  for(int i=0;i<ndata_;i++)
  {
    alphafiltered_[i] = RMISSING;
    betafiltered_[i]  = RMISSING;
    rhofiltered_[i]   = RMISSING;
    alphablock_[i]    = RMISSING;
    betablock_[i]     = RMISSING;
    rhoblock_[i]      = RMISSING;
    facieslog_[i]     = IMISSING;
  }

  for (int w=0 ; w<nWells ; w++)
  {
    if(wells[w]->getUseForFaciesProbabilities())
    { 
      BlockedLogs * bw = wells[w]->getBlockedLogsConstThick();
      bw->getVerticalTrend(bw->getFacies(),vtFacies,random);

      //
      // Set facies log MISSING unless all parameters are defined
      //
      for(int i=0 ; i<nz ; i++)
        if(vtAlpha[w][i] == RMISSING || vtBeta[w][i] == RMISSING || vtRho[w][i] == RMISSING)
          vtFacies[i] = IMISSING;

      //
      // Fill block
      //     
      for(int i=0 ; i<nz ; i++)
      {
        alphafiltered_[i+w*nz] = vtAlphaFiltered[w][i];
        betafiltered_[i+w*nz]  = vtBetaFiltered[w][i];
        rhofiltered_[i+w*nz]   = vtRhoFiltered[w][i];
        alphablock_[i+w*nz]    = vtAlpha[w][i];
        betablock_[i+w*nz]     = vtBeta[w][i];
        rhoblock_[i+w*nz]      = vtRho[w][i];  
        facieslog_[i+w*nz]     = vtFacies[i];
      }
    }
  } 
  delete [] vtFacies;
}
void FaciesProb::setNeededLogsSpatial(SpatialWellFilter * filteredLogs,
                               WellData      ** wells,
                               int              nWells)
{
  int ndata = filteredLogs->getNdata();
  float *alphafiltered = new float[ndata];
  float *betafiltered  = new float[ndata];
  float *rhofiltered   = new float[ndata];
  float *alphablock    = new float[ndata];
  float *betablock     = new float[ndata];
  float *rhoblock      = new float[ndata];
  int *facieslog    = new int[ndata];

  for(int i=0;i<ndata;i++)
  {
    alphafiltered[i] = RMISSING;
    betafiltered[i]  = RMISSING;
    rhofiltered[i]   = RMISSING;
    alphablock[i]    = RMISSING;
    betablock[i]     = RMISSING;
    rhoblock[i]      = RMISSING;
    facieslog[i]     = IMISSING;
  }
 int n, lastn;
 lastn = 0;
  for (int w=0 ; w<nWells ; w++)
  {
    if(wells[w]->getUseForFaciesProbabilities())
    { 
      BlockedLogs * bw = wells[w]->getBlockedLogsOrigThick();
      n = bw->getNumberOfBlocks();
      for(int i=0;i<n;i++)
      {
        alphablock[i+lastn] = bw->getAlpha()[i];
        betablock[i+lastn] = bw->getBeta()[i];
        rhoblock[i+lastn] = bw->getRho()[i];
        if(alphablock[i+lastn]==RMISSING || betablock[i+lastn]==RMISSING || rhoblock[i+lastn]==RMISSING)
          facieslog[i+lastn] = IMISSING;
        else
          facieslog[i+lastn] = bw->getFacies()[i];
      }
      lastn += n;
    }
  }
  ndata_ = lastn;
  alphafiltered_ = new float[ndata_];
  betafiltered_  = new float[ndata_];
  rhofiltered_   = new float[ndata_];
  alphablock_    = new float[ndata_];
  betablock_     = new float[ndata_];
  rhoblock_      = new float[ndata_];
  facieslog_     = new int[ndata_];

  for(int i=0;i<ndata_;i++)
  {
    alphafiltered_[i] = RMISSING;
    betafiltered_[i]  = RMISSING;
    rhofiltered_[i]   = RMISSING;
    alphablock_[i]    = RMISSING;
    betablock_[i]     = RMISSING;
    rhoblock_[i]      = RMISSING;
    facieslog_[i]     = IMISSING;
  }

  for(int i=0;i<ndata_;i++)
  {
    alphablock_[i] = alphablock[i];
    betablock_[i] = betablock[i];
    rhoblock_[i] = rhoblock[i];
    facieslog_[i] = facieslog[i];
    alphafiltered_[i] = filteredLogs->getAlphaFiltered()[i];
    betafiltered_[i] = filteredLogs->getBetaFiltered()[i];
    rhofiltered_[i] = filteredLogs->getRhoFiltered()[i];
  }
  delete [] alphafiltered;
  delete [] betafiltered;
  delete [] rhofiltered;
  delete [] alphablock;
  delete [] betablock;
  delete [] rhoblock;
  delete [] facieslog;

  // For use in smoother
  sigmae_ = new double *[3];
  for(int i=0;i<3;i++)
    sigmae_[i] = new double[3];

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      sigmae_[i][j] = filteredLogs->getSigmae()[i][j];

}

void FaciesProb::normalizeCubes(FFTGrid **priorFaciesCubes)
{
  int i,j,k,l;
  int rnxp = priorFaciesCubes[0]->getRNxp();
  int nyp  = priorFaciesCubes[0]->getNyp();
  int nzp  = priorFaciesCubes[0]->getNzp();
  float sum, value;
  int negative = 0;
  for(i=0;i<nzp;i++)
  {
    for(j=0;j<nyp;j++)
      for(k=0;k<rnxp;k++)
      {
        sum = 0.0;
        for(l=0;l<nFacies_;l++)
        {
          value = priorFaciesCubes[l]->getNextReal();
          if(value<0.0)
          {
            value = 0.0;
            if(negative==0)
            {
              LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: Alt least one negative prior facies probability detected. The value is set to 0.0.\n");
              negative = 1;
            }
          }
          sum+= value;
        }
        if(sum>0.0)
        {
          for(l=0;l<nFacies_;l++)
          { 
            value =priorFaciesCubes[l]->getRealValue(k,j,i)/sum; 
            priorFaciesCubes[l]->setRealValue(k,j,i,value);
          }
        }
      }
  }
}
