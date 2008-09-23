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
#include "src/crava.h"
#include "src/simbox.h"
#include "src/blockedlogs.h"
#include "src/model.h"
#include "src/definitions.h"

FaciesProb::FaciesProb(ModelSettings * modelSettings,
                       int filegrid, const float ** sigma0, float *corrprior, Simbox *simbox, const Simbox& osimbox, int nzp, int nz, 
                       FFTGrid* bgAlpha, FFTGrid* bgBeta, FFTGrid* bgRho, RandomGen *random, float p_undef, float *priorFacies)
  : origSimbox_(osimbox)
{
  modelSettings_ = modelSettings;
  fileGrid_      = filegrid;
  simbox_        = simbox;
  nzp_           = nzp;
  nz_            = nz;
  random_        = random;
  p_undefined_   = p_undef;
  priorFacies_ = priorFacies;

  bgAlpha_ = new FFTGrid(bgAlpha);
  bgBeta_  = new FFTGrid(bgBeta);
  bgRho_   = new FFTGrid(bgRho);

  int i, j;
  sigma0_ = new float*[3];
  sigmaPost_= new float*[3];
  for(i=0;i<3;i++)
  {
    sigma0_[i] = new float[3];
    sigmaPost_[i]=new float[3];
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      sigma0_[i][j] = sigma0[i][j];
      sigmaPost_[i][j] = RMISSING;
    }

  int rnzp = 2*(nzp_/2+1);
  corrprior_ = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
  for(i=0;i<nzp_;i++)
    corrprior_[i] =corrprior[i];
}

FaciesProb::~FaciesProb()
{
  for(int i=0;i<3;i++)
  {
    delete [] sigma0_[i];
    delete [] sigmaPost_[i];
  }
  delete [] sigma0_;
  delete [] sigmaPost_;

  delete [] priorFacies_;
  for(int i=0;i<nFacies_;i++)
  {
    delete density_[i]; 
    delete faciesProb_[i];
  }
  delete [] density_;
  delete [] faciesProb_;
  delete [] nData_;
  delete [] alphafiltered_;
  delete [] betafiltered_;
  delete [] rhofiltered_;
  delete [] alphablock_;
  delete [] betablock_;
  delete [] rhoblock_;
  delete [] facieslog_;

  delete bgAlpha_;
  delete bgBeta_;
  delete bgRho_;

  fftw_free(corrprior_);

  delete simbox_;
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
        smooth[j+k*nbins_+l*nbins_*nbins_] = exp(-0.5f*(jj*dalpha_*jj*dalpha_/(kstda*kstda)+kk*dbeta_*kk*dbeta_/(kstdb*kstdb)+ll*drho_*ll*drho_/(kstdr*kstdr)));
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
      sigma0Inv[i][j]=0.0;
    }
     sigma0Inv[i][i]      = 1.0;
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
  // Block wells
  //
  // NBNB-PAL: De blokkede loggene burde kanskje ha vært mellomlagret ettersom de benyttes to ganger...
  //
  BlockedLogs ** bw = new BlockedLogs * [nWells];  
  int totBlocks = 0;
  int count = 0;
  for (int i = 0 ; i < nWells ; i++)
  {
    if(wells[i]->getUseForFaciesProbabilities())
    { 
      bw[count] = new BlockedLogs(wells[i], simbox_, random_) ;
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
      LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings_->getFaciesName(f));
    LogKit::LogFormatted(LogKit::LOW,"\n------------+");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::LOW,"------------",f);
    LogKit::LogFormatted(LogKit::LOW,"\n");
    for(int f1=0 ; f1 < nFacies_ ; f1++)
    {
      LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings_->getFaciesName(f1));
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
    LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings_->getFaciesName(f));
  LogKit::LogFormatted(LogKit::LOW,"\n------------+");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::LOW,"------------",f);
  LogKit::LogFormatted(LogKit::LOW,"\n");
  for(int f1=0 ; f1 < nFacies_ ; f1++)
  {
    LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings_->getFaciesName(f1));
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

  for(int i=0 ; i < nFacies_ ; i++)
    delete [] condFaciesProb[i];
  delete condFaciesProb;

  for (int i = 0 ; i < nonDeviatedWells ; i++)
    delete bw[i];
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
        //     printf("%d %d %d\n",i,j,k);
        alpha = alphagrid->getNextReal();
        beta = betagrid->getNextReal();
        rho = rhogrid->getNextReal();
        // alpha = alphagrid->getRealValue(k,j,i);
        //  beta = betagrid->getRealValue(k,j,i);
        //   rho = rhogrid->getRealValue(k,j,i);
        if(k<smallrnxp && j<ny && i<nz)
        {
          
          sum = p_undefined_/(nbins_*nbins_*nbinsr_);
          for(l=0;l<nFacies_;l++)
          {
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

void FaciesProb::filterWellLogs(WellData **wells, int nWells,
                                fftw_real *postcova, fftw_real *postcovb, fftw_real *postcovr,
                                fftw_real *postcrab, fftw_real *postcrar, fftw_real *postcrbr, 
                                float lowCut, float highCut, int relative)
{
  LogKit::LogFormatted(LogKit::LOW,"\nFiltering well logs\n");
  LogKit::LogFormatted(LogKit::MEDIUM,"  Wells available:\n");
  for (int w = 0 ; w < nWells ; w++)
    if(wells[w]->getUseForFaciesProbabilities())
      LogKit::LogFormatted(LogKit::MEDIUM,"%s\n",wells[w]->getWellname());
  
  float domega  = static_cast<float> (1000.0/(nzp_*simbox_->getdz()));  //dz in milliseconds
  int w, w1,i,j;
  ndata_ = nWells*nz_;
  alphafiltered_ = new float[ndata_];
  betafiltered_ = new float[ndata_];
  rhofiltered_ = new float[ndata_];
  alphablock_ = new float[ndata_];
  betablock_ = new float[ndata_];
  rhoblock_ = new float[ndata_];
  facieslog_ = new int[ndata_];

  for(i=0;i<ndata_;i++)
  {
    alphafiltered_[i] = RMISSING;
    betafiltered_[i] = RMISSING;
    rhofiltered_[i] = RMISSING;
    alphablock_[i] = RMISSING;
    betablock_[i] = RMISSING;
    rhoblock_[i] = RMISSING;
    facieslog_[i] = IMISSING;
  }
  
  sigmaPost_[0][0]=postcova[0];
  sigmaPost_[1][1]=postcovb[0];
  sigmaPost_[2][2]=postcovr[0];
  sigmaPost_[0][1]=postcrab[0];
  sigmaPost_[1][0]=postcrab[0];
  sigmaPost_[0][2]=postcrar[0];
  sigmaPost_[2][0]=postcrar[0];
  sigmaPost_[1][2]=postcrbr[0];
  sigmaPost_[2][1]=postcrbr[0];


  // Do fourier transform of covariances
  rfftwnd_plan p1 = rfftwnd_create_plan(1, &nzp_, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  int rnzp = 2*(nzp_/2+1);
  fftw_complex *postcova_cAmp  = reinterpret_cast<fftw_complex*>(postcova);
  fftw_complex *postcovb_cAmp  = reinterpret_cast<fftw_complex*>(postcovb);
  fftw_complex *postcovr_cAmp  = reinterpret_cast<fftw_complex*>(postcovr);
  fftw_complex *postcovar_cAmp = reinterpret_cast<fftw_complex*>(postcrar);
  fftw_complex *postcovab_cAmp = reinterpret_cast<fftw_complex*>(postcrab);
  fftw_complex *postcovbr_cAmp = reinterpret_cast<fftw_complex*>(postcrbr);
  fftw_complex *corrprior_cAmp = reinterpret_cast<fftw_complex*>(corrprior_);

  rfftwnd_one_real_to_complex(p1,corrprior_, corrprior_cAmp);
  rfftwnd_one_real_to_complex(p1,postcova, postcova_cAmp);
  rfftwnd_one_real_to_complex(p1,postcovb, postcovb_cAmp);
  rfftwnd_one_real_to_complex(p1,postcovr, postcovr_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrab, postcovab_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrar, postcovar_cAmp);
  rfftwnd_one_real_to_complex(p1,postcrbr, postcovbr_cAmp);

  fftw_real    * alpha_rAmp = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
  fftw_real    * beta_rAmp  = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
  fftw_real    * rho_rAmp   = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp));
 
  fftw_complex * alpha_cAmp = reinterpret_cast<fftw_complex*>(alpha_rAmp);
  fftw_complex * beta_cAmp  = reinterpret_cast<fftw_complex*>(beta_rAmp);
  fftw_complex * rho_cAmp   = reinterpret_cast<fftw_complex*>(rho_rAmp);

  int insideOrigSimbox;
  double x,y,z;

  float * vtAlpha = new float[nz_];            // vt = vertical trend
  float * vtBeta  = new float[nz_];
  float * vtRho   = new float[nz_];

  float * vtAlphaBg = new float[nz_];          // vt = vertical trend
  float * vtBetaBg  = new float[nz_];
  float * vtRhoBg   = new float[nz_];
  int   * vtFacies  = new int[nz_];

  char * fileName     = new char[MAX_STRING];
  char * fullFileName = new char[MAX_STRING];

  sprintf(fileName,"BW_filteredlogs.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file = fopen(fullFileName,"w");

  sprintf(fileName,"BW_originallogs.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file2 = fopen(fullFileName,"w");

  sprintf(fileName,"BW_background.dat");
  fullFileName = ModelSettings::makeFullFileName(fileName);
  FILE *file3 = fopen(fullFileName,"w");

  delete [] fileName;
  delete [] fullFileName;

  for (w1 = 0 ; w1 < nWells ; w1++)
  {
    if(wells[w1]->getUseForFaciesProbabilities())
    { 
      BlockedLogs * bw = new BlockedLogs(wells[w1], simbox_, random_) ;
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bw->getVerticalTrend(bw->getAlphaBackgroundResolution(), vtAlphaBg); // might be missing data at end because of simbox 
      bw->getVerticalTrend(bw->getBetaBackgroundResolution(), vtBetaBg);
      bw->getVerticalTrend(bw->getRhoBackgroundResolution(), vtRhoBg);
      
      extrapolate(vtAlphaBg,nz_);
      extrapolate(vtBetaBg, nz_);
      extrapolate(vtRhoBg, nz_);
      
      bw->getVerticalTrend(bw->getAlpha(),vtAlpha);
      bw->getVerticalTrend(bw->getBeta(),vtBeta);
      bw->getVerticalTrend(bw->getRho(),vtRho);
      bw->getVerticalTrend(bw->getFacies(),vtFacies,random_);
      const int * ipos = bw->getIpos();
      const int * jpos = bw->getJpos();
      const int * kpos = bw->getKpos();

      for(i=0;i<nz_;i++)
      {
        simbox_->getCoord(ipos[0],jpos[0],kpos[0]+i,x,y,z);
        insideOrigSimbox = origSimbox_.getIndex(x,y,z);
        if(vtAlpha[i]!=RMISSING && vtBeta[i]!=RMISSING && vtRho[i]!=RMISSING && insideOrigSimbox!=IMISSING)
          facieslog_[i+w1*nz_] = vtFacies[i];
        else
          facieslog_[i+w1*nz_] = IMISSING;
      }

      float dz = static_cast<float> (simbox_->getdz()*simbox_->getRelThick(ipos[0],jpos[0]));
      Kriging1D::krigVector(vtAlpha, vtAlphaBg, nz_, dz);
      Kriging1D::krigVector(vtBeta, vtBetaBg, nz_, dz);
      Kriging1D::krigVector(vtRho, vtRhoBg, nz_, dz);
      simbox_->getCoord(ipos[0],jpos[0],kpos[0],x,y,z);

      for(i=0;i<nz_;i++)
      {
        alpha_rAmp[i] = vtAlpha[i] - vtAlphaBg[i];
        beta_rAmp[i] = vtBeta[i] - vtBetaBg[i];
        rho_rAmp[i] = vtRho[i] - vtRhoBg[i];
        
        fprintf(file2,"%d %f %f %f %f %d\n",w1, z+i*dz,vtAlpha[i], vtBeta[i], vtRho[i],facieslog_[i+w1*nz_]);
        fprintf(file3,"%d %f %f %f %f %d\n",w1, z+i*dz,vtAlphaBg[i], vtBetaBg[i], vtRhoBg[i],facieslog_[i+w1*nz_] );
      }
      
      j = nz_-1;
      float diffa = (alpha_rAmp[0]-alpha_rAmp[j])/(nzp_-j-1);
      float diffb = (beta_rAmp[0]-beta_rAmp[j])/(nzp_-j-1);
      float diffr = (rho_rAmp[0]-rho_rAmp[j])/(nzp_-j-1);
      //  for(i=sBW[w1]->GetNumberOfObs();i<nzp_;i++) // do padding
      for(i=nz_;i<nzp_;i++)
      {
        alpha_rAmp[i] = alpha_rAmp[j]+(i-j)*diffa;
        beta_rAmp[i] = beta_rAmp[j]+(i-j)*diffb;
        rho_rAmp[i] = rho_rAmp[j]+(i-j)*diffr;
      }
      // do fourier transform of alphadiff,betadiff,rhodiff
      
      rfftwnd_one_real_to_complex(p1, alpha_rAmp, alpha_cAmp);
      rfftwnd_one_real_to_complex(p1, beta_rAmp, beta_cAmp);
      rfftwnd_one_real_to_complex(p1, rho_rAmp, rho_cAmp);
      
      //loop over frequencies, create 3*3 matrices, and the filter
      fftw_complex **sigmaK = new fftw_complex*[3];
      fftw_complex **sigmaE = new fftw_complex*[3];
      fftw_complex *paramvec = new fftw_complex[3];
      fftw_complex *help = new fftw_complex[3];
      int ok;
      double **F = new double*[3];
      
      for(i=0;i<3;i++){
        sigmaK[i] = new fftw_complex[3];
        sigmaE[i] = new fftw_complex[3];
        F[i] =new double[3];
      }

      for(w=0;w<nzp_/2+1;w++)
      {
        //  if(corrprior_cAmp[w].re<delta)
        //    corrprior_cAmp[w].re = 0.0;
        sigmaK[0][0].re = corrprior_cAmp[w].re*sigma0_[0][0];
        //sigmaK[0][0].im = corrprior_cAmp[w].im*sigma0_[0][0];
        sigmaK[0][0].im = 0.0;
        sigmaK[1][1].re = corrprior_cAmp[w].re*sigma0_[1][1];
        //sigmaK[1][1].im = corrprior_cAmp[w].im*sigma0_[1][1];
        sigmaK[1][1].im = 0.0;
        sigmaK[2][2].re = corrprior_cAmp[w].re*sigma0_[2][2];
        //sigmaK[2][2].im = corrprior_cAmp[w].im*sigma0_[2][2];
        sigmaK[2][2].im = 0.0;
        sigmaK[1][0].re = corrprior_cAmp[w].re*sigma0_[1][0];
        //sigmaK[1][0].im = corrprior_cAmp[w].im*sigma0_[1][0];
        sigmaK[1][0].im = 0.0;
        sigmaK[0][1].re = corrprior_cAmp[w].re*sigma0_[1][0];
        //sigmaK[0][1].im = -corrprior_cAmp[w].im*sigma0_[1][0];
        sigmaK[0][1].im = 0.0;
        sigmaK[2][0].re = corrprior_cAmp[w].re*sigma0_[2][0];
        //sigmaK[2][0].im = corrprior_cAmp[w].im*sigma0_[2][0];
        sigmaK[2][0].im = 0.0;
        sigmaK[0][2].re = corrprior_cAmp[w].re*sigma0_[2][0];
        //sigmaK[0][2].im = -corrprior_cAmp[w].im*sigma0_[2][0];
        sigmaK[0][2].im = 0.0;
        sigmaK[1][2].re = corrprior_cAmp[w].re*sigma0_[1][2];
        //sigmaK[1][2].im = corrprior_cAmp[w].im*sigma0_[1][2];
        sigmaK[1][2].im = 0.0;
        sigmaK[2][1].re = corrprior_cAmp[w].re*sigma0_[1][2];
        //sigmaK[2][1].im = -corrprior_cAmp[w].im*sigma0_[1][2];
        sigmaK[2][1].im = 0.0;

        sigmaE[0][0].re = sigmaK[0][0].re - postcova_cAmp[w].re;
        sigmaE[0][0].im = sigmaK[0][0].im - postcova_cAmp[w].im;
        sigmaE[1][1].re = sigmaK[1][1].re - postcovb_cAmp[w].re;
        sigmaE[1][1].im = sigmaK[1][1].im - postcovb_cAmp[w].im;
        sigmaE[2][2].re = sigmaK[2][2].re - postcovr_cAmp[w].re;
        sigmaE[2][2].im = sigmaK[2][2].im - postcovr_cAmp[w].im;
        sigmaE[0][1].re = sigmaK[0][1].re - postcovab_cAmp[w].re;
        sigmaE[0][1].im = sigmaK[0][1].im - postcovab_cAmp[w].im;
        sigmaE[1][0].re = sigmaE[0][1].re;
        sigmaE[1][0].im = -sigmaE[0][1].im;
        sigmaE[0][2].re = sigmaK[0][2].re - postcovar_cAmp[w].re;
        sigmaE[0][2].im = sigmaK[0][2].im - postcovar_cAmp[w].im;
        sigmaE[2][0].re = sigmaE[0][2].re;
        sigmaE[2][0].im = -sigmaE[0][2].im;
        sigmaE[1][2].re = sigmaK[1][2].re - postcovbr_cAmp[w].re;
        sigmaE[1][2].im = sigmaK[1][2].im - postcovbr_cAmp[w].im;
        sigmaE[2][1].re = sigmaE[1][2].re;
        sigmaE[2][1].im = -sigmaE[1][2].im;
        
        /*   
             sigmaK[0][0].re = 2.0;
             sigmaK[0][1].re = 0.9;
             sigmaK[0][2].re = 0.1;
             sigmaK[1][0].re = 0.9;
             sigmaK[1][1].re = 1.0;
             sigmaK[1][2].re = 0.2;
             sigmaK[2][0].re = 0.1;
             sigmaK[2][1].re = 0.2;
             sigmaK[2][2].re = 1.0;
             
             sigmaE[0][0].re = 1.0;
             sigmaE[0][0].im = 0.0;
             sigmaE[0][1].re = 0.4;
             sigmaE[0][1].im = 0.0;
             sigmaE[0][2].re = 0.05;
             sigmaE[0][2].im = 0.0;
             sigmaE[1][0].re = 0.4;
             sigmaE[1][0].im = 0.0;
             sigmaE[1][1].re = 0.5;
             sigmaE[1][1].im = 0.0;
             sigmaE[1][2].re = 0.1;
             sigmaE[1][2].im = 0.0;
             sigmaE[2][0].re = 0.05;
             sigmaE[2][0].im = 0.0;
             sigmaE[2][1].re = 0.1;
             sigmaE[2][1].im = 0.0;
             sigmaE[2][2].re = 0.5;
             sigmaE[2][2].im = 0.0;
        */
        
        //Do cholesky of sigmaK
        ok = lib_matrCholCpx(3,sigmaK);
        if(ok==0 && (w*domega > lowCut && w*domega < highCut) )
        {
          calcFilter(sigmaK, sigmaE, F);
          //apply filter on alpha, beta, rho
          paramvec[0] = alpha_cAmp[w];
          paramvec[1] = beta_cAmp[w];
          paramvec[2] = rho_cAmp[w];
          lib_matrProdMatRVecCpx(F, paramvec, 3,3, help);
          alpha_cAmp[w] = help[0];
          beta_cAmp[w] = help[1];
          rho_cAmp[w] = help[2];
        }
        else
        {
          alpha_cAmp[w].re = 0.0;
          beta_cAmp[w].re = 0.0;
          rho_cAmp[w].re = 0.0;
          alpha_cAmp[w].im = 0.0;
          beta_cAmp[w].im = 0.0;
          rho_cAmp[w].im = 0.0;
        }
      }   
      
      // do inverse fourier transform, and add background model
      rfftwnd_plan p2 = rfftwnd_create_plan(1, &nzp_, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_complex_to_real(p2, alpha_cAmp, alpha_rAmp);
      rfftwnd_one_complex_to_real(p2, beta_cAmp, beta_rAmp);
      rfftwnd_one_complex_to_real(p2, rho_cAmp, rho_rAmp);
      // Add background model to filtered data if not relative method is used.
      
      simbox_->getCoord(ipos[0],jpos[0],kpos[0],x,y,z);
      
      for(i=0;i<nz_;i++)
      {
        alphafiltered_[i+w1*nz_] = alpha_rAmp[i]/nzp_+(1-relative)*vtAlphaBg[i];
        betafiltered_[i+w1*nz_] = beta_rAmp[i]/nzp_+(1-relative)*vtBetaBg[i];
        rhofiltered_[i+w1*nz_] = rho_rAmp[i]/nzp_+(1-relative)*vtRhoBg[i];
        alphablock_[i+w1*nz_] = vtAlpha[i];
        betablock_[i+w1*nz_] = vtBeta[i];
        rhoblock_[i+w1*nz_] = vtRho[i];
        fprintf(file,"%d %f %f %f %f %d \n",w1, z+dz*i,
                alphafiltered_[i+w1*nz_],
                betafiltered_[i+w1*nz_],
                rhofiltered_[i+w1*nz_],
                facieslog_[i+w1*nz_]);   
      }
      fftwnd_destroy_plan(p2);


      for(i=0;i<3;i++)
      {
        delete [] sigmaK[i];
        delete [] sigmaE[i];
        delete [] F[i];
      }
      delete [] sigmaK;
      delete [] sigmaE;
      delete [] F;
      delete [] help;
      delete [] paramvec;
      delete bw;
    } // end if
  } //end for wells
  // Calculate min and max
  fclose(file);
  fclose(file2);
  fclose(file3);
  fftwnd_destroy_plan(p1);

  delete [] vtAlpha;            // vt = vertical trend
  delete [] vtBeta;
  delete [] vtRho;
  delete [] vtFacies;

  delete [] vtAlphaBg;          // vt = vertical trend
  delete [] vtBetaBg;
  delete [] vtRhoBg;

  fftw_free(alpha_rAmp);
  fftw_free(beta_rAmp);
  fftw_free(rho_rAmp);
}

void
FaciesProb::calcFilter(fftw_complex **sigmaK, fftw_complex **sigmaE, double **F)
{
  int i,j;
  int     * ok1 = new int[1];
  double  * ev  = new double[3];
  double ** M   = new double*[3];
  double ** V   = new double*[3];
  double ** EV  = new double*[3];
  double ** sigmaKreal = new double*[3];
  fftw_complex ** help = new fftw_complex*[3];

  for(i=0;i<3;i++){
    help[i] = new fftw_complex[3];
    M[i] = new double[3];
    V[i] = new double[3];
    EV[i] =new double[3];
    sigmaKreal[i] = new double[3];
  }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmaKreal[i][j] = sigmaK[i][j].re;

  lib_matrLXeqBMatCpx(3, sigmaK, sigmaE, 3);
  lib_matrAdjoint(sigmaE,3,3,help);
  lib_matrLXeqBMatCpx(3, sigmaK, help, 3);
  lib_matrAdjoint(help,3,3,sigmaE); //sigmaE = M
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      M[i][j] = sigmaE[i][j].re;

  lib_matr_eigen(M,3,V,ev,ok1);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      if(i==j && ev[i]>0.0)
        EV[i][j] = sqrt(ev[i]);
      else
        EV[i][j] = 0.0;
    }
  lib_matr_prod(V,EV,3,3,3,M);  // M=V*EV
  lib_matrTranspose(V,3,3,EV);  // EV = V^T
  lib_matr_prod(M,EV,3,3,3,V);  // V = A^K=VD^0.5V^T
  //F = L*V*L^-1
    
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      EV[i][j] = sigmaKreal[i][j];
  EV[0][2] = 0.0;
  EV[0][1] = 0.0;
  EV[1][2] = 0.0;
  lib_matr_prod(EV,V,3,3,3,M); //M=LA
  lib_matrTranspose(M,3,3,V);
  lib_matrLtXeqBR(3, sigmaKreal, V, 3);  //L^T*F^T = V
  lib_matrTranspose(V,3,3,F);     
  
  delete [] ok1;
  delete [] ev;
  for (int i=0 ; i < 3 ; i++)
  {
    delete [] M[i];
    delete [] V[i];
    delete [] EV[i];
    delete [] sigmaKreal[i];
    delete [] help[i];

  }
  delete [] M;
  delete [] V;
  delete [] EV;
  delete [] sigmaKreal;
  delete [] help;
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

void
FaciesProb::extrapolate(float * log,
                        int     nz) 
{
  //
  // Extrapolate log[] in both ends if needed
  //  
  int i=0;
  while (i<nz && log[i]==RMISSING)
    i++;
  if (i < nz - 1) 
  { 
    int first_nonmissing = i;
    i = nz - 1;
    while (i>0 && log[i]==RMISSING)
      i--;
    int last_nonmissing = i;
    
    for(int i=0 ; i < first_nonmissing ; i++) { 
      log[i] = log[first_nonmissing]; 
    }
    for(int i=last_nonmissing + 1 ; i < nz ; i++) { 
      log[i] = log[last_nonmissing]; 
    }
    //  if (first_nonmissing > 0)
    //    LogKit::LogFormatted(LogKit::LOW,"Vertical trend for %s extrapolated first %d cells.\n",
    //                      pName, first_nonmissing + 1);
    //    if (nz - 1 - last_nonmissing > 0)
    //     LogKit::LogFormatted(LogKit::LOW,"Vertical trend for %s extrapolated last %d cells.\n",
    //                      pName, nz - 1 - last_nonmissing);
    //for (int i=0 ; i<nz ; i++) {
    //  printf("i log[i]  %d %.3f\n",i,log[i]);
    //}
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"WARNING: All ... trend values are missing.\n");
  }
}
