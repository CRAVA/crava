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
#include "src/simbox.h"
#include "src/io.h"

FaciesProb::FaciesProb(FFTGrid                      * alpha,
                       FFTGrid                      * beta,
                       FFTGrid                      * rho,
                       int                            nFac,
                       float                          p_undef, 
                       const float                  * priorFacies,
                       FFTGrid                     ** priorFaciesCubes,
                       const std::vector<double **> & sigmaEOrig,
                       const WellData              ** wells,
                       int                            nWells,
                       const std::vector<Surface *> & faciesEstimInterval,
                       const double                   dz,
                       bool                           relative,
                       bool                           noVs)
{
  makeFaciesProb(nFac, alpha, beta, rho, sigmaEOrig, wells, nWells, 
                 faciesEstimInterval, dz, relative, noVs,
                 p_undef, priorFacies, priorFaciesCubes);
}

FaciesProb::~FaciesProb()
{
  for(int i=0;i<nFacies_;i++)
  {
    delete faciesProb_[i];
  }
  delete [] faciesProb_;
  delete faciesProbUndef_;
}

std::vector<FFTGrid *>
FaciesProb::makeFaciesHistAndSetPriorProb(const std::vector<float> & alpha,
                                          const std::vector<float> & beta,
                                          const std::vector<float> & rho,
                                          const std::vector<int>   & faciesLog,
                                          const Simbox             * volume)
{
  std::vector<FFTGrid *> hist;
  hist.resize(nFacies_, NULL);
  int i,j,k,l;
  for(i=0;i<nFacies_;i++)
  {
    hist[i] = new FFTGrid(volume->getnx(), volume->getny(), volume->getnz(),
                          volume->getnx(), volume->getny(), volume->getnz());
    hist[i]->createRealGrid();
    hist[i]->setType(FFTGrid::PARAMETER);
    for(l=0;l<volume->getnz();l++)
    {
      for(k=0;k<volume->getny();k++)
      {
        for(j=0;j<volume->getnx();j++)
          hist[i]->setRealValue(j, k, l, 0.0f);
      }
    }
  }
  
  std::vector<int> nData(nFacies_,0);

  int facies;
  for(i=0;i<static_cast<int>(alpha.size());i++)
  {
    if(faciesLog[i]!=IMISSING)
    {
      volume->getIndexes(alpha[i], beta[i], rho[i], j, k, l);
      facies = faciesLog[i];
      nData[facies]++;
      float value = hist[facies]->getRealValue(j, k, l) + 1.0f;
      hist[facies]->setRealValue(j,k,l,value);
    }
  }
 
  for(i=0;i<nFacies_;i++)
  {
    double nf = 1.0/double(nData[i]);
    for(l=0;l<volume->getnz();l++)
    {
      for(k=0;k<volume->getny();k++)
      {
        for(j=0;j<volume->getnx();j++) {
          double value = static_cast<double>(hist[i]->getRealValue(j, k, l))*nf;
          hist[i]->setRealValue(j, k, l, static_cast<float>(value));
        }
      }
    }
  }
  return hist;
}

void       
FaciesProb::makeFaciesDens(int nfac, 
                           const std::vector<double **> & sigmaEOrig,
                           bool                           noVs,
                           const std::vector<float>     & alphaFiltered,
                           const std::vector<float>     & betaFiltered,
                           const std::vector<float>     & rhoFiltered,
                           const std::vector<int>       & faciesLog,
                           std::vector<FFTGrid *>       & density,
                           Simbox                      ** volume)
{ 
  //Note: If noVs is true, the beta dimension is mainly dummy. Due to the lookup mechanism that maps 
  //      values outside the denisty table to the edge, any values should do in this dimension.
  //      Still, we prefer to create reasonable values.
  int i,j,k,l;
  nFacies_ = nfac;

  float kstda, kstdb, kstdr, hopt;
  int nbinsa = 100;
  int nbinsb = 100;
  if(noVs == true)
    nbinsb = 1; //Ignore what happens for Vs.
  int nbinsr = 50;


  int nobs   = 0;
  float *smooth = new float[nbinsa*nbinsb*nbinsr];
  for(i=0;i<static_cast<int>(alphaFiltered.size());i++)
  {
    if(faciesLog[i]!=IMISSING)
      nobs++;
  }

  // Make bins.
  float varAlpha = 0.0f, varBeta = 0.0f, varRho = 0.0f;
  calculateVariances(alphaFiltered, betaFiltered, rhoFiltered, faciesLog,
                     varAlpha, varBeta, varRho);//sets varAlpha etc....
  if(noVs == true)
    varBeta = 5*varAlpha; //Must be large enough to make surface span possible beta values.

  hopt  = static_cast<float>(pow(4.0/7,1.0/7)*pow(static_cast<double>(nobs),-1.0/7));
  kstda = hopt*sqrt(varAlpha);
  kstdb = hopt*sqrt(varBeta);
  kstdr = hopt*sqrt(varRho);

  double ** sigmae = new double *[3];
  for(i=0;i<3;i++)
    sigmae[i] = new double[3];

  if(noVs == false) {
    for(i=0;i<3;i++) {
      for(j=0;j<3;j++)
        sigmae[i][j] = sigmaEOrig[0][i][j];
    }
  }
  else {
    sigmae[0][0] = sigmaEOrig[0][0][0];
    sigmae[0][1] = 0.0;
    sigmae[0][2] = sigmaEOrig[0][0][1];
    sigmae[1][0] = 0.0;
    sigmae[1][1] = 0.0; //Will be overruled by reasonable value.
    sigmae[1][2] = 0.0;
    sigmae[2][0] = sigmaEOrig[0][1][0];
    sigmae[2][1] = 0.0;
    sigmae[2][2] = sigmaEOrig[0][1][1];
  }

  if(sigmae[0][0]<kstda*kstda)
    sigmae[0][0] = kstda*kstda;
  if(sigmae[1][1]<kstdb*kstdb)
    sigmae[1][1] = kstdb*kstdb;
  if(sigmae[2][2]<kstdr*kstdr)
    sigmae[2][2] = kstdr*kstdr;

  // invert sigmae
  double **sigmaeinv = new double *[3];
  for(i=0;i<3;i++)
    sigmaeinv[i] = new double [3];

  for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        if(i==j)
          sigmaeinv[i][j] = 1.0;
        else
          sigmaeinv[i][j] = 0.0;
  lib_matrCholR(3, sigmae);
  lib_matrAXeqBMatR(3, sigmae, sigmaeinv, 3);

  for(int i=0;i<3;i++)
    delete [] sigmae[i];
  delete [] sigmae;

  float alphaMin = *min_element(alphaFiltered.begin(), alphaFiltered.end()) - 5.0f*kstda;
  float alphaMax = *max_element(alphaFiltered.begin(), alphaFiltered.end()) + 5.0f*kstda;
  float betaMin  = *min_element(betaFiltered.begin(), betaFiltered.end()) - 5.0f*kstdb;
  float betaMax  = *max_element(betaFiltered.begin(), betaFiltered.end()) + 5.0f*kstdb;
  float rhoMin   = *min_element(rhoFiltered.begin(), rhoFiltered.end()) - 5.0f*kstdr;
  float rhoMax   = *max_element(rhoFiltered.begin(), rhoFiltered.end()) + 5.0f*kstdr;

  float dAlpha = (alphaMax-alphaMin)/nbinsa;
  float dBeta  = (betaMax-betaMin)/nbinsb;
  float dRho   = (rhoMax-rhoMin)/nbinsr;

  Surface * rhoMinSurf = new Surface(alphaMin, betaMin, alphaMax-alphaMin, betaMax-betaMin, 2, 2, rhoMin);
  //Associate alpha with x, beta with y and rho with z.
  *volume  = new Simbox(alphaMin, betaMin, rhoMinSurf, alphaMax-alphaMin, betaMax-betaMin, rhoMax-rhoMin, 0, dAlpha, dBeta, dRho); 
  density = makeFaciesHistAndSetPriorProb(alphaFiltered, betaFiltered, rhoFiltered, faciesLog, *volume);

  int jj,jjj,kk,kkk,ll,lll;
  lll=2;
  
  float sum = 0.0f;
  for(l=0;l<nbinsr;l++)
  {
    kkk=2;
    if(l<=nbinsr/2)
      ll = l;
    else
    {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0;k<nbinsb;k++)
    {
      jjj=2;
      if(k<=nbinsb/2)
        kk=k;
      else
      {
        kk = -(k-kkk);
        kkk+=2;
      }
      for(j=0;j<nbinsa;j++)
      {
        if(j<=nbinsa/2)
          jj=j;
        else
        {
          jj = -(j-jjj);
          jjj+=2;
        }
        smooth[j+k*nbinsa+l*nbinsa*nbinsb] = float(exp(-0.5f*(jj*dAlpha*jj*dAlpha*sigmaeinv[0][0]
                                                             +kk*dBeta *kk*dBeta *sigmaeinv[1][1]
                                                             +ll*dRho  *ll*dRho  *sigmaeinv[2][2]
                                                           +2*jj*dAlpha*kk*dBeta *sigmaeinv[1][0]
                                                           +2*jj*dAlpha*ll*dRho  *sigmaeinv[2][0]
                                                           +2*kk*dBeta *ll*dRho  *sigmaeinv[2][1])));
        sum = sum+smooth[j+k*nbinsa+l*nbinsa*nbinsb];
      }
    }
  }
  // normalize smoother
  for(l=0;l<nbinsr;l++)
    for(k=0;k<nbinsb;k++)
      for(j=0;j<nbinsa;j++)
        smooth[j+k*nbinsa+l*nbinsa*nbinsb]/=sum;

  FFTGrid *smoother;
  smoother = new FFTGrid(nbinsa, nbinsb ,nbinsr ,nbinsa ,nbinsb ,nbinsr);
  smoother->fillInFromArray(smooth);

  if(ModelSettings::getDebugLevel() >= 1) {
    std::string baseName = "Smoother" + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    smoother->writeAsciiFile(fileName);
  }

  smoother->fftInPlace();

  for(i=0;i<nFacies_;i++)
  {
    density[i]->fftInPlace();
    density[i]->multiply(smoother);
    density[i]->invFFTInPlace();
    density[i]->multiplyByScalar(float(sqrt(double(nbinsa*nbinsb*nbinsr))));
    if(ModelSettings::getDebugLevel() >= 1)
    {
      std::string baseName = "Dens_" + NRLib::ToString(i) + IO::SuffixAsciiFiles();
      std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
      density[i]->writeAsciiFile(fileName);
    }
  }
  delete smoother;
  delete [] smooth;
  for(i=0;i<3;i++)
    delete [] sigmaeinv[i];
  delete [] sigmaeinv;
}


void FaciesProb::makeFaciesProb(int                            nfac, 
                                FFTGrid                      * postAlpha, 
                                FFTGrid                      * postBeta, 
                                FFTGrid                      * postRho,
                                const std::vector<double **> & sigmaEOrig, 
                                const WellData              ** wells, 
                                int                            nWells,
                                const std::vector<Surface *> & faciesEstimInterval,
                                const double                   dz,
                                bool                           relative,
                                bool                           noVs,
                                float                          p_undef,
                                const float                  * priorFacies,
                                FFTGrid                     ** priorFaciesCubes)
{
  std::vector<float> alphaFiltered;
  std::vector<float> betaFiltered;
  std::vector<float> rhoFiltered;
  std::vector<int>   faciesLog;

  setNeededLogsSpatial(wells, nWells, faciesEstimInterval, dz, relative, noVs,
                       alphaFiltered, betaFiltered, rhoFiltered, faciesLog); //Generate these logs.

  std::vector<FFTGrid *> density;
  Simbox   * volume  = NULL;

  makeFaciesDens(nfac, sigmaEOrig, noVs, alphaFiltered, betaFiltered, rhoFiltered, faciesLog, density, &volume);

  if(priorFaciesCubes != NULL)
    normalizeCubes(priorFaciesCubes);
  calculateFaciesProb(postAlpha, postBeta, postRho, density, volume,
                      p_undef, priorFacies, priorFaciesCubes);
}

float FaciesProb::findDensity(float alpha, float beta, float rho, 
                              FFTGrid * density, const Simbox * volume)
{
  double jFull, kFull, lFull;
  volume->getInterpolationIndexes(alpha, beta, rho, jFull, kFull, lFull);
  int j1,k1,l1;
  int j2,k2,l2;
  float wj,wk,wl;
  j1=k1=l1=j2=k2=l2=0;
  j1 = static_cast<int>(floor(jFull));
  if(j1<0) {
    j1 = 0;
    j2 = 0;
    wj = 0;
  }
  else if(j1>=volume->getnx()-1) {
    j1 = volume->getnx()-1;
    j2 = j1;
    wj = 0;
  }
  else {
    j2 = j1 + 1;
    wj = static_cast<float>(jFull-j1);
  }

  k1 = static_cast<int>(floor(kFull));
  if(k1<0) {
    k1 = 0;
    k2 = 0;
    wk = 0;
  }
  else if(k1>=volume->getny()-1) {
    k1 = volume->getny()-1;
    k2 = k1;
    wk = 0;
  }
  else {
    k2 = k1 + 1;
    wk = static_cast<float>(kFull-k1);
  }

  l1 = static_cast<int>(floor(lFull));
  if(l1<0) {
    l1 = 0;
    l2 = 0;
    wl = 0;
  }
  else if(l1>=volume->getnz()-1) {
    l1 = volume->getnz()-1;
    l2 = l1;
    wl = 0;
  }
  else {
    l2 = l1 + 1;
    wl = static_cast<float>(lFull-l1);
  }

  float value1 = MAXIM(0,density->getRealValue(j1,k1,l1));
  float value2 = MAXIM(0,density->getRealValue(j1,k1,l2));
  float value3 = MAXIM(0,density->getRealValue(j1,k2,l1));
  float value4 = MAXIM(0,density->getRealValue(j1,k2,l2));
  float value5 = MAXIM(0,density->getRealValue(j2,k1,l1));
  float value6 = MAXIM(0,density->getRealValue(j2,k1,l2));
  float value7 = MAXIM(0,density->getRealValue(j2,k2,l1));
  float value8 = MAXIM(0,density->getRealValue(j2,k2,l2));
  
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

void FaciesProb::calculateConditionalFaciesProb(WellData                    ** wells, 
                                                int                            nWells, 
                                                const std::vector<Surface *> & faciesEstimInterval,
                                                const ModelSettings          * modelSettings,
                                                const double                   dz)
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

    if (faciesEstimInterval.size() > 0) {
      const double * xPos  = bw[i]->getXpos();
      const double * yPos  = bw[i]->getYpos();
      const double * zPos  = bw[i]->getZpos();

      for (int b = 0 ; b < nBlocks ; b++) {
        const double zTop  = faciesEstimInterval[0]->GetZ(xPos[b],yPos[b]);
        const double zBase = faciesEstimInterval[1]->GetZ(xPos[b],yPos[b]);
        if ( (zPos[b]-0.5*dz) < zTop || (zPos[b]+0.5*dz) > zBase)
          BWfacies[i][b] = IMISSING;
        else
          BWfacies[i][b] = BWfacies_i[b];
        //LogKit::LogFormatted(LogKit::LOW,"ib, BWfacies[i][b] = %d %d\n",ib,BWfacies[ib]);
      }
    }
    else {
      for (int b = 0 ; b < nBlocks ; b++) {
        BWfacies[i][b] = BWfacies_i[b];
        //LogKit::LogFormatted(LogKit::LOW,"ib, BWfacies[i][b] = %d %d\n",ib,BWfacies[ib]);
      }
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
      LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings->getFaciesName(f).c_str());
    LogKit::LogFormatted(LogKit::LOW,"\n------------+");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::LOW,"------------",f);
    LogKit::LogFormatted(LogKit::LOW,"\n");
    for(int f1=0 ; f1 < nFacies_ ; f1++)
    {
      LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings->getFaciesName(f1).c_str());
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
    LogKit::LogFormatted(LogKit::LOW," %11s",modelSettings->getFaciesName(f).c_str());
  LogKit::LogFormatted(LogKit::LOW,"\n------------+");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::LOW,"------------",f);
  LogKit::LogFormatted(LogKit::LOW,"\n");
  for(int f1=0 ; f1 < nFacies_ ; f1++)
  {
    LogKit::LogFormatted(LogKit::LOW,"%-11s |",modelSettings->getFaciesName(f1).c_str());
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

void FaciesProb::calculateFaciesProb(FFTGrid                      * alphagrid, 
                                     FFTGrid                      * betagrid, 
                                     FFTGrid                      * rhogrid,
                                     const std::vector<FFTGrid *> & density, 
                                     const Simbox                 * volume, 
                                     float                          p_undefined,
                                     const float *                  priorFacies,
                                     FFTGrid                    ** priorFaciesCubes)
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
  if(alphagrid->isFile()==1)
    faciesProbUndef_ = new FFTFileGrid(nx, ny, nz, nx, ny, nz);
  else
    faciesProbUndef_ = new FFTGrid(nx, ny, nz, nx, ny, nz);
  alphagrid->setAccessMode(FFTGrid::READ);
  betagrid->setAccessMode(FFTGrid::READ);
  rhogrid->setAccessMode(FFTGrid::READ);
  for(i=0;i<nFacies_;i++)
  {
    if(alphagrid->isFile()==1)
    {
      faciesProb_[i] = new FFTFileGrid(nx, ny, nz, nx, ny, nz);
    }
    else
    {
      faciesProb_[i] = new FFTGrid(nx, ny, nz, nx, ny, nz);
    }
    faciesProb_[i]->setAccessMode(FFTGrid::WRITE);
    faciesProb_[i]->createRealGrid();
  }
  faciesProbUndef_->setAccessMode(FFTGrid::WRITE);
  faciesProbUndef_->createRealGrid();

  smallrnxp = faciesProb_[0]->getRNxp();

  float help;
  float undefSum = p_undefined/(volume->getnx()*volume->getny()*volume->getnz());
  for(i=0;i<nzp;i++)
  {
    for(j=0;j<nyp;j++)
      for(k=0;k<rnxp;k++)
      {
        alpha = alphagrid->getNextReal();
        beta = betagrid->getNextReal();
        rho = rhogrid->getNextReal();
        if(k<smallrnxp && j<ny && i<nz)
        {

          sum = undefSum;
          for(l=0;l<nFacies_;l++)
          {
            if(priorFaciesCubes!=NULL)
              value[l] = priorFaciesCubes[l]->getRealValue(k,j,i)*findDensity(alpha, beta, rho, density[l], volume);
            else
              value[l] = priorFacies[l]*findDensity(alpha, beta, rho, density[l], volume); 
            sum = sum+value[l];
          }
          for(l=0;l<nFacies_;l++)
          {
            help = value[l]/sum;
            if(k<nx)
            {
              faciesProb_[l]->setNextReal(help);         
            }
            else
            {
              faciesProb_[l]->setNextReal(RMISSING);             
            }
          }
          if(k<nx)
            faciesProbUndef_->setNextReal(undefSum/sum);
          else
            faciesProbUndef_->setNextReal(RMISSING);

        }
      }
  }
  delete [] value;
}

void FaciesProb::calculateFaciesProbGeomodel(const float *                  priorFacies,
                                             FFTGrid                    ** priorFaciesCubes)
{
  int i,j,k,l;
  int nx, ny, nz;
  float value;
  nx   = faciesProb_[0]->getNx();
  ny   = faciesProb_[0]->getNy();
  nz   = faciesProb_[0]->getNz();
  float undef;

  for(i=0;i<nz;i++)
  {
    for(j=0;j<ny;j++)
      for(k=0;k<nx;k++)
      {
        undef = faciesProbUndef_->getRealValue(k,j,i);
        for(l=0;l<nFacies_;l++)
        {
          if(priorFaciesCubes!=NULL)
            value = priorFaciesCubes[l]->getRealValue(k,j,i)*undef+faciesProb_[l]->getRealValue(k,j,i);
          else
            value = priorFacies[l]*undef+faciesProb_[l]->getRealValue(k,j,i); 
          faciesProb_[l]->setRealValue(k,j,i,value);

        }
      }
  }
}

void FaciesProb::calculateVariances(const std::vector<float> & alpha,
                                    const std::vector<float> & beta,
                                    const std::vector<float> & rho,
                                    const std::vector<int>   & facies,
                                    float                    & varAlpha,
                                    float                    & varBeta,
                                    float                    & varRho)
{
  int i;
  //bool validA, validB, validR;
  int nA, nB, nR;
  nA = nB = nR = 0;
  float meanA = 0.0f;
  float meanB = 0.0f;
  float meanR = 0.0f;
  int nData = static_cast<int>(alpha.size());
  for(i=0;i<nData;i++)
  {
    if(facies[i]!=IMISSING)
    {
      meanA +=alpha[i];
      nA++;
      meanB +=beta[i];
      nB++;
      meanR +=rho[i];
      nR++;
    }
  }
  meanA/=nA;
  meanB/=nB;
  meanR/=nR;
  varAlpha = varBeta = varRho = 0.0;

  for(i=0;i<nData;i++)
  {
    if(facies[i]!=IMISSING)
    {
      varAlpha += pow(alpha[i]-meanA,2);
      varBeta  += pow(beta[i]-meanB,2);
      varRho   += pow(rho[i]-meanR,2);
    }
  }
  varAlpha /= nA-1;
  varBeta  /= nB-1;
  varRho   /= nR-1;
}


void FaciesProb::setNeededLogsSpatial(const WellData              ** wells,
                                      int                            nWells,
                                      const std::vector<Surface *> & faciesEstimInterval,
                                      const double                   dz,
                                      bool                           relative,
                                      bool                           noVs,
                                      std::vector<float>           & alphaFiltered,
                                      std::vector<float>           & betaFiltered,
                                      std::vector<float>           & rhoFiltered,
                                      std::vector<int>             & faciesLog)
{
  int nData = 0;
  for(int w=0;w<nWells;w++)
    if(wells[w]->getUseForFaciesProbabilities())
      nData += wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
  alphaFiltered.resize(nData, RMISSING);
  betaFiltered.resize(nData, RMISSING);
  rhoFiltered.resize(nData, RMISSING);
  faciesLog.resize(nData, IMISSING);

  int index = 0;
  for (int w=0 ; w<nWells ; w++)
  {
    if(wells[w]->getUseForFaciesProbabilities())
    { 
      BlockedLogs * bw = wells[w]->getBlockedLogsOrigThick();
      int n = bw->getNumberOfBlocks();
      for(int i=0;i<n;i++)
      {
        if(noVs == false) {
          alphaFiltered[index] = bw->getAlphaSeismicResolution()[i];
          betaFiltered[index]  = bw->getBetaSeismicResolution()[i];
          rhoFiltered[index]   = bw->getRhoSeismicResolution()[i];
        }
        else {
          alphaFiltered[index] = bw->getAlphaForFacies()[i];
          //Beta log here is mainly dummy, but must be at correct level.
          betaFiltered[index]  = bw->getBetaHighCutBackground()[i];
          rhoFiltered[index]   = bw->getRhoForFacies()[i];
        }
        if(relative == true) {
          alphaFiltered[index] -= bw->getAlphaHighCutBackground()[i];
          betaFiltered[index]  -= bw->getBetaHighCutBackground()[i];
          rhoFiltered[index]   -= bw->getRhoHighCutBackground()[i];
        }
        if(bw->getAlpha()[i]==RMISSING || bw->getBeta()[i]==RMISSING || bw->getRho()[i]==RMISSING)
          faciesLog[index] = IMISSING;
        else
          faciesLog[index] = bw->getFacies()[i];

        if (faciesEstimInterval.size() > 0) {
          const double * xPos  = bw->getXpos();
          const double * yPos  = bw->getYpos();
          const double * zPos  = bw->getZpos();
          const double   zTop  = faciesEstimInterval[0]->GetZ(xPos[i],yPos[i]);
          const double   zBase = faciesEstimInterval[1]->GetZ(xPos[i],yPos[i]);
          if ( (zPos[i]-0.5*dz) < zTop || (zPos[i]+0.5*dz) > zBase)
            faciesLog[index] = IMISSING;
        }
        index++;
      }
    }
  }
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
