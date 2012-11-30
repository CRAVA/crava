/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string>

#include "lib/lib_matr.h"
#include "lib/random.h"
#include "lib/kriging1d.h"

#include "nrlib/random/chisquared.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/welldata.h"
#include "src/faciesprob.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/corr.h"
#include "src/crava.h"
#include "src/simbox.h"
#include "src/blockedlogs.h"
#include "src/modelsettings.h"
#include "src/definitions.h"
#include "src/filterwelllogs.h"
#include "src/spatialwellfilter.h"
#include "src/simbox.h"
#include "src/io.h"
#include "src/tasklist.h"

FaciesProb::FaciesProb(FFTGrid                        * alpha,
                       FFTGrid                        * beta,
                       FFTGrid                        * rho,
                       const std::vector<std::string> & faciesNames,
                       int                              nFac,
                       float                            p_undef,
                       const float                    * priorFacies,
                       FFTGrid                       ** priorFaciesCubes,
                       const std::vector<double **>   & sigmaEOrig,
                       bool                             useFilter,
                       const WellData                ** wells,
                       int                              nWells,
                       const std::vector<Surface *>   & faciesEstimInterval,
                       const double                     dz,
                       bool                             relative,
                       bool                             noVs,
                       Crava                          * cravaResult,
                       const std::vector<Grid2D *>    & noiseScale,
                       const ModelSettings *            modelSettings,
                       FFTGrid                        * seismicLH)
{
  makeFaciesProb(nFac, alpha, beta, rho, sigmaEOrig, useFilter, wells, nWells,
                 faciesEstimInterval, dz, relative, noVs,
                 p_undef, priorFacies, priorFaciesCubes, cravaResult, noiseScale,
                 modelSettings, seismicLH);

  checkProbabilities(faciesNames,
                     faciesProb_,
                     nFac);
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

void
FaciesProb::checkProbabilities(const std::vector<std::string>  & faciesNames,
                               FFTGrid                        ** faciesProb,
                               int                               nFacies) const
{
  //
  // Check that facies probabilities are in legal range [0, 1]
  //
  for (int i = 0 ; i < nFacies ; i++) {
    faciesProb[i]->calculateStatistics();
    float min_prob = faciesProb[i]->getMinReal();
    float max_prob = faciesProb[i]->getMaxReal();
    if (min_prob < 0.0f) {
      std::string text = "The minimum probability for facies \'" + faciesNames[i] + "\' is less than zero (" + NRLib::ToString(min_prob, 2) + ")\n";
      LogKit::LogFormatted(LogKit::Warning, "\nERROR: " + text + "       This indicates an error in the calculation.\n");
      TaskList::addTask(text + "    Check the facies probability estimation.\n");
    }
    if (min_prob > 1.0f) {
      std::string text = "The maximum probability for facies \'" + faciesNames[i] + "\' is larger than one (" + NRLib::ToString(max_prob, 2) + ")\n";
      LogKit::LogFormatted(LogKit::Warning, "\nERROR: " + text + "       This indicates an error in the calculation.\n");
      TaskList::addTask(text + "    Check the facies probability estimation.\n");
    }
  }
}

std::vector<FFTGrid *>
FaciesProb::makeFaciesHistAndSetPriorProb(const std::vector<float> & alpha, //
                                          const std::vector<float> & beta,
                                          const std::vector<float> & rho,
                                          const std::vector<int>   & faciesLog,
                                          const Simbox             * volume)
{
  std::vector<FFTGrid *> hist;
  hist.resize(nFacies_, NULL);
  int i,j,k,l;
  int rnxp;
  int nx, ny, nz;
  nx = volume->getnx();
  ny = volume->getny();
  nz = volume->getnz();
  for(i=0;i<nFacies_;i++)
  {
    hist[i] = new FFTGrid(nx, ny, nz, nx , ny, nz);
    hist[i]->createRealGrid(false);
    rnxp = hist[i]->getRNxp();
    hist[i]->setType(FFTGrid::PARAMETER);
    hist[i]->setAccessMode(FFTGrid::WRITE);
    for(l=0;l<nz;l++)
    {
      for(k=0;k<ny;k++)
      {
        for(j=0;j<rnxp;j++)
          hist[i]->setNextReal(0.0f);
      }
    }
    hist[i]->endAccess();
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
      hist[facies]->setAccessMode(FFTGrid::RANDOMACCESS);
      float value = hist[facies]->getRealValue(j, k, l) + 1.0f;
      hist[facies]->setRealValue(j,k,l,value);
      hist[facies]->endAccess();
    }
  }

  for(i=0;i<nFacies_;i++)
  {
    hist[i]->setAccessMode(FFTGrid::READANDWRITE);
    rnxp = hist[i]->getRNxp();
    double nf = 1.0/double(nData[i]);
    for(l=0;l<nz;l++)
    {
      for(k=0;k<ny;k++)
      {
        for(j=0;j<rnxp;j++) {
          double value = static_cast<double>(hist[i]->getNextReal())*nf;
          hist[i]->setNextReal(static_cast<float>(value));
        }
      }
    }
    hist[i]->endAccess();
  }
  return hist;
}

void
FaciesProb::makeFaciesDens(int nfac,
                           const std::vector<double **> & sigmaEOrig,
                           bool                           useFilter,
                           bool                           noVs,
                           const std::vector<float>     & alphaFiltered,
                           const std::vector<float>     & betaFiltered,
                           const std::vector<float>     & rhoFiltered,
                           const std::vector<int>       & faciesLog,
                           std::vector<FFTGrid *>       & density,
                           Simbox                      ** volume,
                           int                            index,
                           double                      ** G,
                           Crava                        * cravaResult,
                           const std::vector<Grid2D *>  & noiseScale)
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
  std::vector<float> alphaFilteredNew(alphaFiltered.size());
  std::vector<float> betaFilteredNew(betaFiltered.size());
  std::vector<float> rhoFilteredNew(rhoFiltered.size());
  if(index>0)
  {
    // Multiply alpha, beta, rho with scale
    double ** H    = new double * [3];
    double ** junk = new double * [3];
    for(i=0;i<3;i++) {
      H[i]    = new double[3];
      junk[i] = new double[3];
    }

    int nAng = static_cast<int>(noiseScale.size());
    // int nAng = int(log(float(sigmaEOrig.size()))/log(2.0));

    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }
    //Compute H matrix
    std::vector<double> scale(nAng);
    int factor = 1;
    for(int angle=0;angle<nAng;angle++) {
      if((index & factor) > 0)
        scale[angle] = maxScale[angle];
      else
        scale[angle] = 1.0;
      factor *= 2;
    }

    cravaResult->newPosteriorCovPointwise(H, G, scale, junk);

    double **help1 = new double*[3];
    help1[0] = new double[1];
    help1[1] = new double[1];
    help1[2] = new double[1];
    double **help2 = new double*[3];
    help2[0] = new double[1];
    help2[1] = new double[1];
    help2[2] = new double[1];
    for(i=0;i<int(alphaFiltered.size());i++)
    {
      help1[0][0] = alphaFiltered[i];
      help1[1][0] = betaFiltered[i];
      help1[2][0] = rhoFiltered[i];
      lib_matr_prod(H,help1,3,3,1,help2);  //
      alphaFilteredNew[i] = float(help2[0][0]);
      betaFilteredNew[i] = float(help2[1][0]);
      rhoFilteredNew[i] = float(help2[2][0]);
    }
    for(i=0;i<3;i++)
    {
      delete [] help1[i];
      delete [] help2[i];
    }
    delete [] help1;
    delete [] help2;
  }
  else
  {
    for(i=0;i<int(alphaFiltered.size());i++)
   {
      alphaFilteredNew[i] = alphaFiltered[i];
      betaFilteredNew[i] = betaFiltered[i];
      rhoFilteredNew[i] = rhoFiltered[i];

    }

  }

  // Make bins.
  float varAlpha = 0.0f, varBeta = 0.0f, varRho = 0.0f;
  calculateVariances(alphaFilteredNew, betaFilteredNew, rhoFilteredNew, faciesLog,
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

  if(useFilter == true) {
    if(noVs == false) {
      for(i=0;i<3;i++) {
        for(j=0;j<3;j++)
          sigmae[i][j] = sigmaEOrig[index][i][j];
      }
    }
    else {
      sigmae[0][0] = sigmaEOrig[index][0][0];
      sigmae[0][1] = 0.0;
      sigmae[0][2] = sigmaEOrig[index][0][1];
      sigmae[1][0] = 0.0;
      sigmae[1][1] = 0.0; //Will be overruled by reasonable value.
      sigmae[1][2] = 0.0;
      sigmae[2][0] = sigmaEOrig[index][1][0];
      sigmae[2][1] = 0.0;
      sigmae[2][2] = sigmaEOrig[index][1][1];
    }
  }
  else {
    for(i=0;i<3;i++) {
      for(j=0;j<3;j++)
        sigmae[i][j] = 0.0;
    }
  }

  if(sigmae[0][0]<kstda*kstda)
    sigmae[0][0] = kstda*kstda;
  if(sigmae[1][1]<kstdb*kstdb)
    sigmae[1][1] = kstdb*kstdb;
  if(sigmae[2][2]<kstdr*kstdr)
    sigmae[2][2] = kstdr*kstdr;

  //Establish limits before we invert sigma.
  double alphaMin = *min_element(alphaFilteredNew.begin(), alphaFilteredNew.end()) - 5.0f*sqrt(sigmae[0][0]);
  double alphaMax = *max_element(alphaFilteredNew.begin(), alphaFilteredNew.end()) + 5.0f*sqrt(sigmae[0][0]);
  double betaMin  = *min_element(betaFilteredNew.begin(), betaFilteredNew.end()) - 5.0f*sqrt(sigmae[1][1]);
  double betaMax  = *max_element(betaFilteredNew.begin(), betaFilteredNew.end()) + 5.0f*sqrt(sigmae[1][1]);
  double rhoMin   = *min_element(rhoFilteredNew.begin(), rhoFilteredNew.end()) - 5.0f*sqrt(sigmae[2][2]);
  double rhoMax   = *max_element(rhoFilteredNew.begin(), rhoFilteredNew.end()) + 5.0f*sqrt(sigmae[2][2]);


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

  double dAlpha = (alphaMax-alphaMin)/nbinsa;
  double dBeta  = (betaMax-betaMin)/nbinsb;
  double dRho   = (rhoMax-rhoMin)/nbinsr;

  Surface rhoMinSurf(alphaMin, betaMin, alphaMax-alphaMin, betaMax-betaMin, 2, 2, rhoMin);
  //Associate alpha with x, beta with y and rho with z.
  *volume  = new Simbox(alphaMin, betaMin, rhoMinSurf, alphaMax-alphaMin, betaMax-betaMin, rhoMax-rhoMin, 0, dAlpha, dBeta, dRho);
  density = makeFaciesHistAndSetPriorProb(alphaFilteredNew, betaFilteredNew, rhoFilteredNew, faciesLog, *volume);

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
    if(ModelSettings::getDebugLevel() >= 1)
    {
      std::string baseName = "Hist_" + NRLib::ToString(i) + IO::SuffixAsciiFiles();
      std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
      density[i]->writeAsciiFile(fileName);
    }
    density[i]->fftInPlace();
    density[i]->multiply(smoother);
    density[i]->invFFTInPlace();
    density[i]->multiplyByScalar(float(sqrt(double(nbinsa*nbinsb*nbinsr))));
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
                                bool                           useFilter,
                                const WellData              ** wells,
                                int                            nWells,
                                const std::vector<Surface *> & faciesEstimInterval,
                                const double                   dz,
                                bool                           relative,
                                bool                           noVs,
                                float                          p_undef,
                                const float                  * priorFacies,
                                FFTGrid                     ** priorFaciesCubes,
                                Crava                        * cravaResult,
                                const std::vector<Grid2D *>  & noiseScale,
                                const ModelSettings          * modelSettings,
                                FFTGrid                      * seismicLH)
{
  std::vector<float> alphaFiltered;
  std::vector<float> betaFiltered;
  std::vector<float> rhoFiltered;
  std::vector<int>   faciesLog;

  setNeededLogsSpatial(wells, nWells, faciesEstimInterval, dz, relative, noVs, useFilter,
                       alphaFiltered, betaFiltered, rhoFiltered, faciesLog); //Generate these logs.

  int densdim = static_cast<int>(sigmaEOrig.size());
  std::vector<std::vector<FFTGrid *> > density(densdim);
  std::vector<Simbox*>                 volume(densdim);

  for(int i = 0;i<densdim;i++)
    volume[i] = NULL;

  //int nAng = int(log(float(densdim))/log(2.0));
  int nAng = static_cast<int>(noiseScale.size());
  double **G = new double*[nAng];
  for(int i=0;i<nAng;i++)
    G[i] = new double[3];

  //cravaResult->computeG(G);
  cravaResult->computeG_new(G);

  for(int i=0;i<densdim;i++) {
    makeFaciesDens(nfac, sigmaEOrig, useFilter, noVs, alphaFiltered, betaFiltered, rhoFiltered,
                   faciesLog, density[i], &volume[i], i, G, cravaResult, noiseScale);
    if(((modelSettings->getOtherOutputFlag() & IO::ROCK_PHYSICS) > 0) && (i == 0 || i == densdim-1)) {
      Simbox * expVol = createExpVol(volume[i]);
      for(int j=0; j<static_cast<int>(density[i].size()); j++) {
        std::string baseName;
        if(densdim > 1) {
          if(i == 0)
            baseName = "Rock_Physics_Min_Noise_";
          else
            baseName = "Rock_Physics_Max_Noise_";
        }
        else
          baseName = "Rock_Physics_";
        baseName = baseName + modelSettings->getFaciesName(j);
        bool writeSurface = (j==0); //writeSurface is true if j == 0.
        resampleAndWriteDensity(density[i][j], baseName, volume[i], expVol, i, writeSurface);
      }
      delete expVol;
    }
  }

  if(priorFaciesCubes != NULL)
    normalizeCubes(priorFaciesCubes);

  calculateFaciesProb(postAlpha, postBeta, postRho, density, volume,
                      p_undef, priorFacies, priorFaciesCubes, noiseScale, seismicLH);

  for(int i = 0 ; i < densdim ; i++) {
    for(int j = 0 ; j < static_cast<int>(density[i].size()) ; j++)
      delete density[i][j];
    delete volume[i];
  }

  for(int i=0 ; i<nAng ; i++)
    delete [] G[i];
  delete [] G;
}

float FaciesProb::findDensity(float                                       alpha,
                              float                                       beta,
                              float                                       rho,
                              const std::vector<std::vector<FFTGrid*> > & density,
                              int                                         facies,
                              const std::vector<Simbox *>               & volume,
                              const std::vector<float>                  & t,
                              int                                         nAng)
{
  double jFull, kFull, lFull;

  int i;
  int dim = static_cast<int>(density.size());
  std::vector<float> value(dim);
  for(i=0;i<dim;i++)
  {
  volume[i]->getInterpolationIndexes(alpha, beta, rho, jFull, kFull, lFull);
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
  else if(j1>=volume[i]->getnx()-1) {
    j1 = volume[i]->getnx()-1;
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
  else if(k1>=volume[i]->getny()-1) {
    k1 = volume[i]->getny()-1;
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
  else if(l1>=volume[i]->getnz()-1) {
    l1 = volume[i]->getnz()-1;
    l2 = l1;
    wl = 0;
  }
  else {
    l2 = l1 + 1;
    wl = static_cast<float>(lFull-l1);
  }


    density[i][facies]->setAccessMode(FFTGrid::RANDOMACCESS);
    float value1 = std::max<float>(0,density[i][facies]->getRealValue(j1,k1,l1));
    float value2 = std::max<float>(0,density[i][facies]->getRealValue(j1,k1,l2));
    float value3 = std::max<float>(0,density[i][facies]->getRealValue(j1,k2,l1));
    float value4 = std::max<float>(0,density[i][facies]->getRealValue(j1,k2,l2));
    float value5 = std::max<float>(0,density[i][facies]->getRealValue(j2,k1,l1));
    float value6 = std::max<float>(0,density[i][facies]->getRealValue(j2,k1,l2));
    float value7 = std::max<float>(0,density[i][facies]->getRealValue(j2,k2,l1));
    float value8 = std::max<float>(0,density[i][facies]->getRealValue(j2,k2,l2));
    density[i][facies]->endAccess();
    //float value=0;
    value[i] = 0;
    value[i] += (1.0f-wj)*(1.0f-wk)*(1.0f-wl)*value1;
    value[i] += (1.0f-wj)*(1.0f-wk)*(     wl)*value2;
    value[i] += (1.0f-wj)*(     wk)*(1.0f-wl)*value3;
    value[i] += (1.0f-wj)*(     wk)*(     wl)*value4;
    value[i] += (     wj)*(1.0f-wk)*(1.0f-wl)*value5;
    value[i] += (     wj)*(1.0f-wk)*(     wl)*value6;
    value[i] += (     wj)*(     wk)*(1.0f-wl)*value7;
    value[i] += (     wj)*(     wk)*(     wl)*value8;
  }

  int factor = 1;
  float valuesum = 0;
  for(i=0;i<dim;i++)
  {
    factor = 1;
    for(int j=0;j<nAng;j++)
    {
      if(j>0)
        factor*=2;
      if((i & factor) > 0)
        value[i]*=t[j];
      else
        value[i]*=(1-t[j]);
    }

    valuesum += value[i];
  }

  if (valuesum > 0.0)
    return(valuesum);
  else
    return 0.0;
}


void
FaciesProb::resampleAndWriteDensity(const FFTGrid     * density,
                                    const std::string & fileName,
                                    const Simbox      * origVol,
                                    Simbox            * volume,
                                    int                 gridNo,
                                    bool                writeSurface)
{
  if(writeSurface == true) {
    int format = IO::STORM;
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixDensity() + NRLib::ToString(gridNo);
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixDensity() + NRLib::ToString(gridNo);
    volume->writeTopBotGrids(topSurf, baseSurf, IO::PathToInversionResults(), format);
    volume->setTopBotName(topSurf, baseSurf, format);
  }

  int nx = density->getNx();
  int ny = density->getNy();
  int nz = density->getNz();

  FFTGrid expDens(nx, ny, nz, nx, ny, nz);
  expDens.createRealGrid();

  double sum = 0;
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) {
        double alpha, beta, rho;
        volume->getCoord(i, j, k, alpha, beta, rho);
        double aInd, bInd, rInd;
        origVol->getInterpolationIndexes(log(alpha), log(beta), log(rho), aInd, bInd, rInd);
        double ti = aInd-floor(aInd);
        double tj = bInd-floor(bInd);
        double tk = rInd-floor(rInd);
        int li  = static_cast<int>(floor(aInd));
        int li2 = (li == nx-1) ? li : li+1;
        int lj = static_cast<int>(floor(bInd));
        int lj2 = (lj == ny-1) ? lj : lj+1;
        int lk = static_cast<int>(floor(rInd));
        int lk2 = (lk == nz-1) ? lk : lk+1;
        double tmpFrontLeft  = density->getRealValue(li, lj, lk)*(1-tk)+density->getRealValue(li, lj, lk2)*tk;
        double tmpFrontRight = density->getRealValue(li2, lj, lk)*(1-tk)+density->getRealValue(li2, lj, lk2)*tk;
        double tmpBackLeft   = density->getRealValue(li, lj2, lk)*(1-tk)+density->getRealValue(li, lj2, lk2)*tk;
        double tmpBackRight  = density->getRealValue(li2, lj2, lk)*(1-tk)+density->getRealValue(li2, lj2, lk2)*tk;
        double tmpLeft  = tmpFrontLeft*(1-tj)+tmpBackLeft*tj;
        double tmpRight = tmpFrontRight*(1-tj)+tmpBackRight*tj;
        double value = (tmpLeft*(1-ti)+tmpRight*ti)/alpha/beta/rho;
        expDens.setRealValue(i, j, k, static_cast<float>(value));
        sum += value;
      }
    }
  }
  expDens.multiplyByScalar(static_cast<float>(1.0/sum));
  expDens.writeFile(fileName, "", volume);
}


Simbox *
FaciesProb::createExpVol(const Simbox * volume)
{
  double minA = exp(volume->getx0());
  double minB = exp(volume->gety0());
  double maxA = exp(volume->getx0()+volume->getlx());
  double maxB = exp(volume->gety0()+volume->getly());
  double minLR, maxLR;
  volume->getMinMaxZ(minLR, maxLR);
  double minR = exp(minLR);
  double maxR = exp(maxLR);

  double dA = (maxA-minA)/volume->getnx();
  double dB = (maxB-minB)/volume->getny();
  double dR = (maxR-minR)/volume->getnz();

  Surface rhoMinSurf(minA, minB, maxA-minA, maxB-minB, 2, 2, minR);
  Simbox  * expVol = new Simbox(minA, minB, rhoMinSurf, maxA-minA, maxB-minB, maxR-minR, 0, dA, dB, dR);
  return(expVol);
}


void FaciesProb::calculateConditionalFaciesProb(WellData                      ** wells,
                                                int                              nWells,
                                                const std::vector<Surface *>   & faciesEstimInterval,
                                                const std::vector<std::string> & faciesNames,
                                                const double                     dz)
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
  int nActiveWells = count;

  //
  // Put all blocked facies logs in one vector
  //
  int ** BWfacies        = new int * [nActiveWells];
  int ** faciesCountWell = new int * [nActiveWells];
  int *  faciesCount     = new int[nFacies_];
  for(int f=0; f < nFacies_; f++)
    faciesCount[f] = 0;
  for (int i = 0 ; i < nActiveWells ; i++)
  {
    const int   nBlocks    = bw[i]->getNumberOfBlocks();
    const int * BWfacies_i = bw[i]->getFacies();
    BWfacies[i] = new int[nBlocks];
    faciesCountWell[i] = new int[nFacies_];
    for(int f=0; f < nFacies_; f++)
      faciesCountWell[i][f] = 0;

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
        //LogKit::LogFormatted(LogKit::Low,"ib, BWfacies[i][b] = %d %d\n",ib,BWfacies[ib]);
      }
    }
    else {
      for (int b = 0 ; b < nBlocks ; b++) {
        BWfacies[i][b] = BWfacies_i[b];
        //LogKit::LogFormatted(LogKit::Low,"ib, BWfacies[i][b] = %d %d\n",ib,BWfacies[ib]);
      }
    }
    for(int b = 0; b < nBlocks; b++) {
      for(int f=0; f < nFacies_; f++){
        if(BWfacies[i][b] == f){
          faciesCountWell[i][f]++;
          faciesCount[f]++;
        }
      }
    }
  }

  //
  // Block facies probabilities and put them in one vector
  //
  float *** BWfaciesProb = new float ** [nFacies_];
  for(int f = 0 ; f < nFacies_ ; f++)
  {
    faciesProb_[f]->setAccessMode(FFTGrid::RANDOMACCESS);
    BWfaciesProb[f] = new float * [nActiveWells];
    for (int i = 0 ; i < nActiveWells ; i++)
    {
      const int   nBlocks = bw[i]->getNumberOfBlocks();
      const int * ipos    = bw[i]->getIpos();
      const int * jpos    = bw[i]->getJpos();
      const int * kpos    = bw[i]->getKpos();
      BWfaciesProb[f][i] = new float[nBlocks];
      for(int b = 0 ; b < nBlocks ; b++)
      {
        BWfaciesProb[f][i][b] = faciesProb_[f]->getRealValue(ipos[b],jpos[b],kpos[b]);
        //LogKit::LogFormatted(LogKit::Low,"f,ib, BWfaciesProb[f][i][b] = %d %d %.5f\n",f,ib,BWfaciesProb[f][ib]);
      }
    }
    faciesProb_[f]->endAccess();
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
      sumProb[f1][f2] = new float[nActiveWells];
      numProb[f1][f2] = new int[nActiveWells];
      for (int i = 0 ; i < nActiveWells ; i++)
      {
        sumProb[f1][f2][i] = 0.0f;
        numProb[f1][f2][i] = 0;
        for(int b = 0 ; b < bw[i]->getNumberOfBlocks() ; b++)
        {
          if (BWfacies[i][b] == f1) {
            //LogKit::LogFormatted(LogKit::Low,"f1,f2 = %d,%d   ib = %d    BWfacies[i][b] = %d   sumProb = %.5f\n",f1,f2,ib,BWfacies[i][b],sumProb);
            sumProb[f1][f2][i] += BWfaciesProb[f2][i][b];
            numProb[f1][f2][i] += 1;
          }
        }
      }
    }
  }
  for (int i = 0 ; i < nActiveWells ; i++)
    delete [] BWfacies[i];
  delete [] BWfacies;

  for(int f = 0 ; f < nFacies_ ; f++) {
    for (int i = 0 ; i < nActiveWells ; i++)
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
  bool low_probabilities = false;
  for (int i = 0 ; i < nActiveWells ; i++)
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

    LogKit::LogFormatted(LogKit::Low,"\nWell: "+bw[i]->getWellname()+"\n");
    LogKit::LogFormatted(LogKit::Low,"\nFacies      |");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::Low," %11s",faciesNames[f].c_str());
    LogKit::LogFormatted(LogKit::Low,"\n------------+");
    for(int f=0 ; f < nFacies_ ; f++)
      LogKit::LogFormatted(LogKit::Low,"------------",f);
    LogKit::LogFormatted(LogKit::Low,"\n");
    for(int f1=0 ; f1 < nFacies_ ; f1++)
    {
      LogKit::LogFormatted(LogKit::Low,"%-11s |",faciesNames[f1].c_str());
      for(int f2=0 ; f2 < nFacies_ ; f2++)
      {
        LogKit::LogFormatted(LogKit::Low," %11.3f",condFaciesProb[f2][f1]);
      }
      LogKit::LogFormatted(LogKit::Low,"\n");
    }
    LogKit::LogFormatted(LogKit::Low,"Total Prob  |");
    for(int f=0 ; f < nFacies_ ; f++)
    {
      LogKit::LogFormatted(LogKit::Low," %11.3f",totProb[f]);
    }
    LogKit::LogFormatted(LogKit::Low,"\n");

    bool low_probabilities_this_well = false;
    checkConditionalProbabilities(condFaciesProb, faciesNames, nFacies_, bw[i]->getWellname(), false, low_probabilities_this_well, faciesCountWell[i]);
    low_probabilities = low_probabilities || low_probabilities_this_well;
  }

  for (int i = 0 ; i < nActiveWells ; i++)
    delete [] faciesCountWell[i];
  delete [] faciesCountWell;

  //if (low_probabilities) {
  //   std::string text;
  //   text  = "Low facies probabilities have been detected for one or more wells. Check the conditional\n";
  //   text += "   facies probabilities for all wells.\n";
  //   TaskList::addTask(text);
  // }

  //
  // Estimate P( facies2 | facies1 )
  //
  LogKit::LogFormatted(LogKit::High,"\nThe table below gives the mean conditional probability of finding one of");
  LogKit::LogFormatted(LogKit::High,"\nthe facies specified in the left column when one of the facies specified");
  LogKit::LogFormatted(LogKit::High,"\nin the top row are observed in well logs ==> P(A|B)\n");
  for(int f1 = 0 ; f1 < nFacies_ ; f1++)
  {
    totProb[f1] = 0.0f;
    for(int f2 = 0 ; f2 < nFacies_ ; f2++)
    {
      float sumFaciesProb = 0.0f;
      int   numFaciesProb =   0;
      for (int i = 0 ; i < nActiveWells ; i++)
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

  LogKit::LogFormatted(LogKit::Low,"\nFor all wells:\n");
  LogKit::LogFormatted(LogKit::Low,"\nFacies      |");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::Low," %11s",faciesNames[f].c_str());
  LogKit::LogFormatted(LogKit::Low,"\n------------+");
  for(int f=0 ; f < nFacies_ ; f++)
    LogKit::LogFormatted(LogKit::Low,"------------",f);
  LogKit::LogFormatted(LogKit::Low,"\n");
  for(int f1=0 ; f1 < nFacies_ ; f1++)
  {
    LogKit::LogFormatted(LogKit::Low,"%-11s |",faciesNames[f1].c_str());
    for(int f2=0 ; f2 < nFacies_ ; f2++)
      LogKit::LogFormatted(LogKit::Low," %11.3f",condFaciesProb[f2][f1]);
    LogKit::LogFormatted(LogKit::Low,"\n");
  }
  LogKit::LogFormatted(LogKit::Low,"Total Prob  |");
  for(int f=0 ; f < nFacies_ ; f++)
  {
    LogKit::LogFormatted(LogKit::Low," %11.3f",totProb[f]);
  }
  LogKit::LogFormatted(LogKit::Low,"\n");

  low_probabilities = false;
  checkConditionalProbabilities(condFaciesProb, faciesNames, nFacies_, "NOT SET", true, low_probabilities, faciesCount);

  if (low_probabilities) {
    std::string text;
    text  = "For one or several facies, the total probability of finding that facies is less than\n";
    text += "   0.95. This indicates a problem with the estimation. Check each well carefully.\n";
    TaskList::addTask(text);
  }

  delete [] totProb;
  delete [] faciesCount;

  for(int i=0 ; i < nFacies_ ; i++)
    delete [] condFaciesProb[i];
  delete [] condFaciesProb;
  // We cannot delete blocked logs here, only the array of blocked logs.
  delete [] bw;
}

void FaciesProb::checkConditionalProbabilities(float                         ** condFaciesProb,
                                               const std::vector<std::string> & faciesNames,
                                               const int                        nFacies,
                                               const std::string              & identifier,
                                               const bool                       accumulative,
                                               bool                           & lowProbs,
                                               int                            * faciesCount)
{
  bool tooSmallDiff = false;
  bool wrongMaximum = false;

  std::vector<float> totProb(nFacies); // Automatically initialized to zero

  for (int f1=0 ; f1 < nFacies ; f1++) {
    int   indMin  = -IMISSING;
    int   indMax  = -IMISSING;
    float minProb = 1.0;
    float maxProb = 0.0;
    for (int f2=0 ; f2 < nFacies ; f2++) { // Facies in well-log
      float prob = condFaciesProb[f2][f1];
      totProb[f2] += prob;
      if (prob > maxProb + 0.001) {
          maxProb = prob;
          indMax  = f2;
      }
      if (prob < minProb - 0.001) {
        minProb = prob;
        indMin  = f2;
      }
    }

    float fraction = maxProb/minProb;

    if (fraction < 1.05f) {
      tooSmallDiff = true;
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: For facies \'"+faciesNames[f1]+"\' the difference between smallest and largest probability\n");
      LogKit::LogFormatted(LogKit::Warning,"         is only %.2f percent.\n",(fraction - 1.f)*100.f);
    }
    if (indMax != f1 && indMax !=-IMISSING) {
      wrongMaximum = true;
      if (accumulative && faciesCount[f1]> 0) {
        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The probability of finding facies \'"+faciesNames[f1]+"\' is largest when the logs\n");
        LogKit::LogFormatted(LogKit::Warning,"         show \'"+faciesNames[indMax]+"\' and not when it shows \'"+faciesNames[f1]+"\'.\n");
      }
      else if(faciesCount[f1] > 0) {
        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The probability of finding facies \'"+faciesNames[f1]+"\' is largest when the well log of\n");
        LogKit::LogFormatted(LogKit::Warning,"         well \'"+identifier+"\' shows \'"+faciesNames[indMax]+"\' and not when it shows \'"+faciesNames[f1]+"\'.\n");
      }
      else if(faciesCount[f1] == 0){
        wrongMaximum = false;
      }
    }
  }

  for (int f1=0 ; f1 < nFacies ; f1++) { // This loop cannot be joined with loop above
    if (totProb[f1] < 0.95) {
      lowProbs = true;
      if (accumulative && faciesCount[f1] > 0) {
        std::string text;
        text += "\nWARNING: The total probability for facies \'"+faciesNames[f1]+"\' is only "+NRLib::ToString(totProb[f1],3)+". This means";
        text += "\n         that too often when \'"+faciesNames[f1]+"\' is observed in wells, the inverted (Vp,Vs,Rho)";
        text += "\n         do not predict any of the "+NRLib::ToString(nFacies)+" observed facies. Check";
        text += "\n         the inverted elastic parameters against well logs.";
        LogKit::LogFormatted(LogKit::Warning,text);
      }
      else if(faciesCount[f1] > 0) {
        std::string text;
        text += "\nWARNING: The total probability for facies \'"+faciesNames[f1]+"\' in well "+identifier+" is only "+NRLib::ToString(totProb[f1],3)+". This";
        text += "\n         means that some of the places where \'"+faciesNames[f1]+"\' are observed in the well, the";
        text += "\n         inverted (Vp,Vs,Rho) do not predict any of the "+NRLib::ToString(nFacies)+" observed facies. Check";
        text += "\n         the inverted elastic parameters against well logs.";
        LogKit::LogFormatted(LogKit::Warning,text);
      }
      else {
        lowProbs = false;
        std::string text;
        if(accumulative)
          text += "\nWARNING: Facies \'"+faciesNames[f1]+"\' is not observed in any well. Facies probability can not be estimated for this facies.\n";
        //else
        //  text += "\nWARNING: Facies \'"+faciesNames[f1]+"\' is not observed in well "+identifier+". Facies probability can not be estimated in this well.\n";
        LogKit::LogFormatted(LogKit::Warning,text);
      }

    }
  }

  if (tooSmallDiff) {
    std::string text;
    if (accumulative) {
      text += "For one or more facies the conditional facies probabilities do not distinguish between\n";
      text += "   facies (probabilities are too similar). This indicates that you are using \'absolute\'\n";
      text += "   elastic logs for the probability estimation and that the trends are too large. Use\n";
      text += "   \'relative\' elastic logs instead.\n";
    }
    else {
      text += "For one or more facies the conditional facies probabilities do not distinguish between\n";
      text += "   facies (probabilities are too similar). Check well "+identifier+".\n";
    }
      TaskList::addTask(text);
  }
  if (wrongMaximum) {
    std::string text;
    if (accumulative) {
      text += "For one or more facies the conditional facies probabilities are incorrect. The\n";
      text += "   probability of finding, say \'sand\', is not largest where \'sand\' is observed\n";
      text += "   in the well logs. It seemes that you have major alignment problems between\n";
      text += "   elastic parameters and facies in the logs. Check the alignments.\n";
    }
    else {
      text += "For well "+identifier+" the conditional facies probabilities are incorrect. The\n";
      text += "   probability of finding a facies is not largest where that facies is observed\n";
      text += "   in the well log. You possibly have alignment problems between elastic\n";
      text += "   parameters and facies in the logs. Check the well.\n";
    }
    TaskList::addTask(text);
  }
}

void FaciesProb::calculateFaciesProb(FFTGrid                                    * alphagrid,
                                     FFTGrid                                    * betagrid,
                                     FFTGrid                                    * rhogrid,
                                     const std::vector<std::vector<FFTGrid *> > & density,
                                     const std::vector<Simbox *>                  volume,
                                     float                                        p_undefined,
                                     const float                                * priorFacies,
                                     FFTGrid                                   ** priorFaciesCubes,
                                     const std::vector<Grid2D *>                & noiseScale,
                                     FFTGrid                                    * seismicLH)
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
    faciesProb_[i]->createRealGrid(false);
  }
  faciesProbUndef_->setAccessMode(FFTGrid::WRITE);
  faciesProbUndef_->createRealGrid(false);
  if(seismicLH != NULL)
    seismicLH->setAccessMode(FFTGrid::WRITE);

  if(priorFaciesCubes!=NULL)
    for(i=0;i<nFacies_;i++)
      priorFaciesCubes[i]->setAccessMode(FFTGrid::READ);

  smallrnxp = faciesProb_[0]->getRNxp();

  int nAng = 0;
  for(i=0;i<int(noiseScale.size());i++)
    if(noiseScale[i]!=NULL)
      nAng++;
  std::vector<float> t(nAng);
  double maxS;
  double minS;
  std::vector<Grid2D *> tgrid(nAng);
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      tgrid[angle] = new Grid2D(noiseScale[0]->GetNI(), noiseScale[0]->GetNJ());
      for(size_t ii=0;ii<tgrid[angle]->GetNI();ii++)
        for(size_t jj=0;jj<tgrid[angle]->GetNJ();jj++)
          if(minS==maxS)
            (*tgrid[angle])(ii,jj) = 0.0;
          else
          (*tgrid[angle])(ii,jj) = ((*noiseScale[angle])(ii,jj)-minS)/(maxS-minS);
    }

  LogKit::LogFormatted(LogKit::Low,"\nBuilding facies probabilities:");
  float monitorSize = std::max(1.0f, static_cast<float>(nzp)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  float help;
  float dens;
  float undefSum = p_undefined/(volume[0]->getnx()*volume[0]->getny()*volume[0]->getnz());
  for(i=0;i<nzp;i++)
  {
    for(j=0;j<nyp;j++)
    {
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
            if(k<nx)
            {
              for(int angle = 0;angle<nAng;angle++)
                t[angle] = float((*tgrid[angle])(k,j));
              dens = findDensity(alpha, beta, rho, density, l,volume, t, nAng);
            }
            else
              dens = 1.0;
            if(priorFaciesCubes!=NULL)
              value[l] = priorFaciesCubes[l]->getNextReal()*dens;
            else
              value[l] = priorFacies[l]*dens;
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
          if(k<nx) {
            faciesProbUndef_->setNextReal(undefSum/sum);
            if(seismicLH != NULL)
              seismicLH->setNextReal(sum);
          }
          else {
            faciesProbUndef_->setNextReal(RMISSING);
            if(seismicLH != NULL)
              seismicLH->setNextReal(RMISSING);
          }
        }
      }
    }
    // Log progress
    if (i+1 >= static_cast<int>(nextMonitor)) {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  std::cout << "\n";

  if(priorFaciesCubes!=NULL)
    for(i=0;i<nFacies_;i++)
    {
      priorFaciesCubes[i]->endAccess();
    }
  for(i=0;i<nFacies_;i++)
    faciesProb_[i]->endAccess();
  alphagrid->endAccess();
  betagrid->endAccess();
  rhogrid->endAccess();
  faciesProbUndef_->endAccess();
  if(seismicLH != NULL)
    seismicLH->endAccess();

  delete [] value;
}

void FaciesProb::calculateFaciesProbGeomodel(const float  * priorFacies,
                                             FFTGrid     ** priorFaciesCubes)
{
  int i,j,k,l;
  int nx, ny, nz;
  float value;
  nx   = faciesProb_[0]->getNx();
  ny   = faciesProb_[0]->getNy();
  nz   = faciesProb_[0]->getNz();
  float undef;
  for(i=0;i<nFacies_;i++)
    faciesProb_[i]->setAccessMode(FFTGrid::READANDWRITE);
  faciesProbUndef_->setAccessMode(FFTGrid::READ);
  if(priorFaciesCubes!=NULL)
    for(i=0;i<nFacies_;i++)
      priorFaciesCubes[i]->setAccessMode(FFTGrid::READ);
  int smallrnxp = faciesProb_[0]->getRNxp();
  for(i=0;i<nz;i++)
  {
    for(j=0;j<ny;j++)
      for(k=0;k<smallrnxp;k++)
      {
        undef = faciesProbUndef_->getNextReal();
        for(l=0;l<nFacies_;l++)
        {
          if(priorFaciesCubes!=NULL)
            value = priorFaciesCubes[l]->getNextReal()*undef+faciesProb_[l]->getNextReal();
          else
            value = priorFacies[l]*undef+faciesProb_[l]->getNextReal();
          if(k<nx)
            faciesProb_[l]->setNextReal(value);
          else
            faciesProb_[l]->setNextReal(RMISSING);

        }
      }
  }
  for(i=0;i<nFacies_;i++)
    faciesProb_[i]->endAccess();
  faciesProbUndef_->endAccess();
  if(priorFaciesCubes!=NULL)
    for(i=0;i<nFacies_;i++)
      priorFaciesCubes[i]->endAccess();

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
                                      bool                           useFilter,
                                      std::vector<float>           & alphaFiltered,
                                      std::vector<float>           & betaFiltered,
                                      std::vector<float>           & rhoFiltered,
                                      std::vector<int>             & faciesLog)
{
  int nData = 0;

  for(int w=0;w<nWells;w++) {
    if(wells[w]->getUseForFaciesProbabilities())
      nData += wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
  }

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
        float a, b, r;
        float aBg, bBg, rBg;
        int f;
        if(noVs == false) {
          if(useFilter == true) {
            a = bw->getAlphaSeismicResolution()[i];
            b = bw->getBetaSeismicResolution()[i];
            r = bw->getRhoSeismicResolution()[i];
          }
          else {
            a = bw->getAlphaPredicted()[i];
            b = bw->getBetaPredicted()[i];
            r = bw->getRhoPredicted()[i];
          }
        }
        else {
          //Beta log here is mainly dummy, but must be at correct level.
          b = bw->getBetaHighCutBackground()[i];
          if(useFilter == true) {
            a = bw->getAlphaForFacies()[i];
            r = bw->getRhoForFacies()[i];
          }
          else {
            a = bw->getAlphaPredicted()[i];
            r = bw->getRhoPredicted()[i];
          }
        }

        f = bw->getFacies()[i];

        if (faciesEstimInterval.size() > 0) {
          const double * xPos  = bw->getXpos();
          const double * yPos  = bw->getYpos();
          const double * zPos  = bw->getZpos();
          const double   zTop  = faciesEstimInterval[0]->GetZ(xPos[i],yPos[i]);
          const double   zBase = faciesEstimInterval[1]->GetZ(xPos[i],yPos[i]);
          if ( (zPos[i]-0.5*dz) < zTop || (zPos[i]+0.5*dz) > zBase)
            f = IMISSING;
        }

        if (relative == true)
        {
          aBg = bw->getAlphaHighCutBackground()[i];
          bBg = bw->getBetaHighCutBackground()[i];
          rBg = bw->getRhoHighCutBackground()[i];

          if (a != RMISSING &&  aBg != RMISSING &&
              b != RMISSING &&  bBg != RMISSING &&
              r != RMISSING &&  rBg != RMISSING &&
              f != IMISSING)
          {
            alphaFiltered[index] = a - aBg;
            betaFiltered[index]  = b - bBg;
            rhoFiltered[index]   = r - rBg;
            faciesLog[index]     = f;
            index++;
          }
        }
        else
        {
          if (a != RMISSING && b != RMISSING && r != RMISSING && f != IMISSING)
          {
            alphaFiltered[index] = a;
            betaFiltered[index]  = b;
            rhoFiltered[index]   = r;
            faciesLog[index]     = f;
            index++;
          }
        }
      }
    }
  }
  alphaFiltered.resize(index);
  betaFiltered.resize(index);
  rhoFiltered.resize(index);
  faciesLog.resize(index);
}

void FaciesProb::normalizeCubes(FFTGrid **priorFaciesCubes)
{
  int   i,j,k,l;
  float sum;
  int   rnxp     = priorFaciesCubes[0]->getRNxp();
  int   nyp      = priorFaciesCubes[0]->getNyp();
  int   nzp      = priorFaciesCubes[0]->getNzp();
  int   nx       = priorFaciesCubes[0]->getNx();
  int   totZero  = 0;
  int   large    = 0;
  int   small    = 0;
  int   total    = 0;
  float negative = 0;
  bool missing = false;
  std::vector<float> value(nFacies_);
  for(i=0;i<nFacies_;i++)
    priorFaciesCubes[i]->setAccessMode(FFTGrid::READANDWRITE);

  for(i=0;i<nzp;i++)
  {
    for (j=0; j<nyp; j++)
    {
      for (k=0; k<rnxp; k++)
      {
        sum = 0.0;
        for (l=0; l<nFacies_; l++)
        {
          value[l] = priorFaciesCubes[l]->getNextReal();
          if(value[l]==RMISSING)
          {
            value[l] = 0.0;
            if(k<nx)
              missing = true;

          }
          else if(value[l]<0.0)
          {
            if (value[l] < negative && k<nx)
              negative = value[l];
            value[l] = 0.0;
          }
          sum += value[l];
        }
        if (k<nx)
        {
          total ++;
          if (sum == 0)
            totZero ++;
          else if (sum>1.001)
            large ++;
          else if( sum>0 && sum<0.999)
            small ++;
        }
        for (l=0; l<nFacies_; l++)
        {
          if(sum>0.0)
            value[l] = value[l]/sum;
          priorFaciesCubes[l]->setNextReal(value[l]);
        }
      }
    }
  }

  if(missing == true)
  {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Prior facies probability is undefined in one or more cells. The value is set to zero.\n");
    TaskList::addTask("Missing prior facies probabilities detected. This is a serious problem. \n Make sure that the prior facies probability cubes have positive probability everywhere.");
  }
  if (negative<0)
  {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Negative facies probabilities have been detected, and set to zero. \n         The probabilities have been rescaled. The most negative value is %4.2f.\n",negative);
    TaskList::addTask("Negative prior facies probabilities detected. This is a serious problem. \n Make sure that the prior facies probability cubes have positive probability everywhere.");
  }
  if (totZero > 0)
  {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The sum over the prior facies probabilities is zero in %d of %d locations.\n",totZero,total);
    TaskList::addTask("In some cells, the prior facies probabilities sum up to zero. This is a serious problem. \n Make sure that the prior facies probability cubes are defined everywhere.");
  }
  if (small > 0)
  {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The sum over the prior facies probabilities is less than one in %d of %d locations. \n         The probabilities have been rescaled.\n",small,total);
    TaskList::addTask("The prior facies probability cubes does not sum to one everywhere. This is a serious problem. ");
  }
  if (large > 0)
  {
    LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The sum over the prior facies probabilities is larger than one in %d of %d locations. \n         The probabilities have been rescaled.\n",large, total);
    TaskList::addTask("The sum over prior facies probabilities is larger than one in some locations. \n Make sure that the prior facies probability cubes sum to one everywhere.");
  }
  for(i=0;i<nFacies_;i++)
    priorFaciesCubes[i]->endAccess();
}

std::vector<double> FaciesProb::calculateChiSquareTest(WellData                    ** wells,
                                                       int                            nWells,
                                                       const std::vector<Surface *> & faciesEstimInterval)
{
  int    i, j, k;
  int    count;
  int    df;
  int    thisFacies = IMISSING;
  int    nActualFacies;
  float  sum;
  double chi_i;
  double chi;

  std::vector<double>      pValue(nWells);
  std::vector<float>       prob(nFacies_);
  std::vector<std::string> fit(nWells);

  for (i=0; i<nWells; i++)
  {
    BlockedLogs  * bw        = wells[i]->getBlockedLogsOrigThick();
    const int      nBlocks   = bw->getNumberOfBlocks();
    const int    * BWfacies  = bw->getFacies();
    const int    * ipos      = bw->getIpos();
    const int    * jpos      = bw->getJpos();
    const int    * kpos      = bw->getKpos();
    const double   dz        = bw->getDz();
    const double * xPos      = bw->getXpos();
    const double * yPos      = bw->getYpos();
    const double * zPos      = bw->getZpos();

    df    = 0;
    chi   = 0;
    count = 0;

    for (j=0; j<nFacies_; j++)
      prob[j] = 0;

    // Find the first facies in well
    for (j=0; j<nBlocks; j++)
    {
      if(BWfacies[j] != IMISSING)
      {
        thisFacies = BWfacies[j];
        break;
      }
    }

    //
    // Block facies probabilities and put them in one vector
    //
    std::vector<std::vector<float> > BWfaciesProb(nBlocks,std::vector<float>(nFacies_));

    for (j = 0 ; j < nFacies_ ; j++)
    {
      faciesProb_[j]->setAccessMode(FFTGrid::RANDOMACCESS);
      for(k = 0 ; k < nBlocks ; k++)
      {
        BWfaciesProb[k][j] = faciesProb_[j]->getRealValue(ipos[k],jpos[k],kpos[k]);
      }
      faciesProb_[j]->endAccess();
    }

    nActualFacies = nFacies_;
    for (j = 0; j<nFacies_; j++)
    {
      sum = 0;
      for (k = 0; k<nBlocks; k++)
        sum += BWfaciesProb[k][j];
      if (sum == 0)
        nActualFacies -= 1;
    }

    for (j=0 ; j < nBlocks ; j++)
    {
      if (faciesEstimInterval.size() > 0)
      {
        const double zTop  = faciesEstimInterval[0]->GetZ(xPos[j],yPos[j]);
        const double zBase = faciesEstimInterval[1]->GetZ(xPos[j],yPos[j]);

        if ( !( (zPos[j]-0.5*dz) < zTop || (zPos[j]+0.5*dz) > zBase || BWfacies[j]==IMISSING) )
        {
          if (BWfacies[j] != thisFacies)
          {
            chi_i = 0;
            for (k=0; k<nFacies_; k++)
            {
              if (prob[k]!=0)
              {
                if (k==thisFacies)
                  chi_i += std::pow(count-prob[k],2)/prob[k];
                else
                  chi_i += std::pow(0-prob[k],2)/prob[k];
              }
            }
            thisFacies = BWfacies[j];
            count = 1;
            chi  += chi_i;
            df   += 1;
            for (k=0; k<nFacies_; k++)
              prob[k] = BWfaciesProb[j][k];
          }
          else
          {
            count += 1;
            for (k=0; k<nFacies_; k++)
              prob[k] += BWfaciesProb[j][k];
          }
        }
      }
      else
      {
        if (BWfacies[j]!=IMISSING)
        {
          if (BWfacies[j] != thisFacies)
          {
            chi_i = 0;
            for (k=0; k<nFacies_; k++)
            {
              if (prob[k]!=0)
              {
                if (k==thisFacies)
                  chi_i += std::pow(count-prob[k],2)/prob[k];
                else
                  chi_i += std::pow(0-prob[k],2)/prob[k];
              }
            }
            thisFacies = BWfacies[j];
            count = 1;
            chi  += chi_i;
            df   += 1;
            for (k=0; k<nFacies_; k++)
              prob[k] = BWfaciesProb[j][k];
          }
          else
          {
            count += 1;
            for (k=0; k<nFacies_; k++)
              prob[k] += BWfaciesProb[j][k];
          }
        }
      }
    }
    // Include last observations
    chi_i = 0;
    for (k=0; k<nFacies_; k++)
    {
      if (prob[k]!=0)
      {
        if (k==thisFacies)
          chi_i += std::pow(count-prob[k],2)/prob[k];
        else
          chi_i += std::pow(0-prob[k],2)/prob[k];
      }
    }
    chi += chi_i;
    chi *= 0.3; //Scale chi to give better fit
    df  += 1;

    if (nActualFacies == 1)
    {
      fit[i] = "not calculated";
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The number of actual facies is one, so estimation makes no sence.\n");
      TaskList::addTask("Check that all facies have probability larger than zero");
    }
    else
    {
      NRLib::ChiSquared chisquared((nActualFacies-1)*df);
      pValue[i] = 1-chisquared.Cdf(chi);

      if (pValue[i] >= 0.05)
        fit[i] = "good";
      if (pValue[i] < 0.05 && pValue[i] >= 0.025)
        fit[i] = "poor";
      if (pValue[i] < 0.025 && pValue[i] >= 0.005)
        fit[i] = "bad";
      if (pValue[i] < 0.005)
        fit[i] = "very bad";
    }
  }

  LogKit::LogFormatted(LogKit::Medium,"\nFit between facies probabilities and facies observed in wells: \n");

  LogKit::LogFormatted(LogKit::Medium,"\nWell                   Fit");
  LogKit::LogFormatted(LogKit::Medium,"\n");
  for (i = 0 ; i < 24+13*nFacies_ ; i++)
    LogKit::LogFormatted(LogKit::Medium,"-");
  LogKit::LogFormatted(LogKit::Medium,"\n");

  for (int w = 0 ; w < nWells ; w++)
  {
    LogKit::LogFormatted(LogKit::Medium,"%-23s",wells[w]->getWellname().c_str());
    LogKit::LogFormatted(LogKit::Medium,fit[w]);
    LogKit::LogFormatted(LogKit::Medium,"\n");
  }
  LogKit::LogFormatted(LogKit::Medium,"\n");

  return pValue;
}

 void FaciesProb::writeBWFaciesProb(WellData ** wells,
                                    int         nWells)
{
  int i, j;
  for (i=0; i<nWells; i++)
  {
    BlockedLogs  * bw = wells[i]->getBlockedLogsOrigThick();
    for (j=0; j<nFacies_; j++)
    {
      faciesProb_[j]->setAccessMode(FFTGrid::RANDOMACCESS);
      bw->setLogFromGrid(faciesProb_[j],j,nFacies_,"FACIES_PROB");
      faciesProb_[j]->endAccess();
    }
   }
}

 void FaciesProb::calculateChiSquareAlternativeTest(WellData                    ** wells,
                                                    int                            nWells,
                                                    const std::vector<Surface *> & faciesEstimInterval,
                                                    const ModelSettings          * modelSettings)
{
  // Alternative approach for evaluating facies fit.
  // Is not being called, but is saved for later by approval from Ragnar.

  int    i, j, k;
  int    df;
  int    bin;
  double pValue;
  double chi;

  std::vector<float>  prob(nFacies_);
  std::vector<int>    nPos(10);
  std::vector<int>    nObs(10);

  std::vector<std::vector<std::string> > fit(nWells,std::vector<std::string>(nFacies_));

  for (i=0; i<nWells; i++)
  {
    BlockedLogs  * bw       = wells[i]->getBlockedLogsOrigThick();
    const int      nBlocks  = bw->getNumberOfBlocks();
    const int    * BWfacies = bw->getFacies();
    const int    * ipos     = bw->getIpos();
    const int    * jpos     = bw->getJpos();
    const int    * kpos     = bw->getKpos();
    const double   dz       = bw->getDz();
    const double * xPos     = bw->getXpos();
    const double * yPos     = bw->getYpos();
    const double * zPos     = bw->getZpos();


    std::vector<std::vector<float> > BWfaciesProb(nBlocks,std::vector<float>(nFacies_));
    for (j = 0 ; j < nFacies_ ; j++)
    {
      faciesProb_[j]->setAccessMode(FFTGrid::RANDOMACCESS);
      for(k = 0 ; k < nBlocks ; k++)
      {
        BWfaciesProb[k][j] = faciesProb_[j]->getRealValue(ipos[k],jpos[k],kpos[k]);
      }
      faciesProb_[j]->endAccess();
    }

    for (k=0; k < nFacies_; k++)
    {
      for (j=0; j<10; j++)
      {
        nPos[j] = 0;
        nObs[j] = 0;
      }

      for (j=0 ; j < nBlocks ; j++)
      {
        if (faciesEstimInterval.size() > 0)
        {
          const double zTop  = faciesEstimInterval[0]->GetZ(xPos[j],yPos[j]);
          const double zBase = faciesEstimInterval[1]->GetZ(xPos[j],yPos[j]);
          if ( !( (zPos[j]-0.5*dz) < zTop || (zPos[j]+0.5*dz) > zBase || BWfacies[j]==IMISSING) )
          {
            bin = static_cast<int>(std::floor(BWfaciesProb[j][k]*10));
            nObs[bin] += 1;
            if (BWfacies[j] == k)
              nPos[bin] += 1;
          }
        }
        else
        {
          if (BWfacies[j]!=IMISSING)
          {
            bin = static_cast<int>(std::floor(BWfaciesProb[j][k]*10));
            nObs[bin] += 1;
            if (BWfacies[j] == k)
              nPos[bin] += 1;
          }
        }
      }
      chi = 0;
      df  = 0;
      for (j=0; j<10; j++)
      {
        if (nObs[j] > 0)
        {
          df  += 1;
          chi += std::pow(nPos[j]        -nObs[j]   *(j*0.1+0.05) ,2)/(nObs[j]*   (j*0.1+0.05));
          chi += std::pow(nObs[j]-nPos[j]-nObs[j]*(1-(j*0.1+0.05)),2)/(nObs[j]*(1-(j*0.1+0.05)));
        }
      }
      NRLib::ChiSquared chisquared(df);
      pValue = 1-chisquared.Cdf(chi);

      if (pValue >= 0.05)
        fit[i][k] = "good";
      if (pValue < 0.05 && pValue >= 0.025)
        fit[i][k] = "poor";
      if (pValue < 0.025 && pValue >= 0.005)
        fit[i][k] = "bad";
      if (pValue < 0.005)
        fit[i][k] = "very bad";
    }
  }

  LogKit::LogFormatted(LogKit::Medium,"\nFit between facies probabilities and facies observed in wells: \n");
  LogKit::LogFormatted(LogKit::Medium,"\nWell                      ");
  for (int i = 0 ; i < nFacies_ ; i++)
    LogKit::LogFormatted(LogKit::Medium,"%12s ",modelSettings->getFaciesName(i).c_str());
  LogKit::LogFormatted(LogKit::Medium,"\n");
  for (int i = 0 ; i < 24+13*nFacies_ ; i++)
    LogKit::LogFormatted(LogKit::Medium,"-");
  LogKit::LogFormatted(LogKit::Medium,"\n");
  for (int w = 0 ; w < nWells ; w++)
  {
    LogKit::LogFormatted(LogKit::Medium,"%-24s ",wells[w]->getWellname().c_str());
    for (int i = 0 ; i < nFacies_ ; i++)
    {
      LogKit::LogFormatted(LogKit::Medium,"        ");
      LogKit::LogFormatted(LogKit::Medium,fit[w][i]);
    }
      LogKit::LogFormatted(LogKit::Medium,"\n");
  }
  LogKit::LogFormatted(LogKit::Medium,"\n");
}


FFTGrid *
FaciesProb::createLHCube(FFTGrid     * likelihood,
                         int           fac,
                         const float * priorFacies,
                         FFTGrid    ** priorFaciesCubes)
{
  int nx = likelihood->getNx();
  int ny = likelihood->getNy();
  int nz = likelihood->getNz();
  FFTGrid * result = new FFTGrid(nx, ny, nz, nx, ny, nz);
  result->createRealGrid(false);
  result->setAccessMode(FFTGrid::WRITE);
  faciesProb_[fac]->setAccessMode(FFTGrid::READ);
  if(priorFaciesCubes != NULL)
    priorFaciesCubes[fac]->setAccessMode(FFTGrid::READ);

  int i,j,k, rnx = likelihood->getRNxp();
  for(i=0;i<nz;i++) {
    for(j=0;j<ny;j++) {
      for(k=0;k<rnx;k++) {
        float prob  = faciesProb_[fac]->getNextReal();
        float prior;
        if(priorFaciesCubes != NULL)
          prior = priorFaciesCubes[fac]->getNextReal();
        else
          prior = priorFacies[fac];
        float lh = likelihood->getNextReal();
        if(k<nx) {
          float condLH = lh*prob/prior;
          result->setNextReal(condLH);
        }
        else
          result->setNextReal(RMISSING);
      }
    }
  }
  result->endAccess();
  likelihood->endAccess();
  faciesProb_[fac]->endAccess();
  if(priorFaciesCubes != NULL)
    priorFaciesCubes[fac]->endAccess();

  return(result);
}
