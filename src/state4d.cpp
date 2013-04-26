#include "src/state4d.h"
#include "src/seismicparametersholder.h"
#include "src/timeevolution.h"
#include "src/simbox.h"
#include "lib/lib_matr.h"
#include "src/correlatedrocksamples.h"
#include "rplib/distributionsrock.h"
#include "src/rockphysicsinversion4d.h"
#include <string>
#include "src/vario.h"

State4D::State4D()
{
}

State4D::~State4D()
{
}

void State4D::setDynamicMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho)
{
  mu_dynamic_.resize(3);

  mu_dynamic_[0] = vp;
  mu_dynamic_[1] = vs;
  mu_dynamic_[2] = rho;
}

void State4D::setStaticMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho)
{
  mu_static_.resize(3);

  mu_static_[0] = vp;
  mu_static_[1] = vs;
  mu_static_[2] = rho;
}

void State4D::setStaticSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho)
{
  sigma_static_static_.resize(6);

  sigma_static_static_[0] = vpvp;
  sigma_static_static_[1] = vpvs;
  sigma_static_static_[2] = vprho;
  sigma_static_static_[3] = vsvs;
  sigma_static_static_[4] = vsrho;
  sigma_static_static_[5] = rhorho;
}

void State4D::setDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho)
{
  sigma_dynamic_dynamic_.resize(6);

  sigma_dynamic_dynamic_[0] = vpvp;
  sigma_dynamic_dynamic_[1] = vpvs;
  sigma_dynamic_dynamic_[2] = vprho;
  sigma_dynamic_dynamic_[3] = vsvs;
  sigma_dynamic_dynamic_[4] = vsrho;
  sigma_dynamic_dynamic_[5] = rhorho;
}

void State4D::setStaticDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvp, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhovp, FFTGrid *rhovs, FFTGrid *rhorho)
{
  sigma_static_dynamic_.resize(9);

  sigma_static_dynamic_[0] = vpvp;
  sigma_static_dynamic_[1] = vpvs;
  sigma_static_dynamic_[2] = vprho;
  sigma_static_dynamic_[3] = vsvp;
  sigma_static_dynamic_[4] = vsvs;
  sigma_static_dynamic_[5] = vsrho;
  sigma_static_dynamic_[6] = rhovp;
  sigma_static_dynamic_[7] = rhovs;
  sigma_static_dynamic_[8] = rhorho;
}

void State4D::deleteCovariances()
{
  for(int i=0;i<6;i++)
    delete sigma_static_static_[i];

  for(int i=0;i<6;i++)
    delete sigma_dynamic_dynamic_[i];

  for(int i=0;i<9;i++)
    delete sigma_static_dynamic_[i];
}


void State4D::merge(SeismicParametersHolder & current_state )
{
  LogKit::WriteHeader("Merging  static and dynamic part into the distribution of the current state");
  // We assume FFT transformed grids
  assert(allGridsAreTransformed());

  //First, merge expectations.
  std::vector<FFTGrid *> mu(3);
  mu[0] =  current_state.GetMuAlpha(); //mu_Alpha
  mu[1] =  current_state.GetMuBeta(); //mu_Beta
  mu[2] =  current_state.GetMuRho(); //mu_Rho

  for(int i = 0; i<3; i++)
  {
    mu[i]->setTransformedStatus(true); //Going to fill it with transformed info.
    mu[i]->setAccessMode(FFTGrid::WRITE);
    mu_static_[i]->setAccessMode(FFTGrid::READ);
    mu_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  int nzp = mu[0]->getNzp();
  int nyp = mu[0]->getNyp();
  int cnxp = mu[0]->getCNxp();

  fftw_complex*  muFullPrior=new fftw_complex[6];
  fftw_complex*  muCurrentPrior=new fftw_complex[3];

  for (int k = 0; k < nzp; k++) {
    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {
         muFullPrior[0] = mu_static_[0]->getNextComplex();
         muFullPrior[1] = mu_static_[1]->getNextComplex();
         muFullPrior[2] = mu_static_[2]->getNextComplex();
         muFullPrior[3] = mu_dynamic_[0]->getNextComplex();
         muFullPrior[4] = mu_dynamic_[1]->getNextComplex();
         muFullPrior[5] = mu_dynamic_[2]->getNextComplex();
         for(int l=0;l<3;l++){
           muCurrentPrior[l].re =muFullPrior[l].re+muFullPrior[l+3].re;
           muCurrentPrior[l].im =muFullPrior[l].im+muFullPrior[l+3].im;
         }

         mu[0]->setNextComplex(muCurrentPrior[0]);
         mu[1]->setNextComplex(muCurrentPrior[1]);
         mu[2]->setNextComplex(muCurrentPrior[2]);
      }
    }
  }

  for(int i = 0; i<3; i++)
  {
    mu[i]->endAccess();
    mu_static_[i]->endAccess();
    mu_dynamic_[i]->endAccess();
  }

  delete [] muFullPrior;
  delete [] muCurrentPrior;

  //Merge covariances
  std::vector<FFTGrid *> sigma(6);
  sigma[0]=current_state.GetCovAlpha();
  sigma[1]=current_state.GetCrCovAlphaBeta();
  sigma[2]=current_state.GetCrCovAlphaRho();
  sigma[3]=current_state.GetCovBeta();
  sigma[4]=current_state.GetCrCovBetaRho();
  sigma[5]=current_state.GetCovRho();

  mergeCov(sigma);
}


void State4D::mergeCov(std::vector<FFTGrid * > & sigma)
{
  assert(sigma.size() == 6);

  for(int i = 0; i<6; i++)
  {
    sigma[i]->setTransformedStatus(true); //Going to fill it with transformed info.
    sigma[i]->setAccessMode(FFTGrid::WRITE);
    sigma_static_static_[i]->setAccessMode(FFTGrid::READ);
    sigma_dynamic_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->setAccessMode(FFTGrid::READ);

  std::vector<std::vector<fftw_complex> > sigmaFullPrior(6);
  for(int i=0;i<6;i++)
    sigmaFullPrior[i].resize(6);

  std::vector<std::vector<fftw_complex> > sigmaCurrentPrior(3);
  for(int i=0;i<3;i++)
    sigmaCurrentPrior[i].resize(3);

  int nzp = sigma[0]->getNzp();
  int nyp = sigma[0]->getNyp();
  int cnxp = sigma[0]->getCNxp();

  for (int k = 0; k < nzp; k++) {
    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {
         sigmaFullPrior[0][0] = sigma_static_static_[0]->getNextComplex();
         sigmaFullPrior[0][1] = sigma_static_static_[1]->getNextComplex();
         sigmaFullPrior[0][2] = sigma_static_static_[2]->getNextComplex();
         sigmaFullPrior[0][3] = sigma_static_dynamic_[0]->getNextComplex();
         sigmaFullPrior[0][4] = sigma_static_dynamic_[1]->getNextComplex();
         sigmaFullPrior[0][5] = sigma_static_dynamic_[2]->getNextComplex();
         sigmaFullPrior[1][1] = sigma_static_static_[3]->getNextComplex();
         sigmaFullPrior[1][2] = sigma_static_static_[4]->getNextComplex();
         sigmaFullPrior[1][3] = sigma_static_dynamic_[3]->getNextComplex();
         sigmaFullPrior[1][4] = sigma_static_dynamic_[4]->getNextComplex();
         sigmaFullPrior[1][5] = sigma_static_dynamic_[5]->getNextComplex();
         sigmaFullPrior[2][2] = sigma_static_static_[5]->getNextComplex();
         sigmaFullPrior[2][3] = sigma_static_dynamic_[6]->getNextComplex();
         sigmaFullPrior[2][4] = sigma_static_dynamic_[7]->getNextComplex();
         sigmaFullPrior[2][5] = sigma_static_dynamic_[8]->getNextComplex();
         sigmaFullPrior[3][3] = sigma_dynamic_dynamic_[0]->getNextComplex();
         sigmaFullPrior[3][4] = sigma_dynamic_dynamic_[1]->getNextComplex();
         sigmaFullPrior[3][5] = sigma_dynamic_dynamic_[2]->getNextComplex();
         sigmaFullPrior[4][4] = sigma_dynamic_dynamic_[3]->getNextComplex();
         sigmaFullPrior[4][5] = sigma_dynamic_dynamic_[4]->getNextComplex();
         sigmaFullPrior[5][5] = sigma_dynamic_dynamic_[5]->getNextComplex();
         for(int l=0;l<6;l++)
           for(int m=l+1;m<6;m++)
           {
             sigmaFullPrior[m][l].re=sigmaFullPrior[l][m].re;
             sigmaFullPrior[m][l].im=-sigmaFullPrior[l][m].im;
           }

         for(int l=0;l<3;l++)
           for(int m=0;m<3;m++){
             sigmaCurrentPrior[l][m].re  = sigmaFullPrior[l  ][m  ].re;
             sigmaCurrentPrior[l][m].re += sigmaFullPrior[l+3][m+3].re;
             sigmaCurrentPrior[l][m].re += sigmaFullPrior[l+3][m  ].re;
             sigmaCurrentPrior[l][m].re += sigmaFullPrior[l  ][m+3].re;
             sigmaCurrentPrior[l][m].im  = sigmaFullPrior[l  ][m  ].im;
             sigmaCurrentPrior[l][m].im += sigmaFullPrior[l+3][m+3].im;
             sigmaCurrentPrior[l][m].im += sigmaFullPrior[l+3][m  ].im;
             sigmaCurrentPrior[l][m].im += sigmaFullPrior[l  ][m+3].im;
           }
         sigma[0]->setNextComplex(sigmaCurrentPrior[0][0]);
         sigma[1]->setNextComplex(sigmaCurrentPrior[0][1]);
         sigma[2]->setNextComplex(sigmaCurrentPrior[0][2]);
         sigma[3]->setNextComplex(sigmaCurrentPrior[1][1]);
         sigma[4]->setNextComplex(sigmaCurrentPrior[1][2]);
         sigma[5]->setNextComplex(sigmaCurrentPrior[2][2]);
      }
    }
  }

  for(int i = 0; i<6; i++)
  {
    sigma[i]->endAccess();
    sigma_static_static_[i]->endAccess();
    sigma_dynamic_dynamic_[i]->endAccess();
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->endAccess();
}

void State4D::split(SeismicParametersHolder & current_state )
{

  LogKit::WriteHeader("Splitting inversion result into static and dynamic part");
  // initializing
  assert(allGridsAreTransformed());

  std::vector<FFTGrid *> mu(3);
  mu[0] =  current_state.GetMuAlpha(); //mu_Alpha
  mu[1] =  current_state.GetMuBeta(); //mu_Beta
  mu[2] =  current_state.GetMuRho(); //mu_Rho

  for(int i = 0; i<3; i++)
  {
    assert(mu[i]->getIsTransformed());
    mu[i]->setAccessMode(FFTGrid::READ);
    mu_static_[i]->setAccessMode(FFTGrid::READANDWRITE);
    mu_dynamic_[i]->setAccessMode(FFTGrid::READANDWRITE);
  }

  std::vector<FFTGrid *> sigma(6);
  sigma[0]=current_state.GetCovAlpha();
  sigma[1]=current_state.GetCrCovAlphaBeta();
  sigma[2]=current_state.GetCrCovAlphaRho();
  sigma[3]=current_state.GetCovBeta();
  sigma[4]=current_state.GetCrCovBetaRho();
  sigma[5]=current_state.GetCovRho();



  for(int i = 0; i<6; i++)
  {
    assert(sigma[i]->getIsTransformed());
    sigma[i]->setAccessMode(FFTGrid::READ);
    sigma_static_static_[i]->setAccessMode(FFTGrid::READANDWRITE);
    sigma_dynamic_dynamic_[i]->setAccessMode(FFTGrid::READANDWRITE);
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->setAccessMode(FFTGrid::READANDWRITE);


  int nzp = mu[0]->getNzp();
  int nyp = mu[0]->getNyp();
  int cnxp = mu[0]->getCNxp();

  fftw_complex*  muFullPrior=new fftw_complex[6];
  fftw_complex*  muFullPosterior=new fftw_complex[6];
  fftw_complex*  muCurrentPrior=new fftw_complex[3];
  fftw_complex*  muCurrentPosterior=new fftw_complex[3];

  fftw_complex** sigmaFullPrior          = new fftw_complex*[6];
  fftw_complex** sigmaFullPosterior      = new fftw_complex*[6];
  fftw_complex** sigmaFullVsCurrentPrior = new fftw_complex*[6];
  fftw_complex** adjointSandwich         = new fftw_complex*[6];

  for(int i=0;i<6;i++)
  {
    sigmaFullPrior[i]          = new fftw_complex[6];
    sigmaFullVsCurrentPrior[i] = new fftw_complex[3];
    sigmaFullPosterior[i]      = new fftw_complex[6];
    adjointSandwich[i]          = new fftw_complex[3];
  }

  fftw_complex** sigmaCurrentPrior       = new fftw_complex*[3];
  fftw_complex** sigmaCurrentPriorChol   = new fftw_complex*[3];
  fftw_complex** sigmaCurrentPosterior   = new fftw_complex*[3];
  fftw_complex** sandwich                = new fftw_complex*[3];
  fftw_complex** helper                  = new fftw_complex*[3];

  for(int i=0;i<3;i++)
  {
    sigmaCurrentPrior[i]     = new fftw_complex[3];
    sigmaCurrentPriorChol[i] = new fftw_complex[3];
    sigmaCurrentPosterior[i] = new fftw_complex[3];
    sandwich[i]              = new fftw_complex[6];
    helper[i]                = new fftw_complex[6];
  }
  int counter =0;
  for (int k = 0; k < nzp; k++) {
    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {
         // reading from grids
         muFullPrior[0] = mu_static_[0]->getNextComplex();
         muFullPrior[1] = mu_static_[1]->getNextComplex();
         muFullPrior[2] = mu_static_[2]->getNextComplex();
         muFullPrior[3] = mu_dynamic_[0]->getNextComplex();
         muFullPrior[4] = mu_dynamic_[1]->getNextComplex();
         muFullPrior[5] = mu_dynamic_[2]->getNextComplex();


         sigmaFullPrior[0][0] = sigma_static_static_[0]->getNextComplex();
         sigmaFullPrior[0][1] = sigma_static_static_[1]->getNextComplex();
         sigmaFullPrior[0][2] = sigma_static_static_[2]->getNextComplex();
         sigmaFullPrior[0][3] = sigma_static_dynamic_[0]->getNextComplex();
         sigmaFullPrior[0][4] = sigma_static_dynamic_[1]->getNextComplex();
         sigmaFullPrior[0][5] = sigma_static_dynamic_[2]->getNextComplex();
         sigmaFullPrior[1][1] = sigma_static_static_[3]->getNextComplex();
         sigmaFullPrior[1][2] = sigma_static_static_[4]->getNextComplex();
         sigmaFullPrior[1][3] = sigma_static_dynamic_[3]->getNextComplex();
         sigmaFullPrior[1][4] = sigma_static_dynamic_[4]->getNextComplex();
         sigmaFullPrior[1][5] = sigma_static_dynamic_[5]->getNextComplex();
         sigmaFullPrior[2][2] = sigma_static_static_[5]->getNextComplex();
         sigmaFullPrior[2][3] = sigma_static_dynamic_[6]->getNextComplex();
         sigmaFullPrior[2][4] = sigma_static_dynamic_[7]->getNextComplex();
         sigmaFullPrior[2][5] = sigma_static_dynamic_[8]->getNextComplex();
         sigmaFullPrior[3][3] = sigma_dynamic_dynamic_[0]->getNextComplex();
         sigmaFullPrior[3][4] = sigma_dynamic_dynamic_[1]->getNextComplex();
         sigmaFullPrior[3][5] = sigma_dynamic_dynamic_[2]->getNextComplex();
         sigmaFullPrior[4][4] = sigma_dynamic_dynamic_[3]->getNextComplex();
         sigmaFullPrior[4][5] = sigma_dynamic_dynamic_[4]->getNextComplex();
         sigmaFullPrior[5][5] = sigma_dynamic_dynamic_[5]->getNextComplex();

         muCurrentPosterior[0]=mu[0]->getNextComplex();
         muCurrentPosterior[1]=mu[1]->getNextComplex();
         muCurrentPosterior[2]=mu[2]->getNextComplex();

         sigmaCurrentPosterior[0][0]=sigma[0]->getNextComplex();
         sigmaCurrentPosterior[0][1]=sigma[1]->getNextComplex();
         sigmaCurrentPosterior[0][2]=sigma[2]->getNextComplex();
         sigmaCurrentPosterior[1][1]=sigma[3]->getNextComplex();
         sigmaCurrentPosterior[1][2]=sigma[4]->getNextComplex();
         sigmaCurrentPosterior[2][2]=sigma[5]->getNextComplex();
         // compleating matrixes
         sigmaCurrentPosterior[1][0]=sigmaCurrentPosterior[0][1];
         sigmaCurrentPosterior[2][0]=sigmaCurrentPosterior[0][2];
         sigmaCurrentPosterior[2][1]=sigmaCurrentPosterior[1][2];

         for(int l=0;l<6;l++)
           for(int m=l+1;m<6;m++)
           {
             sigmaFullPrior[m][l].re = sigmaFullPrior[l][m].re;
             sigmaFullPrior[m][l].im = -sigmaFullPrior[l][m].im;
           }

         // computing derived quantities

         for(int l=0;l<3;l++){
           muCurrentPrior[l].re =muFullPrior[l].re+muFullPrior[l+3].re;
           muCurrentPrior[l].im =muFullPrior[l].im+muFullPrior[l+3].im;
         }
         for(int l=0;l<6;l++)
           for(int m=0;m<3;m++){
             sigmaFullVsCurrentPrior[l][m].re =  sigmaFullPrior[l][m].re+sigmaFullPrior[l][m+3].re;
             sigmaFullVsCurrentPrior[l][m].im =  sigmaFullPrior[l][m].im+sigmaFullPrior[l][m+3].im;
           }
         for(int l=0;l<3;l++)
           for(int m=0;m<3;m++){
             sigmaCurrentPrior[l][m].re =  sigmaFullVsCurrentPrior[l][m].re + sigmaFullVsCurrentPrior[l+3][m].re;
             sigmaCurrentPrior[l][m].im =  sigmaFullVsCurrentPrior[l][m].im + sigmaFullVsCurrentPrior[l+3][m].im;
           }

        // solving the matrixequations see NR-Note: SAND/04/2012 page 6.

        // computing: sandwich= inv(sigmaCurrentPrior)*sigmaCurrentVsFullPrior=inv(sigmaCurrentPrior)*adjoint(sigmaFullVsCurrentPrior);
         lib_matrCopyCpx(sigmaCurrentPrior, 3, 3, sigmaCurrentPriorChol);
         lib_matrAdjoint(sigmaFullVsCurrentPrior, 6, 3,sandwich); // here sandwich = adjoint(sigmaFullVsCurrentPrior);

         int flag=lib_matrCholCpx(3, sigmaCurrentPriorChol);                        // these two lines  returns

         if(flag==0){
           lib_matrAXeqBMatCpx(3, sigmaCurrentPriorChol, sandwich, 6);       // sandwich= inv(sigmaCurrentPrior)*adjoint(sigmaFullVsCurrentPrior);

           // computing: sigmaFullPosterior = sigmaFullPrior + adjoint(sandwich)*(sigmaCurrentPosterior-sigmaCurrentPrior)*(sandwich);

           lib_matrSubtMatCpx(sigmaCurrentPrior, 3, 3,sigmaCurrentPosterior);// sigmaCurrentPosterior contains the difference to sigmaCurrentPrior

           lib_matrProdCpx(sigmaCurrentPosterior, sandwich, 3, 3, 6, helper);  // helper= (sigmaCurrentPosterior-sigmaCurrentPrior)*(sandwich);
           lib_matrAdjoint(sandwich,3,6,adjointSandwich);
           lib_matrProdCpx(adjointSandwich, helper, 6, 3, 6, sigmaFullPosterior); // here: sigmaFullPosterior =adjoint(sandwich*)(sigmaCurrentPosterior-sigmaCurrentPrior)*(sandwich);
           lib_matrAddMatCpx(sigmaFullPrior, 6, 6, sigmaFullPosterior); // Final computation

           // computing: muFullPosterior = muFullPrior + adjoint(sandwich)*(muCurrentPosterior-muCurrentPrior)
           lib_matrSubtVecCpx(muCurrentPrior, 3, muCurrentPosterior);// muCurrentPosterior contains: (muCurrentPosterior-muCurrentPrior)
           lib_matrProdMatVecCpx(adjointSandwich, muCurrentPosterior, 6, 3, muFullPosterior); //muFullPosterior=sandwich*(muCurrentPosterior-muCurrentPrior)
           lib_matrAddVecCpx( muFullPrior, 6, muFullPosterior);
         }else
         {
          // /*
          counter++;
          if(counter==100)
          {
            lib_matrDumpCpx("priorFull", sigmaFullPrior, 6,6);
            lib_matrDumpCpx("priorCurrent", sigmaCurrentPrior, 3,3);
            lib_matrDumpCpx("posteriorCurrent", sigmaCurrentPosterior, 3,3); // */
          }
           lib_matrCopyCpx(sigmaFullPrior, 6, 6, sigmaFullPosterior);
           for(int l=0;l<6;l++)
             muFullPosterior[l]= muFullPrior[l];

         }

        // writing to grids
         sigma_static_static_[0]->setNextComplex(  sigmaFullPosterior[0][0]);
         sigma_static_static_[1]->setNextComplex(  sigmaFullPosterior[0][1]);
         sigma_static_static_[2]->setNextComplex(  sigmaFullPosterior[0][2]);
         sigma_static_dynamic_[0]->setNextComplex( sigmaFullPosterior[0][3]);
         sigma_static_dynamic_[1]->setNextComplex( sigmaFullPosterior[0][4]);
         sigma_static_dynamic_[2]->setNextComplex( sigmaFullPosterior[0][5]);
         sigma_static_static_[3]->setNextComplex(  sigmaFullPosterior[1][1]);
         sigma_static_static_[4]->setNextComplex(  sigmaFullPosterior[1][2]);
         sigma_static_dynamic_[3]->setNextComplex( sigmaFullPosterior[1][3]);
         sigma_static_dynamic_[4]->setNextComplex( sigmaFullPosterior[1][4]);
         sigma_static_dynamic_[5]->setNextComplex( sigmaFullPosterior[1][5]);
         sigma_static_static_[5]->setNextComplex(  sigmaFullPosterior[2][2]);
         sigma_static_dynamic_[6]->setNextComplex( sigmaFullPosterior[2][3]);
         sigma_static_dynamic_[7]->setNextComplex( sigmaFullPosterior[2][4]);
         sigma_static_dynamic_[8]->setNextComplex( sigmaFullPosterior[2][5]);
         sigma_dynamic_dynamic_[0]->setNextComplex(sigmaFullPosterior[3][3]);
         sigma_dynamic_dynamic_[1]->setNextComplex(sigmaFullPosterior[3][4]);
         sigma_dynamic_dynamic_[2]->setNextComplex(sigmaFullPosterior[3][5]);
         sigma_dynamic_dynamic_[3]->setNextComplex(sigmaFullPosterior[4][4]);
         sigma_dynamic_dynamic_[4]->setNextComplex(sigmaFullPosterior[4][5]);
         sigma_dynamic_dynamic_[5]->setNextComplex(sigmaFullPosterior[5][5]);

         mu_static_[0]->setNextComplex( muFullPosterior[0]);
         mu_static_[1]->setNextComplex( muFullPosterior[1]);
         mu_static_[2]->setNextComplex( muFullPosterior[2]);
         mu_dynamic_[0]->setNextComplex(muFullPosterior[3]);
         mu_dynamic_[1]->setNextComplex(muFullPosterior[4]);
         mu_dynamic_[2]->setNextComplex(muFullPosterior[5]);
      }
    }
  }

  printf("\n\n #of Shortcuts in split = %d, this is  %f of 100 percent \n",counter, double(counter*100.0)/double(cnxp*nyp*nzp));

  for(int i = 0; i<3; i++)
  {
    mu[i]->endAccess();
    mu_static_[i]->endAccess();
    mu_dynamic_[i]->endAccess();
  }
  for(int i = 0; i<6; i++)
  {
    sigma[i]->endAccess();
    sigma_static_static_[i]->endAccess();
    sigma_dynamic_dynamic_[i]->endAccess();
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->endAccess();


  for(int i=0;i<6;i++)
  {
    delete [] adjointSandwich[i];
    delete [] sigmaFullPrior[i];
    delete [] sigmaFullVsCurrentPrior[i];
    delete [] sigmaFullPosterior[i];
  }
  for(int i=0;i<3;i++)
  {
    delete [] sandwich[i];
    delete [] helper[i];
    delete [] sigmaCurrentPrior[i];
    delete [] sigmaCurrentPriorChol[i];
    delete [] sigmaCurrentPosterior[i];
  }

  delete [] muFullPrior;
  delete [] muFullPosterior;
  delete [] muCurrentPrior;
  delete [] muCurrentPosterior;

  delete [] sandwich;
  delete [] helper;
  delete [] adjointSandwich;
  delete [] sigmaFullPrior;
  delete [] sigmaFullPosterior;
  delete [] sigmaFullVsCurrentPrior;
  delete [] sigmaCurrentPrior;
  delete [] sigmaCurrentPriorChol;
  delete [] sigmaCurrentPosterior;
}

void State4D::evolve(int time_step, const TimeEvolution timeEvolution )
{
  LogKit::WriteHeader("Evolving distribution of dynamic part");
  // Evolution matrix and correction terms from TimeEvolution class
  const NRLib::Matrix evolution_matrix     = timeEvolution.getEvolutionMatrix(time_step);
  const NRLib::Vector mean_correction_term = timeEvolution.getMeanCorrectionTerm(time_step);
  const NRLib::Matrix cov_correction_term  = timeEvolution.getCovarianceCorrectionTerm(time_step);

   // Holders of FFTGrid pointers
  std::vector<FFTGrid *> mu(6);

  std::vector<FFTGrid *> sigma(21);
  // Check the order of the grids here:
  // The order is used later in set-up of symmetric sigma matrix


  mu[0] = getMuVpStatic(); //mu_static_Alpha
  mu[1] = getMuVsStatic(); //mu_static_Beta
  mu[2] = getMuRhoStatic(); //mu_static_Rho
  mu[3] = getMuVpDynamic(); //mu_dynamic_Alpha
  mu[4] = getMuVsDynamic(); //mu_dynamic_Beta
  mu[5] = getMuRhoDynamic(); //mu_dynamic_Rho

  // Note the order of the grids here: The order is used other places in code.
  sigma[0]  = getCovVpVpStaticStatic(); //cov_ss_AlphaStatic AlphaStatic
  sigma[1]  = getCovVpVsStaticStatic();  //cov_ss_AlphaStaticBetaStatic
  sigma[2]  = getCovVpRhoStaticStatic();  //cov_ss_AlphaStaticRhoStatic
  sigma[3]  = getCovVsVsStaticStatic(); //cov_ss_BetaStaticBetaStatic
  sigma[4]  = getCovVsRhoStaticStatic(); //cov_ss_BetaStaticRhoStatic
  sigma[5]  = getCovRhoRhoStaticStatic(); //cov_ss_RhoStaticRhoStatic
  sigma[6]  = getCovVpVpDynamicDynamic(); //cov_dd_AlphaDynamicAlphaDynamic
  sigma[7]  = getCovVpVsDynamicDynamic(); //cov_dd_AlphaDynamicBetaDynamic
  sigma[8]  = getCovVpRhoDynamicDynamic(); //cov_dd_AlphaDynamicRhoDynamic
  sigma[9]  = getCovVsVsDynamicDynamic(); //cov_dd_BetaDynamicBetaDynamic
  sigma[10] = getCovVsRhoDynamicDynamic(); //cov_dd_BetaDynamicRhoDynamic
  sigma[11] = getCovRhoRhoDynamicDynamic(); //cov_dd_RhoDynamicRhoDynamic

  sigma[12] = getCovVpVpStaticDynamic(); //cov_sd_AlphaStaticAlphaDynamic
  sigma[13] = getCovVpVsStaticDynamic();  //cov_sd_AlphaStaticBetaDynamic
  sigma[14] = getCovVpRhoStaticDynamic(); //cov_sd_AlphaStaticRhoDynamic
  sigma[15] = getCovVsVpStaticDynamic();  //cov_sd_BetaStaticAlphaDynamic
  sigma[16] = getCovVsVsStaticDynamic();  //cov_sd_BetaStaticBetaDynamic
  sigma[17] = getCovVsRhoStaticDynamic(); //cov_sd_BetaStaticRhoDynamic
  sigma[18] = getCovRhoVpStaticDynamic(); //cov_sd_RhoStaticAlphaDynamic
  sigma[19] = getCovRhoVsStaticDynamic(); //cov_sd_RhoStaticBetaDynamic
  sigma[20] = getCovRhoRhoStaticDynamic();  //cov_sd_RhoStaticRhoDynamic

  // We assume FFT transformed grids
  for(int i = 0; i<6; i++)
  {
    assert(mu[i]->getIsTransformed());
    mu[i]->setAccessMode(FFTGrid::READANDWRITE);
  }
  for(int i = 0; i<21; i++)
  {
    assert(sigma[i]->getIsTransformed());
    sigma[i]->setAccessMode(FFTGrid::READANDWRITE);
  }

  NRLib::Vector mu_real(6);
  NRLib::Vector mu_imag(6);
  NRLib::Vector mu_real_next(6);
  NRLib::Vector mu_imag_next(6);

  NRLib::Matrix sigma_real(6,6);
  NRLib::Matrix sigma_imag(6,6);
  NRLib::Matrix sigma_real_next(6,6);
  NRLib::Matrix sigma_imag_next(6,6);
  NRLib::Matrix tmp(6,6);

  fftw_complex get_value,return_value;

  int nz   = mu[0]->getNz();
  int ny   = mu[0]->getNy();
  int nx   = mu[0]->getNx();
  int nzp  = mu[0]->getNzp();
  int nyp  = mu[0]->getNyp();
  int nxp  = mu[0]->getNxp();
  int cnxp = mu[0]->getCNxp();


   FFTGrid timeIncSpatialCorr=FFTGrid( nx,  ny,  nz,  nxp,  nyp,  nzp);
   timeIncSpatialCorr.createGrid();
   // NBNB this isa quick fix
   if(time_step==0)
     timeIncSpatialCorr.fillInGenExpCorr(float(nx)/3.0,float(ny)/3.0,50.0,0.0f,0.0f); // longer vertical range in first go
   else
     timeIncSpatialCorr.fillInGenExpCorr(float(nx)/3.0,float(ny)/3.0,5.0,0.0f,0.0f); //
   //NBNB OK  this correlation should come from interface and be more like:
   /* timeIncSpatialCorr.fillInParamCorr(timeEvolution.getPriorCorrXY(time_step),
                                    timeEvolution.getPriorcircCorrT(timeStep),
                                      timeEvolution.gradI(timestep),
                                       timeEvolution.gradJ(timestep));  */
   //timeIncSpatialCorr.writeAsciiRaw("timeIncSpatialCorr.dat");

   timeIncSpatialCorr.fftInPlace();
   timeIncSpatialCorr.setAccessMode(FFTGrid::READ);
  // Iterate through all points in the grid and perform forward transition in time
  for (int k = 0; k < nzp; k++) {
    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {

        fftw_complex ijkLambda = timeIncSpatialCorr.getNextComplex();
        float realTocomplexScaleFactor =  (i==0 && j==0 && k==0 )? float(std::sqrt(double(nxp*nyp*nzp))): 0.0f;  // note add a constant in real domain is
                                                                                                                 // just a value on the 0,0,0 in fft domain
                                                                                                                 // is for the mean what  ijkLambda is for the covariance

        // Set up vectors from the FFT grids
        for (int d = 0; d < 6; d++) {
          get_value  = mu[d]->getNextComplex();
          mu_real(d) = get_value.re;
          mu_imag(d) = get_value.im;
        }

        // Evolve values
        mu_real_next = evolution_matrix*mu_real + mean_correction_term* realTocomplexScaleFactor; // last term is only for (0,0,0 ) coefficient see above
        mu_imag_next = evolution_matrix*mu_imag;

        // Update values in the FFT-grids
        for (int d = 0; d < 6; d++) {
          return_value.re = static_cast<float>(mu_real_next(d));
          return_value.im = static_cast<float>(mu_imag_next(d));
          mu[d]->setNextComplex(return_value);
        }
        // Set up matrices from the FFT-grids.
        // Note: Here we assume a specific order of the elements in the sigma-vector.
        // Static and dynamic parts.
        int counter = 0;
        for (int d1 = 0; d1 < 3; d1++) {
          for (int d2 = d1; d2 < 3; d2++) {
            get_value  = sigma[counter]->getNextComplex();
            sigma_real(d1, d2) = get_value.re;
            sigma_imag(d1, d2) = get_value.im;
            get_value  = sigma[counter+6]->getNextComplex();
            sigma_real(d1+3, d2+3) = get_value.re;  // d1+3 and d2+3 due to block structure of matrix
            sigma_imag(d1+3, d2+3) = get_value.im;

            counter++;

            //Enforcing symmetry of overall matrix.
            if (d1 != d2) {
              sigma_real(d2, d1) = sigma_real(d1, d2);
              sigma_imag(d2 ,d1) = sigma_imag(d1, d2);

              sigma_real(d2+3, d1+3) = sigma_real(d1+3, d2+3);
              sigma_imag(d2+3 ,d1+3) = sigma_imag(d1+3, d2+3);
            }
          }
        }
        // Static-dynamic covariance.
        counter = 12;
        for (int d1 = 0; d1 < 3; d1++) {
          for (int d2 = 3; d2 < 6; d2++) {
            get_value=sigma[counter]->getNextComplex();
            sigma_real(d1, d2) = get_value.re;
            sigma_imag(d1, d2) = get_value.im;
            counter++;

            //Enforcing symmetry of overall matrix.
            sigma_real(d2, d1) = sigma_real(d1, d2);
            sigma_imag(d2 ,d1) = sigma_imag(d1, d2);
          }
        }
        // Evolve values.
        tmp = (evolution_matrix*sigma_real);
        sigma_real_next = tmp*transpose(evolution_matrix) + cov_correction_term *ijkLambda.re;

        tmp = evolution_matrix*sigma_imag;
        sigma_imag_next = tmp*transpose(evolution_matrix);

        counter = 0;
        for (int d1 = 0; d1 < 3; d1++) {// Update values in the FFT-grids.
          for (int d2 = d1; d2 < 3; d2++) {// Static and dynamic parts.
            return_value.re = static_cast<float>(sigma_real_next(d1, d2));
            return_value.im = static_cast<float>(sigma_imag_next(d1, d2));
            sigma[counter]->setNextComplex(return_value);

            return_value.re = static_cast<float>(sigma_real_next(d1+3, d2+3));  //d1+3 and d2+3 due to block structure of matrix
            return_value.im = static_cast<float>(sigma_imag_next(d1+3, d2+3));
            sigma[counter+6]->setNextComplex(return_value);

            counter++;
          }
        }
        // Static-dynamic covariance.
        counter = 12;
        for (int d1 = 0; d1 < 3; d1++) {
          for (int d2 = 3; d2 < 6; d2++) {
            return_value.re = static_cast<float>(sigma_real_next(d1, d2));
            return_value.im = static_cast<float>(sigma_imag_next(d1, d2));
            sigma[counter]->setNextComplex(return_value);
            counter++;
          }
        }
      }
    }
  }

  for(int i = 0; i<6; i++)
    mu[i]->endAccess();
  for(int i = 0; i<21; i++)
    sigma[i]->endAccess();
}

bool
State4D::allGridsAreTransformed()
{
  bool allTransformed=true;

   for(int i = 0; i<3; i++)
     if(!mu_static_[i]->getIsTransformed())
       allTransformed=false;
   for(int i = 0; i<3; i++)
     if(!mu_dynamic_[i]->getIsTransformed())
       allTransformed=false;
   for(int i = 0; i<6; i++)
     if(!sigma_static_static_[i]->getIsTransformed())
       allTransformed=false;
   for(int i = 0; i<6; i++)
     if(!sigma_dynamic_dynamic_[i]->getIsTransformed())
       allTransformed=false;
   for(int i = 0; i<9; i++)
     if(!sigma_static_dynamic_[i]->getIsTransformed())
       allTransformed=false;

   return allTransformed;
}

void
State4D::FFT()
{
  for(int i = 0; i<3; i++)
    mu_static_[i]->fftInPlace();
  for(int i = 0; i<3; i++)
    mu_dynamic_[i]->fftInPlace();
  for(int i = 0; i<6; i++)
    sigma_static_static_[i]->fftInPlace();
  for(int i = 0; i<6; i++)
    sigma_dynamic_dynamic_[i]->fftInPlace();
  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->fftInPlace();
}

void
State4D::iFFT()
{
  for(int i = 0; i<3; i++)
    mu_static_[i]->invFFTInPlace();
  for(int i = 0; i<3; i++)
    mu_dynamic_[i]->invFFTInPlace();
  for(int i = 0; i<6; i++)
    sigma_static_static_[i]->invFFTInPlace();
  for(int i = 0; i<6; i++)
    sigma_dynamic_dynamic_[i]->invFFTInPlace();
  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->invFFTInPlace();
}


void
State4D::iFFTMean()
{
  for(int i = 0; i<3; i++)
    mu_static_[i]->invFFTInPlace();
  for(int i = 0; i<3; i++)
    mu_dynamic_[i]->invFFTInPlace();
}

void
State4D::iFFTCov()
{
  for(int i = 0; i<6; i++)
    sigma_static_static_[i]->invFFTInPlace();
  for(int i = 0; i<6; i++)
    sigma_dynamic_dynamic_[i]->invFFTInPlace();
  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->invFFTInPlace();
}



std::vector<FFTGrid*>
State4D::doRockPhysicsInversion(TimeLine&  time_line, const std::vector<DistributionsRock *> rock_distributions,  TimeEvolution timeEvolution)
{
  LogKit::WriteHeader("Start 4D rock physics inversion");
  bool debug=true; // triggers printouts
  // inverse transform posterior
  iFFT();

  // Note rockSample contains rock sample for all time steps.
  int nSim=10000; // NBNB OK 10000->1000 for speed during debug

  LogKit::LogFormatted(LogKit::Low,"\nSampling rockphysics distribution...");
  std::vector<std::vector<std::vector<double> > > rockSample = timeEvolution.returnCorrelatedSample(nSim,time_line, rock_distributions);
  LogKit::LogFormatted(LogKit::Low,"done\n\n");

  int nTimeSteps=rockSample.size();
  //nSim=rockSample[0].size();
  int nParam=rockSample[0][0].size();
  int nM = 3;  // number of seismic parameters.
  int nRockProperties = nParam-nM; // number of rock parameters



  //write rockSamples to check ok
  if(debug)
  {
    NRLib::Matrix rockSamples(nSim,nParam*nTimeSteps);

    for(int i = 0;i<nSim;i++)
      for(int j=0;j < nTimeSteps;j++ )
        for(int k=0;k<nParam;k++)
        {
          int ind=k+j*nParam;
          rockSamples(i,ind)=rockSample[j][i][k];
        }

    NRLib::WriteMatrixToFile("rockSample.dat", rockSamples);
  }

  std::vector<std::vector<double> > mSamp = makeSeismicParamsFromrockSample(rockSample);

  //write seismic parameters to check ok
  if(debug)
  {
    NRLib::Matrix seisPar(nSim,6);

    for(int i = 0;i<nSim;i++)
      for(int j=0;j < 6;j++ )
        {
          seisPar(i,j)=mSamp[j][i];
        }

    NRLib::WriteMatrixToFile("seisParSample.dat", seisPar);
  }


  NRLib::Vector fullPriorMean    = timeEvolution.computePriorMeanStaticAndDynamicLastTimeStep();
  NRLib::Matrix fullPriorCov     = timeEvolution.computePriorCovStaticAndDynamicLastTimeStep();
  NRLib::Matrix fullPosteriorCov = GetFullCov();

  if(debug)
  {
    NRLib::WriteMatrixToFile("priorCov.dat", fullPriorCov );
    NRLib::WriteMatrixToFile("posteriorCov.dat", fullPosteriorCov);
    NRLib::WriteVectorToFile("priorMean.dat",fullPriorMean );
  }

  LogKit::LogFormatted(LogKit::Low,"\nMaking rock-physics lookup tables, table 1 of %d\n",nRockProperties+1);

  RockPhysicsInversion4D* rockPhysicsInv = new RockPhysicsInversion4D(fullPriorMean,fullPriorCov,fullPosteriorCov,mSamp);
  //LogKit::LogFormatted(LogKit::Low,"done\n\n");

  std::vector<FFTGrid*> prediction(nRockProperties);

  LogKit::LogFormatted(LogKit::Low,"\nDoing rock-physics predictions...");

  for(int i=0;i<nRockProperties;i++)
  {
    LogKit::LogFormatted(LogKit::Low,"\nMaking rock-physics lookup tables, table %d of %d\n",i+2,nRockProperties+1);
    std::vector<double> rSamp = getRockPropertiesFromRockSample(rockSample,i);
    rockPhysicsInv->makeNewPredictionTable(mSamp,rSamp);
    prediction[i] = rockPhysicsInv->makePredictions(mu_static_, mu_dynamic_ );
  }
  LogKit::LogFormatted(LogKit::Low,"done\n\n");

  return prediction;
}

NRLib::Vector
State4D::GetFullMean000()
{
  NRLib::Vector full_mean(6);
  for(int i=0;i<3;i++){
    mu_static_[i]->setAccessMode(FFTGrid::READ);
    mu_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  full_mean(0)=mu_static_[0]->getNextReal();
  full_mean(1)=mu_static_[1]->getNextReal();
  full_mean(2)=mu_static_[2]->getNextReal();
  full_mean(3)=mu_dynamic_[0]->getNextReal();
  full_mean(4)=mu_dynamic_[1]->getNextReal();
  full_mean(5)=mu_dynamic_[2]->getNextReal();

   for(int i=0;i<3;i++){
    mu_static_[i]->endAccess();
    mu_dynamic_[i]->endAccess();
  }

  return full_mean;
}


NRLib::Matrix
State4D::GetFullCov()
{

  NRLib::Matrix covMat(6,6);
  covMat=0.0;

  for(int i=0;i<6;i++)
    sigma_static_static_[i]->setAccessMode(FFTGrid::READ);

  covMat(0,0) = sigma_static_static_[0]->getNextReal();  // [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho (all static)
  covMat(0,1) = sigma_static_static_[1]->getNextReal();
  covMat(0,2) = sigma_static_static_[2]->getNextReal();
  covMat(1,1) = sigma_static_static_[3]->getNextReal();
  covMat(1,2) = sigma_static_static_[4]->getNextReal();
  covMat(2,2) = sigma_static_static_[5]->getNextReal();
  covMat(1,0) = covMat(0,1);
  covMat(2,0) = covMat(0,2);
  covMat(2,1) = covMat(1,2);
  for(int i=0;i<6;i++)
    sigma_static_static_[i]->endAccess();

  for(int i=0;i<6;i++)
    sigma_dynamic_dynamic_[i]->setAccessMode(FFTGrid::READ);

  covMat(3,3) = sigma_dynamic_dynamic_[0]->getNextReal();  // [0] = vp_vp, [1] = vp_vs, [2] = vp_rho ,[3] = vs_vs, [4] = vs_rho, [5] = rho_rho (all static)
  covMat(3,4) = sigma_dynamic_dynamic_[1]->getNextReal();
  covMat(3,5) = sigma_dynamic_dynamic_[2]->getNextReal();
  covMat(4,4) = sigma_dynamic_dynamic_[3]->getNextReal();
  covMat(4,5) = sigma_dynamic_dynamic_[4]->getNextReal();
  covMat(5,5) = sigma_dynamic_dynamic_[5]->getNextReal();
  covMat(4,3) = covMat(3,4);
  covMat(5,3) = covMat(3,5);
  covMat(5,4) = covMat(4,5);
  for(int i=0;i<6;i++)
    sigma_dynamic_dynamic_[i]->endAccess();


  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      int counter=j+3*i;
      sigma_static_dynamic_[counter]->setAccessMode(FFTGrid::READ);
      double value = double(sigma_static_dynamic_[counter]->getNextReal()); // [0] = vpStat_vpDyn, [1] = vpStat_vsDyn, [2] = vpStat_rhoDyn ,[3] = vsStat_vpDyn,
      sigma_static_dynamic_[counter]->endAccess();
      covMat(i,j+3) = value;
      covMat(j+3,i) = value;
    }

  return covMat;
}


std::vector<std::vector<double> >
State4D::makeSeismicParamsFromrockSample(std::vector<std::vector<std::vector<double> > > rS)
{
  int k_max = int(rS.size());// number of surveys
  int i_max = int(rS[0].size());// number of samples
//  int dim   = int(rS[0][0].size());// 3 + number of reservoir variables

  std::vector<std::vector<double> > m(6, std::vector<double>(i_max));

  for(int i=0;i<i_max;i++)
  {
    m[0][i]=rS[0][i][0];
    m[1][i]=rS[0][i][1];
    m[2][i]=rS[0][i][2];
    m[3][i]=rS[k_max-1][i][0]-rS[0][i][0];
    m[4][i]=rS[k_max-1][i][1]-rS[0][i][1];
    m[5][i]=rS[k_max-1][i][2]-rS[0][i][2];
  }

  return m;
}

std::vector<double>
State4D::getRockPropertiesFromRockSample(std::vector<std::vector<std::vector<double> > > rS,int varNumber)
{

  int k_max = int(rS.size());// number of surveys
  int i_max = int(rS[0].size());// number of samples
  //int dim   = int(rS[0][0].size());// 3 + number of reservoir variables

  std::vector<double>  r(i_max);
  for(int i=0;i<i_max;i++)
    r[i]=rS[k_max-1][i][3+varNumber];

  return r;
}
