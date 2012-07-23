#include "src/state4d.h"
#include "src/seismicparametersholder.h"
#include "src/timeevolution.h"
#include "lib/lib_matr.h"

#include "nrlib/flens/nrlib_flens.cpp"


State4D::State4D(){

  // Initilizing the 27 grids. All grids initialize to NULL.
  mu_static_.resize(3, NULL);
  mu_dynamic_.resize(3, NULL);
  sigma_static_static_.resize(6, NULL);
  sigma_dynamic_dynamic_.resize(6, NULL);
  sigma_static_dynamic_.resize(9, NULL);
}

State4D::~State4D(){
}

void State4D::setDynamicMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho){
  mu_dynamic_[0] = vp;
  mu_dynamic_[1] = vs;
  mu_dynamic_[2] = rho;
}

void State4D::setStaticMu(FFTGrid *vp, FFTGrid *vs, FFTGrid *rho){
  mu_static_[0] = vp;
  mu_static_[1] = vs;
  mu_static_[2] = rho;
}

void State4D::setStaticSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho){
  sigma_static_static_[0] = vpvp;
  sigma_static_static_[1] = vpvs;
  sigma_static_static_[2] = vprho;
  sigma_static_static_[3] = vsvs;
  sigma_static_static_[4] = vsrho;
  sigma_static_static_[5] = rhorho;
}

void State4D::setDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhorho){
  sigma_dynamic_dynamic_[0] = vpvp;
  sigma_dynamic_dynamic_[1] = vpvs;
  sigma_dynamic_dynamic_[2] = vprho;
  sigma_dynamic_dynamic_[3] = vsvs;
  sigma_dynamic_dynamic_[4] = vsrho;
  sigma_dynamic_dynamic_[5] = rhorho;
}
void State4D::setStaticDynamicSigma(FFTGrid *vpvp, FFTGrid *vpvs, FFTGrid *vprho, FFTGrid *vsvp, FFTGrid *vsvs, FFTGrid *vsrho, FFTGrid *rhovp, FFTGrid *rhovs, FFTGrid *rhorho){
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

void State4D::merge(SeismicParametersHolder & current_state )
{
  bool IsTested=false;
  assert(IsTested);

  // We assume FFT transformed grids
  assert(allGridsAreTransformed());

  std::vector<FFTGrid *> mu(3);
  mu[0] =  current_state.GetMuAlpha(); //mu_Alpha
  mu[1] =  current_state.GetMuBeta(); //mu_Beta
  mu[2] =  current_state.GetMuRho(); //mu_Rho

  for(int i = 0; i<3; i++)
  {
    assert(mu[i]->getIsTransformed());
    mu[i]->setAccessMode(FFTGrid::WRITE);
    mu_static_[i]->setAccessMode(FFTGrid::READ);
    mu_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  std::vector<FFTGrid *> sigma(6);
  sigma[0]=current_state.GetCovAlpha();
  sigma[1]=current_state.GetCovBeta();
  sigma[2]=current_state.GetCovRho();
  sigma[3]=current_state.GetCrCovAlphaBeta();
  sigma[4]=current_state.GetCrCovAlphaRho();
  sigma[5]=current_state.GetCrCovBetaRho();

  for(int i = 0; i<6; i++)
  {
    assert(sigma[i]->getIsTransformed());
    sigma[i]->setAccessMode(FFTGrid::WRITE);
    sigma_static_static_[i]->setAccessMode(FFTGrid::READ);
    sigma_dynamic_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->setAccessMode(FFTGrid::READ);

  int nzp = mu[0]->getNzp();
  int nyp = mu[0]->getNyp();
  int cnxp = mu[0]->getCNxp();

  fftw_complex*  muFullPrior=new fftw_complex[6];
  fftw_complex*  muCurrentPrior=new fftw_complex[3];

  fftw_complex** sigmaFullPrior=new fftw_complex*[6];
  for(int i=0;i<6;i++)
    sigmaFullPrior[i]=new fftw_complex[6];

  fftw_complex** sigmaCurrentPrior=new fftw_complex*[3];
  for(int i=0;i<3;i++)
    sigmaCurrentPrior[i]=new fftw_complex[3];

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
         for(int l=1;l<6;l++)
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
  for(int i = 0; i<6; i++)
  {
    sigma[i]->endAccess();
    sigma_static_static_[i]->endAccess();
    sigma_dynamic_dynamic_[i]->endAccess();
  }

  for(int i = 0; i<9; i++)
    sigma_static_dynamic_[i]->endAccess();


  for(int i=0;i<6;i++)
    delete [] sigmaFullPrior[i];

  for(int i=0;i<3;i++)
    delete [] sigmaCurrentPrior[i];

  delete [] muFullPrior;
  delete [] muCurrentPrior;

  delete [] sigmaFullPrior;
  delete [] sigmaCurrentPrior;
}

void State4D::split(SeismicParametersHolder current_state )
{
  bool IsTested=false;
  assert(IsTested);

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
  sigma[1]=current_state.GetCovBeta();
  sigma[2]=current_state.GetCovRho();
  sigma[3]=current_state.GetCrCovAlphaBeta();
  sigma[4]=current_state.GetCrCovAlphaRho();
  sigma[5]=current_state.GetCrCovBetaRho();

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


  fftw_complex** sigmaFullPrior=new fftw_complex*[6];
  fftw_complex** sigmaFullPosterior=new fftw_complex*[6];
  fftw_complex** sigmaFullVsCurrentPrior=new fftw_complex*[6];
  fftw_complex** sandwich=new fftw_complex*[6];


  for(int i=0;i<6;i++)
  {
    sigmaFullPrior[i]=new fftw_complex[6];
    sigmaFullVsCurrentPrior[i]=new fftw_complex[3];
    sandwich[i]=new fftw_complex[3];
    sigmaFullPosterior[i]=new fftw_complex[6];
  }
  fftw_complex** sigmaCurrentPrior=new fftw_complex*[3];
  fftw_complex** sigmaCurrentPriorChol=new fftw_complex*[3];
  fftw_complex** sigmaCurrentPosterior=new fftw_complex*[3];
  for(int i=0;i<3;i++)
  {
    sigmaCurrentPrior[i]=new fftw_complex[3];
    sigmaCurrentPosterior[i]=new fftw_complex[3];
  }



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

         for(int l=1;l<6;l++)
           for(int m=l+1;m<6;m++)
           {
             sigmaFullPrior[m][l].re=sigmaFullPrior[l][m].re;
             sigmaFullPrior[m][l].im=-sigmaFullPrior[l][m].im;
           }

         // computing derived quantities

         for(int l=0;l<3;l++){
           muCurrentPrior[l].re =muFullPrior[l].re+muFullPrior[l+3].re;
           muCurrentPrior[l].im =muFullPrior[l].im+muFullPrior[l+3].im;
         }
         for(int l=0;l<6;l++)
           for(int m=0;m<3;m++)
           {
             sigmaFullVsCurrentPrior[l][m].re =  sigmaFullPrior[l][m].re+sigmaFullPrior[l][m+3].re;
             sigmaFullVsCurrentPrior[l][m].im =  sigmaFullPrior[l][m].im+sigmaFullPrior[l][m+3].im;
           }

         for(int l=0;l<3;l++)
           for(int m=0;m<3;m++)
           {
             sigmaCurrentPrior[l][m].re =  sigmaFullVsCurrentPrior[l][m].re + sigmaFullVsCurrentPrior[l+3][m].re;
             sigmaCurrentPrior[l][m].im =  sigmaFullVsCurrentPrior[l][m].im + sigmaFullVsCurrentPrior[l+3][m].im;
           }

        // solving the matrixequations see NR-Note: SAND/04/2012 page 6.

        // computing: sandwich=sigmaFullCurrentPrior*inv(sigmaCurrentPrior);
         lib_matrCopyCpx(sigmaCurrentPrior, 3, 3, sigmaCurrentPriorChol);
         lib_matrCopyCpx(sigmaFullVsCurrentPrior, 6, 3,sandwich);
         lib_matrCholCpx(3, sigmaCurrentPriorChol);
         lib_matrAXeqBMatCpx(3, sigmaCurrentPriorChol, sandwich, 6);

         // computing: sigmaFullPosterior = sigmaFullPrior + sandwich*(sigmaCurrentPosterior-sigmaCurrentPrior)*adjoint(sandwich);
         lib_matrSubtMatCpx( sigmaCurrentPosterior, 3, 3,sigmaCurrentPrior);// sigmaCurrentPosterior contains the difference to sigmaCurrentPrior
         lib_matrProdAdjointCpx(sigmaCurrentPosterior, sandwich, 3, 3, 6, sigmaFullVsCurrentPrior);
         lib_matrProdCpx(sandwich, sigmaFullVsCurrentPrior, 6, 3, 6, sigmaFullPosterior); // here: sigmaFullPosterior =sandwich*(sigmaCurrentPosterior-sigmaCurrentPrior)*adjoint(sandwich);
         lib_matrAddMatCpx(sigmaFullPrior, 6, 6, sigmaFullPosterior);

         // computing: muFullPosterior = muFullPrior + sandwich*(muCurrentPosterior-muCurrentPrior)
         lib_matrSubtVecCpx(muCurrentPrior, 3, muCurrentPosterior);// muCurrentPosterior contains: (muCurrentPosterior-muCurrentPrior)
         lib_matrProdMatVecCpx(sandwich, muCurrentPosterior, 6, 3, muFullPosterior); //muFullPosterior=sandwich*(muCurrentPosterior-muCurrentPrior)
         lib_matrAddVecCpx( muFullPrior, 6, muFullPosterior);

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
    delete [] sandwich[i];
    delete [] sigmaFullPrior[i];
    delete [] sigmaFullVsCurrentPrior;
    delete [] sigmaFullPosterior[i];
  }
  for(int i=0;i<3;i++)
  {
    delete [] sigmaCurrentPrior[i];
    delete [] sigmaCurrentPriorChol[i];
    delete [] sigmaCurrentPosterior[i];
  }

  delete [] muFullPrior;
  delete [] muFullPosterior;
  delete [] muCurrentPrior;
  delete [] muCurrentPosterior;

  delete [] sandwich;
  delete [] sigmaFullPrior;
  delete [] sigmaFullPosterior;
  delete [] sigmaFullVsCurrentPrior;
  delete [] sigmaCurrentPrior;
  delete [] sigmaCurrentPriorChol;
  delete [] sigmaCurrentPosterior;
}

void State4D::evolve(int time_step, const TimeEvolution timeEvolution )
{
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
    assert(mu[i]->getIsTransformed());
  for(int i = 0; i<21; i++)
    assert(sigma[i]->getIsTransformed());

  NRLib::Vector mu_real(6);
  NRLib::Vector mu_imag(6);
  NRLib::Vector mu_real_next(6);
  NRLib::Vector mu_imag_next(6);

  NRLib::Matrix sigma_real(6,6);
  NRLib::Matrix sigma_imag(6,6);
  NRLib::Matrix sigma_real_next(6,6);
  NRLib::Matrix sigma_imag_next(6,6);

  fftw_complex get_value,return_value;

  int nzp_ = mu[0]->getNzp();
  int nyp_ = mu[0]->getNyp();
  int cnxp = mu[0]->getCNxp();


  // Iterate through all points in the grid and perform forward transition in time
  for (int k = 0; k < nzp_; k++) {
    for (int j = 0; j < nyp_; j++) {
      for (int i = 0; i < cnxp; i++) {

        // Set up vectors from the FFT grids
        for (int d = 0; d < 6; d++) {
          get_value  = mu[d]->getNextComplex();
          mu_real(d) = get_value.re;
          mu_imag(d) = get_value.im;
        }

        // Evolve values
        mu_real_next = evolution_matrix*mu_real + mean_correction_term;
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
        sigma_real_next = (evolution_matrix*sigma_real);
        sigma_real_next = sigma_real_next*transpose(evolution_matrix) + cov_correction_term;

        sigma_imag_next = evolution_matrix*sigma_imag;
        sigma_imag_next = sigma_imag_next*transpose(evolution_matrix);

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

