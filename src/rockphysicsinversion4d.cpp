#include "src/rockphysicsinversion4d.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "lib/lib_matr.h"
#include <vector>

RockPhysicsInversion4D::RockPhysicsInversion4D()
{

}

RockPhysicsInversion4D::~RockPhysicsInversion4D()
{
 for(int i =0;i<4;i++)
   fftw_free(smoothingFilter_[i]);

 for(int i=0; i<2; i++){
    for (int j=0; j<nf_[0]; j++){
      delete meanRockPrediction_(i,j);
    }
  }

 fftwnd_destroy_plan(fftplan1_);
 fftwnd_destroy_plan(fftplan2_);
}

RockPhysicsInversion4D::RockPhysicsInversion4D(NRLib::Vector                      priorMean,
                                               NRLib::Matrix                      priorCov,
                                               NRLib::Matrix                      posteriorCov,
                                               std::vector<std::vector<double> >  mSamp)
{
  nf_.resize(4);
  nf_[0] = 60;
  nf_[1] = 60;
  nf_[2] = 60;
  nf_[3] = 60;
  nfp_= 135;

  fftplan1_ = rfftwnd_create_plan(1, &nfp_, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftplan2_ = rfftwnd_create_plan(1, &nfp_, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  v_.resize(4,6);
  SolveGEVProblem(priorCov,posteriorCov, v_);
  NRLib::Matrix tmp;
  NRLib::Matrix priorCovF;
  NRLib::Matrix posteriorCovF;

  meanf_        = priorMean*v_;
  tmp           = priorCov*v_;
  priorCovF     = transpose(v_)*tmp;
  tmp           = posteriorCov*v_;
  posteriorCovF = transpose(v_)*tmp;

  int nSamp = mSamp[0].size();

  NRLib::Vector m(6);
  for(int k=0;k<6;k++)
      m(k)=mSamp[k][0];

  maxf_=m*v_;
  minf_=m*v_;

  for(int i=0;i<nSamp;i++)
  {
    NRLib::Vector f;
    for(int k=0;k<6;k++)
      m(k)=mSamp[k][i];
    f=m*v_;
    for(int k=0;k<4;k++)
    {
      if(minf_(k) > f(k))
        minf_(k)=f(k);
      if(maxf_(k) < f(k))
        maxf_(k)=f(k);
    }
  }

  for(int k=0;k<4;k++)
  {
    minf_(k)-=3*sqrt(posteriorCovF(k,k));
    maxf_(k)+=3*sqrt(posteriorCovF(k,k));
  }

  smoothingFilter_.resize(4);
  priorDistribution_.resize(4);

  for(int i=0;i<4;i++)
  {
    double df=(maxf_(i)-minf_(i))/double(nf_[i]);
    priorDistribution_[i]=MakeGaussKernel(meanf_(i),priorCovF(i,i),minf_(i),df,nfp_);
    smoothingFilter_[i]=MakeSmoothingFilter(posteriorCovF(i,i),df);
  }

  allocatePredictionTables( );

  std::vector<double> dummy(nSamp);
  for(int i=0;i<nSamp;i++)
    dummy[i]=1.0;

  fillInTable( mSamp,dummy,0);
  DivideAndSmoothTable(0,priorDistribution_,smoothingFilter_);
}

void
RockPhysicsInversion4D::makeNewPredictionTable(std::vector<std::vector<double> >  mSamp,std::vector<double>   rSamp)
{
  ClearContentInPredictionTable( );
  fillInTable( mSamp,rSamp,1);
  smoothAllDirectionsAndNormalize();
}

void
RockPhysicsInversion4D::allocatePredictionTables( )
{
   meanRockPrediction_.Resize(2,nf_[0],NULL);

  for(int i=0; i<2; i++){
    for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(i,j) = new FFTGrid(nf_[1], nf_[2], nf_[3], nf_[1], nf_[2], nf_[3]);
      meanRockPrediction_(i,j)->setType(FFTGrid::PARAMETER);
      meanRockPrediction_(i,j)->createRealGrid(false);
      meanRockPrediction_(i,j)->setAccessMode(FFTGrid::WRITE);
      meanRockPrediction_(i,j)->fillInConstant(0.0);
      meanRockPrediction_(i,j)->endAccess();
    }
  }

}

void
RockPhysicsInversion4D::ClearContentInPredictionTable( )
{
  for (int j=0; j<nf_[0]; j++){
    meanRockPrediction_(1,j)->setAccessMode(FFTGrid::WRITE);
    meanRockPrediction_(1,j)->fillInConstant(0.0);
    meanRockPrediction_(1,j)->endAccess();
  }

}

FFTGrid *
RockPhysicsInversion4D::makePredictions(std::vector<FFTGrid *> mu_static_,
                         std::vector<FFTGrid *> mu_dynamic_ )
{
  int nx,ny,nz,rnxp,nxp,nyp,nzp;
  nx=mu_static_[0]->getNx();
  ny=mu_static_[0]->getNy();
  nz=mu_static_[0]->getNz();

  rnxp=mu_static_[0]->getCNxp()*2;
  nxp=mu_static_[0]->getNxp();
  nyp=mu_static_[0]->getNyp();
  nzp=mu_static_[0]->getNzp();


  FFTGrid* prediction= new FFTGrid(nx,ny,nz,nxp,nyp,nzp);
  prediction->createRealGrid();

  for(int i=0;i<3;i++)
  {
    mu_static_[i]->setAccessMode(FFTGrid::READ);
    mu_dynamic_[i]->setAccessMode(FFTGrid::READ);
  }

  prediction->setAccessMode(FFTGrid::WRITE);

  NRLib::Vector m(6);

  for(int k=0;k<nzp;k++)
    for(int j=0;j<nyp;j++)
      for(int i=0;i<rnxp;i++)
      {
        m(0)=mu_static_[0]->getNextReal();
        m(1)=mu_static_[1]->getNextReal();
        m(2)=mu_static_[2]->getNextReal();
        m(3)=mu_dynamic_[0]->getNextReal();
        m(4)=mu_dynamic_[1]->getNextReal();
        m(5)=mu_dynamic_[2]->getNextReal();
        float value;
        NRLib::Vector f;
        f=m*v_;
        value = float( getPredictedValue(f) );
        prediction->setNextReal(value);
      }

  for(int i=0;i<3;i++)
  {
    mu_static_[i]->endAccess();
    mu_dynamic_[i]->endAccess();
  }

  prediction->endAccess();

  return prediction;
}


void
RockPhysicsInversion4D::GetLowerIndexAndW(double minValue,double maxValue,int nValue,double valueIn,int& index, double& w)
{
  double dx    = (maxValue-minValue)/float(nValue);
  double value=valueIn+dx/2; // value of cell center NBNB OK check

  index = int(floor((value-minValue)/dx));
  w = (value-minValue)/dx- index;


  if(index<0)
  {
    index=0;
    w=0.0;
  }

  if(index> nValue-1)
  {
    index=nValue-1;
    w=1.0;
  }
}

int
RockPhysicsInversion4D::GetLowerIndex(double minValue,double maxValue,int nValue,double value)
{
  double dx    = (maxValue-minValue)/float(nValue);
  int index = int(floor((value-minValue)/dx)); // bin number cell center

  if(index<0)
    index=0;

  if(index> nValue-1)
    index=nValue-1;

  return index;
}


double
RockPhysicsInversion4D::getPredictedValue(NRLib::Vector f)
{
  //interploates in a 4D table
  double value;
  std::vector<std::vector<int> > indLoHi(2,std::vector<int>(4));
  std::vector<std::vector<double> > wLoHi(2,std::vector<double>(4));
  double w;
  int index;
  GetLowerIndexAndW(minf_(0),maxf_(0),nf_[0],f(0),index, w);
  wLoHi[0][0]=1-w;
  wLoHi[1][0]=w;
  indLoHi[0][0]=index;
  indLoHi[1][0]=std::min(nf_[0]-1,index+1);

  GetLowerIndexAndW(minf_(1),maxf_(1),nf_[1],f(1),index, w);
  wLoHi[0][1]=1-w;
  wLoHi[1][1]=w;
  indLoHi[0][1]=index;
  indLoHi[1][1]=std::min(nf_[1]-1,index+1);

  GetLowerIndexAndW(minf_(2),maxf_(2),nf_[2],f(2),index, w);
  wLoHi[0][2]=1-w;
  wLoHi[1][2]=w;
  indLoHi[0][2]=index;
  indLoHi[1][2]=std::min(nf_[2]-1,index+1);

  GetLowerIndexAndW(minf_(3),maxf_(3),nf_[3],f(3),index, w);
  wLoHi[0][3]=1-w;
  wLoHi[1][3]=w;
  indLoHi[0][3]=index;
  indLoHi[1][3]=std::min(nf_[3]-1,index+1);

  value=0.0;
  for(int i0=0;i0<2;i0++)
    for(int i1=0;i1<2;i1++)
      for(int i2=0;i2<2;i2++)
        for(int i3=0;i3<2;i3++)
        {
          w=wLoHi[i0][0]*wLoHi[i1][1]*wLoHi[i2][2]*wLoHi[i3][3];
          value+=w*GetGridValue(1,indLoHi[i0][0],indLoHi[i1][1],indLoHi[i2][2],indLoHi[i3][3]);
        }
return value;
}

double
RockPhysicsInversion4D::GetGridValue(int TableNr,int i0,int i1,int i2,int i3)
{
  return meanRockPrediction_(TableNr,i0)->getRealValue(i1,i2,i3);
}

void
RockPhysicsInversion4D::SetGridValue(int TableNr,int i0,int i1,int i2,int i3,double value)
{
  meanRockPrediction_(TableNr,i0)->setRealValue(i1,i2,i3,float(value));
}

void
RockPhysicsInversion4D::AddToGridValue(int TableNr,int i0,int i1,int i2,int i3,double value)
{
  double newValue=meanRockPrediction_(TableNr,i0)->getRealValue(i1,i2,i3)+value;
  meanRockPrediction_(TableNr,i0)->setRealValue(i1,i2,i3,float(newValue));
}


void
RockPhysicsInversion4D::SolveGEVProblem(NRLib::Matrix sigma_prior,
                                        NRLib::Matrix sigma_posterior,
                                        NRLib::Matrix & vOut){
  //Compute transforms v1, v2, v3 and v4  ----------------------------------------

  // computing filter
  NRLib::SymmetricMatrix  sigma_priorInv(6);
  for(int i=0;i<6;i++)
    for(int j=0;j<=i;j++)
      sigma_priorInv(j,i)=sigma_prior(i,j);

  NRLib::CholeskyInvert(sigma_priorInv);
  NRLib::Matrix eye = NRLib::IdentityMatrix(6);
  NRLib::Matrix filter = eye- sigma_priorInv*sigma_posterior;
  // computing components
  NRLib::Vector eigen_values(6);
  NRLib::Matrix eigen_vectors(6,6);
  NRLib::ComputeEigenVectors( filter ,eigen_values,eigen_vectors);

  // extracting four best components
  std::vector<int> current_index(6);
  current_index[0] = 0;current_index[1] = 1;
  current_index[2] = 2;current_index[3] = 3;
  current_index[4] = 4;current_index[5] = 5;

  std::vector<int> best_index(4);
  std::vector<int> keep_index(6);
  keep_index = current_index;

  NRLib::Vector current_value(6);
  NRLib::Vector keep_value(6);
  keep_value  = eigen_values;

  for(int i=0;i<4;i++){
    current_index=keep_index;
    current_value =keep_value;
    double maxVal=-999; // (the theoretical max value is less than 1 and larger than zero)
    for(int j=0;j<6-i;j++){
      if(current_value(j) > maxVal){
        maxVal=current_value(j);
        best_index[i]=current_index[j];
      }
    }
    int counter=0;
    for(int j=0;j<6-i;j++){
      if(current_value(j) != maxVal){
        keep_value(counter)=current_value(j);
        keep_index[counter]=current_index[j];
        counter++;
      }
    }
  }
  vOut.resize(6,4);
  for(int i=0;i<4;i++)
    for(int j=0;j<6;j++)
      vOut(j,i)=eigen_vectors(j,best_index[i]);
}


std::vector<double>
RockPhysicsInversion4D::MakeGaussKernel(double mean, double variance, double minf, double  df,int nf)
{
  std::vector<double> gaussKernel(nf);
  for(int i= 0;i<nf;i++)
  {
    double dist = (minf+(0.5+i)*df -mean);
    gaussKernel[i]=exp(-0.5*dist*dist/variance);
  }

  double sum=0.0;
  for(int i= 0;i<nf;i++)
    sum+=gaussKernel[i];

  for(int i= 0;i<nf;i++)
    gaussKernel[i]/=sum;
  return gaussKernel;
}

fftw_complex*
RockPhysicsInversion4D::MakeSmoothingFilter(double posteriorVariance,double  df)
{
  int cnfp      = nfp_/2 + 1;
  int rnfp      = 2*cnfp;

  fftw_real * gaussKernel = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnfp));
  fftw_complex* smoothingFilter  = reinterpret_cast<fftw_complex*>(gaussKernel );

  for(int i=0;i<cnfp;i++)
    gaussKernel[i] = float(exp(-0.5*(double(i)*df)*(double(i)*df)/ posteriorVariance)) ;

  for(int i=1;i<cnfp;i++) // note 1 is intentional
    gaussKernel[nfp_-i] = float(exp(-0.5*(double(i)*df)*(double(i)*df)/ posteriorVariance)) ;

  double sum=0.0;

  for(int i=0;i<nfp_;i++)
    sum+=gaussKernel[i];

  for(int i=0;i<nfp_;i++)
    gaussKernel[i]/=float(sum);

  rfftwnd_one_real_to_complex(fftplan1_,gaussKernel, smoothingFilter);

  return smoothingFilter;
}


void
RockPhysicsInversion4D::fillInTable( std::vector<std::vector<double> >  mSamp,std::vector<double>   rSamp,int tableInd)
{
  int nSamp = mSamp[0].size();

    for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(tableInd,j)->setAccessMode(FFTGrid::RANDOMACCESS);
    }

  for(int i=0;i<nSamp;i++)
  {
    NRLib::Vector m(6);
    NRLib::Vector f;
    for(int k=0;k<6;k++)
      m(k)=mSamp[k][i];
    f=m*v_;
    int i0 = GetLowerIndex(minf_(0),maxf_(0),nf_[0],f(0));
    int i1 = GetLowerIndex(minf_(1),maxf_(1),nf_[1],f(1));
    int i2 = GetLowerIndex(minf_(2),maxf_(2),nf_[2],f(2));
    int i3 = GetLowerIndex(minf_(3),maxf_(3),nf_[3],f(3));
    AddToGridValue(tableInd,i0,i1,i2,i3,rSamp[i]);
  }


  for (int j=0; j<nf_[0]; j++){
    meanRockPrediction_(tableInd,j)->endAccess();
  }
}

void
RockPhysicsInversion4D::smoothAllDirectionsAndNormalize()
{
  DivideAndSmoothTable(1,priorDistribution_,smoothingFilter_);

  for(int i=0;i<2;i++)
  for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(i,j)->setAccessMode(FFTGrid::RANDOMACCESS);
   }

  for(int i0=0;i0<nf_[0];i0++)
     for(int i1=0;i1<nf_[1];i1++)
        for(int i2=0;i2<nf_[2];i2++)
           for(int i3=0;i3<nf_[3];i3++)
           {
             double number = GetGridValue(1,i0,i1,i2,i3);
             double normalizing =GetGridValue(0,i0,i1,i2,i3);
             double value;
             if(normalizing<1e-5)
               value=RMISSING;
             else;
               value=number/normalizing;

             SetGridValue(1,i0,i1,i2,i3,value);
           }


  for(int i=0;i<2;i++)
    for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(i,j)->endAccess();
    }

}

void
RockPhysicsInversion4D::DivideAndSmoothTable(int tableInd,std::vector<std::vector<double> > priorDistribution, std::vector<fftw_complex*> smoothingFilter)
{
  for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(tableInd,j)->setAccessMode(FFTGrid::RANDOMACCESS);
   }

  int cnfp=nfp_/2+1;
  int rnfp=2*cnfp;

  fftw_real*    rTemp= static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnfp));
  fftw_complex* cTemp  = reinterpret_cast<fftw_complex*>(rTemp);

  double minDivisor = 1e-3;

  LogKit::LogFormatted(LogKit::Low,"\n Smoothing direction 1 of 4\n");
  float monitorSize = std::max(1.0f, static_cast<float>(nf_[0]*nf_[1]*nf_[2]*nf_[3])*0.02f);
  float nextMonitor = monitorSize;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");

  // divide and smooth direction1
  for(int i1=0;i1<nf_[1];i1++)
    for(int i2=0;i2<nf_[2];i2++)
      for(int i3=0;i3<nf_[3];i3++)
      {
        for(int i0=0;i0<nf_[0];i0++)
        {
          double divisor=std::max(minDivisor,priorDistribution_[0][i0]);
          rTemp[i0]=float(GetGridValue(tableInd,i0,i1,i2,i3)/ divisor);
        }
        for(int i0=nf_[0];i0<nfp_;i0++)
          rTemp[i0]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i0=0;i0<cnfp;i0++)
        {
          cTemp[i0].re=cTemp[i0].re*smoothingFilter[0][i0].re;
          cTemp[i0].im=cTemp[i0].im*smoothingFilter[0][i0].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i0=0;i0<nf_[0];i0++)
        {
          SetGridValue(tableInd,i0,i1,i2,i3, rTemp[i0]);
          if ( i1*nf_[2]*nf_[3]*nf_[0] +i2*nf_[3]*nf_[0]+ i3*nf_[0] + i0 + 1 >= static_cast<int>(nextMonitor)) {
            nextMonitor += monitorSize;
            printf("^");
          }
        }
      }

  LogKit::LogFormatted(LogKit::Low,"\n\n Smoothing direction 2 of 4\n");
  monitorSize = std::max(1.0f, static_cast<float>(nf_[0]*nf_[1]*nf_[2]*nf_[3])*0.02f);
  nextMonitor = monitorSize;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");

  // divide and smoothdirection2
  for(int i0=0;i0<nf_[0];i0++)
    for(int i2=0;i2<nf_[2];i2++)
      for(int i3=0;i3<nf_[3];i3++)
      {
        for(int i1=0;i1<nf_[1];i1++)
        {
          double divisor=std::max(minDivisor,priorDistribution_[1][i1]);
          rTemp[i1]=float(GetGridValue(tableInd,i0,i1,i2,i3)/divisor);
        }
        for(int i1=nf_[1];i1<nfp_;i1++)
          rTemp[i1]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i1=0;i1<cnfp;i1++)
        {
          cTemp[i1].re=cTemp[i1].re*smoothingFilter[1][i1].re;
          cTemp[i1].im=cTemp[i1].im*smoothingFilter[1][i1].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i1=0;i1<nf_[1];i1++)
        {
          SetGridValue(tableInd,i0,i1,i2,i3, rTemp[i1]);
          if ( i0*nf_[2]*nf_[3]*nf_[1] +i2*nf_[3]*nf_[1]+ i3*nf_[1] + i1 + 1 >= static_cast<int>(nextMonitor)) {
            nextMonitor += monitorSize;
            printf("^");
          }
        }
      }

  LogKit::LogFormatted(LogKit::Low,"\n\n Smoothing direction 3 of 4\n");
  monitorSize = std::max(1.0f, static_cast<float>(nf_[0]*nf_[1]*nf_[2]*nf_[3])*0.02f);
  nextMonitor = monitorSize;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");

  //  divide and smoothdirection 3
  for(int i0=0;i0<nf_[0];i0++)
    for(int i1=0;i1<nf_[1];i1++)
      for(int i3=0;i3<nf_[3];i3++)
      {
        for(int i2=0;i2<nf_[2];i2++)
        {
          double divisor=std::max(minDivisor,priorDistribution_[2][i2]);
          rTemp[i2]=float(GetGridValue(tableInd,i0,i1,i2,i3)/divisor);

        }
        for(int i2=nf_[2];i2<nfp_;i2++)
          rTemp[i2]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i2=0;i2<cnfp;i2++)
        {
          cTemp[i2].re=cTemp[i2].re*smoothingFilter[2][i2].re;
          cTemp[i2].im=cTemp[i2].im*smoothingFilter[2][i2].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i2=0;i2<nf_[2];i2++)
        {
          SetGridValue(tableInd,i0,i1,i2,i3, rTemp[i2]);
          if ( i0*nf_[1]*nf_[3]*nf_[2] +i1*nf_[3]*nf_[2]+ i3*nf_[2] + i2 + 1 >= static_cast<int>(nextMonitor)) {
            nextMonitor += monitorSize;
            printf("^");
          }
        }
      }

  LogKit::LogFormatted(LogKit::Low,"\n Smoothing last direction \n");
  monitorSize = std::max(1.0f, static_cast<float>(nf_[0]*nf_[1]*nf_[2]*nf_[3])*0.02f);
  nextMonitor = monitorSize;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");
 //  divide and smoothdirection 4
   for(int i0=0;i0<nf_[0];i0++)
    for(int i1=0;i1<nf_[1];i1++)
      for(int i2=0;i2<nf_[2];i2++)
      {
        for(int i3=0;i3<nf_[3];i3++)
        {
          double divisor=std::max(minDivisor,priorDistribution_[3][i3]);
          rTemp[i3]=float(GetGridValue(tableInd,i0,i1,i2,i3)/divisor);
        }
        for(int i3=nf_[3];i3<nfp_;i3++)
          rTemp[i3]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i3=0;i3<cnfp;i3++)
        {
          cTemp[i3].re=cTemp[i3].re*smoothingFilter[3][i3].re;
          cTemp[i3].im=cTemp[i3].im*smoothingFilter[3][i3].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i3=0;i3<nf_[3];i3++)
        {
          SetGridValue(tableInd,i0,i1,i2,i3, rTemp[i3]);
          if ( i0*nf_[1]*nf_[2]*nf_[3] +i1*nf_[2]*nf_[3]+ i2*nf_[3] + i3 + 1 >= static_cast<int>(nextMonitor)) {
            nextMonitor += monitorSize;
            printf("^");
          }
        }
      }

  for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(tableInd,j)->endAccess();
  }

}

void
RockPhysicsInversion4D::writeTableTofile(std::string fileName)
{

  for(int i0=0;i0<nf_[0];i0++)
    for(int i1=0;i1<nf_[1];i1++)
      for(int i2=0;i2<nf_[2];i2++)
        for(int i3=0;i3<nf_[3];i3++)
        {

        }
}