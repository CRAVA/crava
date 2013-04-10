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
  nf_.resize(4,0);
  nf_[0] = 100;
  nf_[1] = 90;
  nf_[2] = 80;
  nf_[3] = 64;
  nfp_= 135;
  fftplan1_ = rfftwnd_create_plan(1, &nfp_, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftplan2_ = rfftwnd_create_plan(1, &nfp_, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  v_.resize(4,6); 
  SolveGEVProblem(priorCov,posteriorCov, v_);
  NRLib::Matrix tmp;
  NRLib::Matrix priorCovF;
  NRLib::Matrix posteriorCovF;
  NRLib::Vector priorMeanF; 

  priorMeanF    = v_*priorMean;  
  tmp           = v_*priorCov;
  priorCovF     = tmp*transpose(v_);
  tmp           = v_*posteriorCov;
  posteriorCovF = tmp*transpose(v_);

  int nSamp = mSamp[0].size(); // NBNB OK check 

  maxf_(4);
  minf_(4);

  for(int i=0;i<nSamp;i++)
  {
    NRLib::Vector m(6);
    NRLib::Vector f;
    for(int k=0;k<6;k++)
      m(k)=mSamp[k][i];
    f=v_*m;
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
    priorDistribution_[i]=MakeGaussKernel(priorMeanF(i),priorCovF(i,i),minf_(i),maxf_(i),nfp_);
    smoothingFilter_[i]=MakeSmoothingFilter(posteriorCovF(i,i),df);  
  }

  std::vector<double> dummy(1,0.0);

  allocatePredictionTables( );
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
  int nx,ny,nz,rnxp,nyp,nzp;
  nx=mu_static_[0]->getNx();
  ny=mu_static_[0]->getNy();
  nz=mu_static_[0]->getNy();

  rnxp=mu_static_[0]->getCNxp()*2;
  nyp=mu_static_[0]->getNyp();
  nzp=mu_static_[0]->getNyp();

  FFTGrid* prediction= new FFTGrid(nx,ny,nz,nx,ny,nz); // no need for padding
  
  for(int i=0;i<3;i++)
  {
    mu_static_[i]->setAccessMode(FFTGrid::READ);
    mu_static_[i]->setAccessMode(FFTGrid::READ);
  }

  prediction->setAccessMode(FFTGrid::RANDOMACCESS);

  NRLib::Vector m(6);

  for(int i=0;i<rnxp;i++)
    for(int j=0;j<nyp;j++)
      for(int k=0;k<nzp;k++)
      {
        m(0)=mu_static_[0]->getNextReal();
        m(1)=mu_static_[1]->getNextReal();
        m(2)=mu_static_[2]->getNextReal();
        m(3)=mu_dynamic_[0]->getNextReal();
        m(4)=mu_dynamic_[1]->getNextReal();
        m(5)=mu_dynamic_[2]->getNextReal();

        if(i<nx && j<ny && k<nz)
        {
          float value;
          NRLib::Vector f;
          f=v_*m;
          value = float( getPredictedValue(f) );
          prediction->setRealValue(i,j,k,value);
        }
      }

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
  for(int i1=0;i1<1;i1++)
    for(int i2=0;i2<1;i2++)
      for(int i3=0;i3<1;i3++)
        for(int i4=0;i4<1;i4++)
        {
          w=wLoHi[i1][0]*wLoHi[i2][1]*wLoHi[i3][2]*wLoHi[i4][3];
          value+=w*GetGridValue(1,indLoHi[i1][0],indLoHi[i2][1],indLoHi[i3][2],indLoHi[i4][3]);
        }
return value;
}

double
RockPhysicsInversion4D::GetGridValue(int TableNr,int i1,int i2,int i3,int i4)
{
  return meanRockPrediction_(TableNr,i1)->getRealValue(i2,i3,i4);
}

void
RockPhysicsInversion4D::SetGridValue(int TableNr,int i1,int i2,int i3,int i4,double value)
{
  meanRockPrediction_(TableNr,i1)->setRealValue(i2,i3,i4,value);
}

void
RockPhysicsInversion4D::AddToGridValue(int TableNr,int i1,int i2,int i3,int i4,double value)
{
  double newValue=meanRockPrediction_(TableNr,i1)->getRealValue(i2,i3,i4)+value;
  meanRockPrediction_(TableNr,i1)->setRealValue(i2,i3,i4,newValue);
}


void
RockPhysicsInversion4D::SolveGEVProblem(NRLib::Matrix sigma_priorIn,
                                        NRLib::Matrix sigma_posteriorIn,
                                        NRLib::Matrix & vOut){
  //Compute transforms v1, v2, v3 and v4  ----------------------------------------

  // Matrix inversion of sigma_prior
  double **sigma_prior_inv = new double *[6];
  double **sigma_prior = new double *[6];
  double **sigma_posterior = new double *[6];

  for(int i=0; i<6; i++)
  {
    sigma_prior_inv[i] = new double [6];
    sigma_prior[i] = new double [6];
    sigma_posterior[i] = new double [6];
  }

  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
    {
      if(i==j)
        sigma_prior_inv[i][j] = 1.0;
      else
        sigma_prior_inv[i][j] = 0.0;
      
      sigma_prior[i][j]     =sigma_priorIn(i,j);
      sigma_posterior[i][j] =sigma_posteriorIn(i,j);
    }

  lib_matrCholR(6, sigma_prior);
  lib_matrAXeqBMatR(6, sigma_prior, sigma_prior_inv, 6);

  // compute product sigma_prior_inv * sigma_posterior

  double **product_mat = new double *[6];

  for (int i=0; i<6; i++)
    product_mat[i] = new double[6];

  lib_matr_prod(sigma_prior_inv,sigma_posterior,6,6,6,product_mat);

  // compute eigenvalues of the product sigma_prior_inv * sigma_posterior

  int     * error       = new int[1];
  double  * eigval      = new double[6];
  double ** eigvec      = new double *[6];

  for(int i=0;i<6;i++)
  {
    eigvec[i]      = new double[6];
  }

  lib_matr_eigen(product_mat,6,eigvec,eigval,error);

  std::vector<int> index_keep;


  // Find index of four largest eigenvalues
  
  double  previous_max=2.0; // maxeigenvalue is 1.0;

  for(int k=0;k<4;k++)
  {
     int max_index = 0;
     double max_eigvalue = 0.0;
     
     for(int i=0; i<6; i++){
       if (max_eigvalue < eigval[i] &&  eigval[i]<previous_max ){
         max_index = i;
         max_eigvalue = eigval[i];
       }
     }
     index_keep.push_back(max_index);
     previous_max= max_eigvalue;
  }

  // The vector index_keep should contain 4 integers
  assert (index_keep.size() == 4);


  // fetch the corresponding eigenvectors into v
  for(int i = 0; i<6; i++){
    vOut(0,i) = eigvec[index_keep[0]][i];
    vOut(1,i) = eigvec[index_keep[1]][i];
    vOut(2,i) = eigvec[index_keep[2]][i];
    vOut(3,i) = eigvec[index_keep[3]][i];
  }


  delete [] error;
  delete [] eigval;
  for(int i=0;i<6;i++){
    delete [] product_mat[i];
    delete [] eigvec[i];
    delete [] sigma_prior_inv[i];
    delete [] sigma_prior[i];
    delete [] sigma_posterior[i];

  }
  delete [] product_mat;
  delete [] sigma_prior_inv;
  delete [] sigma_prior;
  delete [] sigma_posterior;
  delete [] eigvec;
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
    gaussKernel[i] = exp(-0.5*(double(i)*df)*(double(i)*df)/ posteriorVariance) ;

  for(int i=1;i<cnfp;i++) // note 1 is intentional 
    gaussKernel[nfp_-i] = exp(-0.5*(double(i)*df)*(double(i)*df)/ posteriorVariance) ;
 
  double sum=0.0;
  
  for(int i=0;i<nfp_;i++)
    sum+=gaussKernel[i];

  for(int i=0;i<nfp_;i++)
    gaussKernel[i]/=sum;
  
  rfftwnd_one_real_to_complex(fftplan1_,gaussKernel, smoothingFilter);
  
  return smoothingFilter;
}


void
RockPhysicsInversion4D::fillInTable( std::vector<std::vector<double> >  mSamp,std::vector<double>   rSamp,int tableInd)
{
  int nSamp = mSamp[0].size(); // NBNB OK check 

    for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(tableInd,j)->setAccessMode(FFTGrid::RANDOMACCESS);
    }

  for(int i=0;i<nSamp;i++)
  {
    NRLib::Vector m(6);
    NRLib::Vector f;
    for(int k=0;k<6;k++)
      m(k)=mSamp[k][i];
    f=v_*m;
    int i1 = GetLowerIndex(minf_(0),maxf_(0),nf_[0],f(0));
    int i2 = GetLowerIndex(minf_(1),maxf_(1),nf_[0],f(1));
    int i3 = GetLowerIndex(minf_(2),maxf_(2),nf_[0],f(2));
    int i4 = GetLowerIndex(minf_(3),maxf_(3),nf_[0],f(3));
    if(tableInd>0)
      AddToGridValue(1,i1,i2,i3,i4,rSamp[i]);
    else
      AddToGridValue(0,i1,i2,i3,i4,1.0);
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

  for(int i1=0;i1<nf_[0];i1++)
     for(int i2=0;i2<nf_[0];i2++)
        for(int i3=0;i3<nf_[0];i3++)
           for(int i4=0;i4<nf_[0];i4++)
           {
             double number = GetGridValue(1,i1,i2,i3,i4);
             double normalizing =GetGridValue(0,i1,i2,i3,i4);
             SetGridValue(1,i1,i2,i3,i4,number/normalizing);
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

  float maxFFTNumber=1000.0f;


  // divide and smooth direction1
  for(int i2=0;i2<nf_[1];i2++)
    for(int i3=0;i3<nf_[2];i3++)
      for(int i4=0;i4<nf_[3];i4++)
      {
        for(int i1=0;i1<nf_[0];i1++)
        {
          rTemp[i1]=GetGridValue(tableInd,i1,i2,i3,i4)/priorDistribution_[0][i1];
          rTemp[i1] = std::min(maxFFTNumber,rTemp[i1]);
          rTemp[i1] = std::max(-maxFFTNumber,rTemp[i1]);
        }
        for(int i1=nf_[0];i1<nfp_;i1++)
          rTemp[i1]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i1=0;i1<cnfp;i1++)
        {
          cTemp[i1].re=cTemp[i1].re*smoothingFilter[0][i1].re;
          cTemp[i1].im=cTemp[i1].im*smoothingFilter[0][i1].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i1=0;i1<nf_[0];i1++)
          SetGridValue(tableInd,i1,i2,i3,i4, rTemp[i1]);

      }

  
  // divide and smoothdirection2
  for(int i1=0;i1<nf_[0];i1++)
    for(int i3=0;i3<nf_[2];i3++)
      for(int i4=0;i4<nf_[3];i4++)
      {
        for(int i2=0;i2<nf_[1];i2++)
        {
          rTemp[i2]=GetGridValue(tableInd,i1,i2,i3,i4)/priorDistribution_[1][i2];
          rTemp[i2] = std::min(maxFFTNumber,rTemp[i2]);
          rTemp[i2] = std::max(-maxFFTNumber,rTemp[i2]);
        }
        for(int i2=nf_[1];i2<nfp_;i2++)
          rTemp[i2]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i2=0;i2<cnfp;i2++)
        {
          cTemp[i2].re=cTemp[i2].re*smoothingFilter[1][i2].re;
          cTemp[i2].im=cTemp[i2].im*smoothingFilter[1][i2].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i2=0;i2<nf_[1];i2++)
          SetGridValue(tableInd,i1,i2,i3,i4, rTemp[i2]);

      }

  //  divide and smoothdirection 3
    for(int i1=0;i1<nf_[0];i1++)
      for(int i2=0;i2<nf_[1];i2++)
        for(int i4=0;i4<nf_[3];i4++)
      {
        for(int i3=0;i3<nf_[2];i3++)
        {
          rTemp[i3]=GetGridValue(tableInd,i1,i2,i3,i4)/priorDistribution_[2][i3];
          rTemp[i3] = std::min(maxFFTNumber,rTemp[i3]);
          rTemp[i3] = std::max(-maxFFTNumber,rTemp[i3]);
        }
        for(int i3=nf_[1];i3<nfp_;i3++)
          rTemp[i3]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i3=0;i3<cnfp;i3++)
        {
          cTemp[i3].re=cTemp[i3].re*smoothingFilter[2][i3].re;
          cTemp[i3].im=cTemp[i3].im*smoothingFilter[2][i3].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i3=0;i2<nf_[2];i3++)
          SetGridValue(tableInd,i1,i2,i3,i4, rTemp[i3]);

      }   
  
 //  divide and smoothdirection 4
   for(int i1=0;i1<nf_[0];i1++)
    for(int i2=0;i2<nf_[1];i2++)
      for(int i3=0;i3<nf_[2];i3++)
      {
        for(int i4=0;i4<nf_[3];i4++)
        {
          rTemp[i4]=GetGridValue(tableInd,i1,i2,i3,i4)/priorDistribution_[3][i4];
          rTemp[i4] = std::min(maxFFTNumber,rTemp[i4]);
          rTemp[i4] = std::max(-maxFFTNumber,rTemp[i4]);
        }
        for(int i4=nf_[3];i4<nfp_;i4++)
          rTemp[i4]=0.0f;

        rfftwnd_one_real_to_complex(fftplan1_,rTemp,cTemp);

        for(int i4=0;i4<cnfp;i4++)
        {
          cTemp[i4].re=cTemp[i4].re*smoothingFilter[3][i4].re;
          cTemp[i4].im=cTemp[i4].im*smoothingFilter[3][i4].re;
        }

        rfftwnd_one_complex_to_real(fftplan2_,cTemp,rTemp);
        for(int i4=0;i4<nf_[3];i4++)
          SetGridValue(tableInd,i1,i2,i3,i4, rTemp[i4]);
      }

  for (int j=0; j<nf_[0]; j++){
      meanRockPrediction_(tableInd,j)->endAccess();
  }

}
