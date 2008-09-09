#include "src/crava.h"

#include "fft/include/rfftw.h"

#include "src/wavelet1D.h"
#include "src/corr.h"
#include "src/model.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/vario.h"
#include "src/welldata.h"
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/krigingadmin.h"
#include "src/faciesprob.h"
#include "src/definitions.h"
#include "lib/random.h"
#include "lib/lib_matr.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"

#include <assert.h>
#include <time.h>

Crava::Crava(Model * model)
{	
  int   i;
  model_ = model;
  nx_   = model->getBackAlpha()->getNx();
  ny_   = model->getBackAlpha()->getNy();
  nz_   = model->getBackAlpha()->getNz(); 
  nxp_  = model->getBackAlpha()->getNxp();
  nyp_  = model->getBackAlpha()->getNyp();
  nzp_  = model->getBackAlpha()->getNzp();

  lowCut_         = model->getModelSettings()->getLowCut();
  highCut_        = model->getModelSettings()->getHighCut();
  wnc_            = model->getModelSettings()->getWNC();     // white noise component see crava.h
  energyTreshold_ = model->getModelSettings()->getEnergyThreshold();
  theta_          = model->getModelSettings()->getAngle();
  ntheta_         = model->getModelSettings()->getNumberOfAngles();
  fileGrid_       = model->getModelSettings()->getFileGrid();
  outputFlag_     = model->getModelSettings()->getOutputFlag();
  krigingParams_  = model->getModelSettings()->getKrigingParameters();
  nWells_         = model->getModelSettings()->getNumberOfWells();
  nSim_           = model->getModelSettings()->getNumberOfSimulations();
  wells_          = model->getWells();
  simbox_         = model->getTimeSimbox();
  depthSimbox_    = model->getDepthSimbox();
  meanAlpha_      = model->getBackAlpha(); 
  meanBeta_       = model->getBackBeta();
  meanRho_        = model->getBackRho();
  
  Corr * corr     = model->getPriorCorrelations();

  //meanAlpha_->writeAsciiFile("priorA_");
  //meanBeta_->writeAsciiFile("priorB_");
  //meanRho_->writeAsciiFile("priorR_");
  random_ = model->getRandomGen();

  fftw_real * corrT; // =  fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)); 

  if(!model->getModelSettings()->getGenerateSeismic())
  {
    seisData_ = model->getSeisCubes();
    model->releaseGrids(); 
    parPointCov_ = corr->getVar0(); 
    parSpatialCorr_ = createFFTGrid();
    // NBNB   nzp_*0.001*corr->getdt() = T    lowCut = lowIntCut*domega = lowIntCut/T
    int lowIntCut = int(floor(lowCut_*(nzp_*0.001*corr->getdt()))); 
    // computes the integer whis corresponds to the low cut frequency.
    float corrGradI, corrGradJ;
    model->getCorrGradIJ(corrGradI, corrGradJ);
    corrT = parSpatialCorr_-> fillInParamCorr(corr,lowIntCut,corrGradI, corrGradJ);
    parSpatialCorr_->writeStormFile("corrdump",simbox_);
    errCorrUnsmooth_ = createFFTGrid();
    errCorrUnsmooth_->fillInErrCorr(corr,corrGradI,corrGradJ); 
  }
  else
  {
    model->releaseGrids();
    parSpatialCorr_ = NULL;
    corrT = NULL;
    errCorrUnsmooth_ = NULL;
  }

  if((outputFlag_ & (ModelSettings::PRIORCORRELATIONS + ModelSettings::CORRELATION)) > 0)
    dumpCorrT(corrT,corr->getdt());

  fprob_            = NULL;
  pKriging_         = NULL;
  covAlpha_         = NULL; 
  covBeta_          = NULL;
  covRho_           = NULL;
  covCrAlphaBeta_   = NULL;
  covCrAlphaRho_    = NULL;
  covCrBetaRho_     = NULL;
  kd_               = NULL;  
  empSNRatio_       = new float[ntheta_];
  theoSNRatio_      = new float[ntheta_];
  modelVariance_    = new float[ntheta_];
  signalVariance_   = new float[ntheta_];
  errorVariance_    = new float[ntheta_];
  dataVariance_     = new float[ntheta_];
  seisWavelet_      = model->getWavelets();
  scaleWarning_     = 0;
  scaleWarningText_ = new char[12*MAX_STRING*ntheta_]; 
  errThetaCov_      = new float*[ntheta_];  
  for(i=0;i<ntheta_;i++)
    errThetaCov_[i] = new float[ntheta_]; 
 
  // reality check: all dimensions involved match
  assert(meanBeta_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));
  assert(meanRho_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));

  for(i = 0; i< ntheta_; i++)
  {	
    if(!model->getModelSettings()->getGenerateSeismic())
      assert(seisData_[i]->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));  
    assert(seisWavelet_[i]->consistentSize(nzp_, nyp_, nxp_));
  }

  // all  parameters to be estimated are created.
  postAlpha_          = meanAlpha_;         //  The first five write over the input
  postBeta_           = meanBeta_;          //  to save memory
  postRho_            = meanRho_;           // 
  postCovAlpha_       = parSpatialCorr_;    // 
  postCovBeta_        = errCorrUnsmooth_;   // 
  postCovRho_         = createFFTGrid();    // grid not allocated
  postCrCovAlphaBeta_ = createFFTGrid();    // grid not allocate
  postCrCovAlphaRho_  = createFFTGrid();    // grid not allocated
  postCrCovBetaRho_   = createFFTGrid();    // grid not allocated
  postCrCovAlphaBeta_-> setType(FFTGrid::COVARIANCE);
  postCrCovAlphaRho_ -> setType(FFTGrid::COVARIANCE);
  postCrCovBetaRho_  -> setType(FFTGrid::COVARIANCE);
  postCovRho_        -> setType(FFTGrid::COVARIANCE);
  if((outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
  {
    meanAlpha2_ = copyFFTGrid(meanAlpha_);
    meanBeta2_ = copyFFTGrid(meanBeta_);
    meanRho2_ = copyFFTGrid(meanRho_);
  }
 
  if((outputFlag_ & ModelSettings::FACIESPROB) >0 || (outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
  {
    float * corrprior = new float[nzp_];
    //   FILE *test = fopen("test.dat","w");
    int refk;
    for(i = 0 ; i < nzp_ ; i++ )
    {
      if( i < nzp_/2+1)
        refk = i;
      else
        refk = nzp_ - i;
      if(refk < nz_)
        corrprior[i] = parSpatialCorr_->getRealValue(0,0,refk);
      else
        corrprior[i] = 0.0;
      
      //   fprintf(test,"%f\n",corrprior[i]);
    }
    //  fclose(test);

    Simbox *regularSimbox = new Simbox(simbox_);
    assert(typeid(simbox_->GetTopSurface()) == typeid(Surface));
    Surface * tsurf = new Surface(dynamic_cast<const NRLib2::RegularSurface<double> &>
                                 (simbox_->GetTopSurface()));

    regularSimbox->setDepth(tsurf, 0, simbox_->getlz(), simbox_->getdz());
    fprob_ = new FaciesProb(model->getModelSettings(),
                            fileGrid_, parPointCov_,corrprior, regularSimbox, *(const Simbox*)simbox_, 
                            nzp_, nz_, meanAlpha_, meanBeta_, meanRho_, random_, 
                            model->getModelSettings()->getPundef(), model->getPriorFacies());
    delete [] corrprior;
  }

  //meanAlpha_->writeAsciiRaw("alpha");
  //meanBeta_->writeAsciiRaw("beta");
  //meanRho_->writeAsciiRaw("rho");
  meanAlpha_ ->fftInPlace();
  meanBeta_  ->fftInPlace();
  meanRho_   ->fftInPlace();

  // parSpatialCorr_->writeFile("parSpatialCorr",simbox_);
  // parSpatialCorr_->writeAsciiFile("SpatialCorr");         //for debug
  // errCorrUnsmooth_->writeFile("errCorrUnsmooth",simbox_);
  // errCorrUnsmooth_->writeAsciiFile("ErrCorr");            //for debug

  A_ = model->getAMatrix();

  if(!model->getModelSettings()->getGenerateSeismic())
  {
    computeVariances(corrT,model);

    scaleWarning_ = checkScale();  // fills in scaleWarningText_ if needed.

    fftw_free(corrT);
  }

  

  if(!model->getModelSettings()->getGenerateSeismic())
  {
    if(simbox_->getIsConstantThick() == false)
    {
      divideDataByScaleWavelet();
    }
    parSpatialCorr_  ->fftInPlace();
    errCorrUnsmooth_ ->fftInPlace(); 

    for(i = 0; i < ntheta_ ; i++)
    {
      // sprintf(rawSName,"rawSeismic_divided_%i",int(seisData_[i]->getTheta()*180.0/PI+0.5));
      seisData_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
      // seisData_[i]->writeFile(rawSName,simbox_);
      seisData_[i]->fftInPlace();
      seisData_[i]->endAccess();
    }
  }

 }

Crava::~Crava()
{
  delete [] empSNRatio_;
  delete [] theoSNRatio_;
  delete [] modelVariance_;
  delete [] signalVariance_;
  delete [] errorVariance_;
  delete [] dataVariance_;
  if(fprob_!=NULL) delete fprob_;

  int i;
  for(i = 0; i < ntheta_; i++)
    delete [] A_[i] ;
  delete [] A_ ;

  for( i = 0;i<ntheta_;i++)
    delete[] errThetaCov_[i];
  delete [] errThetaCov_; 

  if (pKriging_ != NULL)       delete pKriging_;
  if (covAlpha_)               delete covAlpha_;
  if (covBeta_)                delete covBeta_;
  if (covRho_)                 delete covRho_;
  if (covCrAlphaBeta_)         delete covCrAlphaBeta_;
  if (covCrAlphaRho_)          delete covCrAlphaRho_;
  if (covCrBetaRho_)           delete covCrBetaRho_;
  if (kd_)                     delete kd_;
  
  if(postAlpha_!=NULL)         delete  postAlpha_ ;
  if(postBeta_!=NULL)          delete  postBeta_;
  if(postRho_!=NULL)           delete  postRho_ ;
  if(postCovAlpha_!=NULL)      delete  postCovAlpha_ ;
  if(postCovBeta_!=NULL)       delete  postCovBeta_ ;
  if(postCovRho_!=NULL)        delete  postCovRho_ ;
  if(postCrCovAlphaBeta_!=NULL)delete  postCrCovAlphaBeta_ ;
  if(postCrCovAlphaRho_!=NULL) delete  postCrCovAlphaRho_ ;
  if(postCrCovBetaRho_!=NULL)  delete  postCrCovBetaRho_;
  delete [] scaleWarningText_;
}


void
Crava::computeDataVariance(void)
{
  //
  // Compute variation in raw seismic
  //
  int rnxp = 2*(nxp_/2+1);
  for(int l=0 ; l < ntheta_ ; l++)
  {
    double  totvar = 0;
    long int ndata = 0;
    dataVariance_[l]=0.0;
    seisData_[l]->setAccessMode(FFTGrid::READ);
    for(int k=0 ; k<nzp_ ; k++) 
    {
      double tmpvar1 = 0;
      for(int j=0;j<nyp_;j++) 
      {
        double tmpvar2 = 0;
        for(int i=0; i <rnxp; i++)
        {
          float tmp=seisData_[l]->getNextReal();
          if(k < nz_ && j < ny_ &&  i < nx_ && tmp != 0.0)
          {
            tmpvar2 += double(tmp*tmp);
            ndata++;
          }
        }
        tmpvar1 += tmpvar2;
      }
      totvar += tmpvar1;
    }
    seisData_[l]->endAccess();
    dataVariance_[l] = float(totvar)/float(ndata);
  }
}

void
Crava::setupErrorCorrelation(Model * model)
{
  //
  //  Setup error correlation matrix
  //
  float * e = model->getModelSettings()->getNoiseEnergy() ;       

  for(int l=0 ; l < ntheta_ ; l++)
  {
    if (model->hasSignalToNoiseRatio()) 
    {
      empSNRatio_[l]    = e[l];
      errorVariance_[l] = dataVariance_[l]/empSNRatio_[l];
    }
  else
    {
      errorVariance_[l] = e[l]*e[l];
      empSNRatio_[l]    = float(dataVariance_[l]/errorVariance_[l]);
    }
    if (empSNRatio_[l] < 1.1f) 
    {
      LogKit::LogFormatted(LogKit::LOW,"\nThe empirical signal-to-noise ratio for angle stack %d is %7.1e. Ratios smaller than\n",l+1,empSNRatio_[l]);
      LogKit::LogFormatted(LogKit::LOW," 1 are illegal and CRAVA has to stop. CRAVA was for some reason not able to estimate\n");
      LogKit::LogFormatted(LogKit::LOW," this ratio reliably, and you must give it as input to the model file\n\n");
      exit(1);
    }
  }

  Vario * angularCorr = model->getModelSettings()->getAngularCorr();
  for(int i = 0; i < ntheta_; i++)
    for(int j = 0; j < ntheta_; j++) 
      errThetaCov_[i][j] = float(sqrt(errorVariance_[i])*sqrt(errorVariance_[j])
                                 *angularCorr->corr(theta_[i]-theta_[j],0));
}

void
Crava::computeVariances(fftw_real* corrT,
                        Model * model)
{
  computeDataVariance();
    
  setupErrorCorrelation(model);

  char fileName[MAX_STRING];
  Wavelet ** errorSmooth = new Wavelet*[ntheta_];
  float    * paramVar    = new float[ntheta_] ;
  float    * WDCorrMVar  = new float[ntheta_];

  for(int i=0 ; i < ntheta_ ; i++)
  {
    errorSmooth[i] = new Wavelet1D(seisWavelet_[i],Wavelet::FIRSTORDERFORWARDDIFF,1);
    sprintf(fileName,"Wavediff_%d.dat",i);
    errorSmooth[i]->printToFile(fileName);
  }

  // Compute variation in parameters 
  for(int i=0 ; i < ntheta_ ; i++)
  {
    paramVar[i]=0.0; 
    for(int l=0; l<3 ; l++) 
      for(int m=0 ; m<3 ; m++) 
        paramVar[i] += parPointCov_[l][m]*A_[i][l]*A_[i][m];
  }

  // Compute variation in wavelet
  for(int l=0 ; l < ntheta_ ; l++)
  {
    WDCorrMVar[l] = computeWDCorrMVar(errorSmooth[l],corrT);
  }

  // Compute signal and model variance and theoretical signal-to-noise-ratio
  for(int l=0 ; l < ntheta_ ; l++)
  {
    modelVariance_[l]  = WDCorrMVar[l]*paramVar[l];
    signalVariance_[l] = errorVariance_[l] + modelVariance_[l];
  }

  for(int l=0 ; l < ntheta_ ; l++)
  {
    if (model->getModelSettings()->getMatchEnergies()[l])
    {
      LogKit::LogFormatted(LogKit::LOW,"Matching syntethic and empirical energies:\n");
      float gain = sqrt((errorVariance_[l]/modelVariance_[l])*(empSNRatio_[l] - 1.0f));
      seisWavelet_[l]->scale(gain);
      if((outputFlag_ & ModelSettings::WAVELETS) > 0) 
      {
        sprintf(fileName,"Wavelet_EnergyMatched");
        seisWavelet_[l]->writeWaveletToFile(fileName, 1.0); // dt_max = 1.0;
      }
      modelVariance_[l] *= gain*gain;
      signalVariance_[l] = errorVariance_[l] + modelVariance_[l];
    }
    theoSNRatio_[l] = signalVariance_[l]/errorVariance_[l];
  }

  delete [] paramVar ;
  delete [] WDCorrMVar;
  for(int i=0;i<ntheta_;i++)
    delete errorSmooth[i];
  delete [] errorSmooth;
}


int
Crava::checkScale(void)
{
  char scaleWarning1[240];
  char scaleWarning2[240];
  char scaleWarning3[240];
  char scaleWarning4[240];
  strcpy(scaleWarning1,"The observed variability in seismic data is much larger than in the model.\n   Variance of: \n");
  strcpy(scaleWarning2,"The observed variability in seismic data is much less than in the model.\n   Variance of: \n");
  strcpy(scaleWarning3,"Small signal to noise ratio detected");
  strcpy(scaleWarning4,"Large signal to noise ratio detected");

  bool thisThetaIsOk;
  int  isOk=0;

  for(int l=0 ; l < ntheta_ ; l++)
  {
    thisThetaIsOk=true;

    if(( dataVariance_[l] > 4 * signalVariance_[l]) && thisThetaIsOk) //1 var 4
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 1;
        sprintf(scaleWarningText_,"Model inconsistency in angle %i for seismic data\n%s \n",
          int(theta_[l]/PI*180.0+0.5),scaleWarning1);
      }
      else
      {
        isOk = 1;
        sprintf(scaleWarningText_,"%sModel inconsistency in angle %i for seismic data\n%s \n",
          scaleWarningText_,int(theta_[l]/PI*180.0+0.5),scaleWarning1);
      }
    }
    if( (dataVariance_[l] < 0.1 * signalVariance_[l]) && thisThetaIsOk) //1 var 0.1
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 2;
        sprintf(scaleWarningText_,"Model inconsistency in angle %i for seismic data\n%s\n",
          int(theta_[l]/PI*180.0+0.5),scaleWarning2);
      }
      else
      {
        isOk = 2;
        sprintf(scaleWarningText_,"%sModel inconsistency in angle %i for seismic data\n%s\n",
          scaleWarningText_,int(theta_[l]/PI*180.0+0.5),scaleWarning2);
      }
    }
    if( (modelVariance_[l] < 0.02 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 3;
        sprintf(scaleWarningText_,"%s for angle %i.\n",
          scaleWarning3,int(theta_[l]/PI*180.0+0.5));
      }
      else
      {
        isOk = 3;
        sprintf(scaleWarningText_,"%s%s for angle %i.\n",
          scaleWarningText_,scaleWarning3,int(theta_[l]/PI*180.0+0.5));
      }
    }
    if( (modelVariance_[l] > 50.0 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 4;
        sprintf(scaleWarningText_,"%s for angle %i.\n",
          scaleWarning4,int(theta_[l]/PI*180.0+0.5) );
      }
      else
      {
        isOk = 4;
        sprintf(scaleWarningText_,"%s%s for angle %i.\n",
          scaleWarningText_,scaleWarning4,int(theta_[l]/PI*180.0+0.5)); 
      }
    }
  }
  return isOk;
}

void
Crava:: divideDataByScaleWavelet()
{
  int i,j,k,l,flag;
  float modW, modScaleW;

  fftw_real*    rData;
  fftw_real     tmp;
  fftw_complex* cData ;

  fftw_complex scaleWVal;
  rfftwnd_plan plan1,plan2; 

  rData  = (fftw_real*) fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)); 
  cData  = (fftw_complex* ) rData;
  Wavelet* localWavelet ;


  flag   = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan1  = rfftwnd_create_plan(1,&nzp_,FFTW_REAL_TO_COMPLEX,flag);
  plan2  = rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);

  // FILE * file1 = fopen("wavedumporig.txt","w");
  // FILE * file2 = fopen("wavedumpadj.txt","w");
  // FILE * file3 = fopen("seisdump.txt","w");
  // FILE * file4 = fopen("refldump.txt","w");
  for(l=0 ; l< ntheta_ ; l++ )
  {
    seisWavelet_[l]->fft1DInPlace();
    modW = seisWavelet_[l]->getNorm();
    modW *= modW/float(nzp_);

    seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
    for(i=0; i < nx_; i++)
      for(j=0; j< ny_; j++)
      {
        localWavelet = seisWavelet_[l]->getLocalWavelet(i,j);

        for(k=0;k<nzp_;k++)
        {
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/float(sqrt((float)nzp_));

          if(k > nz_)
          {
            float dist      = seisData_[l]->getDistToBoundary( k, nz_, nzp_);  
            rData[k] *= MAXIM(1-dist*dist,0);
          }

          // if(i == 0 && l == 0)
          // fprintf(file3,"%f ",rData[k]);
        }

        rfftwnd_one_real_to_complex(plan1,rData ,cData);

        double sf     = simbox_->getRelThick(i,j);
        double deltaF = static_cast<double>(nz_)*1000.0/(sf*simbox_->getlz()*static_cast<double>(nzp_));
        //double deltaF= double( nz_*1000.0)/ double(simbox_->getlz()*nzp_);

        for(k=0;k < (nzp_/2 +1);k++) // all complex values
        {
          scaleWVal =  localWavelet->getCAmp(k,static_cast<float>(sf));
          modScaleW =  scaleWVal.re * scaleWVal.re + scaleWVal.im * scaleWVal.im;
          /*         // note scaleWVal is acctually the value of the complex conjugate
          // (see definition of getCAmp)         
          wVal         =  seisWavelet_[l]->getCAmp(k);
          modW           =  wVal.re * wVal.re + wVal.im * wVal.im;
          //  Here we need only the modulus
          // ( see definition of getCAmp)         
          */
          if((modW > 0) && (deltaF*k < highCut_ ) && (deltaF*k > lowCut_ )) //NBNB frequency cleaning
          {
            float tolFac= 0.10f;
            if(modScaleW <  tolFac * modW)
              modScaleW =   float(sqrt(modScaleW*tolFac * modW));
            if(modScaleW == 0)
              modScaleW = 1;
            tmp           = cData[k].re * (scaleWVal.re/modScaleW) 
              -cData[k].im * (scaleWVal.im/modScaleW);
            cData[k].im   = cData[k].im * (scaleWVal.re/modScaleW) 
              +cData[k].re * (scaleWVal.im/modScaleW);
            cData[k].re   = tmp;
            //   if(i==0 && l ==0)
            //   {
            //     fprintf(file1,"%f ",scaleWVal.re);
            //    fprintf(file2,"%f ",sqrt(modScaleW));
            //  }
          }
          else
          {
            cData[k].im = 0.0f;
            cData[k].re = 0.0f;
          }
        }
        delete localWavelet;
        rfftwnd_one_complex_to_real(plan2 ,cData ,rData);
        for(k=0;k<nzp_;k++)
        {
          seisData_[l]->setRealValue(i,j,k,rData[k]/float(sqrt((float)nzp_)),true);
          //  if(i==0 && l==0)
          //    fprintf(file4,"%f ",rData[k]);
        }
      }
      char fName[200];
      if(ModelSettings::getDebugLevel() > 0)
      {
        sprintf(fName,"refl%d",l);
        seisData_[l]->writeFile(fName, simbox_);
      }
      LogKit::LogFormatted(LogKit::LOW,"Interpolating reflections in cube %d: ",l);
      seisData_[l]->interpolateSeismic(energyTreshold_);
      if(ModelSettings::getDebugLevel() > 0)
      {
        sprintf(fName,"reflInterpolated%d",l);
        seisData_[l]->writeFile(fName, simbox_);
      }
      seisData_[l]->endAccess();
  }
  //fclose(file1);
  //fclose(file2);
  //fclose(file3);
  //fclose(file4);

  fftw_free(rData);
  fftwnd_destroy_plan(plan1); 
  fftwnd_destroy_plan(plan2); 
}


void
Crava::multiplyDataByScaleWaveletAndWriteToFile(const char* typeName)
{
  int i,j,k,l,flag;

  fftw_real*    rData;
  fftw_real     tmp;
  fftw_complex* cData;
  fftw_complex scaleWVal;
  rfftwnd_plan plan1,plan2; 

  rData  = (fftw_real*) fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)); 
  cData  = (fftw_complex* ) rData;

  flag   = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan1  = rfftwnd_create_plan(1, &nzp_ ,FFTW_REAL_TO_COMPLEX,flag);
  plan2  = rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);

  Wavelet* localWavelet;

  for(l=0 ; l< ntheta_ ; l++ )
  {
    seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
    seisData_[l]->invFFTInPlace();

    for(i=0; i < nx_; i++)
      for(j=0; j< ny_; j++)
      {
        float sf = static_cast<float>(simbox_->getRelThick(i,j));

        for(k=0;k<nzp_;k++)
        {
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/float(sqrt((float)nzp_));
        }

        rfftwnd_one_real_to_complex(plan1,rData ,cData);
        localWavelet = seisWavelet_[l]->getLocalWavelet(i,j);

        for(k=0;k < (nzp_/2 +1);k++) // all complex values
        {
          scaleWVal    =  localWavelet->getCAmp(k,sf);
          // note scaleWVal is acctually the value of the complex conjugate
          // (see definition of getCAmp)         
          tmp           = cData[k].re * scaleWVal.re + cData[k].im * scaleWVal.im;
          cData[k].im   = cData[k].im * scaleWVal.re - cData[k].re * scaleWVal.im;
          cData[k].re   = tmp;
        }
        delete localWavelet;
        rfftwnd_one_complex_to_real(plan2 ,cData ,rData);
        for(k=0;k<nzp_;k++)
        {
          seisData_[l]->setRealValue(i,j,k,rData[k]/static_cast<float>(sqrt(static_cast<double>(nzp_))),true);
        }

      }
      char fName[200];
      int thetaDeg = int ( theta_[l]/PI*180 + 0.5 );
      sprintf(fName,"%s_%d",typeName,thetaDeg);
      seisData_[l]->writeFile(fName, simbox_);
      seisData_[l]->endAccess();
  }

  fftw_free(rData);
  fftwnd_destroy_plan(plan1); 
  fftwnd_destroy_plan(plan2); 
}



int 
Crava::computePostMeanResidAndFFTCov()
{
  int i,j,k,l,m;

  fftw_complex * kW          = new fftw_complex[ntheta_];

  fftw_complex * errMult1    = new fftw_complex[ntheta_];
  fftw_complex * errMult2    = new fftw_complex[ntheta_];
  fftw_complex * errMult3    = new fftw_complex[ntheta_];

  fftw_complex * ijkData     = new fftw_complex[ntheta_];
  fftw_complex * ijkDataMean = new fftw_complex[ntheta_];
  fftw_complex * ijkRes      = new fftw_complex[ntheta_];
  fftw_complex * ijkMean     = new fftw_complex[3];
  fftw_complex * ijkAns      = new fftw_complex[3]; 
  fftw_complex   kD,kD3;
  fftw_complex   ijkParLam;
  fftw_complex   ijkErrLam;
  fftw_complex   ijkTmp;

  fftw_complex**  K  = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    K[i] = new fftw_complex[3];

  fftw_complex**  KS  = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    KS[i] = new fftw_complex[3];

  fftw_complex**  KScc  = new fftw_complex*[3]; // cc - complex conjugate (and transposed)
  for(i = 0; i < 3; i++)
    KScc[i] = new fftw_complex[ntheta_];  

  fftw_complex**  parVar = new fftw_complex*[3];
  for(i = 0; i < 3; i++)
    parVar[i] = new fftw_complex[3];

  fftw_complex**  margVar = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    margVar[i] = new fftw_complex[ntheta_];

  fftw_complex**  errVar = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    errVar[i] = new fftw_complex[ntheta_];

  fftw_complex** reduceVar = new fftw_complex*[3];
  for(i = 0; i < 3; i++)
    reduceVar[i]= new fftw_complex[3];  

  Wavelet1D * diff1Operator = new Wavelet1D(Wavelet::FIRSTORDERFORWARDDIFF,nz_,nzp_,1);
  Wavelet1D * diff2Operator = new Wavelet1D(diff1Operator,Wavelet::FIRSTORDERBACKWARDDIFF,1);
  Wavelet1D * diff3Operator = new Wavelet1D(diff2Operator,Wavelet::FIRSTORDERCENTRALDIFF,1);

  diff1Operator->fft1DInPlace();
  diff2Operator->fft1DInPlace();
  diff3Operator->fft1DInPlace();

  Wavelet ** errorSmooth  = new Wavelet*[ntheta_];
  Wavelet ** errorSmooth2 = new Wavelet*[ntheta_];
  Wavelet ** errorSmooth3 = new Wavelet*[ntheta_];

  char* fileName = new char[2400] ;
  char* fileNameW = new char[2400] ;
  int thetaDeg;

  for(i = 0; i < ntheta_ ; i++)
  {  
    thetaDeg = int ( theta_[i]/PI*180 + 0.5 );
    errorSmooth[i]  = new Wavelet1D(seisWavelet_[i],Wavelet::FIRSTORDERFORWARDDIFF,1);
    errorSmooth2[i] = new Wavelet1D(errorSmooth[i], Wavelet::FIRSTORDERBACKWARDDIFF,1);
    errorSmooth3[i] = new Wavelet1D(errorSmooth2[i],Wavelet::FIRSTORDERCENTRALDIFF,1); 
    sprintf(fileName,"ErrorSmooth_%i",thetaDeg);
    errorSmooth3[i]->printToFile(fileName);
    errorSmooth[i]->fft1DInPlace();
    errorSmooth2[i]->fft1DInPlace();
    errorSmooth3[i]->fft1DInPlace();

    sprintf(fileName,"Wavelet_%i",thetaDeg);
    seisWavelet_[i]->printToFile(fileName);
    seisWavelet_[i]->fft1DInPlace();
    sprintf(fileNameW,"FourierWavelet_%i",thetaDeg);
    errorSmooth[i]->printToFile(fileNameW);
    seisData_[i]->setAccessMode(FFTGrid::READANDWRITE);
  }

  delete [] fileName;
  delete [] fileNameW;

  meanAlpha_       ->setAccessMode(FFTGrid::READANDWRITE);  //   Note
  meanBeta_        ->setAccessMode(FFTGrid::READANDWRITE);  //   the top five are over written
  meanRho_         ->setAccessMode(FFTGrid::READANDWRITE);  //   does not have the initial meaning.
  parSpatialCorr_  ->setAccessMode(FFTGrid::READANDWRITE);  //   after the prosessing
  errCorrUnsmooth_ ->setAccessMode(FFTGrid::READANDWRITE);  // 

  postCovRho_        ->createComplexGrid();
  postCrCovAlphaBeta_->createComplexGrid();
  postCrCovAlphaRho_ ->createComplexGrid();
  postCrCovBetaRho_  ->createComplexGrid();


  postCovRho_        ->setAccessMode(FFTGrid::WRITE);  
  postCrCovAlphaBeta_->setAccessMode(FFTGrid::WRITE);  
  postCrCovAlphaRho_ ->setAccessMode(FFTGrid::WRITE);  
  postCrCovBetaRho_  ->setAccessMode(FFTGrid::WRITE);  

  // Computes the posterior mean first  below the covariance is computed 
  // To avoid to many grids in mind at the same time
  double priorVarVp,justfactor;
  int cnxp  = nxp_/2+1;
  int cholFlag;
  //   long int timestart, timeend;
  //   time(&timestart);
  float realFrequency;

  for(k = 0; k < nzp_; k++)
  {  
    realFrequency = static_cast<float>((nz_*1000.0f)/(simbox_->getlz()*nzp_)*MINIM(k,nzp_-k)); // the physical frequency
    kD = diff1Operator->getCAmp(k);                   // defines content of kD  
    if (seisWavelet_[0]->getDim() == 1) { //1D-wavelet
      if( simbox_->getIsConstantThick() == true)
      {
        // defines content of K=WDA  
        fillkW(k,errMult1);                                    // errMult1 used as dummy
        lib_matrProdScalVecCpx(kD, errMult1, ntheta_);         // errMult1 used as dummy     
        lib_matrProdDiagCpxR(errMult1, A_, ntheta_, 3, K);     // defines content of (WDA)     K  

        // defines error-term multipliers  
        fillkWNorm(k,errMult1,seisWavelet_);               // defines input of  (kWNorm) errMult1
        fillkWNorm(k,errMult2,errorSmooth3);               // defines input of  (kWD3Norm) errMult2
        lib_matrFillOnesVecCpx(errMult3,ntheta_);          // defines content of errMult3   

      }
      else //simbox_->getIsConstantThick() == false
      {
        kD3 = diff3Operator->getCAmp(k);          // defines  kD3  

        // defines content of K = DA  
        lib_matrFillValueVecCpx(kD, errMult1, ntheta_);      // errMult1 used as dummy     
        lib_matrProdDiagCpxR(errMult1, A_, ntheta_, 3, K); // defines content of ( K = DA )   


        // defines error-term multipliers  
        lib_matrFillOnesVecCpx(errMult1,ntheta_);        // defines content of errMult1
        for(l=0; l < ntheta_; l++)
          errMult1[l].re /= seisWavelet_[l]->getNorm();  // defines content of errMult1

        lib_matrFillValueVecCpx(kD3,errMult2,ntheta_);   // defines content of errMult2
        for(l=0; l < ntheta_; l++) 
        {
          errMult2[l].re  /= errorSmooth3[l]->getNorm(); // defines content of errMult2
          errMult2[l].im  /= errorSmooth3[l]->getNorm(); // defines content of errMult2
        }
        fillInverseAbskWRobust(k,errMult3);              // defines content of errMult3
      } //simbox_->getIsConstantThick() 
    }
    else { //3D-wavelet
      if( simbox_->getIsConstantThick() == true) {
        for (j=0; j<nyp_; j++) {
          for (i=0; i<cnxp; i++) {
            for(l = 0; l < ntheta_; l++) {
              errMult1[l].re  = static_cast<float>( seisWavelet_[l]->getCAmp(k,j,i).re ); 
              errMult1[l].im  = static_cast<float>( seisWavelet_[l]->getCAmp(k,j,i).im ); 
            }
            lib_matrProdScalVecCpx(kD, errMult1, ntheta_); 
            lib_matrProdDiagCpxR(errMult1, A_, ntheta_, 3, K);  
            // defines error-term multipliers  
            for(l = 0; l < ntheta_; l++) {
              errMult1[l].re   = static_cast<float>(seisWavelet_[l]->getCAmp(k,j,i).re / seisWavelet_[l]->getNorm());
              errMult1[l].im   = static_cast<float>(seisWavelet_[l]->getCAmp(k,j,i).im / seisWavelet_[l]->getNorm());
            }       
          }
        }
        for(l = 0; l < ntheta_; l++) {
          errMult2[l].re   = static_cast<float>(errorSmooth3[l]->getCAmp(k).re / errorSmooth3[l]->getNorm());
          errMult2[l].im   = static_cast<float>(errorSmooth3[l]->getCAmp(k).im / errorSmooth3[l]->getNorm());
        }
        lib_matrFillOnesVecCpx(errMult3,ntheta_);
      }
      else {
        LogKit::LogFormatted(LogKit::LOW,"\nERROR: Not implemented inversion with 3D wavelet for non-constant simbox thickness\n");
        exit(1);
      }
    } //3D-wavelet
    for( j = 0; j < nyp_; j++)
    {
      for( i = 0; i < cnxp; i++)
      { 
        ijkMean[0] = meanAlpha_->getNextComplex();
        ijkMean[1] = meanBeta_ ->getNextComplex();
        ijkMean[2] = meanRho_  ->getNextComplex(); 

        for(l = 0; l < ntheta_; l++ )
        {
          ijkData[l] = seisData_[l]->getNextComplex();
          ijkRes[l]  = ijkData[l];
        }  

        ijkTmp       = parSpatialCorr_->getNextComplex();  
        ijkParLam.re = float ( sqrt(ijkTmp.re * ijkTmp.re));
        ijkParLam.im = 0.0;

        for(l = 0; l < 3; l++ )
          for(m = 0; m < 3; m++ )
          {        
            parVar[l][m].re = parPointCov_[l][m] * ijkParLam.re;       
            parVar[l][m].im = 0.0;
            // if(l!=m)
            //   parVar[l][m].re *= 0.75;  //NBNB OK DEBUG TEST
          }

          priorVarVp = parVar[0][0].re; 
          ijkTmp       = errCorrUnsmooth_->getNextComplex();  
          ijkErrLam.re = float( sqrt(ijkTmp.re * ijkTmp.re));
          ijkErrLam.im = 0.0;

          if(realFrequency > lowCut_*simbox_->getMinRelThick() &&  realFrequency < highCut_) // inverting only relevant frequencies
          {
            for(l = 0; l < ntheta_; l++ )
              for(m = 0; m < ntheta_; m++ )
              {        // Note we multiply kWNorm[l] and comp.conj(kWNorm[m]) hence the + and not a minus as in pure multiplication
                errVar[l][m].re  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult1[l].re *  errMult1[m].re +  errMult1[l].im *  errMult1[m].im)); 
                errVar[l][m].re += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult2[l].re *  errMult2[m].re +  errMult2[l].im *  errMult2[m].im)); 

                errVar[l][m].im  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult1[l].re * errMult1[m].im + errMult1[l].im * errMult1[m].re));
                errVar[l][m].im += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult2[l].re * errMult2[m].im + errMult2[l].im * errMult2[m].re));
                if(l==m)
                {
                  errVar[l][m].re += float( wnc_*errThetaCov_[l][m] * errMult3[l].re  * errMult3[l].re);
                  errVar[l][m].im   = 0.0;             
                }        
              }

              lib_matrProdCpx(K, parVar , ntheta_, 3 ,3, KS);              //  KS is defined here  
              lib_matrProdAdjointCpx(KS, K, ntheta_, 3 ,ntheta_, margVar); // margVar = (K)S(K)' is defined here      
              lib_matrAddMatCpx(errVar, ntheta_,ntheta_, margVar);         // errVar  is added to margVar = (WDA)S(WDA)'  + errVar 

              cholFlag=lib_matrCholCpx(ntheta_,margVar);                   // Choleskey factor of margVar is Defined

              if(cholFlag==0)
              { // then it is ok else posterior is identical to prior              

                lib_matrAdjoint(KS,ntheta_,3,KScc);                        //  WDAScc is adjoint of WDAS   
                lib_matrAXeqBMatCpx(ntheta_, margVar, KS, 3);              // redefines WDAS                 
                lib_matrProdCpx(KScc,KS,3,ntheta_,3,reduceVar);            // defines reduceVar
                //double hj=1000000.0;
                //if(reduceVar[0][0].im!=0)
                // hj = MAXIM(reduceVar[0][0].re/reduceVar[0][0].im,-reduceVar[0][0].re/reduceVar[0][0].im); //NBNB DEBUG
                lib_matrSubtMatCpx(reduceVar,3,3,parVar);                  // redefines parVar as the posterior solution

                lib_matrProdMatVecCpx(K,ijkMean, ntheta_, 3, ijkDataMean); //  defines content of ijkDataMean      
                lib_matrSubtVecCpx(ijkDataMean, ntheta_, ijkData);         //  redefines content of ijkData   

                lib_matrProdAdjointMatVecCpx(KS,ijkData,3,ntheta_,ijkAns); // defines ijkAns

                lib_matrAddVecCpx(ijkAns, 3,ijkMean);                      // redefines ijkMean 
                lib_matrProdMatVecCpx(K,ijkMean, ntheta_, 3, ijkData);     // redefines ijkData
                lib_matrSubtVecCpx(ijkData, ntheta_,ijkRes);               // redefines ijkRes      
              } 

              // quality control DEBUG
              if(priorVarVp*4 < ijkAns[0].re*ijkAns[0].re + ijkAns[0].re*ijkAns[0].re)
              {
                justfactor = sqrt(ijkAns[0].re*ijkAns[0].re + ijkAns[0].re*ijkAns[0].re)/sqrt(priorVarVp);
              }
          }
          postAlpha_->setNextComplex(ijkMean[0]);
          postBeta_ ->setNextComplex(ijkMean[1]);
          postRho_  ->setNextComplex(ijkMean[2]);
          postCovAlpha_->setNextComplex(parVar[0][0]);
          postCovBeta_ ->setNextComplex(parVar[1][1]);
          postCovRho_  ->setNextComplex(parVar[2][2]);
          postCrCovAlphaBeta_->setNextComplex(parVar[0][1]);
          postCrCovAlphaRho_ ->setNextComplex(parVar[0][2]);
          postCrCovBetaRho_  ->setNextComplex(parVar[1][2]); 

          for(l=0;l<ntheta_;l++)
            seisData_[l]->setNextComplex(ijkRes[l]);
      }
    }
  }

  //postCovAlpha_->writeAsciiRaw("covalphadump.txt");

  //  time(&timeend);
  // LogKit::LogFormatted(LogKit::LOW,"\n Core inversion finished after %ld seconds ***\n",timeend-timestart);
  // these does not have the initial meaning
  meanAlpha_ = NULL; //the content is taken care of by  postAlpha_
  meanBeta_  = NULL; //the content is taken care of by  postBeta_
  meanRho_   = NULL; //the content is taken care of by  postRho_
  parSpatialCorr_= NULL;   // the content is taken care of by  postCovAlpha_      
  errCorrUnsmooth_=NULL;   // the content is taken care of by  postCovBeta_

  postAlpha_   ->endAccess();  //   
  postBeta_    ->endAccess();  //   
  postRho_     ->endAccess();  //
  postCovAlpha_->endAccess();  //         
  postCovBeta_ ->endAccess();  // 
  postCovRho_        ->endAccess();
  postCrCovAlphaBeta_->endAccess();
  postCrCovAlphaRho_ ->endAccess();
  postCrCovBetaRho_  ->endAccess();  

  assert(postAlpha_->getIsTransformed());
  postAlpha_->invFFTInPlace();

  assert(postBeta_->getIsTransformed());
  postBeta_->invFFTInPlace();

  assert(postRho_->getIsTransformed());
  postRho_->invFFTInPlace();

  if((outputFlag_ & ModelSettings::PREDICTION) > 0)
  {
    if(krigingParams_ != 0) { 
      if (!pKriging_)
        initPostKriging();
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
      LogKit::LogFormatted(LogKit::LOW,"\n***                     Conditioning to wells                       ***"); 
      LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");
      pKriging_->KrigAll(*postAlpha_, *postBeta_, *postRho_);
    }
    writePars(postAlpha_, postBeta_, postRho_, -1);
  }

  char* fileNameS= new char[MAX_STRING];
  for(l=0;l<ntheta_;l++)
    seisData_[l]->endAccess();

  if((outputFlag_ & ModelSettings::RESIDUAL) > 0)
  {
    if(simbox_->getIsConstantThick() != true)
      multiplyDataByScaleWaveletAndWriteToFile("residuals");
    else
    {
      for(l=0;l<ntheta_;l++)
      {
        thetaDeg = int ( theta_[l]/PI*180 + 0.5 );
        sprintf(fileNameS,"residuals_%i",thetaDeg);
        seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
        seisData_[l]->invFFTInPlace();
        seisData_[l]->writeFile(fileNameS,simbox_);
        seisData_[l]->endAccess();
      }
    }
  }
  for(l=0;l<ntheta_;l++)
    delete seisData_[l];

  delete [] fileNameS;

  delete [] seisData_;

  delete [] kW;
  delete [] errMult1;
  delete [] errMult2;
  delete [] errMult3;
  delete [] ijkData;
  delete [] ijkDataMean; 
  delete [] ijkRes;
  delete [] ijkMean ;
  delete [] ijkAns;
  delete    diff1Operator;
  delete    diff2Operator;
  delete    diff3Operator;

  for(i = 0; i < ntheta_; i++)
  {
    delete[]  K[i];
    delete[]  KS[i];
    delete[]  margVar[i] ;
    delete[] errVar[i] ;
    delete errorSmooth[i];
    delete errorSmooth2[i];
    delete errorSmooth3[i];   
  }
  delete[] K;
  delete[] KS;
  delete[] margVar;
  delete[] errVar  ;
  delete[] errorSmooth;
  delete[] errorSmooth2;
  delete[] errorSmooth3;


  for(i = 0; i < 3; i++)
  {
    delete[] KScc[i];
    delete[] parVar[i] ;
    delete[] reduceVar[i];
  }
  delete[] KScc;
  delete[] parVar;  
  delete[] reduceVar;

  return(0);
}


int
Crava::simulate( RandomGen * randomGen)
{   
  if(nSim_>0)
  {
    assert( postCovAlpha_->getIsTransformed()  );
    assert( postCovBeta_->getIsTransformed()  );
    assert( postCovRho_->getIsTransformed()  );     ;
    assert( postCrCovAlphaBeta_->getIsTransformed()  );
    assert( postCrCovAlphaRho_ ->getIsTransformed()  );
    assert( postCrCovBetaRho_ ->getIsTransformed()  );

    int             simNr,i,j,k,l;
    fftw_complex ** ijkPostCov;
    fftw_complex *  ijkSeed;
    FFTGrid *       seed0;
    FFTGrid *       seed1;
    FFTGrid *       seed2;

    ijkPostCov = new fftw_complex*[3];
    for(l=0;l<3;l++)
      ijkPostCov[l]=new fftw_complex[3]; 

    ijkSeed = new fftw_complex[3];


    seed0 =  createFFTGrid();
    seed1 =  createFFTGrid();
    seed2 =  createFFTGrid();
    seed0->createComplexGrid();
    seed1->createComplexGrid();
    seed2->createComplexGrid();

    // long int timestart, timeend;

    for(simNr = 0; simNr < nSim_;  simNr++)
    {
      // time(&timestart);

      seed0->fillInComplexNoise(randomGen);    
      seed1->fillInComplexNoise(randomGen);      
      seed2->fillInComplexNoise(randomGen);

      postCovAlpha_      ->setAccessMode(FFTGrid::READ);  
      postCovBeta_       ->setAccessMode(FFTGrid::READ);  
      postCovRho_        ->setAccessMode(FFTGrid::READ);  
      postCrCovAlphaBeta_->setAccessMode(FFTGrid::READ);  
      postCrCovAlphaRho_ ->setAccessMode(FFTGrid::READ);  
      postCrCovBetaRho_  ->setAccessMode(FFTGrid::READ);  
      seed0 ->setAccessMode(FFTGrid::READANDWRITE); 
      seed1 ->setAccessMode(FFTGrid::READANDWRITE); 
      seed2 ->setAccessMode(FFTGrid::READANDWRITE); 

      int cnxp=nxp_/2+1;
      int cholFlag;
      for(k = 0; k < nzp_; k++)
        for(j = 0; j < nyp_; j++)
          for(i = 0; i < cnxp; i++)
          {
            ijkPostCov[0][0] = postCovAlpha_->getNextComplex();
            ijkPostCov[1][1] = postCovBeta_ ->getNextComplex();
            ijkPostCov[2][2] = postCovRho_  ->getNextComplex();
            ijkPostCov[0][1] = postCrCovAlphaBeta_ ->getNextComplex();
            ijkPostCov[0][2] = postCrCovAlphaRho_ ->getNextComplex(); 
            ijkPostCov[1][2] = postCrCovBetaRho_ ->getNextComplex();

            ijkPostCov[1][0].re =  ijkPostCov[0][1].re;
            ijkPostCov[1][0].im = -ijkPostCov[0][1].im;
            ijkPostCov[2][0].re =  ijkPostCov[0][2].re;
            ijkPostCov[2][0].im = -ijkPostCov[0][2].im;
            ijkPostCov[2][1].re =  ijkPostCov[1][2].re;
            ijkPostCov[2][1].im = -ijkPostCov[1][2].im;

            ijkSeed[0]=seed0->getNextComplex();
            ijkSeed[1]=seed1->getNextComplex();
            ijkSeed[2]=seed2->getNextComplex();

            cholFlag = lib_matrCholCpx(3,ijkPostCov);  // Choleskey factor of posterior covariance write over ijkPostCov
            if(cholFlag == 0)
            {
              lib_matrProdCholVec(3,ijkPostCov,ijkSeed); // write over ijkSeed
            }
            else
            {  
              for(l=0; l< 3;l++)
              {
                ijkSeed[l].re =0.0; 
                ijkSeed[l].im = 0.0;
              }

            }
            seed0->setNextComplex(ijkSeed[0]);
            seed1->setNextComplex(ijkSeed[1]);
            seed2->setNextComplex(ijkSeed[2]);      
          }


          postCovAlpha_->endAccess();  //         
          postCovBeta_ ->endAccess();  // 
          postCovRho_        ->endAccess();
          postCrCovAlphaBeta_->endAccess();
          postCrCovAlphaRho_ ->endAccess();
          postCrCovBetaRho_  ->endAccess(); 
          seed0->endAccess();
          seed1->endAccess();
          seed2->endAccess();

          // time(&timeend);
          // printf("Simulation in FFT domain in %ld seconds \n",timeend-timestart);
          // time(&timestart);

          /*    char* alphaName= new char[MAX_STRING];
          char* betaName = new char[MAX_STRING];
          char* rhoName  = new char[MAX_STRING];
          sprintf(alphaName,"sim_Vp_%i",simNr+1);
          sprintf(betaName,"sim_Vs_%i",simNr+1);
          sprintf(rhoName,"sim_D_%i",simNr+1);  
          */

          seed0->setAccessMode(FFTGrid::RANDOMACCESS);
          seed0->invFFTInPlace();
          seed0->add(postAlpha_);
          seed0->endAccess();

          seed1->setAccessMode(FFTGrid::RANDOMACCESS);
          seed1->invFFTInPlace(); 
          seed1->add(postBeta_);
          seed1->endAccess();

          seed2->setAccessMode(FFTGrid::RANDOMACCESS);
          seed2->invFFTInPlace();   
          seed2->add(postRho_);
          seed2->endAccess();

          //NBNB Bjørn: Sett inn krigingkall for seed0, seed1, seed2 her
          if(krigingParams_ != 0) { 
            if (!pKriging_)
              initPostKriging();
            LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************");
            LogKit::LogFormatted(LogKit::LOW,"\n***                     Conditioning to wells                       ***"); 
            LogKit::LogFormatted(LogKit::LOW,"\n***********************************************************************\n");
            pKriging_->KrigAll(*seed0, *seed1, *seed2);
          }
          writePars(seed0, seed1, seed2, simNr);  

          // time(&timeend);
          // printf("Back transform and write of simulation in %ld seconds \n",timeend-timestart);
    } 

    delete seed0;
    delete seed1;
    delete seed2;

    for(l=0;l<3;l++)
      delete  [] ijkPostCov[l];
    delete [] ijkPostCov;
    delete [] ijkSeed;

  }
  return(0);
}

int
Crava::computePostCov()
{
  assert(postCovAlpha_->getIsTransformed());
  postCovAlpha_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovAlpha_ ->invFFTInPlace();
  int i,j,k;
  float**  postCov  = new float*[3];
  for(i = 0; i < 3; i++)
    postCov[i] = new float[3];

  float*  postCorAlpha  = new float[nz_];
  float c0 = postCovAlpha_->getRealValue(0,0,0);
  for(k=0; k< nz_;k++)
  {
    postCorAlpha[k] = postCovAlpha_->getRealValue(0,0,k)/c0;
  }

  postCov[0][0]=postCovAlpha_->getRealValue(0,0,0);
  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCovAlpha_ ->writeFile("postCovVp",simbox_);
  postCovAlpha_ ->endAccess();
  delete postCovAlpha_;
  postCovAlpha_=NULL;

  assert(postCovBeta_->getIsTransformed());
  postCovBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovBeta_ ->invFFTInPlace();

  float*  postCorBeta  = new float[nz_];
  c0 = postCovBeta_->getRealValue(0,0,0);
  for(k=0; k< nz_;k++)
  {
    postCorBeta[k] = postCovBeta_->getRealValue(0,0,k)/c0;
  }

  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCovBeta_ ->writeFile("postCovVs",simbox_);
  postCov[1][1]=postCovBeta_->getRealValue(0,0,0);
  postCovBeta_ ->endAccess();
  delete postCovBeta_;
  postCovBeta_=NULL;

  assert(postCovRho_->getIsTransformed());
  postCovRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCovRho_ ->invFFTInPlace();

  float*  postCorRho  = new float[nz_];
  c0 = postCovRho_->getRealValue(0,0,0);
  for(k=0; k< nz_;k++)
  {
    postCorRho[k] = postCovRho_->getRealValue(0,0,k)/c0;
  }

  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCovRho_ ->writeFile("postCovD",simbox_);
  postCov[2][2]=postCovRho_->getRealValue(0,0,0);
  postCovRho_ ->endAccess();
  delete postCovRho_;
  postCovRho_=NULL;

  assert(postCrCovAlphaBeta_->getIsTransformed());
  postCrCovAlphaBeta_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaBeta_ ->invFFTInPlace();
  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCrCovAlphaBeta_ ->writeFile("postCrCovVpVs",simbox_);
  postCov[0][1]=postCrCovAlphaBeta_->getRealValue(0,0,0);
  postCov[1][0]=postCrCovAlphaBeta_->getRealValue(0,0,0);
  postCrCovAlphaBeta_ ->endAccess();
  delete postCrCovAlphaBeta_;
  postCrCovAlphaBeta_=NULL;

  assert(postCrCovAlphaRho_->getIsTransformed());
  postCrCovAlphaRho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovAlphaRho_ ->invFFTInPlace();
  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCrCovAlphaRho_->writeFile("postCrCovVpD",simbox_);
  postCov[2][0]=postCrCovAlphaRho_->getRealValue(0,0,0);
  postCov[0][2]=postCrCovAlphaRho_->getRealValue(0,0,0);
  postCrCovAlphaRho_->endAccess();
  delete postCrCovAlphaRho_;
  postCrCovAlphaRho_=NULL;

  assert(postCrCovBetaRho_->getIsTransformed());
  postCrCovBetaRho_->setAccessMode(FFTGrid::RANDOMACCESS);
  postCrCovBetaRho_->invFFTInPlace();
  if((outputFlag_ & ModelSettings::CORRELATION) > 0)
    postCrCovBetaRho_->writeFile("postCrCovVsD",simbox_);
  postCov[2][1]=postCrCovBetaRho_->getRealValue(0,0,0);
  postCov[1][2]=postCrCovBetaRho_->getRealValue(0,0,0);
  postCrCovBetaRho_->endAccess();
  delete postCrCovBetaRho_;
  postCrCovBetaRho_=NULL;

  char * fName = ModelSettings::makeFullFileName("PosteriorVar0",".dat");
  FILE* file=fopen(fName,"w");
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      fprintf(file,"%3e ",postCov[i][j]);
    }
    fprintf(file,"\n");
  }
  fclose(file);
  delete [] fName;

  fName = ModelSettings::makeFullFileName("PosteriorCorrTVp",".dat");
  file=fopen(fName,"w");
  for(k=0;k<nz_;k++)
  {
    fprintf(file,"%3e \n",postCorAlpha[k]);
  }
  fclose(file);
  delete [] fName;

  fName = ModelSettings::makeFullFileName("PosteriorCorrTVs",".dat");
  file=fopen(fName,"w");
  for(k=0;k<nz_;k++)
  {
    fprintf(file,"%3e \n",postCorBeta[k]);
  }
  fclose(file);
  delete [] fName;

  fName = ModelSettings::makeFullFileName("PosteriorCorrTRho",".dat");
  file=fopen(fName,"w");
  for(k=0;k<nz_;k++)
  {
    fprintf(file,"%3e \n",postCorRho[k]);
  }
  fclose(file);
  delete [] fName;

  for(i=0;i<3;i++)
    delete[] postCov[i];
  delete [] postCov;

  delete [] postCorAlpha;
  delete [] postCorBeta;
  delete [] postCorRho;
  return(0);
}


int 
Crava::computeSyntSeismic(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho)
{
  if(!Alpha->getIsTransformed()) Alpha->fftInPlace();
  if(!Beta->getIsTransformed()) Beta->fftInPlace();
  if(!Rho->getIsTransformed()) Rho->fftInPlace();
  int i,j,k,l;

  for(i=0;i<ntheta_;i++)
  {
//    seisWavelet_[i]->flipUpDown();
    if(seisWavelet_[i]->getIsReal() == true)
      seisWavelet_[i]->fft1DInPlace();
  }

  fftw_complex   kD;
  fftw_complex * kWD         = new fftw_complex[ntheta_];
  fftw_complex * ijkSeis     = new fftw_complex[ntheta_];
  fftw_complex * ijkParam    = new fftw_complex[3];
  fftw_complex** WDA         = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    WDA[i] = new fftw_complex[3];

  Wavelet * diffOperator     = new Wavelet1D(Wavelet::FIRSTORDERFORWARDDIFF,nz_,nzp_);
  diffOperator->fft1DInPlace();

  FFTGrid** seisData = new FFTGrid*[ntheta_];
  for(l=0;l<ntheta_;l++)
  {    
    seisData[l]= createFFTGrid();
    seisData[l]->setType(FFTGrid::DATA);
    seisData[l]->createComplexGrid();
    seisData[l]->setAccessMode(FFTGrid::WRITE);
  }

  int cnxp  = nxp_/2+1;
  
  if (seisWavelet_[0]->getDim() == 1) {
    Alpha->setAccessMode(FFTGrid::READ);
    Beta->setAccessMode(FFTGrid::READ);
    Rho->setAccessMode(FFTGrid::READ);
    for(k = 0; k < nzp_; k++)
    {
      kD = diffOperator->getCAmp(k);                   // defines content of kWD    
//        fillkW(k,kWD);                                   // defines content of kWD
      //FRODE - substituted prev. call with following loop - may08
      for(l = 0; l < ntheta_; l++)
      {
        kWD[l].re  = float( seisWavelet_[l]->getCAmp(k).re );// 
        kWD[l].im  = float( -seisWavelet_[l]->getCAmp(k).im );//
      } 
      //FRODE 
      lib_matrProdScalVecCpx(kD, kWD, ntheta_);        // defines content of kWD
      lib_matrProdDiagCpxR(kWD, A_, ntheta_, 3, WDA);  // defines content of WDA  
      for( j = 0; j < nyp_; j++)
      {
        for( i = 0; i < cnxp; i++)
        {
          ijkParam[0]=Alpha->getNextComplex();
          ijkParam[1]=Beta ->getNextComplex();
          ijkParam[2]=Rho  ->getNextComplex();

          lib_matrProdMatVecCpx(WDA,ijkParam,ntheta_,3,ijkSeis); //defines  ijkSeis 

          for(l=0;l<ntheta_;l++)
          {
            seisData[l]->setNextComplex(ijkSeis[l]);
          }
        }
      }
    }
  }
  else { //dim > 1
    Alpha->setAccessMode(FFTGrid::RANDOMACCESS);
    Beta->setAccessMode(FFTGrid::RANDOMACCESS);
    Rho->setAccessMode(FFTGrid::RANDOMACCESS);
    for (k=0; k<nzp_; k++) {
      kD = diffOperator->getCAmp(k);
      for (j=0; j<nyp_; j++) {
        for (i=0; i<cnxp; i++) {
           for(l = 0; l < ntheta_; l++) {
            kWD[l].re  = float( seisWavelet_[l]->getCAmp(k,j,i).re ); 
            kWD[l].im  = float( -seisWavelet_[l]->getCAmp(k,j,i).im );
          }
          lib_matrProdScalVecCpx(kD, kWD, ntheta_);        // defines content of kWD  
          lib_matrProdDiagCpxR(kWD, A_, ntheta_, 3, WDA);  // defines content of WDA  
          ijkParam[0]=Alpha->getComplexValue(i,j,k,true);
          ijkParam[1]=Beta ->getComplexValue(i,j,k,true);
          ijkParam[2]=Rho  ->getComplexValue(i,j,k,true);
          lib_matrProdMatVecCpx(WDA,ijkParam,ntheta_,3,ijkSeis);
          for(l=0;l<ntheta_;l++)
            seisData[l]->setComplexValue(i, j, k, ijkSeis[l], true);
        }
      }
    }
  }
  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  char* fileNameS= new char[MAX_STRING];
  int  thetaDeg;

  for(l=0;l<ntheta_;l++)
  { 
    seisData[l]->endAccess();
    seisData[l]->invFFTInPlace();

    thetaDeg = int (( theta_[l]/PI*180.0 + 0.5) );
    sprintf(fileNameS,"synt_seis_%i",thetaDeg);
    seisData[l]->writeFile(fileNameS,simbox_);
    delete seisData[l];
  }

  delete [] seisData;
  delete diffOperator;

  delete [] kWD;
  delete [] ijkSeis;
  delete [] ijkParam;
  for(i = 0; i < ntheta_; i++)
    delete [] WDA[i];
  delete [] WDA;

  delete [] fileNameS;
  return(0);
}

int 
Crava::computeAcousticImpedance(FFTGrid * Alpha, FFTGrid * Rho ,char * fileName, char * fileName2)
{   
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* prImpedance; 
  prImpedance  = createFFTGrid();
  prImpedance->setType(FFTGrid::PARAMETER);
  prImpedance->createRealGrid();
  prImpedance->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  prImpedance->getrsize(); 
  double ijkA, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkA + ijkR);
    prImpedance->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Rho->endAccess();

  prImpedance->endAccess();
  // prImpedance->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, prImpedance);
  delete prImpedance;

  return(0);
}

int 
Crava::computeShearImpedance(FFTGrid * Beta, FFTGrid * Rho ,char * fileName, char * fileName2)
{

  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);


  FFTGrid* shImpedance; 
  shImpedance  = createFFTGrid();
  shImpedance->setType(FFTGrid::PARAMETER);
  shImpedance->createRealGrid();
  shImpedance->setAccessMode(FFTGrid::WRITE);
  int i;
  int rSize =  shImpedance->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkB + ijkR);
    shImpedance->setNextReal(float( compVal));
  }

  Beta->endAccess();
  Rho->endAccess();

  shImpedance->endAccess(); 
  //shImpedance->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, shImpedance);
  delete shImpedance;

  return(0);
}


int  
Crava::computeVpVsRatio(FFTGrid * Alpha, FFTGrid * Beta,char * fileName, char * fileName2)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace(); 
  if(Beta->getIsTransformed())  Beta->invFFTInPlace();


  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);

  FFTGrid* ratioVpVs; 
  ratioVpVs = createFFTGrid();
  ratioVpVs->setType(FFTGrid::PARAMETER);
  ratioVpVs->createRealGrid();
  ratioVpVs->setAccessMode(FFTGrid::WRITE);
  int i;
  int rSize =  ratioVpVs->getrsize(); 
  double ijkA, ijkB, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    compVal = exp(ijkA - ijkB);
    ratioVpVs->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();

  ratioVpVs->endAccess();
  //ratioVpVs->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, ratioVpVs);
  delete ratioVpVs;

  return(0);
}

int 
Crava::computePoissonRatio(FFTGrid * Alpha, FFTGrid * Beta,char * fileName, char * fileName2)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);


  FFTGrid* poiRat; 
  poiRat  = createFFTGrid();
  poiRat->setType(FFTGrid::PARAMETER);
  poiRat->createRealGrid();
  poiRat->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  poiRat->getrsize(); 
  double ijkA, ijkB, compVal, vRatioSq;
  for(i=0; i  <  rSize; i++)
  {
    ijkA      = Alpha->getNextReal();
    ijkB      = Beta->getNextReal();
    vRatioSq  = exp(2*(ijkA-ijkB));
    compVal   = 0.5*(vRatioSq - 2)/(vRatioSq - 1);
    poiRat->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();

  poiRat->endAccess();
  //poiRat->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, poiRat);
  delete poiRat;

  return(0);
}


int  
Crava::computeLameMu(FFTGrid * Beta, FFTGrid * Rho,char * fileName, char * fileName2 )
{

  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);


  FFTGrid* mu; 
  mu  = createFFTGrid();
  mu->setType(FFTGrid::PARAMETER);
  mu->createRealGrid();
  mu->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  mu->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkR+2*ijkB-13.81551); // -13.81551 in the exponent divides by 1 000 000
    mu->setNextReal(float( compVal));
  }

  Beta->endAccess();
  Rho->endAccess();

  mu->endAccess();
  // mu->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, mu);

  delete mu;

  return(0);
}

int  
Crava::computeLameLambda(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho ,char * fileName, char * fileName2)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);


  FFTGrid* lambda; 
  lambda  = createFFTGrid();
  lambda->setType(FFTGrid::PARAMETER);
  lambda->createRealGrid();
  lambda->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  lambda->getrsize(); 
  double ijkA, ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkR)*(exp(2*ijkA-13.81551)-exp(2*ijkB-13.81551)); // -13.81551 in the exponent divides by 1 000 000
    lambda->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  lambda->endAccess();
  //lambda->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, lambda);

  delete lambda;

  return(0);
}

int  
Crava::computeLambdaRho(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho ,char * fileName, char * fileName2)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);


  FFTGrid* lambdaRho; 
  lambdaRho  = createFFTGrid();
  lambdaRho->setType(FFTGrid::PARAMETER);
  lambdaRho->createRealGrid();
  lambdaRho->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  lambdaRho->getrsize(); 
  double ijkA, ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(2.0*(ijkA +ijkR)-13.81551)-2.0*exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
    lambdaRho->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  lambdaRho->endAccess();
  //lambdaRho->writeFile(fileName,simbox_);
 
  writeToFile(fileName, fileName2, lambdaRho);
  delete lambdaRho;

  return(0);
}


int  
Crava::computeMuRho(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho ,char * fileName, char * fileName2)
{
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);


  FFTGrid* muRho; 
  muRho  = createFFTGrid();
  muRho->setType(FFTGrid::PARAMETER);
  muRho->createRealGrid();
  muRho->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  muRho->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
    muRho->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  muRho->endAccess();
  // muRho->writeFile(fileName,simbox_);
  writeToFile(fileName, fileName2, muRho);

  delete muRho;

  return(0);
}


float  
Crava::computeWDCorrMVar (Wavelet* WD ,fftw_real* corrT)
{
  float var = 0.0;
  int i,j,corrInd;

  for(i=0;i<nzp_;i++)
    for(j=0;j<nzp_;j++)
    {
      corrInd = MAXIM(i-j,j-i);
      var += WD->getRAmp(i)*corrT[corrInd]*WD->getRAmp(j);
    }

    // PAL060213: Temporary printing of CorrT. To be removed
    //
    //  FILE *file;
    //file=fopen("priorCorT.dat","w");
    //for(i=0;i<nzp_;i++)
    //{
    //fprintf(file,"%3e \n",corrT[i]);
    //}
    //fclose(file);

    return var;
}


void
Crava::fillkW(int k, fftw_complex* kW )
{
  int l;
  for(l = 0; l < ntheta_; l++)
  {
    kW[l].re  = float( seisWavelet_[l]->getCAmp(k).re );// 
    kW[l].im  = float( seisWavelet_[l]->getCAmp(k).im );// 
  }

}

void
Crava::fillkWNorm(int k, fftw_complex* kWNorm,Wavelet** wavelet )
{
  int l;
  for(l = 0; l < ntheta_; l++)
  {
    kWNorm[l].re   = float( wavelet[l]->getCAmp(k).re/wavelet[l]->getNorm());
    kWNorm[l].im   = float( wavelet[l]->getCAmp(k).im/wavelet[l]->getNorm());
  }
}

void
Crava::fillInverseAbskWRobust(int k, fftw_complex* invkW )
{
  int l;
  float modulus,modulusFine,maxMod;
  fftw_complex value;
  fftw_complex valueFine;
  for(l = 0; l < ntheta_; l++)
  {
    value  = seisWavelet_[l]->getCAmp(k);
    valueFine = seisWavelet_[l]->getCAmp(k,0.999f);// denser sampling of wavelet

    modulus      = value.re*value.re + value.im*value.im;
    modulusFine  = valueFine.re*valueFine.re + valueFine.im*valueFine.im;
    maxMod       = MAXIM(modulus,modulusFine);

    if(maxMod > 0.0)
    {
      invkW[l].re = float( 1.0/sqrt(maxMod) );
      invkW[l].im = 0.0f;    
    }
    else
    {
      invkW[l].re  =  seisWavelet_[l]->getNorm()*nzp_*nzp_*100.0f; // a big number
      invkW[l].im  =  0.0; // a big number
    }
  }
}


FFTGrid*            
Crava::createFFTGrid()
{
  FFTGrid* fftGrid;

  if(fileGrid_)
    fftGrid =  new FFTFileGrid(nx_,ny_,nz_,nxp_,nyp_,nzp_);
  else
    fftGrid =  new FFTGrid(nx_,ny_,nz_,nxp_,nyp_,nzp_);

  return(fftGrid);
}


FFTGrid*            
Crava::copyFFTGrid(FFTGrid * fftGridOld)
{
  FFTGrid* fftGrid;
  fftGrid =  new FFTGrid(fftGridOld);  
  return(fftGrid);
}


FFTFileGrid*            
Crava::copyFFTGrid(FFTFileGrid * fftGridOld)
{
  FFTFileGrid* fftGrid; 
  fftGrid =  new FFTFileGrid(fftGridOld);
  return(fftGrid);
}


void
Crava::writePars(FFTGrid * alpha, FFTGrid * beta, FFTGrid * rho, int simNum)
{
  char prefix[8];
  char postfix[8];
  if(simNum >= 0)
  {
    sprintf(prefix,"sim_");
    sprintf(postfix,"_%i",simNum+1);
  }
  else
  {
    sprintf(prefix,"pred_");
    sprintf(postfix,"%c",'\0');
  }
  char fileName[MAX_STRING];
  char fileName2[MAX_STRING];

  if((outputFlag_ & ModelSettings::MURHO) > 0)
  {
    sprintf(fileName,"%sMuRho%s",prefix, postfix);
    sprintf(fileName2,"%sMuRho_Depth%s",prefix, postfix);
    computeMuRho(alpha, beta, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::LAMBDARHO) > 0)
  {
    sprintf(fileName,"%sLambdaRho%s",prefix, postfix);
    sprintf(fileName2,"%sLambdaRho_Depth%s",prefix, postfix);
    computeLambdaRho(alpha, beta, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::LAMELAMBDA) > 0)
  {
    sprintf(fileName,"%sLameLambda%s",prefix, postfix);
    sprintf(fileName2,"%sLameLambda_Depth%s",prefix, postfix);
    computeLameLambda(alpha, beta, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::LAMEMU) > 0)
  {
    sprintf(fileName,"%sLameMu%s",prefix, postfix);
   sprintf(fileName2,"%sLameMu_Depth%s",prefix, postfix);
   computeLameMu(beta, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::POISSONRATIO) > 0)
  {
    sprintf(fileName,"%sPoissonRatio%s",prefix, postfix);
    sprintf(fileName2,"%sPoissonRatio_Depth%s",prefix, postfix);
    computePoissonRatio(alpha, beta, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::AI) > 0)
  {
    sprintf(fileName,"%sAI%s",prefix, postfix);
    sprintf(fileName2,"%sAI_Depth%s",prefix, postfix);
    computeAcousticImpedance(alpha, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::SI) > 0)
  {
    sprintf(fileName,"%sSI%s",prefix, postfix);
    sprintf(fileName2,"%sSI_Depth%s",prefix, postfix);
    computeShearImpedance(beta, rho, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::VPVSRATIO) > 0)
  {
    sprintf(fileName,"%sVpVsRatio%s",prefix, postfix);
    sprintf(fileName2,"%sVpVsRatio_Depth%s",prefix, postfix);
    computeVpVsRatio(alpha, beta, fileName, fileName2);
  }
  if((outputFlag_ & ModelSettings::VP) > 0)
  {
    sprintf(fileName,"%sVp%s",prefix, postfix);
    sprintf(fileName2,"%sVp_Depth%s",prefix, postfix);
    alpha->setAccessMode(FFTGrid::RANDOMACCESS);
    alpha->expTransf();
    //alpha->writeFile(fileName,simbox_);
    writeToFile(fileName,fileName2,alpha);
    if(simNum<0) //prediction, need grid unharmed.
      alpha->logTransf();
    alpha->endAccess();
  }
  if((outputFlag_ & ModelSettings::VS) > 0)
  {
    sprintf(fileName,"%sVs%s",prefix, postfix);
    sprintf(fileName2,"%sVs_Depth%s",prefix, postfix);
    beta->setAccessMode(FFTGrid::RANDOMACCESS);
    beta->expTransf();
    //    beta->writeFile(fileName,simbox_);
    writeToFile(fileName, fileName2, beta);
    if(simNum<0) //prediction, need grid unharmed.
      beta->logTransf();
    beta->endAccess();
  }
  if((outputFlag_ & ModelSettings::RHO) > 0)
  {
    sprintf(fileName,"%sRho%s",prefix, postfix);
   sprintf(fileName2,"%sRho_Depth%s",prefix, postfix);
    rho->setAccessMode(FFTGrid::RANDOMACCESS);
    rho->expTransf();
    //rho->writeFile(fileName,simbox_);
    writeToFile(fileName, fileName2, rho);
    if(simNum<0) //prediction, need grid unharmed.
      rho->logTransf();
    rho->endAccess();
  }
}

void
Crava::printEnergyToScreen()
{
  int i;
  LogKit::LogFormatted(LogKit::LOW,"                       ");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"  Seismic %4.1f ",theta_[i]/PI*180);
  LogKit::LogFormatted(LogKit::LOW,"\n----------------------");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"---------------");
  LogKit::LogFormatted(LogKit::LOW,"\nObserved data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",dataVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nModelled data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",signalVariance_[i]);
  //LogKit::LogFormatted(LogKit::LOW,"\nModel variance         :");
  //for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",modelVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nError variance         :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",errorVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nWavelet scale          :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %2.3e  ",seisWavelet_[i]->getScale());
  LogKit::LogFormatted(LogKit::LOW,"\nGiven S/N              :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %5.3f      ",empSNRatio_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nModelled S/N           :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %5.3f      ",theoSNRatio_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\n");
}

void
Crava::dumpCorrT(float* corrT,float dt)
{
  char * filename= ModelSettings::makeFullFileName("PriorCorrT",".dat");
  int i;
  FILE *file = fopen(filename, "w");

  fprintf(file,"%f\n",dt);
  for(i=0;i<nzp_;i++)
  {
    fprintf(file,"%f\n",corrT[i]);
  }
  fclose(file);
  delete [] filename;

}

void Crava::initPostKriging() {
  covAlpha_       = new CovGridSeparated(*getpostCovAlpha());
  covBeta_        = new CovGridSeparated(*getpostCovBeta());
  covRho_         = new CovGridSeparated(*getpostCovRho()); 
  covCrAlphaBeta_ = new CovGridSeparated(*getpostCrCovAlphaBeta());
  covCrAlphaRho_  = new CovGridSeparated(*getpostCrCovAlphaRho());
  covCrBetaRho_   = new CovGridSeparated(*getpostCrCovBetaRho());

  int type = 1; // 1 = full resolution
  kd_ = new KrigingData3D(wells_, nWells_, type);
  kd_->writeToFile("Raw");

  pKriging_ = new CKrigingAdmin(*getSimbox(), 
                                kd_->getData(),
                                kd_->getNumberOfData(),
                                *covAlpha_, *covBeta_, *covRho_, 
                                *covCrAlphaBeta_, *covCrAlphaRho_, *covCrBetaRho_, 
                                int(krigingParams_[0]));
}

void Crava::writeToFile(char * fileName1, char * fileName2, FFTGrid * grid) {

  if(!((outputFlag_ & ModelSettings::NOTIME)>0))
    grid->writeFile(fileName1,simbox_,1, model_->getModelSettings()->getSegyOffset());
  if(depthSimbox_!=NULL)
  {
    if(model_->getVelocity()!=NULL)
    {
      StormContGrid *mapping = model_->getMapping();
      writeDepthStormCube(model_->getVelocity(), mapping, fileName2);
    }
    else
      grid->writeFile(fileName2,depthSimbox_,0);
  }

}

void Crava::computeFaciesProb()
{
  if((outputFlag_ & ModelSettings::FACIESPROB) >0 || (outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
  {
    LogKit::LogFormatted(LogKit::LOW,"Start computing facies probability cubes ...\n");
    int relative;
    if((outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
      relative = 1;
    else relative = 0;

    if (simbox_->getdz() > 4.0f) { // Require this density for estimation of facies probabilities
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: The minimum sampling density is lower than 4.0. The FACIES PROBABILITIES\n");
      LogKit::LogFormatted(LogKit::LOW,"         generated by CRAVA are not reliable. To get more reliable probabilities    \n");
      LogKit::LogFormatted(LogKit::LOW,"         the number of layers must be increased.                                    \n");
    }

    fftw_real *postcova, *postcovb, *postcovr, *postcrab, *postcrar, *postcrbr;

    int rnzp = 2*(nzp_/2+1);
    postcova = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    postcovb = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    postcovr = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    postcrab = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    postcrar = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    postcrbr = (fftw_real*)  fftw_malloc(sizeof(float)*rnzp);
    
    if(postCovAlpha_->getIsTransformed()==true)
      postCovAlpha_->invFFTInPlace();
    
    if(postCovBeta_->getIsTransformed()==true)
      postCovBeta_->invFFTInPlace();
    if(postCovRho_->getIsTransformed()==true)
      postCovRho_->invFFTInPlace();

    if(postCrCovAlphaBeta_->getIsTransformed()==true)
      postCrCovAlphaBeta_->invFFTInPlace();
    if(postCrCovAlphaRho_->getIsTransformed()==true)
      postCrCovAlphaRho_->invFFTInPlace();
    if(postCrCovBetaRho_->getIsTransformed()==true)
      postCrCovBetaRho_->invFFTInPlace();

    int i;
    for(i=0;i<nzp_;i++)
    {
      int refk;
      if( i < nzp_/2+1)
        refk = i;
      else
        refk = nzp_ - i;
      if(refk < nz_)
      {
        postcova[i] = postCovAlpha_->getRealValue(0,0,refk);
        postcovb[i] = postCovBeta_->getRealValue(0,0,refk);
        postcovr[i] = postCovRho_->getRealValue(0,0,refk);
        postcrab[i] = postCrCovAlphaBeta_->getRealValue(0,0,refk);
        postcrar[i] = postCrCovAlphaRho_->getRealValue(0,0,refk);
        postcrbr[i] = postCrCovBetaRho_->getRealValue(0,0,refk);
      }
      else
      {
        postcova[i] = 0.0;
        postcovb[i] = 0.0;
        postcovr[i] = 0.0;
        postcrab[i] = 0.0;
        postcrar[i] = 0.0;
        postcrbr[i] = 0.0;
      }
    }
    
    WellData** ppWellData = (WellData**)(wells_);
    fprob_->filterWellLogs(ppWellData,nWells_,
                           postcova,postcovb,postcovr,
                           postcrab,postcrar,postcrbr, 
                           lowCut_, highCut_, relative);

    if(postCovAlpha_->getIsTransformed()==false)
      postCovAlpha_->fftInPlace();
    if(postCovBeta_->getIsTransformed()==false)
      postCovBeta_->fftInPlace();
    if(postCovRho_->getIsTransformed()==false)
      postCovRho_->fftInPlace();

    if(postCrCovAlphaBeta_->getIsTransformed()==false)
      postCrCovAlphaBeta_->fftInPlace();
    if(postCrCovAlphaRho_->getIsTransformed()==false)
      postCrCovAlphaRho_->fftInPlace();
    if(postCrCovBetaRho_->getIsTransformed()==false)
      postCrCovBetaRho_->fftInPlace();

    int nfac = model_->getModelSettings()->getNumberOfFacies();
    if(relative==0)
      fprob_->makeFaciesProb(nfac,postAlpha_,postBeta_, postRho_);
    else
    {
      meanAlpha2_->subtract(postAlpha_);
      meanAlpha2_->changeSign();
      meanBeta2_->subtract(postBeta_);
      meanBeta2_->changeSign();
      meanRho2_->subtract(postRho_);
      meanRho2_->changeSign();
      fprob_->makeFaciesProb(nfac,meanAlpha2_,meanBeta2_,meanRho2_);

    }
    fprob_->calculateConditionalFaciesProb(wells_, nWells_);

    FFTGrid *grid;
    char fileName[MAX_STRING];
    char fileName2[MAX_STRING];
    char postfix[20];
    LogKit::LogFormatted(LogKit::LOW,"\nProbability cubes done\n");
    if(relative==0)
    {
    for(i=0;i<nfac;i++)
    {
      grid = fprob_->getFaciesProb(i);
      sprintf(postfix,"_%s",model_->getModelSettings()->getFaciesName(i));
      sprintf(fileName,"FaciesProb%s",postfix);
      sprintf(fileName2,"FaciesProb_Depth%s",postfix);
      writeToFile(fileName,fileName2,grid);
    }
    }
    else
    {
    for(i=0;i<nfac;i++)
    {
      grid = fprob_->getFaciesProb(i);
      sprintf(postfix,"_%s",model_->getModelSettings()->getFaciesName(i));
      sprintf(fileName,"FaciesProbRelative%s",postfix);
      sprintf(fileName2,"FaciesProbRelative_Depth%s",postfix);
      writeToFile(fileName,fileName2,grid);
    }
    delete meanAlpha2_;
    delete meanBeta2_;
    delete meanRho2_;
    }


    fftw_free(postcova);
    fftw_free(postcovb);
    fftw_free(postcovr);
    fftw_free(postcrab);
    fftw_free(postcrar);
    fftw_free(postcrbr);
  }
}
void
Crava::writeDepthStormCube(FFTGrid *grid, StormContGrid *mapping, char * fileName)
{
  int i,j,k;
  int nz = depthSimbox_->getnz();
  float time, kindex;
  double x,y;
  StormContGrid outgrid(*mapping);
  for(i=0;i<nx_;i++)
  {
    x = simbox_->getx0()+i*simbox_->getdx();
    for(j=0;j<ny_;j++)
    {
      y = simbox_->gety0()+j*simbox_->getdy();
      for(k=0;k<nz;k++)
      {
        time = (*mapping)(i,j,k);
        kindex = float((time - simbox_->getTop(x,y))/simbox_->getdz());
        outgrid(i,j,k) = grid->getRealValueInterpolated(i,j,kindex);

      }
    }
  }
  outgrid.WriteToFile(fileName);
  

}
