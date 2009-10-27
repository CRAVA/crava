#include "src/crava.h"

#include "fft/include/rfftw.h"

#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
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
#include "src/gridmapping.h"
#include "src/filterwelllogs.h"
#include "src/parameteroutput.h"
#include "src/timings.h"
#include "src/spatialwellfilter.h"
#include "lib/timekit.hpp"
#include "lib/random.h"
#include "lib/utils.h"
#include "lib/lib_matr.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/grid/grid2d.hpp"


#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>

Crava::Crava(Model * model, SpatialWellFilter * spatwellfilter)
{
  Utils::writeHeader("Building Stochastic Model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  model_             = model;
  nx_                = model->getBackAlpha()->getNx();
  ny_                = model->getBackAlpha()->getNy();
  nz_                = model->getBackAlpha()->getNz(); 
  nxp_               = model->getBackAlpha()->getNxp();
  nyp_               = model->getBackAlpha()->getNyp();
  nzp_               = model->getBackAlpha()->getNzp();
  lowCut_            = model->getModelSettings()->getLowCut();
  highCut_           = model->getModelSettings()->getHighCut();
  wnc_               = model->getModelSettings()->getWNC();     // white noise component see crava.h
  energyTreshold_    = model->getModelSettings()->getEnergyThreshold();
  ntheta_            = model->getModelSettings()->getNumberOfAngles();
  fileGrid_          = model->getModelSettings()->getFileGrid();
  outputFlag_        = model->getModelSettings()->getGridOutputFlag();
  writePrediction_   = model->getModelSettings()->getWritePrediction();
  krigingParameter_  = model->getModelSettings()->getKrigingParameter();
  nWells_            = model->getModelSettings()->getNumberOfWells();
  nSim_              = model->getModelSettings()->getNumberOfSimulations();
  wells_             = model->getWells();
  simbox_            = model->getTimeSimbox(); 
  meanAlpha_         = model->getBackAlpha(); 
  meanBeta_          = model->getBackBeta();
  meanRho_           = model->getBackRho();
  correlations_      = model->getCorrelations();
  random_            = model->getRandomGen();
  seisWavelet_       = model->getWavelets();
  A_                 = model->getAMatrix();
  postAlpha_         = meanAlpha_;         // Write over the input to save memory
  postBeta_          = meanBeta_;          // Write over the input to save memory
  postRho_           = meanRho_;           // Write over the input to save memory
  fprob_             = NULL;
  thetaDeg_          = new float[ntheta_];
  empSNRatio_        = new float[ntheta_];
  theoSNRatio_       = new float[ntheta_];
  modelVariance_     = new float[ntheta_];
  signalVariance_    = new float[ntheta_];
  errorVariance_     = new float[ntheta_];
  dataVariance_      = new float[ntheta_];
  scaleWarning_      = 0;
  scaleWarningText_  = new char[12*MAX_STRING*ntheta_]; 
  errThetaCov_       = new double*[ntheta_]; 
  sigmamdnew_        = NULL;
  for(int i=0;i<ntheta_;i++)
    errThetaCov_[i] = new double[ntheta_]; 
  for(int i=0;i<ntheta_;i++)
    thetaDeg_[i] = static_cast<float>(model->getModelSettings()->getAngle(i)*180.0/M_PI); 
  
  fftw_real * corrT = NULL; // =  fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)); 

  // Double-use grids to save memory
  FFTGrid * parSpatialCorr  = NULL;  // Parameter correlation
  FFTGrid * errCorrUnsmooth = NULL;  // Error correlation

  if(!model->getModelSettings()->getGenerateSeismic())
  {
    seisData_       = model->getSeisCubes();
    model->releaseGrids(); 
    bool useFileGrid = (fileGrid_ == 1);
    correlations_->createPostGrids(nx_,ny_,nz_,nxp_,nyp_,nzp_,useFileGrid);
    parPointCov_    = correlations_->getPriorVar0(); 
    parSpatialCorr  = correlations_->getPostCovAlpha();  // Double-use grids to save memory
    errCorrUnsmooth = correlations_->getPostCovBeta();   // Double-use grids to save memory
    // NBNB   nzp_*0.001*corr->getdt() = T    lowCut = lowIntCut*domega = lowIntCut/T
    int lowIntCut = int(floor(lowCut_*(nzp_*0.001*correlations_->getdt()))); 
    // computes the integer whis corresponds to the low cut frequency.
    float corrGradI, corrGradJ;
    model->getCorrGradIJ(corrGradI, corrGradJ);
    corrT = parSpatialCorr->fillInParamCorr(correlations_,lowIntCut,corrGradI, corrGradJ);
    if(spatwellfilter!=NULL)
    {
      for(int i=0;i<nWells_;i++)
        spatwellfilter->setPriorSpatialCorr(parSpatialCorr, wells_[i], i);
    }
    correlations_->setPriorCorrTFiltered(corrT,nz_,nzp_); // Can has zeros in the middle
    errCorrUnsmooth->fillInErrCorr(correlations_,corrGradI,corrGradJ); 
    if((model->getModelSettings()->getOtherOutputFlag() & ModelSettings::PRIORCORRELATIONS) > 0)
      correlations_->writeFilePriorCorrT(corrT,nzp_);     // No zeros in the middle
    // parSpatialCorr->writeFile("parSpatialCorr",simbox_);
    // parSpatialCorr->writeAsciiFile("SpatialCorr");         //for debug
    // errCorrUnsmooth->writeFile("errCorrUnsmooth",simbox_);
    // errCorrUnsmooth->writeAsciiFile("ErrCorr");            //for debug
  }
  else
  {
    model->releaseGrids();
  }

  // reality check: all dimensions involved match
  assert(meanBeta_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));
  assert(meanRho_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));

  for(int i=0 ; i< ntheta_ ; i++)
  {
    if(!model->getModelSettings()->getGenerateSeismic())
      assert(seisData_[i]->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));  
    assert(seisWavelet_[i]->consistentSize(nzp_, nyp_, nxp_));
  }

  if(!model->getModelSettings()->getGenerateSeismic())
  {
    parSpatialCorr->fftInPlace();
    computeVariances(corrT,model->getModelSettings());
    scaleWarning_ = checkScale();  // fills in scaleWarningText_ if needed.
    fftw_free(corrT);  
    if(simbox_->getIsConstantThick() == false)
      divideDataByScaleWavelet();
    errCorrUnsmooth->fftInPlace(); 
    for(int i = 0 ; i < ntheta_ ; i++)
    {
      // sprintf(rawSName,"rawSeismic_divided_%i",int(seisData_[i]->getTheta()*180.0/PI+0.5));
      seisData_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
      // seisData_[i]->writeFile(rawSName,simbox_);
      seisData_[i]->fftInPlace();
      seisData_[i]->endAccess();
    }
  }


  if ((outputFlag_ & ModelSettings::FACIESPROBRELATIVE) > 0 || model->getModelSettings()->noiseIsScaled()==true)
  {
    meanAlpha2_ = copyFFTGrid(meanAlpha_);
    meanBeta2_  = copyFFTGrid(meanBeta_);
    meanRho2_   = copyFFTGrid(meanRho_);
  }

  meanAlpha_->fftInPlace();
  meanBeta_ ->fftInPlace();
  meanRho_  ->fftInPlace();

  Timings::setTimeStochasticModel(wall,cpu);
}

Crava::~Crava()
{
  delete [] thetaDeg_;
  delete [] empSNRatio_;
  delete [] theoSNRatio_;
  delete [] modelVariance_;
  delete [] signalVariance_;
  delete [] errorVariance_;
  delete [] dataVariance_;
  if(fprob_!=NULL) delete fprob_;

  for(int i = 0;i<ntheta_;i++)
    delete[] errThetaCov_[i];
  delete [] errThetaCov_; 

  if(postAlpha_!=NULL) delete  postAlpha_ ;
  if(postBeta_!=NULL)  delete  postBeta_;
  if(postRho_!=NULL)   delete  postRho_ ;

  delete [] scaleWarningText_;
  if(sigmamdnew_!=NULL)
  {
    for(int i=0;i<nx_;i++)
    {
      for(int j=0;j<ny_;j++)
      {
        if((*sigmamdnew_)(i,j)!=NULL)
        {
          for(int ii=0;ii<3;ii++)
            delete [] (*sigmamdnew_)(i,j)[ii];
          delete [] (*sigmamdnew_)(i,j);
        }
      }
    }
     delete sigmamdnew_;

  }
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
Crava::setupErrorCorrelation(ModelSettings * modelSettings)
{
  //
  //  Setup error correlation matrix
  //
  for(int l=0 ; l < ntheta_ ; l++)
  {
    empSNRatio_[l]    = modelSettings->getSNRatio(l);
    if(modelSettings->noiseIsScaled()==true)
      errorVariance_[l] = float(dataVariance_[l]*modelSettings->getMinimumNoiseScaled(l)/empSNRatio_[l]);//NBNB Anne Randi: Sett denne lik minste ved lokal skalering.
    else
      errorVariance_[l] = dataVariance_[l]/empSNRatio_[l];

    if (empSNRatio_[l] < 1.1f) 
    {
      LogKit::LogFormatted(LogKit::LOW,"\nThe empirical signal-to-noise ratio for angle stack %d is %7.1e. Ratios smaller than\n",l+1,empSNRatio_[l]);
      LogKit::LogFormatted(LogKit::LOW," 1 are illegal and CRAVA has to stop. CRAVA was for some reason not able to estimate\n");
      LogKit::LogFormatted(LogKit::LOW," this ratio reliably, and you must give it as input to the model file\n\n");
      exit(1);
    }
  }

  Vario * angularCorr = modelSettings->getAngularCorr();

  for(int i = 0; i < ntheta_; i++)
    for(int j = 0; j < ntheta_; j++) 
      {
        float dTheta = modelSettings->getAngle(i) - modelSettings->getAngle(j);
        errThetaCov_[i][j] = static_cast<float>(sqrt(errorVariance_[i])
                                                *sqrt(errorVariance_[j])
                                                *angularCorr->corr(dTheta,0));
        //errThetaCov_[i][j] = static_cast<float>(sqrt(errorVariance_[i])*sqrt(errorVariance_[j])
        //                                        *angularCorr->corr(theta_[i]-theta_[j],0));
      }

}

void
Crava::computeVariances(fftw_real     * corrT,
                        ModelSettings * modelSettings)
{
  computeDataVariance();
    
  setupErrorCorrelation(modelSettings);

  char fileName[MAX_STRING];
  Wavelet ** errorSmooth = new Wavelet*[ntheta_];
  float    * paramVar    = new float[ntheta_] ;
  float    * WDCorrMVar  = new float[ntheta_];

  for(int i=0 ; i < ntheta_ ; i++)
  {
    if (seisWavelet_[i]->getDim() == 1) {
      errorSmooth[i] = new Wavelet1D(seisWavelet_[i],Wavelet::FIRSTORDERFORWARDDIFF);
      sprintf(fileName,"Wavediff_%d.dat",i);
      errorSmooth[i]->printToFile(fileName);
    }
    else {
      errorSmooth[i] = new Wavelet3D(seisWavelet_[i],Wavelet::FIRSTORDERBACKWARDDIFF);
      errorSmooth[i]->fft1DInPlace();
    }
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
    if (seisWavelet_[l]->getDim() == 1)
      WDCorrMVar[l] = computeWDCorrMVar(errorSmooth[l],corrT);
    else
      WDCorrMVar[l] = computeWDCorrMVar(errorSmooth[l]);
  }

  // Compute signal and model variance and theoretical signal-to-noise-ratio
  for(int l=0 ; l < ntheta_ ; l++)
  {
    modelVariance_[l]  = WDCorrMVar[l]*paramVar[l];
    signalVariance_[l] = errorVariance_[l] + modelVariance_[l];
  }

  for(int l=0 ; l < ntheta_ ; l++)
  {
    if (modelSettings->getMatchEnergies(l))
    {
      LogKit::LogFormatted(LogKit::LOW,"Matching syntethic and empirical energies:\n");
      float gain = sqrt((errorVariance_[l]/modelVariance_[l])*(empSNRatio_[l] - 1.0f));
      seisWavelet_[l]->scale(gain);
      if((modelSettings->getOtherOutputFlag() & ModelSettings::WAVELETS) > 0) 
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
                int(thetaDeg_[l]+0.5),scaleWarning1);
      }
      else
      {
        isOk = 1;
        sprintf(scaleWarningText_,"%sModel inconsistency in angle %i for seismic data\n%s \n",
                scaleWarningText_,int(thetaDeg_[l]+0.5),scaleWarning1);
      }
    }
    if( (dataVariance_[l] < 0.1 * signalVariance_[l]) && thisThetaIsOk) //1 var 0.1
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 2;
        sprintf(scaleWarningText_,"Model inconsistency in angle %i for seismic data\n%s\n",
                int(thetaDeg_[l]+0.5),scaleWarning2);
      }
      else
      {
        isOk = 2;
        sprintf(scaleWarningText_,"%sModel inconsistency in angle %i for seismic data\n%s\n",
                scaleWarningText_,int(thetaDeg_[l]+0.5),scaleWarning2);
      }
    }
    if( (modelVariance_[l] < 0.02 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
        {
        isOk = 3;
        sprintf(scaleWarningText_,"%s for angle %i.\n",
                scaleWarning3,int(thetaDeg_[l]+0.5));
      }
      else
        {
        isOk = 3;
        sprintf(scaleWarningText_,"%s%s for angle %i.\n",
                scaleWarningText_,scaleWarning3,int(thetaDeg_[l]+0.5));
      }
    }
    if( (modelVariance_[l] > 50.0 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
        {
        isOk = 4;
        sprintf(scaleWarningText_,"%s for angle %i.\n",
                scaleWarning4,int(thetaDeg_[l]+0.5) );
      }
      else
      {
        isOk = 4;
        sprintf(scaleWarningText_,"%s%s for angle %i.\n",
                scaleWarningText_,scaleWarning4,int(thetaDeg_[l]+0.5)); 
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

  rData  = static_cast<fftw_real*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real))); 
  cData  = reinterpret_cast<fftw_complex*>(rData);
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
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/static_cast<float>(sqrt(static_cast<float>(nzp_)));

          if(k > nz_)
          {
            float dist = seisData_[l]->getDistToBoundary( k, nz_, nzp_);  
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
          seisData_[l]->setRealValue(i,j,k,rData[k]/static_cast<float>(sqrt(static_cast<float>(nzp_))),true);
        }
      }
      char fName[200];
      if(ModelSettings::getDebugLevel() > 0)
      {
        sprintf(fName,"refl%d",l);
        std::string sgriLabel("Reflection coefficients for incidence angle ");
        sgriLabel += NRLib::ToString(thetaDeg_[l]+0.5);
        seisData_[l]->writeFile(fName, simbox_, sgriLabel);
      }
      LogKit::LogFormatted(LogKit::MEDIUM,"Interpolating reflections in volume %d: ",l);
      seisData_[l]->interpolateSeismic(energyTreshold_);
      if(ModelSettings::getDebugLevel() > 0)
      {
        sprintf(fName,"reflInterpolated%d",l);
        std::string sgriLabel("Interpolated reflections for incidence angle ");
        sgriLabel += NRLib::ToString(thetaDeg_[l]+0.5);
        seisData_[l]->writeFile(fName, simbox_, sgriLabel);
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

  rData  = static_cast<fftw_real*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real))); 
  cData  = reinterpret_cast<fftw_complex*>(rData);

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
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/static_cast<float>(sqrt(static_cast<float>(nzp_)));
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
      sprintf(fName,"%s_%d",typeName,int(thetaDeg_[l]+0.5));
      std::string sgriLabel(typeName);
      sgriLabel += " for incidence angle ";
      seisData_[l]->writeFile(fName, simbox_, sgriLabel);
      seisData_[l]->endAccess();
  }

  fftw_free(rData);
  fftwnd_destroy_plan(plan1); 
  fftwnd_destroy_plan(plan2); 
}


int 
Crava::computePostMeanResidAndFFTCov()
{
  Utils::writeHeader("Posterior model / Performing Inversion");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
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

  Wavelet1D * diff1Operator = new Wavelet1D(Wavelet::FIRSTORDERFORWARDDIFF,nz_,nzp_);
  Wavelet1D * diff2Operator = new Wavelet1D(diff1Operator,Wavelet::FIRSTORDERBACKWARDDIFF);
  Wavelet1D * diff3Operator = new Wavelet1D(diff2Operator,Wavelet::FIRSTORDERCENTRALDIFF);

  diff1Operator->fft1DInPlace();
  delete diff2Operator;
  diff3Operator->fft1DInPlace();

  Wavelet ** errorSmooth  = new Wavelet*[ntheta_];
  Wavelet ** errorSmooth2 = new Wavelet*[ntheta_];
  Wavelet ** errorSmooth3 = new Wavelet*[ntheta_];

  char* fileName = new char[2400] ;
  char* fileNameW = new char[2400] ;
  int cnxp  = nxp_/2+1;

  for(l = 0; l < ntheta_ ; l++)
  {
    seisData_[l]->setAccessMode(FFTGrid::READANDWRITE);
    if (seisWavelet_[0]->getDim() == 1) {
      errorSmooth[l]  = new Wavelet1D(seisWavelet_[l],Wavelet::FIRSTORDERFORWARDDIFF);
      errorSmooth2[l] = new Wavelet1D(errorSmooth[l], Wavelet::FIRSTORDERBACKWARDDIFF);
      errorSmooth3[l] = new Wavelet1D(errorSmooth2[l],Wavelet::FIRSTORDERCENTRALDIFF); 
      sprintf(fileName,"ErrorSmooth_%i",int(thetaDeg_[l]+0.5));
      errorSmooth3[l]->printToFile(fileName);
      errorSmooth3[l]->fft1DInPlace();

      sprintf(fileName,"Wavelet_%i",int(thetaDeg_[l]+0.5));
      seisWavelet_[l]->printToFile(fileName);
      seisWavelet_[l]->fft1DInPlace();
      sprintf(fileNameW,"FourierWavelet_%i",int(thetaDeg_[l]+0.5));
      seisWavelet_[l]->printToFile(fileNameW);
      delete errorSmooth[l];
      delete errorSmooth2[l];
    }
    else { //3D-wavelet
      errorSmooth3[l] = new Wavelet3D(seisWavelet_[l]);
      errorSmooth3[l]->fft1DInPlace();
      errorSmooth3[l]->multiplyByR(3.0);
      seisWavelet_[l]->fft1DInPlace();

      FFTGrid * errCorrUnsmooth = correlations_->getPostCovBeta(); // Stored in postCovBeta
      errCorrUnsmooth->setAccessMode(FFTGrid::RANDOMACCESS);
      int endX;
      if (nxp_/2 == (nxp_+1)/2) //nxp_ even
        endX = cnxp-1;
      else 
        endX = cnxp;
      float ijkNorm1 = 0.0;
      float ijkNorm2 = 0.0;
      float scale = static_cast<float>( 1.0/(nxp_*nyp_*nzp_));
      for (k=0; k<nzp_; k++) {
        float ijNorm1 = 0.0;
        float ijNorm2 = 0.0;
        for (j=0; j<nyp_; j++) {
          float iNorm1 = 0.0;
          float iNorm2 = 0.0;
          fftw_complex s1 = seisWavelet_[l]->getCAmp(k,j,0);
          fftw_complex s2 = errorSmooth3[l]->getCAmp(k,j,0);
          fftw_complex rho = errCorrUnsmooth->getComplexValue(0,j,k,true);
          float lambda = sqrt(rho.re * rho.re);
          iNorm1 += (s1.re * s1.re + s1.im * s1.im) * lambda;
          iNorm2 += (s2.re * s2.re + s2.im * s2.im) * lambda;
          for (i=1; i<endX; i++) {
            s1 = seisWavelet_[l]->getCAmp(k,j,i);
            s2 = errorSmooth3[l]->getCAmp(k,j,i);
            rho = errCorrUnsmooth->getComplexValue(i,j,k,true);
            lambda = sqrt(rho.re * rho.re);
            iNorm1 += 2.0f * (s1.re * s1.re + s1.im * s1.im) * lambda;
            iNorm2 += 2.0f * (s2.re * s2.re + s2.im * s2.im) * lambda;
          }
          if (endX == cnxp-1) { //nxp_ even, takes the last element once
            s1 = seisWavelet_[l]->getCAmp(k,j,cnxp-1);
            s2 = errorSmooth3[l]->getCAmp(k,j,cnxp-1);
            rho = errCorrUnsmooth->getComplexValue(cnxp-1,j,k,true);
            lambda = sqrt(rho.re * rho.re);
            iNorm1 += (s1.re * s1.re + s1.im * s1.im) * lambda;
            iNorm2 += (s2.re * s2.re + s2.im * s2.im) * lambda;
          }
          ijNorm1 += iNorm1;
          ijNorm2 += iNorm2;
        }
        ijkNorm1 += scale * ijNorm1;
        ijkNorm2 += ijNorm2;
      }
      ijkNorm2 /= (scale);
      seisWavelet_[l]->setNorm(sqrt(ijkNorm1));
      errorSmooth3[l]->setNorm(sqrt(ijkNorm2));
    }
  }
  delete[] errorSmooth;
  delete[] errorSmooth2;
  delete [] fileName;
  delete [] fileNameW;

  meanAlpha_->setAccessMode(FFTGrid::READANDWRITE);  //   Note
  meanBeta_ ->setAccessMode(FFTGrid::READANDWRITE);  //   the top five are over written
  meanRho_  ->setAccessMode(FFTGrid::READANDWRITE);  //   does not have the initial meaning.

  FFTGrid * parSpatialCorr     = correlations_->getPostCovAlpha(); // NB! Note double usage of postCovAlpha
  FFTGrid * errCorrUnsmooth    = correlations_->getPostCovBeta();  // NB! Note double usage of postCovBeta
  FFTGrid * postCovAlpha       = correlations_->getPostCovAlpha();
  FFTGrid * postCovBeta        = correlations_->getPostCovBeta();
  FFTGrid * postCovRho         = correlations_->getPostCovRho();
  FFTGrid * postCrCovAlphaBeta = correlations_->getPostCrCovAlphaBeta();
  FFTGrid * postCrCovAlphaRho  = correlations_->getPostCrCovAlphaRho();
  FFTGrid * postCrCovBetaRho   = correlations_->getPostCrCovBetaRho();
  parSpatialCorr    ->setAccessMode(FFTGrid::READANDWRITE);  //   after the prosessing
  errCorrUnsmooth   ->setAccessMode(FFTGrid::READANDWRITE);  // 
  postCovRho        ->setAccessMode(FFTGrid::WRITE);  
  postCrCovAlphaBeta->setAccessMode(FFTGrid::WRITE);  
  postCrCovAlphaRho ->setAccessMode(FFTGrid::WRITE);  
  postCrCovBetaRho  ->setAccessMode(FFTGrid::WRITE);  

  // Computes the posterior mean first  below the covariance is computed 
  // To avoid to many grids in mind at the same time
  double priorVarVp,justfactor;

  int cholFlag;
  //   long int timestart, timeend;
  //   time(&timestart);
  float realFrequency;

  LogKit::LogFormatted(LogKit::LOW,"\nBuilding posterior distribution:");
  float monitorSize = std::max(1.0f, static_cast<float>(nzp_)*0.02f);
  float nextMonitor = monitorSize;
  std::cout 
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^"; 

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

    for( j = 0; j < nyp_; j++) {
      for( i = 0; i < cnxp; i++) { 
        if (seisWavelet_[0]->getDim() == 3) { //3D-wavelet
          if( simbox_->getIsConstantThick() == true) {
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
            for(l = 0; l < ntheta_; l++) {
              errMult2[l].re   = static_cast<float>(errorSmooth3[l]->getCAmp(k,j,i).re / errorSmooth3[l]->getNorm());
              errMult2[l].im   = static_cast<float>(errorSmooth3[l]->getCAmp(k,j,i).im / errorSmooth3[l]->getNorm());
            }
            lib_matrFillOnesVecCpx(errMult3,ntheta_);
          }
          else {
            LogKit::LogFormatted(LogKit::LOW,"\nERROR: Not implemented inversion with 3D wavelet for non-constant simbox thickness\n");
            exit(1);
          }
        } //3D-wavelet
        ijkMean[0] = meanAlpha_->getNextComplex();
        ijkMean[1] = meanBeta_ ->getNextComplex();
        ijkMean[2] = meanRho_  ->getNextComplex(); 

        for(l = 0; l < ntheta_; l++ )
        {
          ijkData[l] = seisData_[l]->getNextComplex();
          ijkRes[l]  = ijkData[l];
        }  

        ijkTmp       = parSpatialCorr->getNextComplex();  
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
          ijkTmp       = errCorrUnsmooth->getNextComplex();  
          ijkErrLam.re = float( sqrt(ijkTmp.re * ijkTmp.re));
          ijkErrLam.im = 0.0;

          if(realFrequency > lowCut_*simbox_->getMinRelThick() &&  realFrequency < highCut_) // inverting only relevant frequencies
          {
            for(l = 0; l < ntheta_; l++ )
              for(m = 0; m < ntheta_; m++ )
              {        // Note we multiply kWNorm[l] and comp.conj(kWNorm[m]) hence the + and not a minus as in pure multiplication
                errVar[l][m].re  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult1[l].re *  errMult1[m].re +  errMult1[l].im *  errMult1[m].im)); 
                errVar[l][m].re += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult2[l].re *  errMult2[m].re +  errMult2[l].im *  errMult2[m].im)); 
                if(l==m) {
                  errVar[l][m].re += float( wnc_*errThetaCov_[l][m] * errMult3[l].re  * errMult3[l].re);
                  errVar[l][m].im   = 0.0;             
                }   
                else {
                  errVar[l][m].im  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult1[l].re * errMult1[m].im + errMult1[l].im * errMult1[m].re));
                  errVar[l][m].im += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult2[l].re * errMult2[m].im + errMult2[l].im * errMult2[m].re));
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
          postCovAlpha->setNextComplex(parVar[0][0]);
          postCovBeta ->setNextComplex(parVar[1][1]);
          postCovRho  ->setNextComplex(parVar[2][2]);
          postCrCovAlphaBeta->setNextComplex(parVar[0][1]);
          postCrCovAlphaRho ->setNextComplex(parVar[0][2]);
          postCrCovBetaRho  ->setNextComplex(parVar[1][2]); 

          for(l=0;l<ntheta_;l++)
            seisData_[l]->setNextComplex(ijkRes[l]);
      }
    }
    // Log progress
    if (static_cast<float>(k+1) >= nextMonitor) 
    { 
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  std::cout << "\n";

  //  time(&timeend);
  // LogKit::LogFormatted(LogKit::LOW,"\n Core inversion finished after %ld seconds ***\n",timeend-timestart);
  // these does not have the initial meaning
  meanAlpha_      = NULL; // the content is taken care of by  postAlpha_
  meanBeta_       = NULL; // the content is taken care of by  postBeta_
  meanRho_        = NULL; // the content is taken care of by  postRho_
  parSpatialCorr  = NULL; // the content is taken care of by  postCovAlpha      
  errCorrUnsmooth = NULL; // the content is taken care of by  postCovBeta

  postAlpha_->endAccess();
  postBeta_->endAccess(); 
  postRho_->endAccess();  

  postCovAlpha->endAccess();
  postCovBeta->endAccess(); 
  postCovRho->endAccess();
  postCrCovAlphaBeta->endAccess();
  postCrCovAlphaRho->endAccess();
  postCrCovBetaRho->endAccess();  

  postAlpha_->invFFTInPlace();
  postBeta_->invFFTInPlace();
  postRho_->invFFTInPlace();

  if(model_->getVelocityFromInversion() == true) { //Conversion undefined until prediction ready. Complete it.
    postAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);
    postAlpha_->expTransf();
    GridMapping * tdMap = model_->getTimeDepthMapping();
    const GridMapping * dcMap = model_->getTimeCutMapping();
    const Simbox * timeSimbox = simbox_;
    if(dcMap != NULL)
      timeSimbox = dcMap->getSimbox();

    tdMap->setMappingFromVelocity(postAlpha_, timeSimbox);
    postAlpha_->logTransf();
    postAlpha_->endAccess();
  }

  //NBNB Anne Randi: Skaler traser ihht notat fra Hugo
  if(model_->getModelSettings()->noiseIsScaled()==true)
  {
  
  correctAlphaBetaRho(model_->getModelSettings());
  
  }

  if(writePrediction_ == true)
  {
    doPostKriging(*postAlpha_, *postBeta_, *postRho_); 
    ParameterOutput::writeParameters(simbox_, model_, postAlpha_, postBeta_, postRho_, 
                                     outputFlag_, fileGrid_, -1);
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
        sprintf(fileNameS,"residuals_%i",int(thetaDeg_[l]+0.5));
        std::string sgriLabel("Residuals for incidence angle");
        sgriLabel += NRLib::ToString(int(thetaDeg_[l]+0.5));
        seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
        seisData_[l]->invFFTInPlace();
        seisData_[l]->writeFile(fileNameS,simbox_, sgriLabel);
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
  delete    diff3Operator;

  for(i = 0; i < ntheta_; i++)
  {
    delete[]  K[i];
    delete[]  KS[i];
    delete[]  margVar[i] ;
    delete[] errVar[i] ;
    delete errorSmooth3[i];   
  }
  delete[] K;
  delete[] KS;
  delete[] margVar;
  delete[] errVar  ;
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

  Timings::setTimeInversion(wall,cpu);
  return(0);
}


int
Crava::simulate(RandomGen * randomGen)
{   
  Utils::writeHeader("Simulating from posterior model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  if(nSim_>0)
  {
    FFTGrid * postCovAlpha       = correlations_->getPostCovAlpha();
    FFTGrid * postCovBeta        = correlations_->getPostCovBeta();
    FFTGrid * postCovRho         = correlations_->getPostCovRho();
    FFTGrid * postCrCovAlphaBeta = correlations_->getPostCrCovAlphaBeta();
    FFTGrid * postCrCovAlphaRho  = correlations_->getPostCrCovAlphaRho();
    FFTGrid * postCrCovBetaRho   = correlations_->getPostCrCovBetaRho();

    assert( postCovAlpha->getIsTransformed() );
    assert( postCovBeta->getIsTransformed() );
    assert( postCovRho->getIsTransformed() );
    assert( postCrCovAlphaBeta->getIsTransformed() );
    assert( postCrCovAlphaRho->getIsTransformed() );
    assert( postCrCovBetaRho->getIsTransformed() );

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

      postCovAlpha      ->setAccessMode(FFTGrid::READ);  
      postCovBeta       ->setAccessMode(FFTGrid::READ);  
      postCovRho        ->setAccessMode(FFTGrid::READ);  
      postCrCovAlphaBeta->setAccessMode(FFTGrid::READ);  
      postCrCovAlphaRho ->setAccessMode(FFTGrid::READ);  
      postCrCovBetaRho  ->setAccessMode(FFTGrid::READ);  
      seed0 ->setAccessMode(FFTGrid::READANDWRITE); 
      seed1 ->setAccessMode(FFTGrid::READANDWRITE); 
      seed2 ->setAccessMode(FFTGrid::READANDWRITE); 

      int cnxp=nxp_/2+1;
      int cholFlag;
      for(k = 0; k < nzp_; k++)
        for(j = 0; j < nyp_; j++)
          for(i = 0; i < cnxp; i++)
          {
            ijkPostCov[0][0] = postCovAlpha      ->getNextComplex();
            ijkPostCov[1][1] = postCovBeta       ->getNextComplex();
            ijkPostCov[2][2] = postCovRho        ->getNextComplex();
            ijkPostCov[0][1] = postCrCovAlphaBeta->getNextComplex();
            ijkPostCov[0][2] = postCrCovAlphaRho ->getNextComplex(); 
            ijkPostCov[1][2] = postCrCovBetaRho  ->getNextComplex();

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

          postCovAlpha->endAccess();  //         
          postCovBeta->endAccess();   // 
          postCovRho->endAccess();
          postCrCovAlphaBeta->endAccess();
          postCrCovAlphaRho->endAccess();
          postCrCovBetaRho->endAccess(); 
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
         
          seed1->setAccessMode(FFTGrid::RANDOMACCESS);
          seed1->invFFTInPlace(); 
         

          seed2->setAccessMode(FFTGrid::RANDOMACCESS);
          seed2->invFFTInPlace(); 

          if(model_->getModelSettings()->noiseIsScaled()==true)
          {
            float alpha,beta, rho;
            float alphanew, betanew, rhonew;

            for(j=0;j<ny_;j++)
              for(i=0;i<nx_;i++)
                for(k=0;k<nz_;k++)
                {
                  alpha = seed0->getRealValue(i,j,k);
                  beta = seed1->getRealValue(i,j,k);
                  rho = seed2->getRealValue(i,j,k);
                  alphanew = float((*sigmamdnew_)(i,j)[0][0]*alpha+ (*sigmamdnew_)(i,j)[0][1]*beta+(*sigmamdnew_)(i,j)[0][2]*rho);
                  betanew = float((*sigmamdnew_)(i,j)[1][0]*alpha+ (*sigmamdnew_)(i,j)[1][1]*beta+(*sigmamdnew_)(i,j)[1][2]*rho);
                  rhonew = float((*sigmamdnew_)(i,j)[2][0]*alpha+ (*sigmamdnew_)(i,j)[2][1]*beta+(*sigmamdnew_)(i,j)[2][2]*rho);
                  seed0->setRealValue(i,j,k,alphanew);
                  seed1->setRealValue(i,j,k,betanew);
                  seed2->setRealValue(i,j,k,rhonew);
                }
          }

          seed0->add(postAlpha_);
          seed0->endAccess();
          seed1->add(postBeta_);
          seed1->endAccess();
          seed2->add(postRho_);
          seed2->endAccess();



          doPostKriging(*seed0, *seed1, *seed2); 
          ParameterOutput::writeParameters(simbox_, model_, seed0, seed1, seed2,
                                           outputFlag_, fileGrid_, simNr);
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
  Timings::setTimeSimulation(wall,cpu);
  return(0);
}

void 
Crava::doPostKriging(FFTGrid & postAlpha, 
                     FFTGrid & postBeta, 
                     FFTGrid & postRho) 
{
  if(krigingParameter_ > 0) { 
    Utils::writeHeader("Conditioning to wells");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    CovGridSeparated covGridAlpha      (*correlations_->getPostCovAlpha()      );
    CovGridSeparated covGridBeta       (*correlations_->getPostCovBeta()       );
    CovGridSeparated covGridRho        (*correlations_->getPostCovRho()        ); 
    CovGridSeparated covGridCrAlphaBeta(*correlations_->getPostCrCovAlphaBeta());
    CovGridSeparated covGridCrAlphaRho (*correlations_->getPostCrCovAlphaRho() );
    CovGridSeparated covGridCrBetaRho  (*correlations_->getPostCrCovBetaRho()  );

    KrigingData3D kd(wells_, nWells_, 1); // 1 = full resolution logs
    kd.writeToFile("Raw");
    
    CKrigingAdmin pKriging(*simbox_, 
                           kd.getData(), kd.getNumberOfData(),
                           covGridAlpha, covGridBeta, covGridRho, 
                           covGridCrAlphaBeta, covGridCrAlphaRho, covGridCrBetaRho, 
                           krigingParameter_);

    pKriging.KrigAll(postAlpha, postBeta, postRho);
    Timings::setTimeKriging(wall,cpu);
  } 
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

  for(l=0;l<ntheta_;l++)
  { 
    seisData[l]->endAccess();
    seisData[l]->invFFTInPlace();

    sprintf(fileNameS,"synt_seis_%i",int(thetaDeg_[l]+0.5));
    std::string sgriLabel("Synthetic seismic for incidence angle ");
    sgriLabel += NRLib::ToString(int(thetaDeg_[l]+0.5));
    seisData[l]->writeFile(fileNameS,simbox_,sgriLabel);
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
    return var;
}

float  
Crava::computeWDCorrMVar (Wavelet* WD)
{
  float var = 0.0;
  int cnxp = nxp_/2 + 1;
  
  int endX;
  if (nxp_/2 == (nxp_+1)/2) //nxp_ even
    endX = cnxp-1;
  else 
    endX = cnxp;

  FFTGrid * parSpatialCorr = correlations_->getPostCovAlpha(); // Stored in postCovAlpha
  parSpatialCorr->setAccessMode(FFTGrid::READANDWRITE); 
  for(int k=0;k<nzp_;k++) {
    for(int j=0;j<nyp_;j++) {
      fftw_complex lambdaTmp = parSpatialCorr->getNextComplex();
      float ijkLambda = sqrt(lambdaTmp.re * lambdaTmp.re);
      fftw_complex wdTmp = WD->getCAmp(k,j,0);
      float ijkWD = (wdTmp.re * wdTmp.re) + (wdTmp.im * wdTmp.im);
      var += ijkLambda * ijkWD;
      for (int i=1;i<endX;i++) {
        lambdaTmp = parSpatialCorr->getNextComplex();
        ijkLambda = sqrt(lambdaTmp.re * lambdaTmp.re);
        wdTmp = WD->getCAmp(k,j,i);
        ijkWD = (wdTmp.re * wdTmp.re) + (wdTmp.im * wdTmp.im);
        var += 2.0f * ijkLambda * ijkWD;
      }
      if (endX == cnxp-1) { //nxp_ even, takes the last element once
        lambdaTmp = parSpatialCorr->getNextComplex();
        ijkLambda = sqrt(lambdaTmp.re * lambdaTmp.re);
        wdTmp = WD->getCAmp(k,j,cnxp-1);
        ijkWD = (wdTmp.re * wdTmp.re) + (wdTmp.im * wdTmp.im);
        var += ijkLambda * ijkWD;
      }
    }
  }
  var /= static_cast<float>(nxp_*nyp_*nzp_);
  parSpatialCorr->endAccess();
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
Crava::printEnergyToScreen()
{
  int i;
  LogKit::LogFormatted(LogKit::LOW,"\n\n                       ");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"  Seismic %4.1f ",thetaDeg_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\n----------------------");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"---------------");
  LogKit::LogFormatted(LogKit::LOW,"\nObserved data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",dataVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nModelled data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",signalVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nError variance         :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %1.3e  ",errorVariance_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nWavelet scale          :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %2.3e  ",seisWavelet_[i]->getScale());
  LogKit::LogFormatted(LogKit::LOW,"\nEmpirical S/N          :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %5.2f      ",empSNRatio_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\nModelled S/N           :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::LOW,"    %5.2f      ",theoSNRatio_[i]);
  LogKit::LogFormatted(LogKit::LOW,"\n");
}


void 
Crava::computeFaciesProb(SpatialWellFilter *filteredlogs)
{
  if((outputFlag_ & ModelSettings::FACIESPROB) >0 || (outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
  {
    Utils::writeHeader("Facies probability volumes");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    if (simbox_->getdz() > 4.0f) { // Require this density for estimation of facies probabilities
      LogKit::LogFormatted(LogKit::LOW,"\nWARNING: The minimum sampling density is lower than 4.0. The FACIES PROBABILITIES\n");
      LogKit::LogFormatted(LogKit::LOW,"         generated by CRAVA are not reliable. To get more reliable probabilities    \n");
      LogKit::LogFormatted(LogKit::LOW,"         the number of layers must be increased.                                    \n");
    }

    LogKit::LogFormatted(LogKit::LOW,"\n");
    LogKit::LogFormatted(LogKit::LOW,"Well                    Use    SyntheticVs    Deviated\n");
    LogKit::LogFormatted(LogKit::LOW,"------------------------------------------------------\n");
    for(int i=0 ; i<nWells_ ; i++) {
      LogKit::LogFormatted(LogKit::LOW,"%-23s %3s        %3s          %3s\n",
                           wells_[i]->getWellname(),
                           ( wells_[i]->getUseForFaciesProbabilities() ? "yes" : " no" ),
                           ( wells_[i]->hasSyntheticVsLog()            ? "yes" : " no" ),
                           ( wells_[i]->isDeviated()                   ? "yes" : " no" ));
    }

    int nfac = model_->getModelSettings()->getNumberOfFacies();

    std::string baseName;
    if((outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
    {
      meanAlpha2_->subtract(postAlpha_);
      meanAlpha2_->changeSign();
      meanBeta2_->subtract(postBeta_);
      meanBeta2_->changeSign();
      meanRho2_->subtract(postRho_);
      meanRho2_->changeSign();
      fprob_ = new FaciesProb(meanAlpha2_,
                              meanBeta2_,
                              meanRho2_,
                              nfac,
                              model_->getModelSettings()->getPundef(), 
                              model_->getPriorFacies(), 
                              model_->getPriorFaciesCubes(),
                              const_cast<const double **>(filteredlogs->getSigmae()),
                              const_cast<const WellData **>(wells_), 
                              nWells_,
                              true);
      delete meanAlpha2_;
      delete meanBeta2_;
      delete meanRho2_;

      baseName = "FaciesProbRelative_";
    }
    else
    {
      fprob_ = new FaciesProb(postAlpha_,
                              postBeta_,
                              postRho_,
                              nfac,
                              model_->getModelSettings()->getPundef(), 
                              model_->getPriorFacies(), 
                              model_->getPriorFaciesCubes(),
                              const_cast<const double **>(filteredlogs->getSigmae()),
                              const_cast<const WellData **>(wells_), 
                              nWells_,
                              false);
      baseName = "FaciesProb_";
    }
    fprob_->calculateConditionalFaciesProb(wells_, nWells_, model_->getModelSettings());
    LogKit::LogFormatted(LogKit::LOW,"\nProbability cubes done\n");
    for(int i=0;i<nfac;i++)
    {
      FFTGrid * grid = fprob_->getFaciesProb(i);
      std::string fileName = baseName + model_->getModelSettings()->getFaciesName(i);
      ParameterOutput::writeToFile(simbox_,model_,grid,fileName,"");
    }
    Timings::setTimeFaciesProb(wall,cpu);
  }
}

void 
Crava::filterLogs(Simbox          * timeSimboxConstThick,
                  FilterWellLogs *& filterlogs)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
  int relative;
  if((outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
    relative = 1;
  else 
    relative = 0;

  filterlogs = new FilterWellLogs(timeSimboxConstThick, 
                                  simbox_,
                                  correlations_,
                                  nzp_, nz_, 
                                  wells_, nWells_, 
                                  lowCut_, highCut_, 
                                  relative);
  Timings::setTimeFiltering(wall,cpu);
}
int Crava::getRelative()
{
int relative;
  if((outputFlag_ & ModelSettings::FACIESPROBRELATIVE)>0)
    relative = 1;
  else 
    relative = 0;

  return relative;
}

void Crava::computeG(double **G)
{
  correlations_->invFFT();
  correlations_->getPostVariances();
  correlations_->FFT();
 // correlations_->writeFilePostVariances();
 // correlations_->writeFilePriorVariances();
  double **sigmam = new double*[3];
  int i,j;
  for(i=0;i<3;i++)
  {
    sigmam[i] = new double[3];
    for(j=0;j<3;j++)
      sigmam[i][j] = double(correlations_->getPriorVar0()[i][j]);
  }
  float **sigmamd = correlations_->getPostVar0();
  double **sigmadelta = new double*[3];
  for(i=0;i<3;i++)
    sigmadelta[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmadelta[i][j] = double(sigmam[i][j]-sigmamd[i][j]);
 
  double  * eigval = new double[3];
  int     * error  = new int[1];
  double ** eigvec = new double *[3];
  double ** eigvalmat   = new double *[3];
  double ** help        = new double *[3];
  double ** eigvectrans = new double *[3];
  for(i=0;i<3;i++)
  {
    eigvec[i] = new double[3];
    eigvalmat[i]   = new double[3];
    help[i]        = new double[3];
    eigvectrans[i] = new double [3];
  }

  lib_matr_eigen(sigmadelta,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j] = eigval[i];
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmadelta);


  lib_matr_eigen(sigmam,3,eigvec,eigval,error);
  
  double ** sigmaminv   = new double *[3];
  double ** A           = new double *[3];
 
  for(i=0;i<3;i++)
  { 
    sigmaminv[i]   = new double[3];
    A[i]           = new double[3]; 
  }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0000001)
        eigvalmat[i][j] = 1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmaminv);
  lib_matr_prod(sigmaminv,sigmadelta,3,3,3,help);
  lib_matr_prod(help,sigmaminv,3,3,3,A);

  lib_matr_eigen(A,3,eigvec,eigval,error);
  lib_matr_sort3x3(eigval,eigvec);
  lib_matrTranspose(eigvec,3,3,eigvectrans);// V_A transponert

  double ** lambdag = new double*[ntheta_];
  for(i=0;i<ntheta_;i++)
    lambdag[i] = new double[3];
 
    for(i=0;i<ntheta_;i++)
      for(j=0;j<3;j++)
        if(i==j)
          lambdag[i][j] = sqrt(eigval[j]/(1.0-eigval[j]));
        else
          lambdag[i][j] = 0.0;
 

  double  * eigvale       = new double[ntheta_];
  double ** eigvece       = new double *[ntheta_];
  double ** eigvalmate    = new double*[ntheta_];
  double ** eigvecetrans  = new double *[ntheta_];
  double ** helpe         = new double *[ntheta_];
  

  for(i=0;i<ntheta_;i++)
  {
    eigvece[i] = new double[ntheta_];
    eigvalmate[i] = new double[ntheta_];
    eigvecetrans[i] = new double[ntheta_];
    helpe[i]        = new double[ntheta_];
  }
  lib_matr_eigen(errThetaCov_,ntheta_,eigvece,eigvale,error);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j && eigvale[i]>0.0)
        eigvalmate[i][j] = sqrt(eigvale[i]);
      else
        eigvalmate[i][j] = 0.0;
  
  lib_matrTranspose(eigvece,ntheta_,ntheta_,eigvecetrans);
  lib_matr_prod(eigvece,eigvalmate,ntheta_,ntheta_,ntheta_,helpe);
  lib_matr_prod(helpe,eigvecetrans,ntheta_,ntheta_,ntheta_,eigvece);
  double **help1= new double*[ntheta_];
  double **help2 = new double*[ntheta_];
  
  for(i=0;i<ntheta_;i++)
  {
    help1[i] = new double[3];
    help2[i] =new double[3];
  }
  lib_matr_prod(eigvece,lambdag,ntheta_,ntheta_,3,help1);
  lib_matr_prod(help1,eigvectrans,ntheta_,3,3,help2); //
  lib_matr_prod(help2,sigmaminv,ntheta_,3,3,G);
  
  for(i=0;i<3;i++)
  {
    delete [] eigvec[i];
    delete [] help[i];
    delete [] eigvectrans[i];
    delete [] eigvalmat[i];
    delete [] A[i];
    delete [] sigmadelta[i];
    delete [] sigmaminv[i];
  }
  delete [] eigval;
  delete [] eigvec;
  delete [] help;
  delete [] eigvectrans;
  delete [] eigvalmat;
  delete [] A;
  delete [] sigmadelta;
  delete [] sigmaminv;

  for(i=0;i<ntheta_;i++)
  {
    delete [] help1[i];
    delete [] help2[i];
    delete [] helpe[i];
    delete [] eigvecetrans[i];
    delete [] eigvece[i];
    delete [] eigvalmate[i];
  }

  delete [] eigvale;
  delete [] help1;
  delete [] help2;
  delete [] helpe;
  delete [] eigvecetrans;
  delete [] eigvece;
  delete [] eigvalmate;
  delete [] error;
  
}
void Crava::newPosteriorCovPointwise(double ** sigmanew, double **G, int igrid, int jgrid, ModelSettings *modelSettings, double **sigmamdnew)
{
  double **sigmaenew = new double*[ntheta_];
  double **D         = new double*[ntheta_];
  double **help      = new double*[ntheta_];
  int i,j;
  for(i=0;i<ntheta_;i++)
  {
    sigmaenew[i] = new double[ntheta_];
    D[i]         = new double[ntheta_];
    help[i]      = new double[ntheta_];
  }
  
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j)
        D[i][j] = sqrt(modelSettings->getNoiseScaled(i,igrid,jgrid)/modelSettings->getMinimumNoiseScaled(i));
      else
        D[i][j] = 0.0;

  lib_matr_prod(D,errThetaCov_,ntheta_,ntheta_,ntheta_,help);
  lib_matr_prod(help,D,ntheta_,ntheta_,ntheta_,sigmaenew);

  double **GT    = new double *[3];
  double **help1 = new double *[ntheta_];
  double **help2 = new double *[3];
  double **help4 = new double *[3];
  for(i=0;i<3;i++)
  {
    GT[i] = new double[ntheta_];
    help2[i] = new double[ntheta_];
    help4[i] = new double[ntheta_];
  }
  for(i=0;i<ntheta_;i++)
    help1[i] = new double[3];
  lib_matrTranspose(G,ntheta_,3,GT);
  double **sigmam = new double*[3];
  for(i=0;i<3;i++)
  {
    sigmam[i] = new double[3];
    for(j=0;j<3;j++)
      sigmam[i][j] = double(correlations_->getPriorVar0()[i][j]);
  }
  lib_matr_prod(G,sigmam,ntheta_,3,3,help1);
  lib_matr_prod(help1,GT,ntheta_,3,ntheta_,help);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      help[i][j]+=sigmaenew[i][j];


  int     * error         = new int[1];
  double  * eigvale       = new double[ntheta_];
  double ** eigvece       = new double *[ntheta_];
  double ** eigvalmate    = new double*[ntheta_];
  double ** eigvecetrans  = new double *[ntheta_];
  
  for(i=0;i<ntheta_;i++)
  {
    eigvece[i] = new double[ntheta_];
    eigvalmate[i] = new double[ntheta_];
    eigvecetrans[i] = new double[ntheta_];  
  }
  lib_matr_eigen(help,ntheta_,eigvece,eigvale,error);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j && eigvale[i]>0.00000001)
        eigvalmate[i][j] = 1.0/eigvale[i];
      else
        eigvalmate[i][j] = 0.0;

  lib_matrTranspose(eigvece,ntheta_,ntheta_,eigvecetrans);
  lib_matr_prod(eigvece,eigvalmate,ntheta_,ntheta_,ntheta_,help);
  lib_matr_prod(help,eigvecetrans,ntheta_,ntheta_,ntheta_,eigvece); // eigvece er inversmatrisen

  double ** help3 = new double *[3];
  for(i=0;i<3;i++)
    help3[i] = new double[3];
  double **deltanew = new double *[3];
  for(i=0;i<3;i++)
    deltanew[i] = new double[3];
  lib_matr_prod(sigmam,GT,3,3,ntheta_,help4);
  lib_matr_prod(help4,eigvece,3,ntheta_,ntheta_,help2);
  lib_matr_prod(help2,G,3,ntheta_,3,help3);
  lib_matr_prod(help3,sigmam,3,3,3,deltanew);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      sigmamdnew[i][j] = -deltanew[i][j];
      sigmamdnew[i][j]+=sigmam[i][j];
    }

  double  * eigval = new double[3];
  double     **eigvalmat  = new double*[3];
  double ** eigvec = new double *[3];
  double ** eigvectrans = new double *[3];
  for(i=0;i<3;i++)
  {
    eigvec[i] = new double[3];
    eigvalmat[i] = new double[3];
    eigvectrans[i] = new double[3];
  }

  lib_matr_eigen(sigmamdnew,3,eigvec,eigval,error); // take square root of sigmamdnew, because this is what is needed to save for later use in simulation.
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,sigmamdnew);

  lib_matr_eigen(deltanew,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,deltanew);

  float **sigmamd = correlations_->getPostVar0();
  double **sigmadelta = new double*[3];
  
  for(i=0;i<3;i++)
    sigmadelta[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmadelta[i][j] = double(sigmam[i][j]-sigmamd[i][j]);

  lib_matr_eigen(sigmadelta,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j&& eigval[i]>0.0000001)
        eigvalmat[i][j]=1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,sigmadelta);

  lib_matr_prod(deltanew,sigmadelta,3,3,3,sigmanew);

  for(i=0;i<ntheta_;i++)
  {
    delete [] D[i];
    delete [] sigmaenew[i];
    delete [] help[i];
    delete [] eigvece[i];
    delete [] eigvalmate[i];
    delete [] eigvecetrans[i];
    delete [] help1[i];
    
  }
  delete [] help;
  delete [] D;
  delete [] sigmaenew;
  delete [] eigvece;
  delete [] eigvalmate;
  delete [] eigvecetrans;
  delete [] help1;

  for(i=0;i<3;i++)
  {
    delete [] eigvalmat[i];
    delete [] eigvectrans[i];
    delete [] sigmadelta[i];
    delete [] sigmam[i];
    delete [] eigvec[i];
    delete [] help3[i];
    delete [] GT[i];
    delete [] help2[i];
    delete [] help4[i];
    
  }
  delete [] eigvalmat;
  delete [] eigvectrans;
  delete [] sigmadelta;
  delete [] sigmam;
  delete [] eigvec;
  delete [] help3;
  delete [] GT;
  delete [] help2;
  delete [] help4;
  delete [] error;


}

void Crava::correctAlphaBetaRho(ModelSettings *modelSettings)
{
  int i,j,k;
  double **G = new double*[ntheta_];
  for(i=0;i<ntheta_;i++)
    G[i] = new double[3];
  computeG(G);
 
  double **sigmanew = new double *[3];
  double **sigmamd = new double *[3];
  for(i=0;i<3;i++)
  {
    sigmanew[i] = new double[3];
    sigmamd[i] = new double[3];
  }
  double **sigmamdold = new double*[3];
  
  for(i=0;i<3;i++)
    sigmamdold[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmamdold[i][j] = correlations_->getPostVar0()[i][j];

  double  * eigval       = new double[3];
  double ** eigvalmat    = new double*[3];
  double ** eigvec       = new double *[3];
  double ** eigvectrans  = new double *[3];
  int     * error        = new int[1];
  double ** help         = new double*[3];
  for(i=0;i<3;i++)
  {
    eigvec[i] = new double[3];
    eigvalmat[i] = new double[3];
    eigvectrans[i] = new double[3];
    help[i] = new double[3];
  }

  lib_matr_eigen(sigmamdold,3,eigvec,eigval,error); 
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmamdold);

  


  float *alpha = NULL;
  float *beta = NULL;
  float *rho = NULL;
  float *meanalpha = NULL;
  float *meanbeta = NULL;
  float *meanrho = NULL;
  float alphadiff, betadiff, rhodiff;

  if(modelSettings->getNumberOfSimulations()>0)
    sigmamdnew_ = new NRLib::Grid2D<double **>(nx_,ny_,NULL);
  else
    sigmamdnew_ = NULL;
  for(i=0;i<nx_;i++)
    for(j=0;j<ny_;j++)
    {
      newPosteriorCovPointwise(sigmanew,G, i, j, modelSettings, sigmamd);
      lib_matr_prod(sigmamd,sigmamdold,3,3,3,eigvec); // store product in eigvec
      if(sigmamdnew_!=NULL)
      {
        (*sigmamdnew_)(i,j) = new double*[3];
        for(int ii=0;ii<3;ii++)
        {
          (*sigmamdnew_)(i,j)[ii] =new double[3];
          for(int jj=0;jj<3;jj++)
            (*sigmamdnew_)(i,j)[ii][jj] = eigvec[ii][jj];
        }

      }
      alpha = postAlpha_->getRealTrace(i,j);
      beta = postBeta_->getRealTrace(i,j);
      rho = postRho_->getRealTrace(i,j);
      meanalpha = meanAlpha2_->getRealTrace(i,j);
      meanbeta = meanBeta2_->getRealTrace(i,j);
      meanrho = meanRho2_->getRealTrace(i,j);
      for(k=0;k<nz_;k++)
      {
        alphadiff = alpha[k]-meanalpha[k];
        betadiff = beta[k]-meanbeta[k];
        rhodiff = rho[k]-meanrho[k];
        alpha[k] = float(meanalpha[k]+sigmanew[0][0]*alphadiff+sigmanew[0][1]*betadiff+sigmanew[0][2]*rhodiff);
        beta[k] = float(meanbeta[k]+sigmanew[1][0]*alphadiff+sigmanew[1][1]*betadiff+sigmanew[1][2]*rhodiff);
        rho[k] = float(meanrho[k]+sigmanew[2][0]*alphadiff+sigmanew[2][1]*betadiff+sigmanew[2][2]*rhodiff);
      }
      postAlpha_->setRealTrace(i,j, alpha);
      postBeta_->setRealTrace(i,j,alpha);
      postRho_->setRealTrace(i,j,alpha);

    }
    for(i=0;i<ntheta_;i++)
      delete [] G[i];

    for(i=0;i<3;i++)
    {
      delete [] eigvalmat[i];
      delete [] eigvectrans[i];
      delete [] eigvec[i];
      delete [] sigmamdold[i];
      delete [] help[i];
    }
    delete [] eigval;
    delete [] eigvalmat;
    delete [] eigvec;
    delete [] sigmamdold;
    delete [] error;
    delete [] help;

    delete [] G;
    delete [] alpha;
    delete [] beta;
    delete [] rho;
    delete [] meanalpha;
    delete [] meanbeta;
    delete [] meanrho;

}
