#include "wavelet1D.h"

#include <string.h>
#include <assert.h>
#include <math.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/irapgrid.h"
#include "lib/log.h"
#include "lib/sgri.h"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet1D::Wavelet1D(Simbox         * simbox,
                     FFTGrid        * seisCube,
                     WellData      ** wells,
                     ModelSettings  * modelSettings,
                     float          * coeff,
                     int              dim)
  : Wavelet (dim)
{
  LogKit::writeLog("\n  Estimating 1D wavelet from seismic data and (nonfiltered) blocked wells\n");
  readtype_ = ESTIMATE;
  // initialization

  int   nWells        = modelSettings->getNumberOfWells();
  float waveletLength = modelSettings->getWaveletTaperingL();
  maxShift_           = modelSettings->getMaxWaveletShift();
  minRelativeAmp_     = modelSettings->getMinRelWaveletAmp();

  scale_=1.0f; 
  gridNI_=0;   
  gridNJ_=0;
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  errCode_=0;
  dz_        = static_cast<float>(simbox->getdz());
  nz_        = simbox->getnz();
  nzp_       = seisCube->getNzp();
  cnzp_      = nzp_/2+1;
  rnzp_	     = 2*cnzp_;
  theta_     = seisCube->getTheta();

  inFFTorder_= true;
  isReal_    = true; 
  cz_        = 0;

  float* dz = new float[nWells];
  float*  wellWeight=new float[nWells];
  int i, j, k;
  for(i=0;i<3;i++)
    coeff_[i] = coeff[i];

  int   nz   = simbox->getnz();
  float dz0  = static_cast<float>(simbox->getdz());
  int   nzp  = seisCube->getNzp();
  int   cnzp = (nzp/2+1);
  int   rnzp = 2*cnzp;

  float * alpha    = new float[nz];
  float * beta     = new float[nz];
  float * rho      = new float[nz];
  float * seisData = new float[nz];
  bool  * hasData  = new bool[nz];

  //Wavelet estimation
  fftw_real    ** cpp_r = new fftw_real*[nWells];
  fftw_complex ** cpp_c = (fftw_complex** ) cpp_r;
  
  fftw_real    ** seis_r = new fftw_real*[nWells];
  fftw_complex ** seis_c = (fftw_complex** ) seis_r;

  fftw_real    ** synt_seis_r = new fftw_real*[nWells];
  fftw_complex ** synt_seis_c = (fftw_complex** ) synt_seis_r;

  fftw_real    ** cor_cpp_r = new fftw_real*[nWells];
  fftw_complex ** cor_cpp_c = (fftw_complex** ) cor_cpp_r;
  
  fftw_real    ** ccor_seis_cpp_r = new fftw_real*[nWells];
  fftw_complex ** ccor_seis_cpp_c = (fftw_complex** ) ccor_seis_cpp_r;

  fftw_real    ** wavelet_r = new fftw_real*[nWells];
  fftw_complex ** wavelet_c = (fftw_complex** ) wavelet_r; 

  int maxBlocks = 0;
  for(i=0;i<nWells;i++)
  {
    cpp_r[i]           = new fftw_real[rnzp];
    seis_r[i]          = new fftw_real[rnzp];
    synt_seis_r[i]     = new fftw_real[rnzp];
    cor_cpp_r[i]       = new fftw_real[rnzp];
    ccor_seis_cpp_r[i] = new fftw_real[rnzp];
    wavelet_r[i]       = new fftw_real[rnzp];
    wellWeight[i]      = 0;
    const int * ipos = wells[i]->getBlockedLogsPropThick()->getIpos();
    const int * jpos = wells[i]->getBlockedLogsPropThick()->getJpos();
    dz[i] = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])*simbox->getdz());
    int nBlocks = wells[i]->getBlockedLogsPropThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) 
  {
    if (!wells[w]->isDeviated())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bl->getVerticalTrend(bl->getAlpha(), alpha);
      bl->getVerticalTrend(bl->getBeta(), beta);
      bl->getVerticalTrend(bl->getRho(), rho);
      bl->getVerticalTrend(seisLog, seisData);

      for (k = 0 ; k < nz ; k++) {
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      }

      int start,length;
      findContiniousPartOfData(hasData,nz,start,length);

      if(length*dz0 > waveletLength  ) // must have enough data
      {
        fillInCpp(alpha,beta,rho,start,length,cpp_r[w],nzp);
        fft(cpp_r[w],cpp_c[w],nzp);
        fillInSeismic(seisData,start,length,seis_r[w],nzp);
        fft(seis_r[w],seis_c[w],nzp);
        estimateCor(cpp_c[w],cpp_c[w],cor_cpp_c[w],cnzp);
        fftInv(cor_cpp_c[w],cor_cpp_r[w],nzp);
        estimateCor(seis_c[w],cpp_c[w],ccor_seis_cpp_c[w],cnzp);
        fftInv(ccor_seis_cpp_c[w],ccor_seis_cpp_r[w],nzp);
        fftInv(cpp_c[w],cpp_r[w],nzp);
        fftInv(seis_c[w],seis_r[w],nzp);
        wellWeight[w] = length*dz[w]*(cor_cpp_r[w][0]+cor_cpp_r[w][1]);// Gives most weight to long datasets with  
                                                                       // large reflection coefficients
      }
      else
        wellWeight[w] = 0; 
    }
  }
  delete [] seisLog;
  float* shiftWell = new float[nWells];
  float shiftAvg = shiftOptimal(ccor_seis_cpp_r,wellWeight,dz,nWells,nzp,shiftWell);

  multiplyPapolouis(ccor_seis_cpp_r,dz,nWells,nzp, waveletLength);
  multiplyPapolouis(cor_cpp_r,dz,nWells,nzp, waveletLength);
  getWavelet(ccor_seis_cpp_r,cor_cpp_r,wavelet_r,wellWeight, nWells, nzp);
  rAmp_      = averageWavelets(wavelet_r,nWells,nzp,wellWeight,dz,dz0); // wavelet centered
  cAmp_      = (fftw_complex*) rAmp_;
  char* fileName = new char[MAX_STRING];
  //sprintf(fileName,"wavelet");
  //printToFile(fileName,rAmp_, nzp);

  for(int w=0;w<nWells;w++) // gets syntetic seismic with estimated wavelet
  {
    fillInnWavelet(wavelet_r[w],nzp,dz[w]);
    shiftReal(shiftWell[w]/dz[w],wavelet_r[w],nzp);
    //sprintf(fileName,"waveletShift");
    //printToFile(fileName,wavelet_r[w], nzp);
    fft(wavelet_r[w],wavelet_c[w],nzp);
    //sprintf(fileName,"cpp");
    //printToFile(fileName,cpp_r[w], nzp);
    fft(cpp_r[w],cpp_c[w],nzp);
    convolve(wavelet_c[w],cpp_c[w],synt_seis_c[w],cnzp);
    fftInv(synt_seis_c[w],synt_seis_r[w],nzp); // 
    //sprintf(fileName,"syntSeis");
    //printToFile(fileName,synt_seis_r[w], nzp);
    //sprintf(fileName,"seis");
    //printToFile(fileName,seis_r[w], nzp);
    //printf("Test run\n");
  }

  float* scaleOptWell    = new float[nWells];
  float* errWellOptScale = new float[nWells];
  float* errWell = new float[nWells];
  float err;
  float scaleOpt = findOptimalWaveletScale(synt_seis_r,seis_r,nWells,nzp,wellWeight,err,errWell,scaleOptWell,errWellOptScale);

  delete [] scaleOptWell ;
  delete [] errWellOptScale;
  delete [] errWell ;
  
//  shiftReal(shiftAvg/dz0,rAmp_,nzp);
  
  shiftAndScale(shiftAvg,scaleOpt);//shifts wavelet average from wells
  invFFT1DInPlace();
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( LogKit::getDebugLevel() > 0 )
  {
    //flipUpDown();
    sprintf(fileName,"estimated_wavelet");
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
    //flipUpDown();
  }
  double norm2=0.0;
  for(i=0; i < nzp_; i++ )
      norm2 += rAmp_[i]*rAmp_[i];
  norm_= float( sqrt(norm2));
  //char* fileName=new char[MAX_STRING];


  if(LogKit::getDebugLevel() > 1) {
    sprintf(fileName,"estimated_wavelet_full_%d.txt",int(ceil((theta_*180/PI)-0.5)) );
    FILE* fid = fopen(fileName,"w");
    for(i=0;i<nzp_;i++)
      fprintf(fid,"%1.3e\n",rAmp_[i]);
    fclose(fid);
    
    for(j=0;j<nWells;j++)
    {
      sprintf(fileName,"seis_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",seis_r[0][i]);
      fclose(fid);
        
      sprintf(fileName,"cor_cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",cor_cpp_r[0][i]);
      fclose(fid);
      
      sprintf(fileName,"ccor_seis_cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",ccor_seis_cpp_r[0][i]);
      fclose(fid);

      sprintf(fileName,"cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",cpp_r[0][i]);
      fclose(fid);
    }
  }
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  for(i=0;i<nWells;i++)
  {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    delete [] cor_cpp_r[i] ;
    delete [] ccor_seis_cpp_r[i] ;
    delete [] wavelet_r[i];
  }
  flipUpDown(); //NB ODD temporary fix - FRODE
}

Wavelet1D::Wavelet1D(char * fileName, ModelSettings * modelSettings, int fileFormat, int dim)
  : Wavelet(modelSettings, dim)
{
  switch (fileFormat)
	{
	case OLD: WaveletReadOld(fileName);
		break;
	case JASON: WaveletReadJason(fileName);
		break;
  }
  
	for(int i=0; i < rnzp_ ;i++)
	{  
		if(i < nzp_)
		{
			rAmp_[i]*=scale_;
		}
		else
		{
			rAmp_[i]=RMISSING;
		}// endif
	}//end for i
}

Wavelet1D::Wavelet1D(Wavelet * wavelet, int difftype, int dim)
  : Wavelet(wavelet, difftype, dim)
{
}

Wavelet1D::Wavelet1D(int difftype, int nz, int nzp, int dim)
  : Wavelet(difftype, nz, nzp, dim)
{
}

Wavelet1D::Wavelet1D(Wavelet * wavelet, int dim)
  : Wavelet(wavelet, dim)
{
}

void
Wavelet1D::WaveletReadJason(char * fileName)
{
  readtype_=JASON;
  FILE* file = fopen(fileName,"r");
  bool lineIsComment = true; 
  char* dummyStr= new char[MAX_STRING];
  while( lineIsComment ==true)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      readToEOL(file);
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
      return;
    } // endif
    else
    {
      readToEOL(file);
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
      {
        lineIsComment = false;
      }
    }
  }  // end while
  float shift=float(atof(dummyStr));
  if(fscanf(file,"%s",dummyStr) == EOF)
  {
    readToEOL(file);
    if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
    else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
    errCode_=2; 
    return;
  } // endif
  dz_=float(atof(dummyStr));
  if(fscanf(file,"%s",dummyStr) == EOF)
  {
    readToEOL(file);
    if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
    else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
    errCode_=2; 
    return;
  } // endif
  nz_=atoi(dummyStr);

  cz_ =  int(floor((fabs(shift/dz_))+0.5));
  nzp_         = nz_;
  cnzp_        = nzp_/2+1;
  rnzp_        = 2*cnzp_; 
  rAmp_        = (fftw_real* ) fftw_malloc(sizeof(float)*rnzp_);
  cAmp_        = (fftw_complex* ) rAmp_;
  norm_        = RMISSING;

  for(int i=0; i<nz_;i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      readToEOL(file);
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
      return;
    } // endif
    else
    {
      rAmp_[i] = (fftw_real) atof(dummyStr);
    }
  }
  fclose(file);
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("\n  Estimated wavelet length:  %.1fms.",waveletLength_);
  delete [] dummyStr;
}

void
Wavelet1D::WaveletReadOld(char * fileName)
{
  readtype_=OLD;
  FILE* file = fopen(fileName,"r");
 
  int  maxWaveletL=10000;

  
  int   i,pos,shift,nSamples;
  char  headStr[MAX_STRING];
  char  tmpStr[MAX_STRING];
  char  targetString[MAX_STRING];
  char  number[MAX_STRING];
  float dz,ampMult;

  for(i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",headStr) == EOF)
    {
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
    } // endif
  }  // end for i


  strcpy(targetString,"CMX");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    if (errCode_ == 0)  sprintf(errText_,"Error when reading wavelet amplitude from file  %s.\n", fileName);
    else sprintf(errText_,"%sError when reading wavelet amplitude from file  %s. \n", errText_,fileName);
    errCode_=3; 
    ampMult = RMISSING; // Dummy setting to avoid g++ warning 
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    ampMult= float ( atof(number) );
  } //endif


  strcpy(targetString,"SI");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    if (errCode_ == 0)  sprintf(errText_,"Error when reading sampling interval from file  %s.\n", fileName);
    else sprintf(errText_,"%sError when reading sampling interval from file  %s.\n",errText_,fileName);
    errCode_=4; 
    dz = RMISSING; // Dummy setting to avoid g++ warning
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    dz= float ( atof(number) );
  } //endif

  float * tempWave= (float*) fftw_malloc(sizeof(float)* maxWaveletL );

  nSamples=0;

  while(fscanf(file,"%s",tmpStr) != EOF)
  {      
    if( maxWaveletL  > nSamples )
    {
      char * target = new char[2];
      strcpy(target,"F");
      pos = findEnd(tmpStr, 0, target);
      delete [] target;
      if(pos == -1)
      {
        tempWave[nSamples]=float( atof(tmpStr) );
      }
      else     
      {
        nSamples--;
      }
      nSamples++;
    }
    else
    {
      if (errCode_ == 0)  sprintf(errText_,"Error in memory use when reading wavelet from file  %s.\n", fileName);
      else  sprintf(errText_,"%sError in memory use when reading wavelet from file  %s.\n",errText_,fileName);
      errCode_=5;
    }//endif
  }//endwhile

  fclose(file);

  strcpy(targetString,"SHIFT");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    shift=nSamples/2; // integer division
    if(shift*2 == nSamples)
    {
      if (errCode_ == 0)  
        sprintf(errText_,"Error when reading wavelet shift from file  %s.\n    --> No SHIFT and even number of data.\n", fileName);
      else  
        sprintf(errText_,"%sError when reading wavelet shift from file  %s.\n    --> No SHIFT and even number of data.\n",errText_,fileName); 
      errCode_=6; 
      //
      // NBNB-PAL: Temporary? hack since the errorCode is never returned from constructor
      //
      LogKit::writeLog("\nError when reading wavelet shift from file  %s.",fileName);
      LogKit::writeLog("\n    --> No SHIFT and even number of data.\n");
      exit(1);
    }
    //cz_ = shift;   // case no flip
    shift =-shift-1; // case flip
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    shift= int ( atof(number) );

    //cz_=nSamples+shift;  // case no flip
  }//endif

  if(errCode_ == 0) 
  {
    cz_          = -shift-1; // case flip
    theta_       = RMISSING;
    dz_          = dz;
    nz_          = nSamples;
    nzp_         = nSamples;
    cnzp_        = nzp_/2+1;
    rnzp_        = 2*cnzp_; 
    inFFTorder_  = false;
    isReal_      = true;
    rAmp_        = (fftw_real* ) fftw_malloc(sizeof(float)*rnzp_);
    cAmp_        = (fftw_complex* ) rAmp_;
    norm_        = RMISSING;

    // Note the wavelet is fliped left to right  
    // since it is used for a convolution 

    for(i=0; i < rnzp_ ;i++)
    {  
      if(i < nzp_)
      {
        rAmp_[i]=ampMult*tempWave[nzp_-i-1];
      }
      else
      {
        rAmp_[i]=RMISSING;
      }// endif
    }//end for i 
  }

  fftw_free(tempWave);
  // LogKit::writeLog("\nReading wavelet file %s  ... done.\n",fileName);

  //
  // Estimate wavelet length
  //
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("\n  Estimated wavelet length:  %.1fms.\n",waveletLength_);
}

void
Wavelet1D::resample(float dz, int nz, float pz, float theta) 
{
  theta_=theta;
  if(errCode_==0){
    //LogKit::writeLog("  Resampling wavelet\n");
    assert(isReal_);
    assert(!inFFTorder_);
    int nzp,cnzp,rnzp,k;
    float z;
    fftw_real* wlet;

    nzp   =  findClosestFactorableNumber( (int) ceil( nz*(1.0f+pz) ) );
    cnzp  =  nzp/2 + 1;
    rnzp  =  2*cnzp;

    wlet  = (fftw_real *) fftw_malloc( sizeof(fftw_real)*rnzp );

    for(k=0; k < rnzp; k++)
    {
      if(k < nzp)
      {
        if(k < nzp/2+1)
        {
          z = float( dz*k );
        }
        else
        {
          z= float( dz*(k-nzp) );
        }
        wlet[k] = getWaveletValue(z, rAmp_ , cz_, nz_, dz_);
      }
      else
      {
        wlet[k] =RMISSING;
      }
    }
    fftw_free( rAmp_);

    float norm2 = 0.0; 
    for(k=0; k < nzp; k++) norm2 +=wlet[k]*wlet[k];

    rAmp_       = (fftw_real *)    wlet; // rAmp_ is not allocated 
    cAmp_       = (fftw_complex* ) rAmp_;
    nzp_        = nzp;
    rnzp_       = rnzp;
    cnzp_       = cnzp;
    cz_         = 0;
    nz_         = nz;
    dz_         = dz;
    norm_       = float( sqrt( norm2) );
    inFFTorder_ = true;
  }
  if(readtype_ == OLD) //FRODE
     flipUpDown();
  if( LogKit::getDebugLevel() > 0 )
  {
    //flipUpDown();// ODD temporary debugfix
    char* fileName = new char[MAX_STRING];
    sprintf(fileName,"resampled_wavelet");
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
    //flipUpDown();// ODD temporary debugfix
    delete [] fileName;
  }
}

bool           
Wavelet1D::consistentSize(int nzp, int, int) const
{ 
  if (nzp!=nzp_) 
    printf("nzp=%d  nzp_wavelet1D=%d\n",nzp,nzp_); 
  return (nzp==nzp_);
}

void           
Wavelet1D::flipUpDown()
{
  if(isReal_==true)
  {
    float tmp;
    for(int i=1;i<nzp_/2;i++)
    {
      tmp=rAmp_[i];
      rAmp_[i] = rAmp_[nzp_-i];
      rAmp_[nzp_-i] =tmp;
    }
  }
  else
  {
    for(int i=0;i<cnzp_;i++)
    {
      cAmp_[i].im *=-1.0;
    }
  }
}

float 
Wavelet1D::getWaveletValue(float z, float *Wavelet, int center, int nz, float dz)
{
  // returns the value of the vavelet in the location z. Wavelet have the length nz 
  // and the center value is Wavelet[center]
  // uses kriging with ricker 20Hz wavelet as correlation function.
  float value;
  // double a,b,c;
  //int    im1,i0,ip1,
  int    k,l;
  int*   ind=new int[6];// iL1,iL2,iL3,iR1,iR2,iR3;
  double* val=new double[6];//vL1,vL2,vL3,vR1,vR2,vR3;

  ind[2]= int( floor( (z/dz) ) );
  for(k=0;k<6;k++)
    ind[k]=  ind[2]+k-2;

  //i0    =  int( ceil( (z/dz) - 0.5) );
  //ip1   =  i0 + 1;
  //im1   =  i0 - 1;
  //hz    =  z - i0*dz;

  for(k=0;k<6;k++)
    val[k]=  getArrayValueOrZero(ind[k]+center , Wavelet,  nz); 

  double** Cov = new double*[6];
  double*  cov = new double[6];
  double   nu = 20;
  double   deltaT;

  for(k=0;k<6;k++)
  { 
    Cov[k] = new double[6];
    deltaT = (dz*ind[k]-z)*0.001;
    cov[k] = (1-2*nu*nu*PI*PI*(deltaT)*(deltaT))*exp(-nu*nu*PI*PI*(deltaT)*(deltaT));
    for(l=0;l<6;l++)
    {
      deltaT =(dz*ind[k]-dz*ind[l])*0.001;
      Cov[k][l] = (1-2*nu*nu*PI*PI*deltaT * deltaT )*exp(-nu*nu*PI*PI*deltaT*deltaT);
    }
  }
  //OK not very intellegent implementation since chol is done for each time step.
  lib_matrCholR(6,  Cov);
  lib_matrAxeqbR(6, Cov, cov); // cov contains now the kriging weigths; 

  value = 0.0;
  for(k=0;k<6;k++)
  {
    value+= float(val[k]*cov[k]);
    delete [] Cov[k];
  }
  delete [] Cov;
  delete [] cov;


  //v0    =  getArrayValueOrZero(i0+center  , Wavelet,  nz);    
  //vp1   =  getArrayValueOrZero(ip1+center , Wavelet,  nz);  
  //vm1   =  getArrayValueOrZero(im1+center , Wavelet,  nz);  

  //c     =  v0;
  //b     =  (vp1-vm1)/(2.0*dz);
  //a     =  (vp1+vm1-2.0*v0)/(2.0*dz*dz);
  //value2 =  float(a*hz*hz + b*hz + c);

  delete [] val;
  delete [] ind;
  return value;
}

int Wavelet1D::getWaveletLengthI()
{
  bool trans=false;
  if(isReal_==false)
  {
    invFFT1DInPlace();
    trans=true;
  }
  float maxAmp =  fabs(rAmp_[0]); // gets max amp 
  for(int i=1;i <nzp_;i++)
    if(fabs(rAmp_[i]) > maxAmp)
      maxAmp = fabs(rAmp_[i]);

  float minAmp= maxAmp*minRelativeAmp_; // minimum relevant amplitude

  int wLength=nzp_;

  for(int i=nzp_/2;i>0;i--)
  {
    if(fabs(rAmp_[i]) >minAmp)
    {
      wLength= (i*2+1);// adds both sides 
      break;
    }
    if(fabs(rAmp_[nzp_-i]) > minAmp)
    {
      wLength= (2*i+1);// adds both sides 
      break;
    }
  }
  wLength =MINIM(wLength,2*(nzp_/2)-1); // allways even number
  if(trans==true)
    fft1DInPlace();

  return wLength;
}
float
Wavelet1D::getWaveletLengthF()
{
  return dz_*float( getWaveletLengthI() );
}

float
Wavelet1D::shiftOptimal(fftw_real** ccor_seis_cpp_r,float* wellWeight,float* dz,int nWells,int nzp,float* shiftWell)
{
  float shift=0.0f;
  float totalWeight=0;
  float sum=0;
  int w,i,polarity;
  // if the sum from -maxShift_ to maxShift_ ms is 
  // positive then polarity is positive
  for(w=0;w<nWells;w++)
  {
    shiftWell[w]=0;
    if(wellWeight[w] > 0) // only include those estimated
    {
      for(i=0;i<ceil(maxShift_/dz[w]);i++)//zero included
        sum+=ccor_seis_cpp_r[w][i];
      for(i=0;i<floor(maxShift_/dz[w]);i++)
        sum+=ccor_seis_cpp_r[w][nzp-i-1];
    }
  }

  polarity=-1;
  if(sum > 0)
    polarity=1;

  // gets optimal shift
  float maxValue;
  float shiftF;
  int shiftI;
  float f1,f2,f3;

  for(w=0;w<nWells;w++)
  {
    if(wellWeight[w]>0)
    {
      maxValue = 0.0f;
      shiftI=0;

      for(i=0;i<ceil(maxShift_/dz[w]);i++)
      {
        if(ccor_seis_cpp_r[w][i]*polarity > maxValue)
        {
          maxValue = ccor_seis_cpp_r[w][i]*polarity;
          shiftI = i;
        }
      }
      for(i=0;i<floor(maxShift_/dz[w]);i++)
      {
        if(ccor_seis_cpp_r[w][nzp-1-i]*polarity > maxValue)
        {
          maxValue = ccor_seis_cpp_r[w][nzp-1-i]*polarity;
          shiftI = -1-i;
        }
      }
      if(shiftI < 0)
      {
        if(ccor_seis_cpp_r[w][nzp+shiftI-1]*polarity < maxValue) //then local max
        {
          f1 = ccor_seis_cpp_r[w][nzp+shiftI-1];
          f2 = ccor_seis_cpp_r[w][nzp+shiftI];
          int ind3;
          if(shiftI==-1)
            ind3 = 0;
          else
            ind3=nzp+shiftI+1;
          f3 = ccor_seis_cpp_r[w][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      else //positive or zero shift
      {
        if(ccor_seis_cpp_r[w][shiftI+1]*polarity < maxValue) //then local max
        {
          f3 = ccor_seis_cpp_r[w][shiftI+1];
          f2 = ccor_seis_cpp_r[w][shiftI];
          int ind1;
          if(shiftI==0)
            ind1 = nzp-1;
          else
            ind1=shiftI-1;
          f1 = ccor_seis_cpp_r[w][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      shiftWell[w] = shiftF*dz[w];
      shiftReal(-shiftF, ccor_seis_cpp_r[w],nzp);// 
      shift += wellWeight[w]*shiftF*dz[w];//weigthing shift according to wellWeight
      totalWeight += wellWeight[w];
    }
  }
  shift/=totalWeight;
  return shift;
}

void
Wavelet1D::multiplyPapolouis(fftw_real** vec, float* dz,int nWells,int nzp, float waveletLength) const
{
  int i,w;
  float wHL=float( waveletLength/2.0);
  float weight,dist;
  for(w=0;w<nWells;w++)
    for(i=1;i<nzp;i++)
    {
      dist =MINIM(i,nzp-i)*dz[w];
      if(dist < wHL)
      {
        weight  = float(1.0/PI*fabs(sin(PI*(dist)/wHL))); 
        weight += float((1-fabs(dist)/wHL)*cos(PI*dist/wHL));
      }
      else
        weight=0;

      vec[w][i]*=weight;
    }
}

void
Wavelet1D::getWavelet(fftw_real** ccor_seis_cpp_r,fftw_real** cor_cpp_r,fftw_real** wavelet_r,float* wellWeight,int nWells,int nt)
{
  fftw_complex* c_sc,*c_cc,*wav;

  int cnzp = nt/2+1;
  for(int w=0;w<nWells;w++)
  {
    if(wellWeight[w] > 0)
    {
      c_sc   = (fftw_complex*)ccor_seis_cpp_r[w];
      fft(ccor_seis_cpp_r[w],c_sc,nt);
      c_cc   = (fftw_complex*)cor_cpp_r[w];
      fft(cor_cpp_r[w],c_cc,nt);
      wav    = (fftw_complex*) wavelet_r[w];

      for(int i=0;i<cnzp;i++)
      {
        wav[i].re=c_sc[i].re/c_cc[i].re;//note c_cc[i].im =0
        wav[i].im=c_sc[i].im/c_cc[i].re;
      }
      fftInv(c_sc,ccor_seis_cpp_r[w],nt);
      fftInv(c_cc,cor_cpp_r[w],nt);
      fftInv(wav,wavelet_r[w],nt);
    }
  }
}

fftw_real* 
Wavelet1D::averageWavelets(fftw_real** wavelet_r,int nWells,int nzp,float* wellWeight,float* dz,float dzOut) const
{
  // assumes dz[w] < dzOut for all w
  fftw_real* wave= (fftw_real*) fftw_malloc(rnzp_*sizeof(fftw_real));  

  int w,i;

  float* weight= new float[nWells];// weight is length of data interval in well
  float sum = 0;
  for(i=0;i<nWells;i++)
    sum+=wellWeight[i];

  for(i=0;i<nWells;i++)
    weight[i]=wellWeight[i]/sum;

  
  for(i=0;i<nzp;i++)
    wave[i] = 0.0; // initialize
  float ww;

  for( w=0;w<nWells;w++)
  {
    if(wellWeight[w]>0)
    {
      wave[0] += weight[w]* wavelet_r[w][0]; // central

      for(i=1;i <nzp/2-1;i++) // positive time
      {
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];  // fraction of distance to ind*dz[w]
        if(ind<nzp/2)
            ww= (1-t)* wavelet_r[w][ind]+t*(wavelet_r[w][ind+1]);
        else
          ww=0;
        wave[i]    += weight[w]*ww;

      }
      for(i=1;i < nzp/2-1;i++)// negative time Note dz[w] <= dzOut for all w
      {
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];
        if(ind<nzp/2)
          ww    = (1-t)* wavelet_r[w][nzp-ind]+t*(wavelet_r[w][nzp-(ind+1)]);
        else
          ww=0;

        wave[nzp-i] += weight[w]*ww;
      }
    }
  }
  
  char* fileName = new char[MAX_STRING];
  sprintf(fileName,"wavelet_%d_fftOrder_noshift",int(floor(theta_/PI*180+0.5)));
  printToFile(fileName,wave,nzp_);
  
  delete [] weight;
  return wave;
}

float 
Wavelet1D::getArrayValueOrZero(int i  ,float * Wavelet, int nz) const
{
  float value;

  if(i > -1 && i < nz)
    value = Wavelet[i];
  else
    value = 0.0; 
  return value;
}

