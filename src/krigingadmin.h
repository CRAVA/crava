
#ifndef CKRIGINGADMIN_H
#define CKRIGINGADMIN_H

class FFTGrid;
class CBWellPt;
class Model;
class Simbox;
class CovGridSeparated;
#include "src/box.h"

class CKrigingAdmin
{
private: enum DataBoxSize { DBS_TOO_SMALL, DBS_TOO_BIG, DBS_RIGHT};
public:
  CKrigingAdmin(const Simbox& simbox, CBWellPt ** pBWellPt, int noData,              
                CovGridSeparated& covAlpha, 
                CovGridSeparated& covBeta, 
                CovGridSeparated& covRho,   
                CovGridSeparated& covCrAlphaBeta, 
                CovGridSeparated& covCrAlphaRho, 
                CovGridSeparated& covCrBetaRho,
                int  dataTarget = 200,
                bool backgroundModel = false);
  ~CKrigingAdmin(void);
  enum Gamma {ALPHA_KRIG, BETA_KRIG, RHO_KRIG};
  void KrigAll(FFTGrid& trendAlpha, FFTGrid& trendBeta, FFTGrid& trendRho, 
               bool trendsAlreadySubtracted = false);

private:
  void            Init();
  void            KrigAll(Gamma gamma);
  void            KrigBlock(Gamma gamma);
  /* Finds the data by using the following rule: Cokriging 3 variables X,Y,Z.  
  If you are doing kriging on X. Then for each well obs: if you have info on X use it and 
  ignore the two others Y,Z. Else use info on Y and Z. 
  */
  void            SubtractTrends(FFTGrid& trend_alpha, FFTGrid& trend_beta, FFTGrid& trend_rho);
  void            FindDataInDataBlockLoop(Gamma gamma);
  DataBoxSize     FindDataInDataBlock(Gamma gamma, const CBox & dataBlock);
  int             NBlocks(int dBlocks, int lSBox) const;
  void            AllocateSpaceForMatrixEq();
  void            DeAllocateSpaceForMatrixEq();
  void            SetMatrix(Gamma gamma);
  void            SetKrigVector(Gamma gamma);
  void            SolveKrigingEq(Gamma gamma);
  void            DoCholesky();
  void            EstimateSizeOfBlock();
  void            EstimateSizeOfBlock2();
  float           CalcCPUTime(float dxBlock, float dyBlockExt, float& nd, bool& rapidInc);
  const FFTGrid & CreateAndFillFFTGridWithCov(int i);    // needed in kriging hack modus, remember delete
  const FFTGrid & CreateAndFillFFTGridWithCovRot(int i); // needed in kriging hack modus, remember delete
  const FFTGrid & CreateAndFillFFTGridWithCrCov();       // needed in kriging hack modus, remember delete
  double       ** CopyMatrix(double** inMatrix, double** outMatrix, int size);
  void            RotateVec(float& rx, float& ry, float& rz, const float mat[][3]);
  void            SmoothKrigedResult(Gamma gamma);
  void            CalcSmoothWeights(Gamma gamma, int direction); // direction X(1), Y(2), Z(3)
  int             GetSmoothBlockNx() { return 2 * (dxSmoothBlock_ + 1); }
  int             GetSmoothBlockNy() { return 2 * (dySmoothBlock_ + 1); }
  int             GetSmoothBlockNz() { return 2 * (dzSmoothBlock_ + 1); }
  void            WriteDebugOutput() const;
  void            WriteDebugOutput2() const;
  void            Require(bool test, const std::string msg = "Requirement failed") const;
  FFTGrid       * CreateValidGrid() const;

  const Simbox  & simbox_;
  FFTGrid       * trendAlpha_, *trendBeta_, *trendRho_; // pointers to current kriging variables 
  CovGridSeparated &covAlpha_, &covBeta_, &covRho_, &covCrAlphaBeta_, &covCrAlphaRho_, &covCrBetaRho_;
  FFTGrid       * pBWellGrid_; // a "bool" grid that says "true" (1.0f), (or NOT -1.0f) if there is at least one blocked valid well data in the cell 
  CBWellPt     ** pBWellPt_;
  CBox            currDataBox_, currBlock_;                  // current data neightbourhood and kriging area
  int             dxBlock_, dyBlock_, dzBlock_;              // number of cells to define a kriging block 
  int             dxBlockExt_, dyBlockExt_, dzBlockExt_;     // number of additional cells to reach data neighbourhood
  int             i_, j_, k_;                                // current kriging indexes
  int           * pIndexAlpha_, *pIndexBeta_, *pIndexRho_;  // holds an array of indexes into pBWells_
  int             sizeAlpha_, sizeBeta_, sizeRho_;           // current sizes
  int             noValidAlpha_, noValidBeta_, noValidRho_;  // number of valid a, b og r data
  int             noValid_;                                  // total number of valid data
  int             noData_;                                   // number kriging data (blocks)
  int             dataTarget_;                               // target values
  int             noKrigedCells_;                            // number of actually kriged cells
  int             noEmptyDataBlocks_;                        // number of empty datablocks
  int             noKrigedVariables_;                        // number of kriged variables so far
  int             noCholeskyDecomp_;                         // number of cholesky decompositions
  int             monitorSize_;                              // for progress monitor

  double       ** krigMatrix_, **krigMatrix2_, *krigVector_; 
  float         * krigDataVector_;
  int             rangeAlphaX_, rangeAlphaY_, rangeAlphaZ_;
  int             rangeBetaX_, rangeBetaY_, rangeBetaZ_;
  int             rangeRhoX_, rangeRhoY_, rangeRhoZ_;        // ranges estimated from covariance cubes
  float           rangeX_, rangeY_, rangeZ_;                 // max ranges
  enum            { maxDataTolerance_         = 10,          // in % of dataTarget_ 
                    maxDataBlockLoopCounter_  =  7,          // ca. max number of attempts to find a right data block neighbourhood, not used
                    maxCholeskyLoopCounter_   = 20,          // max number of attempts to cholesky decomposition
                    switchFailed_             =  1};         // assert flag

  int             totalNoDataInCurrKrigBlock_;               // total number of data in current kriging block
  int             noSolvedMatrixEq_;                         // total number of times we have actually solved the matrix eq, for debug
  int              noRMissing_;                               // total number of times we have missing real values
  bool            failed2EstimateRange_, failed2EstimateDefaultDataBoxAndBlock_;             // bool flags if we failed 2 estimate true
  bool            backgroundModel_;
  int             dxSmoothBlock_, dySmoothBlock_, dzSmoothBlock_;                            // normal value is 2, data size is 2*n + 2
  double       ** ppKrigSmoothWeightsX_, **ppKrigSmoothWeightsY_, **ppKrigSmoothWeightsZ_; // first index is kriged point, second is data 
};

#endif
