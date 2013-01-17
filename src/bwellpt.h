/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef CBWELLPT_H
#define CBWELLPT_H

// a simple helping class: a blocked well is build up by a array of these blocked well points.

class CBWellPt
{
public:
  enum Gamma {ALPHA, BETA, RHO, ERROR};
  CBWellPt(int i = 0, int j = 0, int k = 0);
  ~CBWellPt(void);
  void  AddLog(float alpha, float beta, float rho); // Add log(data) if valid, increase count
  void  SubtractOnly(Gamma type, float gamma); // Subtract if data is valid, do NOT increase count
  void  Divide(); // divide by number of valid observations
  void  GetIJK(int& i, int& j, int& k) const {i = i_; j = j_; k = k_;}
  int   GetI() const { return i_; }
  int   GetJ() const { return j_; }
  int   GetK() const { return k_; }
  void  GetAlphaBetaRho(float& alpha, float& beta, float& rho) const {alpha = alpha_; beta = beta_; rho = rho_;}
  float GetAlpha() const { return alpha_ ;}
  float GetBeta()  const { return beta_  ;}
  float GetRho()   const { return rho_   ;}
  bool  IsIndexesEqual(int i, int j, int k) const
  {
    return (i_ == i && j_ == j && k_ == k);
  }
  void IsValidObs(bool& alpha, bool& beta, bool& rho) const
  {
    alpha = (noValidObsAlpha_ > 0);
    beta = (noValidObsBeta_ > 0);
    rho = (noValidObsRho_ > 0);
  }

private:
  int i_, j_, k_;
  float alpha_, beta_, rho_; // modified values
  float alphaOrig_, betaOrig_, rhoOrig_; // original values from well
  int noValidObsAlpha_, noValidObsBeta_, noValidObsRho_;
};

#endif
