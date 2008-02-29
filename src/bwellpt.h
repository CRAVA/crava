#ifndef CBWELLPT_H
#define CBWELLPT_H

// a simple helping class: a blocked well is build up by a array of these blocked well points.

class CBWellPt
{
public:
  enum Gamma {ALPHA, BETA, RHO, ERROR};
  CBWellPt(int i = 0, int j = 0, int k = 0);
  ~CBWellPt(void);
  void AddLog(float alpha, float beta, float rho); // Add log(data) if valid, increase count
  void AddLog(Gamma gamma, float val);
  void SubtractOnly(float alpha, float beta, float rho); // Subtract if data is valid, do NOT increase count 
  void SubtractOnly(Gamma type, float gamma); // Subtract if data is valid, do NOT increase count 
  bool IsIndexesEqual(int i, int j, int k) const {
    return (i_ == i && j_ == j && k_ == k);
  }
  void Divide(); // divide by number of valid observations
  void GetIJK(int& i, int& j, int& k) const {i = i_; j = j_; k = k_;}
  int GetI() const { return i_; }
  int GetJ() const { return j_; }
  int GetK() const { return k_; }
  void GetAlphaBetaRho(float& alpha, float& beta, float& rho) const {alpha = alpha_; beta = beta_; rho = rho_;}
  float GetAlpha() const { return alpha_; }
  float GetBeta() const { return beta_; }
  float GetRho() const { return rho_;}
  /*
  void GetNumberOfValidObs(int& noAlpha, int& noBeta, int& noRho) const { 
  noAlpha = noValidObsAlpha_; noBeta = noValidObsBeta_; noRho = noValidObsRho_;}
  */
  void IsValidObs(bool& alpha, bool& beta, bool& rho) const {
    alpha = (noValidObsAlpha_ > 0); beta = (noValidObsBeta_ > 0); rho = (noValidObsRho_ > 0); }
  bool IsValidObsAlpha() const { return (noValidObsAlpha_ > 0); }
  bool IsValidObsBeta() const { return (noValidObsBeta_ > 0); }
  bool IsValidObsRho() const { return (noValidObsRho_ > 0); }
  void SetFaciesParam(int nfac, int *fnr);
  void FindLargestFaciesCount(double unif01);
  void AddFaciesLog(int val);
  bool FaciesDefined() const;
  int GetFacies() const {return facies_;}
  bool FaciesNotMissing() const;

private:
  int i_, j_, k_;
  float alpha_, beta_, rho_; // modified values
  float alphaOrig_, betaOrig_, rhoOrig_; // original values from well
  int noValidObsAlpha_, noValidObsBeta_, noValidObsRho_;
  int *faciesnr_;
  int *faciescount_;
  int facies_;
  //int faciesOrig_;
  int nFacies_;
  
};

#endif
