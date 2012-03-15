#include "rplib/orddiffeqsolver.h"

#include "nrlib/exception/exception.hpp"

#include <cmath>


OrdDiffEqSolver::OrdDiffEqSolver() {

}

OrdDiffEqSolver::~OrdDiffEqSolver()
{

}

void
OrdDiffEqSolver::
Ode45(std::vector<double>                 (*func_ptr)(std::vector<double>&, double),
      double                               t0,
      double                               tfinal,
      std::vector<double>&                 y0,
      std::vector<double>&                 tout,
      std::vector< std::vector<double> >&  yout,
      double                               tol) {


  // constant matrices initialization
  static std::vector<double> alpha(5);
  alpha[0] = 1.0/4.0;
  alpha[1] = 3.0/8.0;
  alpha[2] = 12.0/13.0;
  alpha[3] = 1.0;
  alpha[4] = 1.0/2.0;

  static std::vector< std::vector<double> > beta(5);

  beta[0].resize(6, 0.0);
  beta[0][0] = 1.0/4.0;

  beta[1].resize(6, 0.0);
  beta[1][0] = 3.0/32.0;
  beta[1][1] = 9.0/32.0;

  beta[2].resize(6, 0.0);
  beta[2][0] = 1932.0/2197.0;
  beta[2][1] = -7200.0/2197.0;
  beta[2][2] = 7296.0/2197.0;

  beta[3].resize(6, 0.0);
  beta[3][0] = 8341.0/4104.0;
  beta[3][1] = -32832.0/4104.0;
  beta[3][2] = 29440.0/4104.0;
  beta[3][3] = -845.0/4104.0;

  beta[4].resize(6, 0.0);
  beta[4][0] = -6080.0/20520.0;
  beta[4][1] = 41040.0/20520.0;
  beta[4][2] = -28352.0/20520.0;
  beta[4][3] = 9295.0/20520.0;
  beta[4][4] = -5643.0/20520.0;

  static std::vector< std::vector<double> > gamma(2);

  gamma[0].resize(6, 0.0);
  gamma[0][0] = 902880.0/7618050.0;

  gamma[0][2] = 3953664.0/7618050.0;
  gamma[0][3] = 3855735.0/7618050.0;
  gamma[0][4] = -1371249.0/7618050.0;
  gamma[0][5] = 277020.0/7618050.0;

  gamma[1].resize(6, 0.0);

  gamma[1][0] = -2090.0/752400.0;

  gamma[1][2] = 22528.0/752400.0;
  gamma[1][3] = 21970.0/752400.0;
  gamma[1][4] = -15048.0/752400.0;
  gamma[1][5] = -27360.0/752400.0;

  // Other initialization
  double t = t0;
  double hmax = (tfinal - t)/16.0;
  double h = hmax/8.0;
  double power = 1.0/5.0;

  std::vector<double> y = y0;
  std::vector< std::vector<double> > f; // 6 x 2 matrix
  f.resize(6);
  for (size_t i = 0; i < f.size(); i++)
    f[i].resize(2, 0.0);

  unsigned int chunk = 128;
  tout.reserve(chunk);
  yout.reserve(chunk);
  for (size_t i = 0; i < yout.size(); i++)
    yout[i].reserve(y.size());

  tout.push_back(t0);
  yout.push_back(y0);

  while (t < tfinal && (t + h) > t) {
    if (t+h > tfinal)
      h = tfinal - t; //NBNB fjellvoll is this correct in c++

    //Compute the slopes
    std::vector<double> temp = (*func_ptr)(y, t);
    f[0] = temp;

    for (unsigned int j = 0; j < 5; j++) {
      double t1 = t + alpha[j]*h;
      std::vector<double> y1(y);
      CalcVector(beta, f, h, j, y1);
      temp = (*func_ptr)(y1, t1);
      f[j+1] = temp;
    } // end loop j

    //estimate error and acceptable error
    std::vector<double> d(y.size(), 0.0);
    CalcVector(gamma, f, h, 1, d);

    double delta = std::abs(d[0]);
    if (std::abs(d[1]) > delta)
      delta = std::abs(d[1]);

    double tau = std::abs(y[0]);
    if (std::abs(y[1]) > tau)
      tau = std::abs(y[1]);

    if (1.0 > tau)
      tau = 1.0;

    tau *= tol;

    if (delta <= tau) {
      t += h;
      CalcVector(gamma, f, h, 0, y);
      tout.push_back(t);
      yout.push_back(y);
    }

    if (delta != 0.0) {
      h = 0.8*h*pow(tau/delta, power);
      if (hmax < h)
        h = hmax;
    }
  } // end while

  if (t < tfinal)
    throw NRLib::Exception("DEM: Singularity likely.");

}


void
OrdDiffEqSolver::
CalcVector(const std::vector< std::vector<double> >& matrix,
           const std::vector< std::vector<double> >& f,
           double                                    h,
           size_t                                    row,
           std::vector<double>&                      y) {

  for (unsigned int i1 = 0; i1 < y.size(); i1++)
    for (unsigned int j1 = 0; j1 < 6; j1++)
      y[i1] += h*matrix[row][j1]*f[j1][i1];
}
