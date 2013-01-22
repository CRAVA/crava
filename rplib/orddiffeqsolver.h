#ifndef RPLIB_ORDDIFFEQSOLVER_H
#define RPLIB_ORDDIFFEQSOLVER_H

#include <cstring>
#include <vector>

class OrdDiffEqSolver {
 public:
   OrdDiffEqSolver();
   ~OrdDiffEqSolver();

 //ODE45 integrates a system of ordinary differential equations using
 //4th and 5th order Runge-Kutta formulas.
 static void Ode45(std::vector<double>                 (func_ptr)(std::vector<double>&, double),
                   double                               t0,
                   double                               tfinal,
                   std::vector<double>&                 y0,
                   std::vector<double>&                 tout,
                   std::vector< std::vector<double> >&  yout,
                   double                               tol = 1.e-6);

private:
static void CalcVector(const std::vector< std::vector<double> >& matrix,
                       const std::vector< std::vector<double> >& f,
                       double                                    h,
                       size_t                                    row,
                       std::vector<double>&                      y);

};
#endif
