#include "src/fftgrid.h"
#include <math.h>
#include <algorithm>
#include "lib/lib_matr.h"
#include <src/posteriorelasticpdf.h>
#include <src/simbox.h>



void PosteriorElasticPDF::CalculateVariance2D(double                      ** sigma_smooth,
                                              double                      ** sigma_2d,
                                              const std::vector<double>      v1,
                                              const std::vector<double>      v2)
{
  sigma_2d[0][0] = v1[0]*v1[0]*sigma_smooth[0][0] + v1[1]*v1[1]*sigma_smooth[1][1] + v1[2]*v1[2]*sigma_smooth[2][2]
                  + 2*v1[0]*v1[1]*sigma_smooth[0][1]+ 2*v1[0]*v1[2]*sigma_smooth[0][2] + 2*v1[1]*v1[2]*sigma_smooth[1][2];
  sigma_2d[0][1] = v1[0]*v2[0]*sigma_smooth[0][0] + v1[0]*v2[1]*sigma_smooth[0][1] + v1[0]*v2[2]*sigma_smooth[0][2]
                +v1[1]*v2[0]*sigma_smooth[0][1] + v1[1]*v2[1]*sigma_smooth[1][1] + v1[1]*v2[2]*sigma_smooth[1][2]
                +v1[2]*v2[0]*sigma_smooth[0][2] + v1[2]*v2[1]*sigma_smooth[1][2] + v1[2]*v2[2]*sigma_smooth[2][2];
  sigma_2d[1][0] = sigma_2d[0][1];
  sigma_2d[1][1] = v2[0]*v2[0]*sigma_smooth[0][0] + v2[1]*v2[1]*sigma_smooth[1][1] + v2[2]*v2[2]*sigma_smooth[2][2]
                + 2*v2[0]*v2[1]*sigma_smooth[0][1]+ 2*v2[0]*v2[2]*sigma_smooth[0][2] + 2*v2[1]*v2[2]*sigma_smooth[1][2];

}

void PosteriorElasticPDF::CalculateTransform2D(const std::vector<double> d1,
                            const std::vector<double> d2,
                            const std::vector<double> d3,
                            std::vector<double> x,
                            std::vector<double> y,
                            const std::vector<double> v1,
                            const std::vector<double> v2)
{
  int dim = static_cast<int>(d1.size());
  for (int i=0; i<dim; i++){
    x[i] = v1[i]*d1[i] + v1[i]*d2[i] + v1[i]*d3[i];
    y[i] = v2[i]*d1[i] + v2[i]*d2[i] + v2[i]*d3[i];
  }
}

void PosteriorElasticPDF::InvertSquareMatrix(double   ** matrix,
                                             double   ** inv_matrix,
                                             int         n)
{
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(i==j)
        inv_matrix[i][j] = 1.0;
      else
        inv_matrix[i][j] = 0.0;

  lib_matrCholR(n, matrix);
  lib_matrAXeqBMatR(n, matrix, inv_matrix, n);
}


void PosteriorElasticPDF::SolveGEVProblem(double **sigma_prior,
                       double **sigma_posterior,
                       std::vector<double> & v1,
                       std::vector<double> & v2)
{
  //Compute transforms v1 and v2 ----------------------------------------

  // Matrix inversion of sigma_prior
  double **sigma_prior_inv = new double *[3];
  for(int i=0; i<3; i++)
    sigma_prior_inv[i] = new double [3];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      if(i==j)
        sigma_prior_inv[i][j] = 1.0;
      else
        sigma_prior_inv[i][j] = 0.0;

  lib_matrCholR(3, sigma_prior);
  lib_matrAXeqBMatR(3, sigma_prior, sigma_prior_inv, 3);

  // compute product sigma_prior_inv * sigma_posterior

  double **product_mat = new double *[3];

  for (int i=0; i<3; i++)
    product_mat[i] = new double[3];

  lib_matr_prod(sigma_prior_inv,sigma_posterior,3,3,3,product_mat);

  // compute eigenvalues of the product sigma_prior_inv * sigma_posterior

  int     * error       = new int[1];
  double  * eigval      = new double[3];
  double ** eigvec      = new double *[3];

  for(int i=0;i<3;i++)
  {
    eigvec[i]      = new double[3];
  }

  lib_matr_eigen(product_mat,3,eigvec,eigval,error);

  std::vector<int> index_keep;


  // Find index of max eigenvalue

  int max_index = 0;
  double max_eigvalue = 0.0;

  for(int i=0; i<3; i++){
    if (max_eigvalue<eigval[i]){
      max_index = i;
      max_eigvalue = eigval[i];
    }
    for(int j=0; j<3; j++){
      if(j!=max_eigvalue)
        index_keep.push_back(j);
    }
  }

  // The vector index_keep should contain two and only two integers
  assert (index_keep.size() == 2);


  // fetch the corresponding eigenvectors into v1 and v2
  for(int i = 0; i<3; i++){
    v1[i] = eigvec[index_keep[0]][i];
    v2[i] = eigvec[index_keep[1]][i];
  }


  delete [] error;
  delete [] eigval;
  for(int i=0;i<3;i++)
    delete [] eigvec[i];
  delete [] eigvec;
}
