#include "src/fftgrid.h"
#include <math.h>
#include <src/posteriorelasticpdf.h>
#include <src/simbox.h>



void PosteriorElasticPDF::CalculateVariance2D(NRLib::Matrix & sigma_smooth,//double                      ** sigma_smooth,
                                              NRLib::Matrix & sigma_2d,//double                      ** sigma_2d,
                                              const std::vector<double>      & v1_std,
                                              const std::vector<double>      & v2_std)
{
  NRLib::Vector v1(3);
  NRLib::Vector v2(3);
  for (int i = 0; i < 3; i++){
    v1(i) = v1_std[i];
    v2(i) = v2_std[i];
  }
  sigma_2d(0,0) = v1(0)*v1(0)*sigma_smooth(0,0) + v1(1)*v1(1)*sigma_smooth(1,1) + v1(2)*v1(2)*sigma_smooth(2,2)
                  + 2*v1(0)*v1(1)*sigma_smooth(0,1)+ 2*v1(0)*v1(2)*sigma_smooth(0,2) + 2*v1(1)*v1(2)*sigma_smooth(1,2);
  sigma_2d(0,1) = v1(0)*v2(0)*sigma_smooth(0,0) + v1(0)*v2(1)*sigma_smooth(0,1) + v1(0)*v2(2)*sigma_smooth(0,2)
                +v1(1)*v2(0)*sigma_smooth(0,1) + v1(1)*v2(1)*sigma_smooth(1,1) + v1(1)*v2(2)*sigma_smooth(1,2)
                +v1(2)*v2(0)*sigma_smooth(0,2) + v1(2)*v2(1)*sigma_smooth(1,2) + v1(2)*v2(2)*sigma_smooth(2,2);
  sigma_2d(1,0) = sigma_2d(0,1);
  sigma_2d(1,1) = v2(0)*v2(0)*sigma_smooth(0,0) + v2(1)*v2(1)*sigma_smooth(1,1) + v2(2)*v2(2)*sigma_smooth(2,2)
                + 2*v2(0)*v2(1)*sigma_smooth(0,1)+ 2*v2(0)*v2(2)*sigma_smooth(0,2) + 2*v2(1)*v2(2)*sigma_smooth(1,2);

}

void PosteriorElasticPDF::CalculateTransform2D(const std::vector<double>                & d1,
                                               const std::vector<double>                & d2,
                                               const std::vector<double>                & d3,
                                               std::vector<std::vector<double> >        & x,
                                               const NRLib::Matrix                      & v)
{
  int dim = static_cast<int>(d1.size());
  int sizeOfV = static_cast<int>(v.numRows());
  for (int i=0; i<sizeOfV; i++){
    for (int j=0; j<dim; j++){
      x[i][j] = v(i,0)*d1[j] + v(i,1)*d2[j] + v(i,2)*d3[j];
    }
  }
}

void PosteriorElasticPDF::InvertSquareMatrix(NRLib::Matrix  & matrix,//double   ** matrix,
                                             NRLib::Matrix  & inv_matrix,//double   ** inv_matrix,
                                             int              n)
{
  assert (matrix.numCols() == n && matrix.numRows() == n);
  //assert (matrix.numCols() == inv_matrix.numCols() && matrix.numRows() == inv_matrix.numRows());

  NRLib::SymmetricMatrix sym_mat(n);
  for (int i= 0; i < n; i++)
    for (int j= i; j < n; j++)
      sym_mat(i,j) = matrix(i,j);

  inv_matrix = NRLib::IdentityMatrix(n);

  /*
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(i==j)
        inv_matrix(i,j) = 1.0;
      else
        inv_matrix(i,j) = 0.0;
  */
  NRLib::CholeskySolve(sym_mat, inv_matrix);
  //lib_matrCholR(n, matrix);
  //lib_matrAXeqBMatR(n, matrix, inv_matrix, n);
}


void PosteriorElasticPDF::SolveGEVProblem(NRLib::Matrix           & sigma_prior,  //double               ** sigma_prior,
                                          NRLib::Matrix           & sigma_post, //double               ** sigma_posterior,
                                          std::vector<double>     & v1,         //std::vector<double>     & v1,
                                          std::vector<double>     & v2)
{
  //Compute transforms v1 and v2 ----------------------------------------

  // Matrix inversion of sigma_prior
  //double **sigma_prior_inv = new double *[3];
  //for(int i=0; i<3; i++)
  //  sigma_prior_inv[i] = new double [3];
  NRLib::Matrix sigma_pri_inv = NRLib::IdentityMatrix(3);

  NRLib::SymmetricMatrix sym_mat;
  for(int i=0; i<3; i++)
    for(int j=i; j<3; j++)
      sym_mat(i,j) = sigma_prior(i,j);



  NRLib::CholeskySolve(sym_mat, sigma_pri_inv);

  //lib_matrCholR(3, sigma_prior);
  //lib_matrAXeqBMatR(3, sigma_prior, sigma_prior_inv, 3);

  // compute product sigma_prior_inv * sigma_posterior

  //double **product_mat = new double *[3];

  //for (int i=0; i<3; i++)
  //  product_mat[i] = new double[3];

  //lib_matr_prod(sigma_prior_inv,sigma_posterior,3,3,3,product_mat);

  NRLib::Matrix product = sigma_pri_inv*sigma_post;

  // compute eigenvalues of the product sigma_prior_inv * sigma_posterior

  //int     * error       = new int[1];
  //double  * eigval      = new double[3];
  NRLib::Vector eigval(3, 0);
  NRLib::Matrix eigvec(3,3,0);
  //double ** eigvec      = new double *[3];

  //for(int i=0;i<3;i++)
  //{
  //  eigvec[i]      = new double[3];
  //}

  NRLib::ComputeEigenVectors(product, eigval, eigvec);
  //lib_matr_eigen(product_mat,3,eigvec,eigval,error);

  std::vector<int> index_keep;


  // Find index of max eigenvalue

  int max_index = 0;
  double max_eigvalue = 0.0;

  for(int i=0; i<3; i++){
    if (max_eigvalue<eigval(i)){
      max_index = i;
      max_eigvalue = eigval(i);
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
    v1[i] = eigvec(index_keep[0],i);
    v2[i] = eigvec(index_keep[1],i);
  }


  //delete [] error;
  //delete [] eigval;
  //for(int i=0;i<3;i++)
  //  delete [] eigvec[i];
  //delete [] eigvec;
}

void PosteriorElasticPDF::SetupSmoothingGaussian2D(FFTGrid                 * smoother,
                                                   const NRLib::Matrix     & sigma_inv,
                                                   int                       n1,
                                                   int                       n2,
                                                   int                       n3,
                                                   double                    dx,
                                                   double                    dy)
{
  float *smooth = new float[n1*n2*n3];
  int j,k,l,jj,jjj,kk,kkk,ll,lll;
  lll=2;

  float sum = 0.0f;
  for(l=0; l<n3; l++) {
    kkk=2;
    if(l<=n3/2)
      ll = l;
    else {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0; k<n2; k++) {
      jjj=2;
      if(k<=n2/2)
        kk=k;
      else {
        kk = -(k-kkk);
        kkk+=2;
      }
      for(j=0; j<n1; j++) {
        if(j<=n1/2)
          jj=j;
        else {
          jj = -(j-jjj);
          jjj+=2;
        }
        smooth[j+k*n1+l*n1*n2] = float(exp(-0.5f*(jj*dx*jj*dx*sigma_inv(0,0)
                                                   +kk*dy*kk*dy*sigma_inv(1,1)
                                                   //+ll*dz_*ll*dz_*sigma_inv[2][2]
                                                   +2*jj*dx*kk*dy*sigma_inv(1,0))));
                                                   //+2*jj*dx_*ll*dz_*sigma_inv[2][0]
                                                   //+2*kk*dy_*ll*dz_*sigma_inv[2][1])));
        sum += smooth[j+k*n1+l*n1*n2];
      }
    }
  }

  // normalize smoother
  for(l=0;l<n3;l++)
    for(k=0;k<n2;k++)
      for(j=0;j<n1;j++)
        smooth[j+k*n1+l*n1*n2]/=sum;

  smoother->fillInFromArray(smooth); //No mode/randomaccess
  //normalizing constant for the smoother
  //smoother->multiplyByScalar(static_cast<float>(1.0/sum)); //No mode required

  delete [] smooth;
}
