#include "rplib/pdf2dempiricalwithtrend.h"
#include "rplib/pdf3dempirical.h"
#include "src/definitions.h"

Pdf2DEmpiricalWithTrend::Pdf2DEmpiricalWithTrend(const std::vector<double> & d1,   // first dimension of data points
                                                 const std::vector<double> & d2,                     // second dimension of data points
                                                 const std::vector<double> & d3,                     // third dimension of data points
                                                 const std::vector<double> & t1,                     // first dimension of trend values
                                                 const std::vector<double> & t2,                     // second dimension of trend values
                                                 const NRLib::Trend2D & a1,                          // Coefficients for first dimension reduction
                                                 const NRLib::Trend2D & b1,
                                                 const NRLib::Trend2D & c1,
                                                 const NRLib::Trend2D & a2,                          // Coefficients for second dimesion reduction
                                                 const NRLib::Trend2D & b2,
                                                 const NRLib::Trend2D & c2,
                                                 int n1,
                                                 int n2,
                                                 int n3,
                                                 double smooth_var1,
                                                 double smooth_var2,
                                                 double smooth_var3,
                                                 double smooth_corr12,
                                                 double smooth_corr13,
                                                 double smooth_corr23)
                                                 :
a1_(a1),
b1_(b1),
c1_(c1),
a2_(a2),
b2_(b2),
c2_(c2)
{

  int data_size = static_cast<int>(t1.size());
  std::vector<double> linear_combination_dim1(data_size);
  std::vector<double> linear_combination_dim2(data_size);

  int dummy = -1;

  // Loop through all data points
  for(int i = 0; i<data_size; i++){
    // Finding linear combinations of the datapoints given the dimension reduction coefficients and the trend values
    linear_combination_dim1[i] = a1_.GetValue(t1[i], t2[i], dummy)*d1[i] +
                                 b1_.GetValue(t1[i], t2[i], dummy)*d2[i] +
                                 c1_.GetValue(t1[i], t2[i], dummy)*d3[i];
    linear_combination_dim2[i] = a2_.GetValue(t1[i], t2[i], dummy)*d1[i] +
                                 b2_.GetValue(t1[i], t2[i], dummy)*d2[i] +
                                 c2_.GetValue(t1[i], t2[i], dummy)*d3[i];
  }

  // Creating probability density functions based on the linear combination of data points and the trend values
  density_pdf1_ = new Pdf3DEmpirical(linear_combination_dim1, t1, t2, n1, n2, n3, smooth_var1, smooth_var2, smooth_var3, smooth_corr12, smooth_corr13, smooth_corr23);
  density_pdf2_ = new Pdf3DEmpirical(linear_combination_dim2, t1, t2, n1, n2, n3, smooth_var1, smooth_var2, smooth_var3, smooth_corr12, smooth_corr13, smooth_corr23);


  // A note on the normalization of the probability density function:
  // Each of the probability density functions density_pdf1_ and density_pdf2_ are normalized (with respect to "themselves"). The probability density function
  // that this class as a total represents, is the pointwise product of these two probability density functions, density_pdf1_*density_pdf2_ given trend values. (By pointwise meaning
  // for each set of trend values). The pointwise product should for each set of trend values be normalized, that is equal to 1. (Ragnar vær så snill og korriger forståelsen her!!)
  // Thus we need to find normalizing coefficients for each set of trend values.

  // Setting parameters regarding the trend vectors
  // These are class variables for later loop-up in scaling_coefficients_ in the density-function.
  t1_min_ = *(std::min_element(t1.begin(), t1.end()));  // Bør det være litt utvidelse i grensene her?
  t1_max_ = *(std::max_element(t1.begin(), t1.end()));
  t2_min_ = *(std::min_element(t2.begin(), t2.end()));
  t2_max_ = *(std::max_element(t2.begin(), t2.end()));
  dt1_ = (t1_max_ - t1_min_)/t1.size();
  dt2_ = (t2_max_ - t2_min_)/t2.size();

  // For storing of the normalizing coefficients
  scaling_coefficients_.Resize(t1.size(), t2.size(), 0);
  int i = 0;  // index in scaling_coefficients_
  int j = 0;

  //Loop through all combinations of trend values to find normalizing coefficient for each set of trend values
  for(double t1_iter = t1_min_; t1_iter<t1_max_; t1_iter+=dt1_){
    i++;
    for(double t2_iter = t2_min_; t2_iter<t2_max_; t2_iter+=dt2_){
      j++;

      double sum_pdf1 = 0;
      double sum_pdf2 = 0;

      for (double x1 = density_pdf1_->getV1Min(); x1<density_pdf1_->getV1Max(); x1+=density_pdf1_->getdV1()){
        sum_pdf1 += density_pdf1_->density(x1, t1_iter, t2_iter);
      }
      for (double x2 = density_pdf2_->getV1Min(); x2<density_pdf2_->getV1Max(); x2+=density_pdf2_->getdV1()) {
        sum_pdf2 += density_pdf2_->density(x2, t1_iter, t2_iter);
      }

      scaling_coefficients_(i,j) = 1/(sum_pdf1*sum_pdf2);
    }
  }

}

Pdf2DEmpiricalWithTrend::~Pdf2DEmpiricalWithTrend()
{
}

double
Pdf2DEmpiricalWithTrend::density(const double & vp,
                                 const double & vs,
                                 const double & rho,
                                 const double & s1,
                                 const double & s2) const
{

  double dummy = -1;

  // Find indices of input values for the trend values to access correct scaling coefficient
  int i = static_cast<int>(floor((s1 - t1_min_)/dt1_));
  int j = static_cast<int>(floor((s2 - t2_min_)/dt2_));

  double coef = scaling_coefficients_(i,j);

  // linear combination of input datapoint
  double x0 = a1_.GetValue(s1, s2, dummy)*vp +
              b1_.GetValue(s1, s2, dummy)*vs +
              c1_.GetValue(s1, s2, dummy)*rho;

  double returnvalue = coef * density_pdf1_->density(x0, s1, s2, dummy, dummy) * density_pdf2_->density(x0, s1, s2, dummy, dummy);

  return (returnvalue);
}
