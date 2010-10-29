
#include "MultivariateGaussian.h"

#include "MersenneTwisterRNG.h"

#include <cmath>

MultivariateGaussian
::MultivariateGaussian()
{
  m_Mean = VectorType(1, 0.0);
  m_Covariance = MatrixType(1, 1, 1.0);

  m_NormalizingCoefficient = 1.0 / (2.0*M_PI);
}

MultivariateGaussian
::~MultivariateGaussian()
{

}

void
MultivariateGaussian
::SetMean(const VectorType& mu)
{
  m_Mean = mu;
}

void
MultivariateGaussian
::SetCovariance(const MatrixType& cov)
{
  m_Covariance = cov;

  m_InverseCovariance = MatrixInverseType(cov);

  double dim = cov.rows();

  // Update normalizing coeff
  //m_NormalizingCoefficient = 1.0 / pow(2.0*M_PI*detcov, 1.0/dim);
  m_ZInv = 1.0 / pow(2.0*M_PI*detcov, 1.0/dim);
  m_LogZInv = log(m_ZInv);
}

double
MultivariateGaussian
::Evaluate(const VectorType& x)
{
  double m = this->EvaluateMahalanobis(x);

  //return exp(-m) * m_NormalizingCoefficient;
  return exp(-0.5*m) * m_ZInv;
}

double
MultivariateGaussian
::EvaluateMahalanobisDistance(const VectorType& x)
{
  MatrixType d = x.transpose() * m_InverseCovariance * x;

  return d(0, 0);
}

VectorType
MultivariateGaussian
::GenerateRandomVariate()
{
  // Generate standard multivariate normal

  // Z = cov^0.5 * X + mu
}
