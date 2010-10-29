
#ifndef _MultivariateGaussian_h
#define _MultivariateGaussian_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class MultivariateGaussian
{

public:

  MultivariateGaussian();
  ~MultivariateGaussian();

  unsigned int GetDimension() const { return m_Dimension; }

  void SetMean(const VectorType& mu);
  void SetCovariance(const MatrixType& cov);

  double Evaluate(const VectorType& x);

  double EvaluateMahalanobisDistance(const VectorType& x);

  VectorType GenerateRandomVariate();

protected:

  void UpdateIntegral();

private:

  VectorType m_Mean;
  MatrixType m_Covariance;

  double m_IntegralValue;

  unsigned int m_Dimension;

};

#endif
