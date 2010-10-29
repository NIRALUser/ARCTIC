
////////////////////////////////////////////////////////////////////////////////
//
// Kernel density estimation using spherical Gaussian kernel
//
// Bandwidth selection:
// Barendgea???
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2007


#ifndef _SphericalKernelDensityEstimator_h
#define _SphericalKernelDensityEstimator_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class SphericalKernelDensityEstimator
{

public:

  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  typedef vnl_matrix<float> SampleMatrixType;
  typedef vnl_vector<float> SampleVectorType;

  SphericalKernelDensityEstimator();
  ~SphericalKernelDensityEstimator();

  void SetMaximumIterations(unsigned int n)
  { m_MaximumIterations = n; }

  void SetBurnInIterations(unsigned int n)
  { m_BurnInIterations = n; }

  void SetMinWidth(double h)
  { m_MinWidth = h; }

  // Automatically determine kernel bandwidth using MCMC
  double ComputeKernelBandwidth(const SampleMatrixType& samples, double h0);

  double ComputeDensity(double h, const SampleVectorType& x);

  // Faster density computation by truncating contributions from far sample
  // points
//TODO:
  //double ComputeTruncatedDensity(double h, const SampleVectorType& x);

protected:

  double ComputeLeaveOneOutDensity(double h, unsigned int iexclude);

  double ComputeCrossValidationPosterior(double h);

private:

  unsigned int m_MaximumIterations;
  unsigned int m_BurnInIterations;

  double m_MinWidth;

  SampleMatrixType m_Samples;

};

#endif
