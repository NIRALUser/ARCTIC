
////////////////////////////////////////////////////////////////////////////////
//
// Kernel density estimation using Gaussian kernel with diagonal covariance
//
// Fast density estimation performed using the truncation method described in:
// Byeungwoo Jeon and David A. Landgrebe, "Fast Parzen Density Estimation Using
// Clustering-Based Branch and Bound," IEEE Transa PAMI, Vol. 16, No. 9,
// pp 950-954, September 1994.
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2007


#ifndef _DiagonalKernelDensityEstimator_h
#define _DiagonalKernelDensityEstimator_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class DiagonalKernelDensityEstimator
{

public:

  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  typedef vnl_matrix<double> SampleMatrixType;
  typedef vnl_vector<double> SampleVectorType;

  typedef vnl_vector<unsigned int> LabelVectorType;
  typedef vnl_vector<unsigned char> MarkVectorType;

  DiagonalKernelDensityEstimator();
  ~DiagonalKernelDensityEstimator();

  DiagonalKernelDensityEstimator&
    operator=(const DiagonalKernelDensityEstimator& d);

  void SetMaximumIterations(unsigned int n)
  { m_MaximumIterations = n; }

  void SetBurnInIterations(unsigned int n)
  { m_BurnInIterations = n; }

  void SetMinWidth(double h)
  { m_MinWidth = h; }

  void SetLambda(double d)
  { m_Lambda = d; }

  void SetTruncationDistance(double d)
  { m_TruncationDistance = d; }

  // Min number of clusters for computing truncated density
  void SetMinClusterCount(unsigned int n)
  { m_MinClusterCount = n; }

  void SetInputSamples(const SampleMatrixType& samples)
  { m_Samples = samples; this->ClusterSamples(); }

  // Automatically determine kernel bandwidth using MCMC, give initial
  // guess and random walk step (std dev)
  VectorType ComputeKernelBandwidth(const VectorType& h_init, double hstep);

  double ComputeDensity(const VectorType& h, const SampleVectorType& x);

//TODO:
  // Faster density computation by truncating contributions from far sample
  // points
  double ComputeTruncatedDensity(
    const VectorType& h, const SampleVectorType& x);

  void ResetTruncationBandwidth();

protected:

  double ComputeLeaveOneOutDensity(const VectorType& h, unsigned int iexclude);

  double ComputeCrossValidationPosterior(const VectorType& h);

  void ClusterSamples();

  void UpdateDMax(const VectorType& h);
  void UpdateScaledSamples(const VectorType& h);

private:

  unsigned int m_MaximumIterations;
  unsigned int m_BurnInIterations;

  double m_MinWidth;

  double m_Lambda;

  unsigned int m_MinClusterCount;

  double m_TruncationDistance;

  SampleMatrixType m_Samples;

  LabelVectorType m_SampleClusterLabels;
  SampleMatrixType m_SampleClusterMeans;

  // Extent of each cluster, for truncation
  VectorType m_SampleClusterDMax;

  MarkVectorType m_SampleTestMarkers;

  // Last bandwidth used for computing truncated density
  VectorType m_LastTruncationBandwidth;

  // Samples scaled by 1/h using cached bandwidth
  SampleMatrixType m_ScaledSamples;

};

#endif
