
#include "SphericalKernelDensityEstimator.h"

#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "muException.h"

#include "vnl/vnl_math.h"

#include <cmath>

#define MIN_PROB 1e-20

SphericalKernelDensityEstimator
::SphericalKernelDensityEstimator()
{
  m_MaximumIterations = 1000;

  m_BurnInIterations = 100;

  m_MinWidth = 0.01;
}

SphericalKernelDensityEstimator
::~SphericalKernelDensityEstimator()
{

}

double
SphericalKernelDensityEstimator
::ComputeKernelBandwidth(const SampleMatrixType& samples, double h0)
{
  if (m_BurnInIterations > m_MaximumIterations)
    muExceptionMacro(<< "Burn in iterations must be < max iterations");

  unsigned int n = samples.rows();
  unsigned int dim = samples.cols();

  if (n < 2)
    muExceptionMacro(<< "Not enough samples");
  
  if (dim == 0)
    return 0;

  double hvar = 1.0;

  m_Samples = samples;

  // Metropolis-Hastings for computing ergodic average
  double sumH = 0;
  double movingAveH = 0;

  double h = h0;
  double L = this->ComputeCrossValidationPosterior(h);

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  unsigned int transitions = 0;

  unsigned int checkIters = m_MaximumIterations / 10;

  unsigned int lastAdj = 1;

  double K = 1.5;

  unsigned int iter = 0;
  while (true)
  {
    iter++;

    double hnext = h + rng->GenerateNormal(0, hvar);

    //if ((iter % 100) == 1)
    //  K += 0.1;

    double Lnext = this->ComputeCrossValidationPosterior(hnext);

    //double acc = Lnext / L;
    //acc = pow(acc, K);
    double acc = exp(K * (Lnext - L));
    if (acc > 1.0)
      acc = 1.0;

    double q = rng->GenerateUniformRealOpenInterval();

    if (q <= acc)
    {
      h = hnext;
      L = Lnext;
//std::cout << "new h = " << h << " new L = " << L << std::endl;
      transitions++;
    }

    // Burn in iterations
    if (iter <= m_BurnInIterations)
      continue;

    if ((iter % checkIters) == 0)
    {
      double acc_ratio = transitions / (double)(iter - lastAdj);
      if (acc_ratio < 0.6)
      {
        hvar /= 2.25;
        transitions = 0;
        lastAdj = iter;
      }
    }

    sumH += h;

    double a = 1.0 / (double)(iter-m_BurnInIterations);
    movingAveH = (1.0-a)*movingAveH + a*h;

    if (iter > m_MaximumIterations)
      break;
  }

  h = sumH / (iter-m_BurnInIterations);

/*
SampleVectorType zerov(dim, 0.0);
std::cout << "p(0) = " << this->ComputeDensity(h, zerov) << std::endl;
SampleVectorType x1 = m_Samples.get_row(0);
x1[0] = 2.0;
x1[1] = 2.0;
std::cout << "x1 = " << x1 << std::endl;
std::cout << "p(x1) = " << this->ComputeDensity(h, x1) << std::endl;

std::cout << "Final h = " << h << " log posterior = " << this->ComputeCrossValidationPosterior(h) << std::endl;
*/

  //return h;
  return movingAveH;
}

double
SphericalKernelDensityEstimator
::ComputeDensity(double h, const SampleVectorType& x)
{
  if (h < m_MinWidth)
    return 0.0;

  double h_sq = h*h;

  unsigned int n = m_Samples.rows();
  unsigned int d = m_Samples.columns();

  double p = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    SampleVectorType dvec = m_Samples.get_row(k) - x;

    p += exp(-0.5 * dvec.squared_magnitude() / h_sq);
  }

  p /= n * pow(2.0*h_sq*vnl_math::pi, d/2.0) + 1e-20;

  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;

  return p;
}

double
SphericalKernelDensityEstimator
::ComputeLeaveOneOutDensity(double h,  unsigned int iexclude)
{
  double h_sq = h*h;

  unsigned int n = m_Samples.rows();
  unsigned int d = m_Samples.columns();

  SampleVectorType x = m_Samples.get_row(iexclude);

  double std_denom = sqrt(2.0*vnl_math::pi);

  double p = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    if (k == iexclude)
      continue;

    SampleVectorType dvec = m_Samples.get_row(k) - x;

    p += exp(-0.5 * dvec.squared_magnitude() / h_sq);
  }

  p /= (n-1) * pow(2*h_sq*vnl_math::pi, d/2.0) + 1e-20;

  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;

  return p;
}

double
SphericalKernelDensityEstimator
::ComputeCrossValidationPosterior(double h)
{
  if (h < m_MinWidth)
    //return 0.0;
    return -vnl_huge_val(1.0);

  unsigned int n = m_Samples.rows();
  unsigned int d = m_Samples.columns();

  // Prior on kernel widths
  double lambda = 0.5;
  double prior = 0;
  for (unsigned int i = 0; i < d; i++)
  {
    double pr_i = 1.0 / (1.0 + lambda*h*h);
    if (pr_i >= 1e-20)
      prior += log(pr_i);
  }

  double likelihood = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    double fh_k = this->ComputeLeaveOneOutDensity(h, k);
    if (fh_k >= 1e-20)
      likelihood += log(fh_k);
  }
  //likelihood /= n;

  double post = prior + likelihood;

  //return exp(post) + 1e-20;
  return post;
}
