
#include "DiagonalKernelDensityEstimator.h"

#include "DynArray.h"
#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "KMeansEstimator.h"
#include "MSTClusteringProcess.h"

#include "muException.h"

#include "vnl/vnl_math.h"

#include <cmath>

#define MIN_PROB 1e-20

//#define LOGP(x) (100.0 * log((x) + 1.0))
//#define EXPP(x) (exp((x) / 100.0) - 1.0)
#define LOGP(x) (log(x))
#define EXPP(x) (exp(x))

DiagonalKernelDensityEstimator
::DiagonalKernelDensityEstimator()
{
  m_MaximumIterations = 150;

  m_BurnInIterations = 10;

  m_Lambda = 0.05;

  m_MinWidth = 0.01;

  m_TruncationDistance = 4.0;

  m_MinClusterCount = 10;

  m_LastTruncationBandwidth = VectorType(0, 0.0);
}

DiagonalKernelDensityEstimator
::~DiagonalKernelDensityEstimator()
{

}

DiagonalKernelDensityEstimator&
DiagonalKernelDensityEstimator
::operator=(const DiagonalKernelDensityEstimator& d)
{
  this->m_MaximumIterations = d.m_MaximumIterations;
  this->m_BurnInIterations = d.m_BurnInIterations;

  this->m_MinWidth = d.m_MinWidth;
  this->m_Lambda = d.m_Lambda;

  this->m_MinClusterCount = d.m_MinClusterCount;
  this->m_TruncationDistance = d.m_TruncationDistance;

  this->m_Samples = d.m_Samples;

  this->m_SampleClusterLabels = d.m_SampleClusterLabels;
  this->m_SampleClusterMeans = d.m_SampleClusterMeans;

  this->m_SampleClusterDMax = d.m_SampleClusterDMax;

  this->m_SampleTestMarkers = d.m_SampleTestMarkers;
  this->m_LastTruncationBandwidth = d.m_LastTruncationBandwidth;
  this->m_ScaledSamples = d.m_ScaledSamples;

  return *this;
}

void
DiagonalKernelDensityEstimator
::ResetTruncationBandwidth()
{
  m_LastTruncationBandwidth = VectorType(0, 0.0);
}

DiagonalKernelDensityEstimator::VectorType
DiagonalKernelDensityEstimator
::ComputeKernelBandwidth(const VectorType& h_init, double hstep)
{
  if (m_BurnInIterations > m_MaximumIterations)
    muExceptionMacro(<< "Burn in iterations must be < max iterations");

  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.cols();

  if (n < 2)
    muExceptionMacro(<< "Not enough samples");

  if (dim == 0)
    return VectorType(0, 0.0);

std::cout << "bandwidth est with " << n << " samples" << std::endl;

  // Metropolis-Hastings for computing the ergodic average of h
  VectorType sumH(dim, 0.0);
  VectorType movingAveH(dim, 0.0);

  VectorType h = h_init;
  double post = this->ComputeCrossValidationPosterior(h);

  VectorType h_opt = h_init;
  double max_post = post;

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  unsigned int transitions = 0;

  unsigned int lastAdjIter = 1;

  unsigned int checkIters = m_MaximumIterations / 10;

  double K = 1.0;

  double hvar = hstep * hstep;

  unsigned int iter = 0;
  while (true)
  {
    iter++;

    VectorType hnext = h;
    for (unsigned int i = 0; i < dim; i++)
      //hnext[i] += rng->GenerateNormal(0, hvar);
      hnext[i] += hstep * (2.0*rng->GenerateUniformRealClosedInterval() - 1.0);

    if ((iter % 100) == 1)
      K += 0.1;

    double postnext = this->ComputeCrossValidationPosterior(hnext);

    //double acc = postnext / post;
    //acc = pow(acc, K);

    double acc = EXPP(K * (postnext - post));
    //double acc = EXPP(postnext - post);

    if (acc > 1.0)
      acc = 1.0;

    double q = rng->GenerateUniformRealOpenInterval();

    if (q <= acc)
    {
      h = hnext;
      post = postnext;
//std::cout << "new h = " << h << " new post = " << post << std::endl;

      if (post > max_post)
      {
        h_opt = h;
        max_post = post;
//std::cout << "new h_opt = " << h_opt << " new max post = " << max_post << std::endl;
      }

      transitions++;
    }

#if 0
    if ((iter % checkIters) == 0)
    {
      K += 0.2;
      double accRate = transitions / (double)(iter-lastAdjIter);
      if (accRate < 0.6)
      {
        hvar /= 2.25;
        transitions = 0;
        lastAdjIter = iter;
      }
    }
#endif

    // Burn in iterations
    if (iter <= m_BurnInIterations)
      continue;

    sumH += h;

    // Update moving average
    double a = 1.0 / (double)(iter - m_BurnInIterations);
    movingAveH = (movingAveH * (1.0-a)) + (h * a);

    if (iter > m_MaximumIterations)
      break;
  }

  //h = sumH / (iter-m_BurnInIterations);

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
  //return h_opt;
}

double
DiagonalKernelDensityEstimator
::ComputeDensity(const VectorType& h, const SampleVectorType& x)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  for (unsigned int i = 0; i < dim; i++)
    if (h[i] < m_MinWidth)
      return 0.0;

  double p = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    SampleVectorType y_k = m_Samples.get_row(k);

    VectorType dvec(dim);
    for (unsigned int i = 0; i < dim; i++)
      dvec[i] = (x[i] - y_k[i]) / h[i];

    p += exp(-0.5 * dvec.squared_magnitude());
  }

  double detH = 1.0;
  for (unsigned int i = 0; i < dim; i++)
    detH *= h[i];

  p /= n * detH * pow(2.0*vnl_math::pi, dim/2.0) + 1e-20;

/*
  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;
*/

  return p;
}

double
DiagonalKernelDensityEstimator
::ComputeLeaveOneOutDensity(const VectorType& h,  unsigned int iexclude)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  SampleVectorType x = m_Samples.get_row(iexclude);

  double p = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    if (k == iexclude)
      continue;

    SampleVectorType y_k = m_Samples.get_row(k);

    VectorType dvec(dim);
    for (unsigned int i = 0; i < dim; i++)
      dvec[i] = (x[i] - y_k[i]) / h[i];

    p += exp(-0.5 * dvec.squared_magnitude());
  }

  double detH = 1.0;
  for (unsigned int i = 0; i < dim; i++)
    detH *= h[i];

  p /= (n-1) * detH * pow(2.0*vnl_math::pi, dim/2.0) + 1e-20;

//std::cout << "PP tot lo = " << p << std::endl;

/*
  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;
*/

  return p;
}

double
DiagonalKernelDensityEstimator
::ComputeCrossValidationPosterior(const VectorType& h)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  for (unsigned int i = 0; i < dim; i++)
    if (h[i] < m_MinWidth)
      //return 0.0;
      return -vnl_huge_val(1.0);

  // Prior on kernel widths
  double prior = 0;
  for (unsigned int i = 0; i < dim; i++)
  {
    double pr_i = 1.0 / (1.0 + m_Lambda * h[i]*h[i]);
    if (pr_i > 1e-20)
      prior += LOGP(pr_i);
//PP
    //prior -= m_Lambda * h[i]*h[i];
    // Should have smaller kernels with more samples ???
    //prior -= m_Lambda * log((double)n) * fabs(h[i]);
  }

#if 0
  double likelihood = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    double fh_k = this->ComputeLeaveOneOutDensity(h, k);
    if (fh_k > 1e-20)
      likelihood += LOGP(fh_k);
  }
#else
  double likelihood = 0;
  unsigned int numTest = 0;
  for (unsigned int k = 0; k < n; k++)
  {
    if (m_SampleTestMarkers[k] == 0)
      continue;
    double fh_k = this->ComputeLeaveOneOutDensity(h, k);
    if (fh_k >= 1e-20)
      likelihood += LOGP(fh_k);
    numTest++;
  }
#endif

// PP TEST
// likelihood is negative????

  double post = prior + likelihood;

  //return exp(post) + 1e-20;
  //return post + 1e-20;
  return post;
}

void
DiagonalKernelDensityEstimator
::ClusterSamples()
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

#if 0
  DynArray<MSTClusteringProcess::VertexType> mstsamples;
  mstsamples.Allocate(n);
  for (unsigned int k = 0; k < n; k++)
  {
    MSTClusteringProcess::VertexType v(dim);
    for (unsigned int i = 0; i < dim; i++)
      v[i] = m_Samples(k, i);
    mstsamples.Append(v);
  }

  // Set up MST clustering, with sorting based on size
  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(mstsamples);
  mstProc.SortOn();

  // Compute cluster labels
  m_SampleClusterLabels = LabelVectorType(n, 0);

  unsigned int numClusters = 0;
  for (double T = 3.0; T >= 1.0; T -= 0.1)
  {
std::cout << "PP KDE MST clust T = " << T << std::endl;
    numClusters = mstProc.GetClusters(m_SampleClusterLabels.data_block(), T);
    if (numClusters > m_MinClusterCount)
      break;
  }
std::cout << "KDE truncation with " << numClusters << " clusters" << std::endl;

  // Compute cluster means
  m_SampleClusterMeans = SampleMatrixType(numClusters, dim, 0.0);

  for (unsigned int c = 0; c < numClusters; c++)
  {
    unsigned int countMembers = 0;
    for (unsigned int k = 0; k < n; k++)
    {
      if (m_SampleClusterLabels[k] != c)
        continue;
      countMembers++;
      for (unsigned int i = 0; i < dim; i++)
        m_SampleClusterMeans(c, i) += m_Samples(k, i);
    }

    if (countMembers > 0)
    {
      for (unsigned int i = 0; i < dim; i++)
        m_SampleClusterMeans(c, i) /= countMembers;
    }
  }

//std::cout << "Cluster means = \n" << m_SampleClusterMeans << std::endl;
#endif

  unsigned int numClusters = m_MinClusterCount;

  if (numClusters > n)
    numClusters = n;

  KMeansEstimator kmeans;
  kmeans.SetNumberOfClusters(numClusters);
  kmeans.SetNumberOfStarts(10);
  kmeans.UseKdTreeOn();
  kmeans.SetMaximumIterations(100);
  kmeans.SetInput(m_Samples);

  m_SampleClusterMeans = kmeans.GetMeans();

  DynArray<unsigned int> labels = kmeans.GetLabels();
  m_SampleClusterLabels = LabelVectorType(n, 0);
  for (unsigned int i = 0; i < n; i++)
    m_SampleClusterLabels[i] = labels[i];

//TODO: not exactly part of clustering

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  DynArray<unsigned int> testInd;
  for (unsigned int c = 0; c < numClusters; c++)
  {
    DynArray<unsigned int> clusterInd;
    for (unsigned int i = 0; i < n; i++)
      if (m_SampleClusterLabels[i] == c)
        clusterInd.Append(i);

    unsigned int numMembers = clusterInd.GetSize();
    if (numMembers == 0)
      continue;
    if (numMembers == 1)
    {
      testInd.Append(clusterInd[0]);
      continue;
    }

    unsigned int numTest = numMembers/5 + 1;

    unsigned int* picks =
      rng->GenerateIntegerSequence(numTest, numMembers-1);
    for (unsigned int k = 0; k < numTest; k++)
      testInd.Append(picks[k]);
    delete [] picks;
  }

  m_SampleTestMarkers = MarkVectorType(n, 0);
  for (unsigned int i = 0; i < testInd.GetSize(); i++)
    m_SampleTestMarkers[testInd[i]] = 1;

/*
  unsigned int* testInd = rng->GenerateIntegerSequence(numTestPts, n-1);

  m_SampleTestMarkers = MarkVectorType(n, 0);
  for (unsigned int i = 0; i < numTestPts; i++)
    m_SampleTestMarkers[testInd[i]] = 1;

  delete [] testInd;
*/

}

void
DiagonalKernelDensityEstimator
::UpdateDMax(const VectorType& h)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  unsigned int numClusters = m_SampleClusterMeans.rows();

  // Compute maximum within-cluster displacements (unscaled wrt kernel width)
  m_SampleClusterDMax = VectorType(numClusters, 0.0);

  for (unsigned int c = 0; c < numClusters; c++)
  {
    SampleVectorType mu_c = m_SampleClusterMeans.get_row(c);

    double dmax = 0.0;
    for (unsigned int k = 0; k < n; k++)
    {
      if (m_SampleClusterLabels[k] != c)
        continue;
      SampleVectorType dvec = m_Samples.get_row(k) - mu_c;
      for (unsigned int i = 0; i < dim; i++)
        dvec[i] /= h[i];
      double d = dvec.magnitude();
      if (d > dmax)
        dmax = d;
    }

    m_SampleClusterDMax[c] = dmax;
  }

//std::cout << "DMax = " << m_SampleClusterDMax << std::endl;

}

void
DiagonalKernelDensityEstimator
::UpdateScaledSamples(const VectorType& h)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  m_ScaledSamples.set_size(n, dim);
  for (unsigned int i = 0; i < n; i++)
  {
    SampleVectorType y = m_Samples.get_row(i);
    for (unsigned int d = 0; d < dim; d++)
      y[d] /= h[d];
    m_ScaledSamples.set_row(i, y);
  }
}

double
DiagonalKernelDensityEstimator
::ComputeTruncatedDensity(const VectorType& h, const SampleVectorType& x)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  for (unsigned int i = 0; i < dim; i++)
    if (h[i] < m_MinWidth)
      return 0.0;

  if (m_LastTruncationBandwidth.size() != h.size())
  {
    this->UpdateDMax(h);
    this->UpdateScaledSamples(h);
  }
  else
  {
    for (unsigned int i = 0; i < dim; i++)
      if (fabs(m_LastTruncationBandwidth[i]-h[i]) > 1e-2)
      {
        this->UpdateDMax(h);
        this->UpdateScaledSamples(h);
        break;
      }
  }

  m_LastTruncationBandwidth = h;

  unsigned int numClusters = m_SampleClusterMeans.rows();

  unsigned char* sampleMask = new unsigned char[n];
  for (unsigned int k = 0; k < n; k++)
    sampleMask[k] = 0;

  SampleVectorType x_scaled(dim);
  for (unsigned int i = 0; i < dim; i++)
    x_scaled[i] = x[i] / h[i]; 

  for (unsigned int c = 0; c < numClusters; c++)
  {
    SampleVectorType mu_c = m_SampleClusterMeans.get_row(c);
    for (unsigned int i = 0; i < dim; i++)
      mu_c[i] /= h[i];

    //SampleVectorType dvec = mu_c - x;
    //for (unsigned int i = 0; i < dim; i++)
    //  dvec[i] /= h[i];
    SampleVectorType dvec = mu_c - x_scaled;

    double d = dvec.squared_magnitude() - m_SampleClusterDMax[c];

    if (d >= m_TruncationDistance)
      continue;

    for (unsigned int k = 0; k < n; k++)
    {
      if (m_SampleClusterLabels[k] == c)
        sampleMask[k] = 1;
    }
  }

  unsigned int n_trunc = 0;
  for (unsigned int k = 0; k < n; k++)
    if (sampleMask[k] != 0)
      n_trunc++;

//std::cout << "  using " << n_trunc << " of " << n << " samples" << std::endl;

  if (n_trunc == 0)
  {
    delete [] sampleMask;
    return 0.0;
  }

  double p = 0.0;
  for (unsigned int k = 0; k < n; k++)
  {
    if (sampleMask[k] == 0)
      continue;

    SampleVectorType y_k = m_ScaledSamples.get_row(k);

    VectorType dvec(dim);
    for (unsigned int i = 0; i < dim; i++)
      dvec[i] = x_scaled[i] - y_k[i];
      //dvec[i] = (x[i] - y_k[i]) / h[i];

    p += exp(-0.5 * dvec.squared_magnitude());
  }

  delete [] sampleMask;

  double detH = 1.0;
  for (unsigned int i = 0; i < dim; i++)
    detH *= h[i];

  p /= n * detH * pow(2.0*vnl_math::pi, dim/2.0) + 1e-20;
  //p /= n_trunc * detH * pow(2.0*vnl_math::pi, dim/2.0) + 1e-20;

/*
  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;
*/

  return p;
}
