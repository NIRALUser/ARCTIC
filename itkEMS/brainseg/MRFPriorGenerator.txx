
#ifndef _MRFPriorGenerator_txx
#define _MRFPriorGenerator_txx

#include "MRFPriorGenerator.h"

#include "itkNeighborhoodIterator.h"

template <class TLabelImage, class TProbabilityImage>
MRFPriorGenerator <TLabelImage, TProbabilityImage>
::MRFPriorGenerator()
{

}

template <class TLabelImage, class TProbabilityImage>
void
MRFPriorGenerator <TLabelImage, TProbabilityImage>
::EstimateInteractionWeights(LabelImagePointer labelImg)
{
  SymmetricMatrixType W(numLabels, numLabels);
  W.set_identity();

  unsigned int numParams = numLabels * (numLabels+1) / 2;

  unsigned int iter = 0;
  while (true)
  {
    ++iter;

    SymmetricMatrixType Wnext = W;

    for (unsigned int r = 0; r < numLabels; r++)
      for (unsigned int c = r; c < numLabels; c++)
      {
        Wnext(r, c) += rng->GenerateNormal(0.0, wvar);
      }

    if (iter <= m_BurnInIterations)
      continue;

    double Enext = this->ComputeTotalEnergy(labelImg, Wnext);

    if (iter > m_MaximumIterations)
      break;
  }
}

template <class TLabelImage, class TProbabilityImage>
double
MRFPriorGenerator <TLabelImage, TProbabilityImage>
::ComputeTotalEnergy(LabelImagePointer labelImg, const SymmetricMatrixType& W)
{

}


#endif
