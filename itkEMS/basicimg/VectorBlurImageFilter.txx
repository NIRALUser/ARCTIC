
#ifndef _VectorBlurImageFilter_txx
#define _VectorBlurImageFilter_txx
#include "VectorBlurImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

#include <vector>

template <class TInputImage, class TOutputImage>
VectorBlurImageFilter<TInputImage, TOutputImage>
::VectorBlurImageFilter()
{
  m_Radius.Fill(1);
  m_Variance = 1.0;
}

template <class TInputImage, class TOutputImage>
void 
VectorBlurImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (itk::InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
  
  if ( !inputPtr || !outputPtr )
  {
    return;
  }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
  {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
  }
  else
  {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
    itk::OStringStream msg;
    msg << static_cast<const char *>(this->GetNameOfClass())
        << "::GenerateInputRequestedRegion()";
    e.SetLocation(msg.str().c_str());
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
  }
}


template< class TInputImage, class TOutputImage>
void
VectorBlurImageFilter< TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  unsigned int i;
  itk::ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

  itk::ConstNeighborhoodIterator<InputImageType> bit;
  itk::ImageRegionIterator<OutputImageType> it;
  
  // Allocate output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename  InputImageType::ConstPointer input  = this->GetInput();

  typename InputImageType::SpacingType spacing = input->GetSpacing();
  
  // Find the data-set boundary "faces"
  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
  faceList = bC(input, outputRegionForThread, m_Radius);

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;

  // support progress methods/callbacks
  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  for (fit=faceList.begin(); fit != faceList.end(); ++fit)
  { 
    bit = itk::ConstNeighborhoodIterator<InputImageType>(m_Radius,
                                                    input, *fit);
    unsigned int neighborhoodSize = bit.Size();
    it = itk::ImageRegionIterator<OutputImageType>(output, *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    std::vector<double> weights(neighborhoodSize, 0.0);

    double sumWeights = 1e-20;

    for (i = 0; i < neighborhoodSize; ++i)
    {
      typename OutputImageType::OffsetType o = bit.GetOffset(i);
      double w = 0;
      for (unsigned int dim = 0; dim < OutputImageType::GetImageDimension(); dim++)
      {
        double step = o[dim] * spacing[dim];
        w += step*step;
      }
      weights[i] = exp(-0.5 * w / m_Variance);
      sumWeights += weights[i];
    }

    for (i = 0; i < neighborhoodSize; ++i)
      weights[i] = weights[i] / sumWeights;
  
    InputPixelType sum;
    sum.Fill(0);

    while ( ! bit.IsAtEnd() )
    {

      sum.Fill(0);
      for (i = 0; i < neighborhoodSize; ++i)
      {
        sum += bit.GetPixel(i) * weights[i];
      }
      
      it.Set(sum);
      
      ++bit;
      ++it;
      progress.CompletedPixel();
    }
  }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
VectorBlurImageFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os, 
  itk::Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius: " << m_Radius << std::endl;

}

#endif
