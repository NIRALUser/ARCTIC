
//
// Generate MRF priors from a label image
//

#ifndef _MRFPriorGenerator_h
#define _MRFPriorGenerator_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkSmartPointer.h"

#include "vnl/vnl_sym_matrix.h"

#include <vector>

template <class TLabelImage, class TProbabilityImage>
class MRFPriorGenerator: public itk::Object
{
public:
  // Standard class typedefs
  typedef MRFPriorGenerator Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // Method for creation through the object factory
  itkNewMacro(Self);

  // The dimension of the image we're working with
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TLabelImage::ImageDimension);

  typedef vnl_sym_matrix<double> SymmetricMatrixType;

  typedef TLabelImage LabelImageType;
  typedef typename LabelImageType::Pointer LabelImagePointer;
  typedef typename LabelImageType::IndexType LabelImageIndexType;
  typedef typename LabelImageType::OffsetType LabelImageOffsetType;
  typedef typename LabelImageType::PixelType LabelImagePixelType;
  typedef typename LabelImageType::RegionType LabelImageRegionType;
  typedef typename LabelImageType::SizeType LabelImageSizeType;
  typedef typename LabelImageType::SpacingType LabelImageSpacingType;

  typedef TProbabilityImage ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::OffsetType ProbabilityImageOffsetType;
  typedef typename ProbabilityImageType::PixelType ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType ProbabilityImageSizeType;
  typedef typename ProbabilityImageType::SpacingType ProbabilityImageSpacingType;

  typedef std::vector<ProbabilityImagePointer> ProbabilityImageList;

  ProbabilityImageList GenerateMRFPriors(LabelImagePointer img);

protected:

  MRFPriorGenerator();
  ~MRFPriorGenerator() { }

  void EstimateInteractionWeights(LabelImagePointer img);

  unsigned int m_NeighborhoodRadius;
  double m_NeighborhoodStdDev;

  SymmetricMatrixType m_InteractionWeights;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "MRFPriorGenerator.txx"
#endif

#endif
