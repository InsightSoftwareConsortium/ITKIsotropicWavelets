/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkFrequencyShrinkImageFilter_h
#define itkFrequencyShrinkImageFilter_h

#include <itkShrinkImageFilter.h>
#include <itkFrequencyBandImageFilter.h>
#include <itkEnableIf.h>

namespace itk
{
/** \class FrequencyShrinkImageFilter
 * \brief Reduce the size of an image in the frequency domain by an integer
 * factor --fixed to 2 at the moment-- in each dimension.
 * This filter discard all the high frequency bins depending on the shrink factor
 * Example: N is even.
 * | 0        1   ... N/2-1   N/2             : N/2+1 ... N-1 |
 * | 0 (DC)   Low ... High    Highest-Nyquist : High  ... Low |
 * Considering a shrink factor of 2.
 * | 0  1 ... N/4 ... N/2 : N/2+1 ... N-1-N/4 ... N-1 |
 * The output involves no interpolation, just chopping off high-frequencies.
 * | 0  1 ... N/4 : N-1-N/4 ... N-1 |
 * Example: N is odd.
 * | 0        1   ... N/2-1   N/2     : N/2+1    ... N-1 |
 * | 0 (DC)   Low ... High    Nyq-pos : Nyq-neg  ... Low |
 * If freq are generated from an FFT of a real image, then the input is hermitian;
 * 0 (DC)
 * I(1) == I(N-1)
 * I(2) == I(N-2)
 * ...
 * I(N/2) if N=even. Unique Nyquist: shared between pos and neg freqs.
 * OR
 * I(N/2) == I((N+1)/2) if N=odd. Nyquist has pos and neg components.
 *
 * Example (Odd):
 * inputSize     = 9
 * shrinkFactors = [2]
 * outputSize    = 4
 *                         f_pos  |  f_neg
 * inputImageIndices   = 0 1 2 3 4 5 6 7 8
 * outputImageIndices  = 0 1      2      3
 * So: input(8) == output(3), etc
 *
 * Example (Even):
 * inputSize     = 8
 * shrinkFactors = [2]
 * outputSize    = 4
 *                         f_pos | f_neg
 * inputImageIndices   = 0 1 2 3 4 5 6 7
 * outputImageIndices  = 0 1     2     3
 * So: input(7) == output(3), etc

 * The output image size in each dimension is given by:
 * outputSize[j] = std::floor(inputSize[j]/shrinkFactor[j]);
 *
 * This code was contributed in the Insight Journal paper:
 * https://hdl.handle.net....
 *
 * \ingroup IsotropicWavelets
 */
template<typename TImageType>
class FrequencyShrinkImageFilter:
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  /** Standard class typedefs. */
  typedef FrequencyShrinkImageFilter                 Self;
  typedef ImageToImageFilter<TImageType, TImageType> Superclass;
  typedef SmartPointer<Self>                         Pointer;
  typedef SmartPointer<const Self>                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FrequencyShrinkImageFilter, ImageToImageFilter);

  /** Typedef to images */
  typedef TImageType                       ImageType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename TImageType::IndexType   IndexType;
  typedef typename TImageType::PixelType   PixelType;

  /** Typedef to describe the output image region type. */
  typedef typename TImageType::RegionType ImageRegionType;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImageType::ImageDimension );
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TImageType::ImageDimension );

  typedef FixedArray< unsigned int, ImageDimension > ShrinkFactorsType;
  typedef FrequencyBandImageFilter<TImageType>       FrequencyBandFilterType;

  /** Set the shrink factors. Values are clamped to
   * a minimum value of 1. Default is 1 for all dimensions. */
  itkSetMacro(ShrinkFactors, ShrinkFactorsType);
  void SetShrinkFactors(unsigned int factor);

  void SetShrinkFactor(unsigned int i, unsigned int factor);

  /** Get the shrink factors. */
  itkGetConstReferenceMacro(ShrinkFactors, ShrinkFactorsType);

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /** FrequencyShrinkImageFilter needs a larger input requested region than the output
   * requested region.  As such, FrequencyShrinkImageFilter needs to provide an
   * implementation for GenerateInputRequestedRegion() in order to inform the
   * pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( ImageTypeHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< typename TImageType::PixelType > ) );
  /** End concept checking */
#endif

  itkGetConstReferenceMacro(ApplyBandFilter, bool);
  itkSetMacro(ApplyBandFilter, bool);
  itkBooleanMacro(ApplyBandFilter);

  itkGetMacro(FrequencyBandFilter, typename FrequencyBandFilterType::Pointer);

protected:
  FrequencyShrinkImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

private:
  FrequencyShrinkImageFilter(const Self&) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;

  ShrinkFactorsType                         m_ShrinkFactors;
  bool                                      m_ApplyBandFilter;
  typename FrequencyBandFilterType::Pointer m_FrequencyBandFilter;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFrequencyShrinkImageFilter.hxx"
#endif

#endif
