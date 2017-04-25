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
#ifndef itkFrequencyShrinkImageFilter_hxx
#define itkFrequencyShrinkImageFilter_hxx

#include <itkFrequencyShrinkImageFilter.h>
#include <itkProgressReporter.h>
#include <numeric>
#include <functional>
#include "itkInd2Sub.h"
#include <itkPasteImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageLinearIteratorWithIndex.h>
// #include <itkGaussianSpatialFunction.h>
// #include <itkFrequencyImageRegionIteratorWithIndex.h>

namespace itk
{
template<class TImageType>
FrequencyShrinkImageFilter< TImageType >
::FrequencyShrinkImageFilter()
: m_ApplyBandFilter(false)
{
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    m_ShrinkFactors[j] = 2;
    }

  this->m_FrequencyBandFilter = FrequencyBandFilterType::New();
  // The band filter only let pass half of the frequencies.
  this->m_FrequencyBandFilter->SetFrequencyThresholdsInRadians( 0.0, Math::pi_over_2 );
  bool lowFreqThresholdPassing  = true;
  bool highFreqThresholdPassing = true;
  this->m_FrequencyBandFilter->SetPassBand( lowFreqThresholdPassing, highFreqThresholdPassing );
  // The band is not radial, but square like.
  this->m_FrequencyBandFilter->SetRadialBand(false);
  // Pass high positive freqs but stop negative high ones, to avoid overlaping.
  this->m_FrequencyBandFilter->SetPassNegativeHighFrequencyThreshold(false);
}

template<class TImageType>
void
FrequencyShrinkImageFilter<TImageType>
::SetShrinkFactors(unsigned int factor)
{
  unsigned int j;

  for ( j = 0; j < ImageDimension; j++ )
    {
    if ( factor != m_ShrinkFactors[j] )
      {
      break;
      }
    }
  if ( j < ImageDimension )
    {
    for ( j = 0; j < ImageDimension; j++ )
      {
      m_ShrinkFactors[j] = factor;
      if ( m_ShrinkFactors[j] < 1 )
        {
        m_ShrinkFactors[j] = 1;
        }
      }
    }
  this->Modified();
}

template<class TImageType>
void
FrequencyShrinkImageFilter<TImageType>
::SetShrinkFactor(unsigned int i, unsigned int factor)
{
  if ( m_ShrinkFactors[i] == factor )
    {
    return;
    }

  m_ShrinkFactors[i] = factor;
  this->Modified();
}

/**
 * Implementation Detail:
 * The implementation calculate the number of different regions in an image,
 * depending on the dimension:
 * numberOfRegions = 2^dim (positive and negative frequencies per dim)
 * then uses function to convert a linear array of regions [0, ..., numberOfRegions - 1]
 * to binary subindices (only two options: positive or negative region)
 * In 3D: numberOfRegions = Nr = 2^3 = 8
 * sizeOfSubindices = [2,2,2]
 * Region = 0       -----> Ind2Sub(   0, [2,2,2]) = [0,0,0]
 * Region = 1       -----> Ind2Sub(   1, [2,2,2]) = [1,0,0]
 * Region = Nr - 1  -----> Ind2Sub(Nr-1, [2,2,2]) = [1,1,1]
 * So, if the result of Ind2Sub is 0 we paste the positive frequencies, if 1, negative freq
 */
template<class TImageType>
void
FrequencyShrinkImageFilter<TImageType>
::GenerateData()
{
  // Get the input and output pointers
  const ImageType *           inputPtr = this->GetInput();
  typename ImageType::Pointer outputPtr = this->GetOutput();
  this->AllocateOutputs();
  outputPtr->FillBuffer(0);

  typename TImageType::SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  // outputSize is Floor(inputSize/ShrinkFactor)
  typename TImageType::SizeType outputSize = outputPtr->GetLargestPossibleRegion().GetSize();
  FixedArray<bool, TImageType::ImageDimension> outputSizeIsEven;
  FixedArray<bool, TImageType::ImageDimension> inputSizeIsEven;
  // sizeInputRegionToCopy is the region to copy shared by all quadrants
  // this excludes the zero bands (that are not shared for all quadrants)
  // also excludes the nyquist band if output size is even.
  // (nyquist bands need to be recalculated anyway)
  // but includes the largest(+region)/smallest(-region) frequency bins if output size is odd.
  typename TImageType::SizeType sizeInputRegionToCopy;
  for ( unsigned int dim=0; dim < TImageType::ImageDimension; ++dim )
  {
    outputSizeIsEven[dim] = (outputSize[dim] % 2 == 0);
    inputSizeIsEven[dim] = (inputSize[dim] % 2 == 0);
    // -1 to exclude the zero band.
    sizeInputRegionToCopy[dim]  = Math::Ceil<SizeValueType>(outputSize[dim]/2.0) - 1;
  }

  const typename TImageType::IndexType indexOrigOut = outputPtr->GetLargestPossibleRegion().GetIndex();
  const typename TImageType::IndexType indexOrigIn= inputPtr->GetLargestPossibleRegion().GetIndex();
  // Manage ImageDimension array linearly:{{{
  FixedArray<unsigned int , ImageDimension> nsizes;
  unsigned int numberOfRegions = 1;
  for (unsigned int dim = 0; dim < ImageDimension; ++dim)
    {
    nsizes[dim] = 2;
    numberOfRegions *= nsizes[dim];
    }
  FixedArray<unsigned int, ImageDimension> subIndices;
  /// }}}

  // Prepare filter to paste the different regions into output.
  typedef itk::PasteImageFilter<ImageType> PasteFilterType;
  typename PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
  pasteFilter->SetSourceImage(inputPtr);
  pasteFilter->SetDestinationImage(outputPtr);
  pasteFilter->InPlaceOn();

  typedef typename ImageType::RegionType RegionType;
  ProgressReporter progress(this, 0, numberOfRegions );

  for (unsigned int n = 0; n < numberOfRegions; ++n)
    {
    subIndices = itk::Ind2Sub<ImageDimension>(n, nsizes);
    RegionType zoneRegion;
    typename ImageType::SizeType zoneSize;
    typename TImageType::IndexType inputIndex  = indexOrigIn;
    typename TImageType::IndexType outputIndex = indexOrigOut;
    for (unsigned int dim = 0; dim < ImageDimension; ++dim)
      {
      if(subIndices[dim] == 0) // positive frequencies
        {
        // + 1 to include the zero band that is in positive regions.
        zoneSize[dim]    = sizeInputRegionToCopy[dim] + 1;
        inputIndex[dim]  = indexOrigIn[dim];
        outputIndex[dim] = indexOrigOut[dim];
        }
      else // negative frequencies
        {
        zoneSize[dim]    = sizeInputRegionToCopy[dim];
        inputIndex[dim]  = indexOrigIn[dim] + inputSize[dim]  - zoneSize[dim];
        outputIndex[dim] = indexOrigOut[dim] + outputSize[dim] - zoneSize[dim];
        }
      // if(subIndices[dim] == 0) // positive frequencies
      //   {
      //   zoneSize[dim]    = sizeInputRegionToCopy[dim];
      //   inputIndex[dim]  = indexOrigOut[dim];
      //   outputIndex[dim] = indexOrigOut[dim];
      //   }
      // else // negative frequencies
      //   {
      //   zoneSize[dim]    = sizeInputRegionToCopy[dim];
      //   inputIndex[dim]  = indexOrigOut[dim] + inputSize[dim]  - zoneSize[dim];
      //   outputIndex[dim] = indexOrigOut[dim] + outputSize[dim] - zoneSize[dim];
      //   }
      }
    zoneRegion.SetIndex(inputIndex);
    zoneRegion.SetSize(zoneSize);
    itkDebugMacro( << "n:" << n << " region: " << zoneRegion);

    pasteFilter->SetSourceRegion(zoneRegion);
    pasteFilter->SetDestinationIndex(outputIndex);
    if (n == numberOfRegions - 1) // Graft the output.
      {
      pasteFilter->GraftOutput(outputPtr);
      pasteFilter->Update();
      this->GraftOutput(pasteFilter->GetOutput());
      }
    else // update output
      {
      pasteFilter->Update();
      outputPtr = pasteFilter->GetOutput();
      }
    progress.CompletedPixel();
    }

  /** Ensure image is hermitian in the Nyquist bands (even)
   * Example: Image 2D size 8, index = [0,...,7]
   * Each quadrant is a region pasted from the original image. The index refers to the input image of size 8. The input image is hermitian, so:
   *   0               0
   * 1 == 7          1 == 3
   * 2 == 6            2      <-2/6 is Nyq in new image.
   * 3 == 5
   *   4  <- Nyq original
   * Hermitian table, using the equivalences above. (assuming imag part zero to avoid working with conjugates) . Note that index 6 is the new Nyquist.
   *    0    1  | 2/6 |  7
   * 0 0,0  1,0 | 2,0 | 1,0
   * 1 0,1  1,1 | 2,1 | 1,1
   * ----------------------
   2/6 0,2  1,2 | 2,2 | 1,2
   * ----------------------
   * 7 0,1  1,1 | 2,1 | 1,1
   *
   * Size = 9x9, [0,...,8], ShrinkedSize = 5x5
   *   0              0
   * 1 == 8        1 == 4
   * 2 == 7        2 == 3
   * 3 == 6
   * 4 == 5
   *    0    1  |  2    7  |  8
   * 0 0,0  1,0 | 2,0  2,0 | 1,0
   * 1 0,1  1,1 | 2,1  2,1 | 1,1
   * ---------------------------
  *  2 0,2  1,2 | 2,2  2,2 | 1,2
   * 7 0,2  1,2 | 2,2  2,2 | 1,2
   * ---------------------------
   * 8 0,1  1,1 | 2,1  2,1 | 1,1
   *
   * Size = 7x7, [0,...,6], ShrinkedSize = 4x4
   *   0              0
   * 1 == 6         1 == 3
   * 2 == 5           2
   * 3 == 4
   *
   *    0    1  | 2/5 |  6
   * 0 0,0  1,0 | 2,0 | 1,0
   * 1 0,1  1,1 | 2,1 | 1,1
   * ----------------------
   2/5 0,2  1,2 | 2,2 | 1,2
   * ----------------------
   * 6 0,1  1,1 | 2,1 | 1,1
   *
   * Size = 10x7, [0,...,9], [0,...6], ShrinkedSize = 5x4
   * {7}  0         0      {10}  0       0
   *   1 == 6    1 == 3       1 == 9  1 == 4
   *   2 == 5       2         2 == 8  2 == 3
   *   3 == 4                 3 == 7
   *                          4 == 6
   *                            5
   *    0    1  | 2/5 |  6
   * 0 0,0  1,0 | 2,0 | 1,0
   * 1 0,1  1,1 | 2,1 | 1,1
   * ----------------------
   * 2 0,2  1,2 | 2,2 | 1,2
   * 8 0,2  1,2 | 2,2 | 1,2
   * ----------------------
   * 9 0,1  1,1 | 2,1 | 1,1
   * /
  */
// Fix Nyquist band. Folding results for hermiticity.
// Nyquist band has to be recalculated for the case that the outputSize is even. Nyquist band is shared between + and - regions, even though only the + region is stored. (The - region value can be calculated with the conjugate from the postiive one)
    {
    // typedef ImageDimension - 1 NyquistDimension1;
    // typedef itk::Size<NyquistDimension1> NyquistSizeType;
    // typedef itk::ImageRegion<NyquistDimension1> NyquistImageRegion;
    for(unsigned int dim = 0; dim < ImageDimension; ++dim)
      {
      // Nyquist band/slice only exist on even dimension.
      if(!outputSizeIsEven[dim])
        {
        continue;
        }

      typename ImageType::SizeType nyquistRegionSize;
      typename ImageType::IndexType nyquistRegionIndex;
      typename ImageType::SizeType positiveRegionSize;
      typename ImageType::IndexType positiveRegionIndex;
      typename ImageType::SizeType negativeRegionSize;
      typename ImageType::IndexType negativeRegionIndex;
      for(unsigned int dimInner = 0; dimInner < ImageDimension; ++dimInner)
        {
        if(dimInner == dim)
          {
          nyquistRegionSize[dimInner] = 1;
          positiveRegionSize[dimInner] = 1;
          negativeRegionSize[dimInner] = 1;
          nyquistRegionIndex[dimInner] = indexOrigOut[dimInner] + outputSize[dimInner]/2;
          positiveRegionIndex[dimInner] = indexOrigIn[dimInner] + (sizeInputRegionToCopy[dimInner] + 1);
          negativeRegionIndex[dimInner] = indexOrigIn[dimInner] +
            inputSize[dimInner] - (sizeInputRegionToCopy[dimInner] + 2);
          }
        else
          {
          nyquistRegionSize[dimInner] = outputSize[dimInner];
          positiveRegionSize[dimInner] = outputSize[dimInner];
          negativeRegionSize[dimInner] = outputSize[dimInner];
          nyquistRegionIndex[dimInner] = indexOrigOut[dimInner];
          positiveRegionIndex[dimInner] = indexOrigIn[dimInner];
          negativeRegionIndex[dimInner] = indexOrigIn[dimInner] +
            inputSize[dimInner] - outputSize[dimInner];
          }
        }

      typedef itk::ImageLinearIteratorWithIndex< ImageType > ImageIterator;
      typedef itk::ImageLinearConstIteratorWithIndex< ImageType > ImageConstIterator;
      unsigned int iteratorLineDirection = (dim + 1) % ImageDimension;
      RegionType nyquistRegion;
      nyquistRegion.SetSize(nyquistRegionSize);
      nyquistRegion.SetIndex(nyquistRegionIndex);
      ImageIterator nyquistIt(outputPtr, nyquistRegion);
      nyquistIt.SetDirection(iteratorLineDirection);
      nyquistIt.GoToBegin();

      RegionType positiveRegion;
      positiveRegion.SetSize(positiveRegionSize);
      positiveRegion.SetIndex(positiveRegionIndex);
      ImageConstIterator positiveIt(inputPtr, positiveRegion);
      positiveIt.SetDirection(iteratorLineDirection);
      positiveIt.GoToBegin();

      RegionType negativeRegion;
      negativeRegion.SetSize(negativeRegionSize);
      negativeRegion.SetIndex(negativeRegionIndex);
      ImageConstIterator negativeIt(inputPtr, negativeRegion);
      negativeIt.SetDirection(iteratorLineDirection);
      negativeIt.GoToBegin();
      negativeIt.GoToReverseBeginOfLine();
      while( !nyquistIt.IsAtEnd() )
        {
        while( !nyquistIt.IsAtEndOfLine() )
          {
          nyquistIt.Set(
            (positiveIt.Get()
             + std::conj(negativeIt.Get())
            ) * static_cast<typename ImageType::PixelType::value_type>(0.5) );
          ++nyquistIt; ++positiveIt;
          --negativeIt;
          }
        nyquistIt.NextLine(); positiveIt.NextLine();
        negativeIt.NextLine();
        negativeIt.GoToReverseBeginOfLine();
        }
      }
    }
}

template<class TImageType>
void
FrequencyShrinkImageFilter<TImageType>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  ImagePointer inputPtr  = const_cast<TImageType *>(this->GetInput() );
  ImagePointer outputPtr = this->GetOutput();

  // The filter chops high frequencys [0 1...H,H-1 H-2...1].
  // We need the whole input image, indepently of the RequestedRegion.
  inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template<class TImageType>
void
FrequencyShrinkImageFilter<TImageType>
::GenerateOutputInformation()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // Get pointers to the input and output
  ImageConstPointer inputPtr  = this->GetInput();
  ImagePointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // Compute the output spacing, the output image size, and the
  // output image start index
  const typename TImageType::SpacingType & inputSpacing =
    inputPtr->GetSpacing();
  const typename TImageType::SizeType & inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TImageType::IndexType & inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();
  const typename TImageType::PointType & inputOrigin =
    inputPtr->GetOrigin();

  // ContinuousIndex<double,ImageDimension> inputIndexOutputOrigin;

  typename TImageType::SpacingType outputSpacing(inputSpacing);
  typename TImageType::SizeType    outputSize;
  typename TImageType::PointType   outputOrigin;
  typename TImageType::IndexType   outputStartIndex;

  // TODO Check if you want to modify metadata in this filter.
  // Reduce Spacing, a frequency shrinker deletes high frequency domain.
  // The spacing is taken into account by FrequencyIterators method GetFrequency().
  for ( unsigned int i = 0; i < TImageType::ImageDimension; i++ )
    {
    outputSpacing[i] = inputSpacing[i] * m_ShrinkFactors[i];
    // inputIndexOutputOrigin[i] = 0.5*(m_ShrinkFactors[i]-1);
    // outputStartIndex[i] =
    //   Math::Ceil<SizeValueType>(inputStartIndex[i]/static_cast<double>( m_ShrinkFactors[i]) );
    // outputSize[i] = Math::Floor<SizeValueType>(
    //     static_cast<double>( inputSize[i] -
    //       outputStartIndex[i]*m_ShrinkFactors[i]+inputStartIndex[i])
    //     / static_cast<double>(m_ShrinkFactors[i])
    //     );
    outputStartIndex[i] = inputStartIndex[i];
    outputSize[i] = Math::Floor<SizeValueType>(
        static_cast<double>( inputSize[i] )
        / static_cast<double>(m_ShrinkFactors[i])
        );

    if ( outputSize[i] < 1 )
      {
      itkExceptionMacro("InputImage is too small! An output pixel does not map to a whole input bin.");
      }
    }

  // inputPtr->TransformContinuousIndexToPhysicalPoint(inputIndexOutputOrigin, outputOrigin);
  outputOrigin = inputOrigin;

  outputPtr->SetSpacing(outputSpacing);
  outputPtr->SetOrigin(outputOrigin);

  // Set region
  typename TImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
}

template<class TImageType>
void
FrequencyShrinkImageFilter< TImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Shrink Factor: ";
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    os << m_ShrinkFactors[j] << " ";
    }
  os << std::endl;
  os << "ApplyBandFilter: " << this->m_ApplyBandFilter << std::endl;

  itkPrintSelfObjectMacro(FrequencyBandFilter);

}
} // end namespace itk

#endif
