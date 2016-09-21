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
#ifndef itkMonogenicPhaseAnalysisEigenValuesImageFilter_hxx
#define itkMonogenicPhaseAnalysisEigenValuesImageFilter_hxx
#include "itkMonogenicPhaseAnalysisEigenValuesImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include <numeric>
// Eigen Calculations
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkGaussianImageSource.h>
#include <itkMatrix.h>
#include <itkSymmetricEigenAnalysis.h>

#include "itkProgressReporter.h"
#include "itkStatisticsImageFilter.h"
#include "itkRieszFrequencyFunction.h"
namespace itk {
template< typename TInputImage, typename TOutputImage >
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::MonogenicPhaseAnalysisEigenValuesImageFilter()
: m_GaussianWindowRadius(2),
  m_GaussianWindowSigma(1.0),
  m_ApplySoftThreshold(true)
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(3);

  for (unsigned int n_output = 0; n_output < 3; ++n_output)
    {
    this->SetNthOutput(n_output, this->MakeOutput(n_output));
    }
}

template< typename TInputImage, typename TOutputImage >
void
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Threshold : " << m_Threshold << std::endl;
  os << indent << "Mean Amplitude : " << m_MeanAmp << std::endl;
  os << indent << "Sigma Amplitude: " << m_SigmaAmp << std::endl;
}

template< typename TInputImage, typename TOutputImage >
void
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  m_NC = this->GetInput()->GetNumberOfComponentsPerPixel();
  if (this->GetInput()->GetNumberOfComponentsPerPixel() != ImageDimension + 1 )
    {
    itkExceptionMacro(<<"Number of components of input image ("<<m_NC<<") is not ImageDimension+1 (" << ImageDimension + 1 << ")");
    }

  ThreadIdType nbOfThreads = this->GetNumberOfThreads();
  m_Barrier1 = Barrier::New();
  m_Barrier1->Initialize(nbOfThreads);
  m_Barrier2 = Barrier::New();
  m_Barrier2->Initialize(nbOfThreads);

  m_MeanAmp = 0;
  m_SigmaAmp= 0;
  m_Threshold = 0;

  m_EigenVectorsImage = EigenVectorsImageType::New();
  m_EigenValuesImage  = EigenValuesImageType::New();
  m_EigenVectorsImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_EigenVectorsImage->Allocate();
  m_EigenValuesImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_EigenValuesImage->Allocate();
}

template< typename TInputImage, typename TOutputImage >
void
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
        ThreadIdType threadId)
{
  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  typename OutputImageType::Pointer outputPtr = this->GetOutput();
  typename OutputImageType::Pointer amplitudePtr = this->GetOutput(1);
  typename OutputImageType::Pointer phasePtr = this->GetOutput(2);

  OutputImageRegionIterator outIt(outputPtr, outputRegionForThread);
  OutputImageRegionIterator ampIt(amplitudePtr, outputRegionForThread);
  OutputImageRegionIterator phaseIt(phasePtr, outputRegionForThread);

  InputImageRegionConstIterator monoIt(this->GetInput(), outputRegionForThread);
  //TODO call it with sensible paramenters
  this->ComputeEigenAnalysisInNeighborhoodWindow(
    this->m_GaussianWindowRadius, this->m_GaussianWindowSigma ,outputRegionForThread);
// Copied from  MaximumResponse function
//
  ImageRegionConstIterator<EigenVectorsImageType> eigenVectorsIt(
    m_EigenVectorsImage, outputRegionForThread);

  // monoProjection = eigenVector * R.
  itk::VariableLengthVector<OutputImagePixelType> monoProjection;
  monoProjection.SetSize(ImageDimension + 1);
  // OutputImagePixelType rieszProjection;
  DirectionType direction;
  for(monoIt.GoToBegin(), eigenVectorsIt.GoToBegin(),
    ampIt.GoToBegin(), phaseIt.GoToBegin();
    !monoIt.IsAtEnd();
    ++monoIt, ++eigenVectorsIt, ++ampIt, ++phaseIt )
    {
    InputImagePixelType monoValue = monoIt.Get();
    monoProjection[0] = monoValue[0];
    for(unsigned int column=0; column < ImageDimension; ++column )
      {
      for (unsigned int i = 0; i < ImageDimension; ++i)
        direction[i] = eigenVectorsIt.Get()[i][column];
      monoProjection[column + 1] = ComputeRieszProjection(monoValue, direction);
      }
    // TODO hack to get only the first component.
    // OutputImagePixelType projectionNorm = monoProjection[1]*monoProjection[1];
    OutputImagePixelType projectionNorm = this->ComputeRieszNormSquare(monoProjection);
    ampIt.Set(this->ComputeAmplitude(monoProjection, projectionNorm));
    phaseIt.Set(cos(this->ComputePhase(monoProjection, projectionNorm)));
    }

  // Wait for mean/variance calculation. Stats only once.
  m_Barrier1->Wait();
  if(threadId == this->GetNumberOfThreads() - 1)
    {
    typedef itk::StatisticsImageFilter<OutputImageType> StatisticsImageFilter;
    typename StatisticsImageFilter::Pointer statsFilter = StatisticsImageFilter::New();
    statsFilter->SetInput(amplitudePtr);
    statsFilter->Update();
    m_MeanAmp = statsFilter->GetMean();
    m_SigmaAmp = sqrt(statsFilter->GetVariance());
    // 2 sigmas, it could be user chosen.
    m_Threshold = m_MeanAmp + 2 * m_SigmaAmp;
    }
  m_Barrier2->Wait();

  for( outIt.GoToBegin(),ampIt.GoToBegin(), phaseIt.GoToBegin();
    !outIt.IsAtEnd();
    ++outIt, ++ampIt, ++phaseIt)
    {
    OutputImagePixelType out_value = cos(phaseIt.Get());

    if(this->GetApplySoftThreshold() == true)
      if (ampIt.Get() < m_Threshold)
        out_value = out_value * ampIt.Get() / m_Threshold;

    outIt.Set(out_value);
    progress.CompletedPixel();
    }
}

template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::OutputImagePixelType
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ComputeRieszNormSquare(const InputImagePixelType & monoPixel) const
{

  OutputImagePixelType out(0);
  for(unsigned int r = 1; r < this->GetNC(); r++)
    {
    out += monoPixel[r]*monoPixel[r];
    }
  return out;
}

template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::OutputImagePixelType
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ComputeAmplitude(const InputImagePixelType & monoPixel,
  const OutputImagePixelType & rieszNormSquare ) const
{
  return sqrt( monoPixel[0]*monoPixel[0] + rieszNormSquare );
}

template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::OutputImagePixelType
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ComputePhase(const InputImagePixelType & monoPixel,
  const OutputImagePixelType & rieszNormSquare ) const
{
  return atan2(sqrt(rieszNormSquare), monoPixel[0]);
}

template< typename TInputImage, typename TOutputImage >
itk::FixedArray<
  typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::OutputImagePixelType,
  MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::ImageDimension - 1 >
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ComputePhaseOrientation(const InputImagePixelType & monoPixel,
  const OutputImagePixelType & rieszNormSquare ) const
{
  // the angles of the polar coordinates of the normed vector:
  // V = (R1*f, ..., Rn*f) / RieszNorm
  itk::FixedArray< OutputImagePixelType, ImageDimension - 1> out(0);
  OutputImagePixelType rNorm = sqrt(rieszNormSquare);
  OutputImagePixelType r1Unitary = monoPixel[1] / rNorm;
  for(unsigned int i = 0; i < ImageDimension - 1; i++)
    {
    out[i] = atan2(monoPixel[i+2]/rNorm, r1Unitary) +
     ( (monoPixel[i+2] >= 0) ? 0 : itk::Math::pi  );
    }
  return out;
}

template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::OutputImagePixelType
MonogenicPhaseAnalysisEigenValuesImageFilter< TInputImage, TOutputImage >
::ComputeRieszProjection(
    const InputImagePixelType & monoPixel,
    const DirectionType & direction ) const
{
  OutputImagePixelType out(0);
  for (unsigned int r = 1; r < this->GetNC(); r++)
    out += direction[r] * monoPixel[r];
  return out;
}

/** For each pixel of eigenOut (size of TInput)
 * For each RieszComponent:
 * Use NeighborhoodIterator in RieszComponents using gaussian_radius
 * Weight value of pixels in the neighborhood using the GaussianImage
 *
 * Compute Matrix (2D,Dimension x Dimension) eigenMatrix,
 * which is the weighted contribution of the RieszComponents.
 * J(x_0)[m][n] =
 * Sum_each_neighbor_pixel_x(gaussian(x) * RieszComponent(x)[m] * RieszComponent(x)[n] )
 * Compute eigenVectors and eigenValues of the matrix.
 * Store them in the output.
 */
template< typename TInputImage, typename TOutputImage >
std::pair<
  typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::EigenVectorsImageType::Pointer,
  typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::EigenValuesImageType::Pointer >
MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
::ComputeEigenAnalysisInNeighborhoodWindow(
    const unsigned int & gaussian_window_radius,
    const FloatType & gaussian_window_sigma,
    const OutputImageRegionType & outputRegionForThread) const
{
  // Set GaussianImageSource
  Size<ImageDimension> radius;
  Size<ImageDimension> domainKernelSize;
  radius.Fill(gaussian_window_radius);
  domainKernelSize.Fill(2 * gaussian_window_radius + 1);
  typedef GaussianImageSource< OutputImageType > GaussianSourceType;
  typename GaussianSourceType::Pointer gaussianImage =
    GaussianSourceType::New();
  typename GaussianSourceType::ArrayType mean;
  typename GaussianSourceType::ArrayType sigma;

  const SpacingType inputSpacing = this->GetInput()->GetSpacing();
  const typename InputImageType::PointType inputOrigin = this->GetInput()->GetOrigin();
  gaussianImage->SetSize( domainKernelSize );
  gaussianImage->SetSpacing(inputSpacing);
  gaussianImage->SetOrigin(inputOrigin);
  gaussianImage->SetScale(1.0);
  gaussianImage->SetNormalized(true);

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    mean[i] = inputSpacing[i] * radius[i] + inputOrigin[i]; // center pixel pos
    sigma[i] = gaussian_window_sigma;
    }
  gaussianImage->SetSigma(sigma);
  gaussianImage->SetMean(mean);
  gaussianImage->Update();

  ImageRegionIterator<EigenVectorsImageType> eigenVectorsIt(
      m_EigenVectorsImage,
      outputRegionForThread);
  ImageRegionIterator<EigenValuesImageType> eigenValuesIt(
      m_EigenValuesImage,
      outputRegionForThread);
  ConstNeighborhoodIterator<InputImageType> monoIt(
      radius, this->GetInput(),
      outputRegionForThread);
  ConstNeighborhoodIterator<typename GaussianSourceType::OutputImageType>
    gaussianIt(
      radius, gaussianImage->GetOutput(),
      gaussianImage->GetOutput()->GetRequestedRegion());

  typedef itk::Matrix<FloatType, ImageDimension,ImageDimension>         EigenMatrixType;
  typedef itk::FixedArray<FloatType, ImageDimension>                    EigenValuesType;
  typedef itk::SymmetricEigenAnalysis<EigenMatrixType, EigenValuesType> SymmetricEigenAnalysisType;
  EigenMatrixType            eigenMatrix;
  EigenMatrixType            eigenVectors;
  EigenValuesType            eigenValues;
  SymmetricEigenAnalysisType eigenSystem(ImageDimension);

  for(monoIt.GoToBegin(), gaussianIt.GoToBegin(),
    eigenVectorsIt.GoToBegin(), eigenValuesIt.GoToBegin();
    !monoIt.IsAtEnd();
    ++monoIt, ++eigenVectorsIt, ++eigenValuesIt )
    {
    //Set the matrix
    eigenMatrix.Fill(0);
    for (unsigned int r = 0; r <= gaussian_window_radius; ++r)
      for (unsigned int m = 0; m < ImageDimension; ++m)
        for (unsigned int n = m; n < ImageDimension; ++n)
          if (r == 0)
            eigenMatrix[m][n] = eigenMatrix[n][m] +=
              gaussianIt.GetCenterPixel() +
              monoIt.GetCenterPixel()[m] + monoIt.GetCenterPixel()[n];
          else
            for (unsigned int axis = 0; axis < ImageDimension; ++axis)
              {
              eigenMatrix[m][n] = eigenMatrix[n][m] +=
                gaussianIt.GetNext(axis, r) +
                monoIt.GetNext(axis,r)[m] + monoIt.GetNext(axis,r)[n];
              eigenMatrix[m][n] +=
                gaussianIt.GetPrevious(axis, r) +
                monoIt.GetPrevious(r)[m] + monoIt.GetPrevious(r)[n];
              }
    eigenSystem.ComputeEigenValuesAndVectors(
      eigenMatrix, eigenValues, eigenVectors );
    // Copy to Output
    eigenVectorsIt.Set(eigenVectors);
    eigenValuesIt.Set(eigenValues);

    } // end NeighborhoodIterator monoIt

  return std::make_pair(m_EigenVectorsImage, m_EigenValuesImage);
}

// template< typename TInputImage, typename TOutputImage >
// typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::RieszComponentsImageType::Pointer
// MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
// ::ComputeRieszComponentsWithMaximumResponse(
//     const typename EigenVectorsImageType::Pointer eigenVectors,
//     const RieszComponentsImageType* rieszComponents ) const
// {
//   typename RieszComponentsImageType::Pointer maxLocalRiesz =
//     RieszComponentsImageType::New();
//   maxLocalRiesz->SetRegions(rieszComponents->GetLargestPossibleRegion());
//   maxLocalRiesz->SetVectorLength(ImageDimension);
//   maxLocalRiesz->Allocate();
//
//   ImageRegionIterator<RieszComponentsImageType> maxLocalRieszIt(
//     maxLocalRiesz, maxLocalRiesz->GetRequestedRegion());
//   ImageRegionConstIterator<RieszComponentsImageType> rieszIt(
//     rieszComponents, rieszComponents->GetRequestedRegion());
//   ImageRegionConstIterator<EigenVectorsImageType> eigenVectorsIt(
//     eigenVectors, eigenVectors->GetRequestedRegion());
//
//   itk::VariableLengthVector<InputImagePixelType> rieszProjection;
//   rieszProjection.SetSize(ImageDimension);
//   DirectionType direction;
//   for(rieszIt.GoToBegin(), eigenVectorsIt.GoToBegin(),
//     maxLocalRieszIt.GoToBegin() ;
//     !rieszIt.IsAtEnd() ;
//     ++rieszIt, ++eigenVectorsIt, ++maxLocalRieszIt )
//     {
//       for (unsigned int column = 0 ; column < ImageDimension ; ++column){
//         for (unsigned int i = 0 ; i < ImageDimension ; ++i)
//           direction[i] = eigenVectorsIt.Get()[i][column];
//         rieszProjection[column] = ComputeLocalRieszProjection(direction, rieszIt);
//       }
//       maxLocalRieszIt.Set( rieszProjection );
//     }
//
//   return maxLocalRiesz;
// }

/************************************
template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::InputImageType::Pointer
MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
::ComputeRieszWeightedNormByEigenValues(
    const RieszComponentsImageType* rieszComponents,
    const typename EigenValuesImageType::Pointer eigenValues) const
{
  typename InputImageType::Pointer normOutput = InputImageType::New();
  normOutput->SetRegions(rieszComponents->GetBufferedRegion());
  normOutput->Allocate();
  normOutput->FillBuffer(NumericTraits<InputImagePixelType>::Zero);

  ImageRegionIterator<InputImageType> normOutputIt(
    normOutput, normOutput->GetRequestedRegion());
  ImageRegionConstIterator<RieszComponentsImageType> rieszIt(
    rieszComponents, rieszComponents->GetRequestedRegion());
  ImageRegionConstIterator<EigenValuesImageType> eigenValuesIt(
    eigenValues, eigenValues->GetRequestedRegion());

  InputImagePixelType localValue;
  itk::VariableLengthVector<InputImagePixelType> localEigenValues;
  localEigenValues.SetSize(ImageDimension);
  DirectionType weights;
  typename DirectionType::ValueType eigen_values_sum ;
  for(rieszIt.GoToBegin(), eigenValuesIt.GoToBegin(),
    normOutputIt.GoToBegin() ;
    !rieszIt.IsAtEnd() ;
    ++rieszIt, ++eigenValuesIt, ++normOutputIt )
    {
      eigen_values_sum = 0;
      localValue = 0;
      // Get the sum of eigenValues first
      for (unsigned int i = 0 ; i < ImageDimension ; ++i)
      {
        localEigenValues[i] = eigenValuesIt.Get()[i];
        eigen_values_sum += localEigenValues[i];
      }
      // std::cout << "Sum eigenValues: " << eigen_values_sum << std::endl;
      // Compute the weighted norm using eigenValues
      for (unsigned int i = 0 ; i < ImageDimension ; ++i)
      {
          weights[i] = localEigenValues[i] / eigen_values_sum;
          localValue += localValue + rieszIt.Get()[i] * rieszIt.Get()[i] *
            weights[i] * weights[i];
          // std::cout << "Weight[" << i << "]: " << weights[i] <<
          //   "localValue (temp)" << localValue << std::endl;
      }

      normOutputIt.Set( sqrt(localValue) );
    }

  return normOutput;
}


template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
::ComputeLocalPhaseInDirection(
    const DirectionType & unitary_direction,
    const InputImageType* rieszReal,
    const RieszComponentsImageType* rieszComponents) const
{
  typename InputImageType::Pointer rieszProjection =
    ComputeRieszProjection(unitary_direction, rieszComponents);
  rieszProjection->DisconnectPipeline();
  typedef itk::DivideImageFilter<TOutputImage,TOutputImage, TOutputImage>
    DivideFilterType;
  typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
  divideFilter->SetInput1(rieszProjection);
  divideFilter->SetInput2(rieszReal);
  divideFilter->Update();

  typedef itk::AtanImageFilter<TOutputImage,TOutputImage> AtanFilterType;
  typename AtanFilterType::Pointer atanFilter = AtanFilterType::New();
  atanFilter->SetInput(divideFilter->GetOutput());
  atanFilter->Update();
  typename InputImageType::Pointer output = atanFilter->GetOutput();
  output->DisconnectPipeline();
  return output;
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
::ComputeTotalRieszProjection(const DirectionType & unitary_direction) const
    // const RieszComponentsImageType* rieszComponents) const
{
  OutputImageConstPointer rieszComponents = this->GetOutput();

  typedef itk::ImageDuplicator< InputImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator= DuplicatorType::New();
  duplicator->SetInputImage(this->GetInput());
  duplicator->Update();
  typename InputImageType::Pointer output = duplicator->GetModifiableOutput();
  OutputImageRegionConstIterator outIt(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

  for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
    {
    // TODO AQUI TAS QUEDAO BRO

    }
  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,
          InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);

  typedef itk::MultiplyImageFilter<TOutputImage,TOutputImage,TOutputImage>
    MultiplyImageFilterType;
  typename MultiplyImageFilterType::Pointer multiplyByConstant =
    MultiplyImageFilterType::New();

  typedef AddImageFilter<InputImageType,InputImageType,InputImageType>
    AddImageFilterType;
  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  for (unsigned int i = 0 ; i<ImageDimension ; ++i){
    castIndex->SetIndex(i);
    castIndex->Update();
    multiplyByConstant->SetInput(castIndex->GetOutput());
    multiplyByConstant->SetConstant(unitary_direction[i]);
    multiplyByConstant->Update();

    addFilter->SetInput1(output);
    addFilter->SetInput2(multiplyByConstant->GetOutput());
    addFilter->Update();
    output = addFilter->GetOutput();
    output->DisconnectPipeline();
  }
  return output;
}

template< typename TInputImage, typename TOutputImage >
typename MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>::InputImageType::Pointer
MonogenicPhaseAnalysisEigenValuesImageFilter<TInputImage, TOutputImage>
::ComputeRieszWeightedNorm(
    const RieszComponentsImageType* rieszComponents,
    const DirectionType & weights) const
{
  typename InputImageType::Pointer output = InputImageType::New();
  output->SetRegions(rieszComponents->GetBufferedRegion());
  output->Allocate();
  output->FillBuffer(NumericTraits<InputImagePixelType>::Zero);

  typedef VectorIndexSelectionCastImageFilter<RieszComponentsImageType,
          InputImageType> CastIndexType;
  typename CastIndexType::Pointer castIndex = CastIndexType::New();
  castIndex->SetInput(rieszComponents);

  typedef SquareImageFilter<InputImageType,InputImageType>
    SquareImageFilterType;
  typename SquareImageFilterType::Pointer squareFilter =
    SquareImageFilterType::New();

  typedef AddImageFilter<InputImageType,InputImageType,InputImageType>
    AddImageFilterType;
  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  typedef itk::MultiplyImageFilter<TOutputImage,TOutputImage,TOutputImage> MultiplyImageFilterType;
  typename MultiplyImageFilterType::Pointer multiplyByConstant = MultiplyImageFilterType::New();
  for (unsigned int i = 0 ; i<ImageDimension ; ++i){
    castIndex->SetIndex(i);
    castIndex->Update();
    multiplyByConstant->SetInput(castIndex->GetOutput());
    multiplyByConstant->SetConstant(weights[i]);
    multiplyByConstant->Update();
    squareFilter->SetInput(multiplyByConstant->GetOutput());
    squareFilter->Update();
    typename InputImageType::Pointer squareResult = squareFilter->GetOutput();
    squareResult->DisconnectPipeline();
    addFilter->SetInput1(squareResult);
    addFilter->SetInput2(output);
    addFilter->Update();
    output = addFilter->GetOutput();
    output->DisconnectPipeline();
  }
  typedef SqrtImageFilter<InputImageType,InputImageType> SqrtImageFilterType;
  typename SqrtImageFilterType::Pointer sqrtFilter = SqrtImageFilterType::New();
  sqrtFilter->SetInput(output);
  sqrtFilter->InPlaceOn();
  sqrtFilter->Update();
  return sqrtFilter->GetOutput();
}
**********************/

} // end namespace itk
#endif
