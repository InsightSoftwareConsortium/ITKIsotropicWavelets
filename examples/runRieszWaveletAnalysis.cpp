/*=========================================================================
 *
 *  Copyright Pablo Hernandez-Cerdan
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

#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyInverse.h"
#include "itkWaveletFrequencyForwardUndecimated.h"
#include "itkWaveletFrequencyInverseUndecimated.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkHeldIsotropicWavelet.h"
#include "itkVowIsotropicWavelet.h"
#include "itkSimoncelliIsotropicWavelet.h"
#include "itkShannonIsotropicWavelet.h"

#include "itkMonogenicSignalFrequencyImageFilter.h"
#include "itkVectorInverseFFTImageFilter.h"
#include "itkPhaseAnalysisSoftThresholdImageFilter.h"
#include "itkZeroDCImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkNumberToString.h"
#include <string>

std::string
AppendToFilenameRiesz(const std::string& filename, const std::string & appendix)
{
  std::size_t foundDot = filename.find_last_of('.');
  return filename.substr( 0, foundDot ) + appendix + filename.substr( foundDot );
}

// 1. Wavelet analysis (forward) on input image.
// 2. Create a Monogenic Signal (from Riesz function ) on each wavelet output..
// 3. Do a PhaseAnalysis on each Monogenic Signal.
// 4. Wavelet reconstruction (inverse) using as coefficients the output of the PhaseAnalysis.
// Without applying reconstruction factors: ApplyReconstructionFactorOff()
// 5. The result of the reconstruction will be an image that uses phase information at each level/band for improving local structure information, and can also work as an equalizator of brightness.
template< unsigned int VDimension, typename TWaveletFunction >
int
runRieszWaveletPhaseAnalysis( const std::string& inputImage,
  const std::string & outputImage,
  const unsigned int& inputLevels,
  const unsigned int& inputBands,
  const bool applySoftThreshold,
  const double thresholdNumOfSigmas = 2.0)
{
  const unsigned int Dimension = VDimension;

  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;
  using ReaderType = itk::ImageFileReader< ImageType >;

  itk::NumberToString< unsigned int > n2s;
  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  reader->Update();

  using BoundaryConditionType = itk::ConstantBoundaryCondition<ImageType>;
  BoundaryConditionType bounds;
  bounds.SetConstant(itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
 
  using FFTPadFilterType = itk::FFTPadImageFilter<ImageType>;
  auto fftPadFilter = FFTPadFilterType::New();
  fftPadFilter->SetInput(reader->GetOutput());
  fftPadFilter->SetBoundaryCondition(&bounds);
  fftPadFilter->Update();

  auto sizeOriginal = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  auto sizeAfterPad = fftPadFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout << "Image Padded" << std::endl;
  std::cout << "Original Size:" << sizeOriginal << std::endl;
  std::cout << "After Pad Size:" << sizeAfterPad << std::endl;

  unsigned int scaleFactor = 2;
  unsigned int maxNumberOfLevels = itk::utils::ComputeMaxNumberOfLevels(sizeAfterPad, scaleFactor);
  std::cout << "MaxNumberOfLevels allowed for the padded size: " << maxNumberOfLevels <<
    ". Recommended: " << static_cast<int>(maxNumberOfLevels) - 3  << std::endl;

  using ZeroDCFilterType = itk::ZeroDCImageFilter< ImageType >;
  auto zeroDCFilter = ZeroDCFilterType::New();
  zeroDCFilter->SetInput( fftPadFilter->GetOutput() );
  zeroDCFilter->Update();

  // Perform FFT on input image.
  using FFTForwardFilterType = itk::ForwardFFTImageFilter< typename ZeroDCFilterType::OutputImageType >;
  auto fftForwardFilter = FFTForwardFilterType::New();
  fftForwardFilter->SetInput( zeroDCFilter->GetOutput() );
  fftForwardFilter->Update();
  using ComplexImageType = typename FFTForwardFilterType::OutputImageType;

  using InverseFFTFilterType = itk::InverseFFTImageFilter< ComplexImageType, ImageType >;

  // Forward Wavelet
  using WaveletFunctionType = TWaveletFunction;
  using WaveletFilterBankType = itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType >;
  using ForwardWaveletType = itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  // using ForwardWaveletType = itk::WaveletFrequencyForwardUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  auto forwardWavelet = ForwardWaveletType::New();
  unsigned int highSubBands = inputBands;
  unsigned int levels = inputLevels;
  forwardWavelet->SetHighPassSubBands( highSubBands );
  forwardWavelet->SetLevels( levels );
  forwardWavelet->SetInput( fftForwardFilter->GetOutput() );
  forwardWavelet->Update();
  typename ForwardWaveletType::OutputsType analysisWavelets =
    forwardWavelet->GetOutputs();

  // Apply Monogenic signal to wavelet results
  using MonogenicSignalFrequencyFilterType = itk::MonogenicSignalFrequencyImageFilter< ComplexImageType >;
  using VectorMonoOutputType = typename MonogenicSignalFrequencyFilterType::OutputImageType;
  using VectorInverseFFTType = itk::VectorInverseFFTImageFilter< VectorMonoOutputType >;
  using PhaseAnalysisFilter = itk::PhaseAnalysisSoftThresholdImageFilter< typename VectorInverseFFTType::OutputImageType >;

  typename ForwardWaveletType::OutputsType modifiedWavelets;
  unsigned int numberOfOutputs = forwardWavelet->GetNumberOfOutputs();
  for ( unsigned int i = 0; i < forwardWavelet->GetNumberOfOutputs(); ++i )
    {
    std::cout << "Output #: " << i << " / " << numberOfOutputs - 1 << std::endl;
    auto monoFilter = MonogenicSignalFrequencyFilterType::New();
    auto vecInverseFFT = VectorInverseFFTType::New();
    auto phaseAnalyzer = PhaseAnalysisFilter::New();
    auto fftForwardPhaseFilter = FFTForwardFilterType::New();

    // Generate a monogenic signal (vector valued)
    monoFilter->SetInput( analysisWavelets[i] );
    monoFilter->Update();

    vecInverseFFT->SetInput( monoFilter->GetOutput() );
    vecInverseFFT->Update();

    phaseAnalyzer->SetInput( vecInverseFFT->GetOutput() );
    phaseAnalyzer->SetApplySoftThreshold( applySoftThreshold );
    if (applySoftThreshold)
      {
      phaseAnalyzer->SetNumOfSigmas(thresholdNumOfSigmas);
      }
    phaseAnalyzer->Update();

    fftForwardPhaseFilter->SetInput( phaseAnalyzer->GetOutputCosPhase() );
    fftForwardPhaseFilter->Update();

    modifiedWavelets.push_back( fftForwardPhaseFilter->GetOutput() );
    modifiedWavelets.back()->DisconnectPipeline();
    }

  using InverseWaveletType = itk::WaveletFrequencyInverse< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  // using InverseWaveletType = itk::WaveletFrequencyInverseUndecimated< ComplexImageType, ComplexImageType, WaveletFilterBankType >;
  auto inverseWavelet = InverseWaveletType::New();
  inverseWavelet->SetHighPassSubBands( highSubBands );
  inverseWavelet->SetLevels( levels );
  inverseWavelet->SetInputs( modifiedWavelets );
  // The coefficients are now phases, do not apply reconstruction factors.
  inverseWavelet->ApplyReconstructionFactorsOff();
  inverseWavelet->Update();

  auto inverseFFT = InverseFFTFilterType::New();
  inverseFFT->SetInput( inverseWavelet->GetOutput() );
  inverseFFT->Update();

  using ChangeInformationFilterType = itk::ChangeInformationImageFilter< ImageType >;
  auto changeInfoFilter = ChangeInformationFilterType::New();
  changeInfoFilter->SetInput( inverseFFT->GetOutput());
  changeInfoFilter->SetReferenceImage(fftPadFilter->GetOutput());
  changeInfoFilter->UseReferenceImageOn();
  changeInfoFilter->ChangeAll();
  changeInfoFilter->Update();

  // Cast To Float for save as tiff.
  using ImageFloatType = itk::Image< float, Dimension >;
  using CastFloatType = itk::CastImageFilter< ImageType, ImageFloatType>;
  auto caster = CastFloatType::New();
  caster->SetInput(changeInfoFilter->GetOutput());
  caster->Update();

  // using WriterType = itk::ImageFileWriter< typename InverseFFTFilterType::OutputImageType >;
  using WriterType = itk::ImageFileWriter< ImageFloatType >;
  auto writer = WriterType::New();
  std::string appendString = "_L" + n2s(inputLevels) + "_B" + n2s(inputBands) + "_S" + n2s(thresholdNumOfSigmas);
  std::string outputFile = AppendToFilenameRiesz(outputImage, appendString);
  writer->SetFileName( outputFile );
  writer->SetInput( caster->GetOutput() );

  writer->Update();
  std::cout << "Output generated in: " << outputFile << std::endl;
  return EXIT_SUCCESS;
}

template< typename WaveletScalarType, unsigned int VDimension >
int runWithChosenWavelet(const std::string & waveletFunction,
  const std::string& inputImage,
  const std::string & outputImage,
  const unsigned int& inputLevels,
  const unsigned int& inputBands,
  const bool applySoftThreshold,
  const double thresholdNumOfSigmas)
{
  constexpr unsigned int ImageDimension = VDimension;
  if ( waveletFunction == "Held" )
    {
    using WaveletType = itk::HeldIsotropicWavelet< WaveletScalarType, ImageDimension>;
    return runRieszWaveletPhaseAnalysis< ImageDimension, WaveletType >(
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else if ( waveletFunction == "Vow" )
    {
    using WaveletType = itk::VowIsotropicWavelet< WaveletScalarType, ImageDimension>;
    return runRieszWaveletPhaseAnalysis< ImageDimension, WaveletType >(
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else if ( waveletFunction == "Simoncelli" )
    {
    using WaveletType = itk::SimoncelliIsotropicWavelet< WaveletScalarType, ImageDimension>;
    return runRieszWaveletPhaseAnalysis< ImageDimension, WaveletType >(
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else if ( waveletFunction == "Shannon" )
    {
    using WaveletType = itk::ShannonIsotropicWavelet< WaveletScalarType, ImageDimension>;
    return runRieszWaveletPhaseAnalysis< ImageDimension, WaveletType >(
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else
    {
    std::cerr << " failed!" << std::endl;
    std::cerr << waveletFunction << " wavelet type not supported." << std::endl;
    return EXIT_FAILURE;
    }
}

int main( int argc, char *argv[] )
{
  if ( argc < 8 || argc > 9 )
    {
    std::cerr << "Usage: " << argv[0]
              << " inputImage outputImage inputLevels inputBands waveletFunction dimension Apply|NoApply [thresholdNumOfSigmas]"
              << std::endl;
    return EXIT_FAILURE;
    }

  const std::string inputImage  = argv[1];
  const std::string outputImage = argv[2];
  const unsigned int inputLevels = atoi( argv[3] );
  const unsigned int inputBands  = atoi( argv[4] );
  const std::string waveletFunction = argv[5];
  const unsigned int dimension = atoi( argv[6] );
  const std::string applySoftThresholdInput = argv[7];
  bool applySoftThreshold = false;
  if ( applySoftThresholdInput == "Apply" )
    {
    applySoftThreshold = true;
    }
  else if ( applySoftThresholdInput == "NoApply" )
    {
    applySoftThreshold = false;
    }
  else
    {
    std::cerr << "Unkown string: " + applySoftThresholdInput + " . Use Apply or NoApply." << std::endl;
    return EXIT_FAILURE;
    }

  double thresholdNumOfSigmas = 2.0;
  if ( argc == 9 )
    {
    thresholdNumOfSigmas = atof(argv[8]);
    }


  using WaveletScalarType = double;
  if ( dimension == 2 )
    {
    constexpr unsigned int ImageDimension = 2;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(waveletFunction,
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else if ( dimension == 3 )
    {
    constexpr unsigned int ImageDimension = 3;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(waveletFunction,
      inputImage, outputImage, inputLevels, inputBands,
      applySoftThreshold, thresholdNumOfSigmas);
    }
  else
    {
    std::cerr << "Failed!" << std::endl;
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
    }
}
