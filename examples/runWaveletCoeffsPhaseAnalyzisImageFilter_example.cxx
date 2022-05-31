/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumberToString.h"
#include <string>

#include "itkWaveletCoeffsPhaseAnalyzisImageFilter.h"

std::string
AppendToFilenameRiesz(const std::string & filename, const std::string & appendix)
{
  std::size_t foundDot = filename.find_last_of('.');
  return filename.substr(0, foundDot) + appendix + filename.substr(foundDot);
}

// 1. Wavelet analysis (forward) on input image.
// 2. Create a Monogenic Signal (from Riesz function ) on each wavelet output..
// 3. Do a PhaseAnalysis on each Monogenic Signal.
// 4. Wavelet reconstruction (inverse) using as coefficients the output of the PhaseAnalysis.
// Without applying reconstruction factors: ApplyReconstructionFactorOff()
// 5. The result of the reconstruction will be an image that uses phase information at each level/band for improving
// local structure information, and can also work as an equalizator of brightness.
template <unsigned int VDimension, typename TWaveletFunction>
int
runWaveletCoeffsPhaseAnalysis(const std::string &  inputImage,
                              const std::string &  outputImage,
                              const unsigned int & inputLevels,
                              const unsigned int & inputBands,
                              const bool           applySoftThreshold,
                              const double         thresholdNumOfSigmas = 2.0)
{
  const unsigned int Dimension = VDimension;

  using PixelType = float;
  using ImageType = itk::Image<PixelType, Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;

  typename ReaderType ::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImage);

  using FilterType = itk::WaveletCoeffsPhaseAnalyzisImageFilter<ImageType, TWaveletFunction>;
  typename FilterType ::Pointer filter = FilterType::New();
  filter->SetInput(reader->GetOutput());

  filter->SetLevels(inputLevels);
  filter->SetHighPassSubBands(inputBands);

  filter->SetApplySoftThreshold(applySoftThreshold);
  filter->SetThresholdNumOfSigmas(thresholdNumOfSigmas);

  itk::NumberToString<unsigned int> n2s;
  using WriterType = itk::ImageFileWriter<ImageType>;
  typename WriterType ::Pointer writer = WriterType::New();
  writer->SetInput(filter->GetOutput());
  std::string appendString = "_L" + n2s(inputLevels) + "_B" + n2s(inputBands) + "_S" + n2s(thresholdNumOfSigmas);
  std::string outputFile = AppendToFilenameRiesz(outputImage, appendString);
  writer->SetFileName(outputFile);
  writer->UseCompressionOn();

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Error : " << e << std::endl;
  }

  return EXIT_SUCCESS;
}

template <typename WaveletScalarType, unsigned int VDimension>
int
runWithChosenWavelet(const std::string &  waveletFunction,
                     const std::string &  inputImage,
                     const std::string &  outputImage,
                     const unsigned int & inputLevels,
                     const unsigned int & inputBands,
                     const bool           applySoftThreshold,
                     const double         thresholdNumOfSigmas)
{
  constexpr unsigned int ImageDimension = VDimension;
  if (waveletFunction == "Held")
  {
    using WaveletType = itk::HeldIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsPhaseAnalysis<ImageDimension, WaveletType>(
      inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else if (waveletFunction == "Vow")
  {
    using WaveletType = itk::VowIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsPhaseAnalysis<ImageDimension, WaveletType>(
      inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else if (waveletFunction == "Simoncelli")
  {
    using WaveletType = itk::SimoncelliIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsPhaseAnalysis<ImageDimension, WaveletType>(
      inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else if (waveletFunction == "Shannon")
  {
    using WaveletType = itk::ShannonIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsPhaseAnalysis<ImageDimension, WaveletType>(
      inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else
  {
    std::cerr << " failed!" << std::endl;
    std::cerr << waveletFunction << " wavelet type not supported." << std::endl;
    return EXIT_FAILURE;
  }
}

int
main(int argc, char * argv[])
{
  if (argc < 8 || argc > 9)
  {
    std::cerr << "Usage : " << std::endl;
    std::cerr << argv[0]
              << " inputImageFile outputImageFile inputLevels inputBands waveletFunction dimension Apply|NoApply "
                 "[thresholdNumOfSigmas]"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string  inputImage = argv[1];
  const std::string  outputImage = argv[2];
  const unsigned int inputLevels = atoi(argv[3]);
  const unsigned int inputBands = atoi(argv[4]);
  const std::string  waveletFunction = argv[5];
  const unsigned int dimension = atoi(argv[6]);
  const std::string  applySoftThresholdInput = argv[7];
  bool               applySoftThreshold = false;
  if (applySoftThresholdInput == "Apply")
  {
    applySoftThreshold = true;
  }
  else if (applySoftThresholdInput == "NoApply")
  {
    applySoftThreshold = false;
  }
  else
  {
    std::cerr << "Unkown string: " + applySoftThresholdInput + " . Use Apply or NoApply." << std::endl;
    return EXIT_FAILURE;
  }
  double thresholdNumOfSigmas = 2.0;
  if (argc == 9)
  {
    thresholdNumOfSigmas = atof(argv[8]);
  }


  using WaveletScalarType = double;
  if (dimension == 2)
  {
    constexpr unsigned int ImageDimension = 2;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(
      waveletFunction, inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else if (dimension == 3)
  {
    constexpr unsigned int ImageDimension = 3;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(
      waveletFunction, inputImage, outputImage, inputLevels, inputBands, applySoftThreshold, thresholdNumOfSigmas);
  }
  else
  {
    std::cerr << "Failed!" << std::endl;
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
  }
}
