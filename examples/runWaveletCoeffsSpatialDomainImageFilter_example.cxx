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

#include "itkWaveletCoeffsSpatialDomainImageFilter.h"

std::string
AppendToFilenameRiesz(const std::string & filename, const std::string & appendix)
{
  std::size_t foundDot = filename.find_last_of('.');
  return filename.substr(0, foundDot) + appendix + filename.substr(foundDot);
}

// 1. Wavelet analysis (forward) on input image.
// 2. Save those each wavelet coefficients on disk after converting them into spatial domain.
template <unsigned int VDimension, typename TWaveletFunction>
int
runWaveletCoeffsSpatialDomain(const std::string &  inputImage,
                              const std::string &  outputImage,
                              const unsigned int & inputLevels,
                              const unsigned int & inputBands)
{
  const unsigned int Dimension = VDimension;

  using PixelType = float;
  using ImageType = itk::Image<PixelType, Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;

  typename ReaderType ::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImage);


  using FilterType = itk::WaveletCoeffsSpatialDomainImageFilter<ImageType, TWaveletFunction>;
  typename FilterType ::Pointer filter = FilterType::New();
  filter->SetInput(reader->GetOutput());

  filter->SetLevels(inputLevels);
  filter->SetHighPassSubBands(inputBands);
  // filter->Print(std::cout); // Debugging line

  unsigned int TotalOutputs = inputLevels * inputBands + 1;
  std::cout << "Number Of Outputs Expected Total: " << TotalOutputs << std::endl;

  unsigned int computedNumberOfOutputs = filter->GetNumberOfOutputs();
  std::cout << "Number Of Outputs Calculated by Constructor : " << computedNumberOfOutputs << std::endl;

  itk::NumberToString<unsigned int> n2s;
  using WriterType = itk::ImageFileWriter<ImageType>;
  typename WriterType ::Pointer writer = WriterType::New();

  filter->Update();

  std::cout << "Number Of Outputs After Filter Built Up : " << filter->GetNumberOfOutputs() << std::endl;

  for (unsigned int i = 0; i < TotalOutputs; ++i)
  {
    writer->SetInput(filter->GetOutput(i));

    std::string appendString = "_L" + n2s(inputLevels) + "_B" + n2s(inputBands) + "_i" + n2s(i);
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
    writer->ResetPipeline();

    // Print the output volume size while retrieving the filter result.
    typename ImageType::Pointer outputWavelets;
    outputWavelets = filter->GetOutput(i);
    outputWavelets->Update();
    typename ImageType::SizeType sizeOfImage = outputWavelets->GetLargestPossibleRegion().GetSize();
    std::cout << "Wavelet Coefficients for index : " << i << " is saved ! " << std::endl;
    std::cout << "Wavelet Coefficients for index : " << i << " Size     : " << sizeOfImage << std::endl;
  }

  return EXIT_SUCCESS;
}

template <typename WaveletScalarType, unsigned int VDimension>
int
runWithChosenWavelet(const std::string &  waveletFunction,
                     const std::string &  inputImage,
                     const std::string &  outputImage,
                     const unsigned int & inputLevels,
                     const unsigned int & inputBands)
{
  constexpr unsigned int ImageDimension = VDimension;
  if (waveletFunction == "Held")
  {
    using WaveletType = itk::HeldIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsSpatialDomain<ImageDimension, WaveletType>(inputImage, outputImage, inputLevels, inputBands);
  }
  else if (waveletFunction == "Vow")
  {
    using WaveletType = itk::VowIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsSpatialDomain<ImageDimension, WaveletType>(inputImage, outputImage, inputLevels, inputBands);
  }
  else if (waveletFunction == "Simoncelli")
  {
    using WaveletType = itk::SimoncelliIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsSpatialDomain<ImageDimension, WaveletType>(inputImage, outputImage, inputLevels, inputBands);
  }
  else if (waveletFunction == "Shannon")
  {
    using WaveletType = itk::ShannonIsotropicWavelet<WaveletScalarType, ImageDimension>;
    return runWaveletCoeffsSpatialDomain<ImageDimension, WaveletType>(inputImage, outputImage, inputLevels, inputBands);
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
  if (argc != 7)
  {
    std::cerr << "Usage : " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile inputLevels inputBands waveletFunction dimension"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string  inputImage = argv[1];
  const std::string  outputImage = argv[2];
  const unsigned int inputLevels = atoi(argv[3]);
  const unsigned int inputBands = atoi(argv[4]);
  const std::string  waveletFunction = argv[5];
  const unsigned int dimension = atoi(argv[6]);
  using WaveletScalarType = double;
  if (dimension == 2)
  {
    constexpr unsigned int ImageDimension = 2;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(
      waveletFunction, inputImage, outputImage, inputLevels, inputBands);
  }
  else if (dimension == 3)
  {
    constexpr unsigned int ImageDimension = 3;
    return runWithChosenWavelet<WaveletScalarType, ImageDimension>(
      waveletFunction, inputImage, outputImage, inputLevels, inputBands);
  }
  else
  {
    std::cerr << "Failed!" << std::endl;
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
  }
}
