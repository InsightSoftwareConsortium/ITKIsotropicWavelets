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
#include <string>
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyInverse.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkHeldIsotropicWavelet.h"
#include "itkVowIsotropicWavelet.h"

#include "itkMonogenicSignalFrequencyImageFilter.h"
#include "itkVectorInverseFFTImageFilter.h"
#include "itkMonogenicPhaseAnalysisEigenValuesImageFilter.h"
#include "itkMonogenicPhaseAnalysisSoftThresholdImageFilter.h"
#include "itkZeroDCImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumberToString.h"

//Visualize for dev/debug purposes. Set in cmake file. Require VTK
#if ITK_VISUALIZE_TESTS != 0
#include "itkViewImage.h"
#endif
using namespace std;
using namespace itk;

template <unsigned int N>
int runRieszWaveletPhaseAnalysisTest( const std::string& inputImage,
    const std::string& outputImage,
    const unsigned int& inputLevels,
    const unsigned int& inputBands)
{
  const unsigned int dimension = N;
  typedef float                            PixelType;
  typedef itk::Image<PixelType, dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImage);
  reader->Update();
  reader->UpdateLargestPossibleRegion();

  itk::NumberToString<float> n2s;
  // Wavelet analysis (forward)
  // Monogenic Analysis on each output. (In-place style)
  // Wavelet reconstruction (inverse) from PhaseAnalysis output
  // Profit.

  // TODO check  differences applying zeroDC image filter.

  typedef itk::ZeroDCImageFilter<ImageType> ZeroDCFilterType;
  typename ZeroDCFilterType::Pointer zeroDCFilter = ZeroDCFilterType::New();
  zeroDCFilter->SetInput(reader->GetOutput());
  zeroDCFilter->Update();

  // Perform FFT on input image.
  typedef itk::ForwardFFTImageFilter<ImageType> FFTForwardFilterType;
  typename FFTForwardFilterType::Pointer fftForwardFilter = FFTForwardFilterType::New();
  fftForwardFilter->SetInput(zeroDCFilter->GetOutput());
  fftForwardFilter->Update();
  typedef typename FFTForwardFilterType::OutputImageType ComplexImageType;

  // Forward Wavelet
  typedef itk::HeldIsotropicWavelet<PixelType> WaveletFunctionType;
  // typedef itk::VowIsotropicWavelet<> WaveletFunctionType;
  typedef itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType> WaveletFilterBankType;
  typedef itk::WaveletFrequencyForward<ComplexImageType, ComplexImageType, WaveletFilterBankType> ForwardWaveletType;
  typename ForwardWaveletType::Pointer forwardWavelet = ForwardWaveletType::New();
  unsigned int high_sub_bands = inputBands;
  unsigned int levels = inputLevels;
  forwardWavelet->SetHighPassSubBands(high_sub_bands);
  forwardWavelet->SetLevels(levels);
  forwardWavelet->SetInput(fftForwardFilter->GetOutput());
  forwardWavelet->Update();

  typedef itk::MonogenicSignalFrequencyImageFilter<ComplexImageType> MonogenicSignalFrequencyFilterType;

  typedef typename MonogenicSignalFrequencyFilterType::OutputImageType VectorMonoOutputType;
  typedef itk::VectorInverseFFTImageFilter<VectorMonoOutputType> VectorInverseFFTType;

  // Input to the PhaseAnalysisEigenValues
  typedef MonogenicPhaseAnalysisSoftThresholdImageFilter<typename VectorInverseFFTType::OutputImageType> PhaseAnalysisFilter;
  // typedef MonogenicPhaseAnalysisEigenValuesImageFilter<typename VectorInverseFFTType::OutputImageType> PhaseAnalysisFilter;

  typename ForwardWaveletType::OutputsType analysisWavelets = forwardWavelet->GetOutputs();
  typename ForwardWaveletType::OutputsType modifiedWavelets;
  unsigned int noutputs = forwardWavelet->GetNumberOfOutputs();
  for (unsigned int nout = 0; nout < forwardWavelet->GetNumberOfOutputs(); ++nout)
    {
    std::cout << "************Nout: "<< nout << " / " <<noutputs << std::endl;
    if(nout == 0) // TODO check this. avoid phase analysis in last approximation (low_pass).
      {
      modifiedWavelets.push_back(analysisWavelets[nout]);
    //TODO remove this visualize
#if ITK_VISUALIZE_TESTS != 0
    typedef itk::InverseFFTImageFilter<ComplexImageType> FFTInverseFilterType;
    typename FFTInverseFilterType::Pointer fftInv = FFTInverseFilterType::New();
    fftInv->SetInput(analysisWavelets[nout]);
    fftInv->Update();
    Testing::ViewImage(fftInv->GetOutput(),  "Wavelet coef 0 (LowPass) Original" );
#endif
      continue;
      }
    typename MonogenicSignalFrequencyFilterType::Pointer monoFilter = MonogenicSignalFrequencyFilterType::New();
    typename VectorInverseFFTType::Pointer vecInverseFFT = VectorInverseFFTType::New();
    typename PhaseAnalysisFilter::Pointer phaseAnalyzer = PhaseAnalysisFilter::New();
    typename FFTForwardFilterType::Pointer fftForwardPhaseFilter = FFTForwardFilterType::New();
    monoFilter->SetInput(analysisWavelets[nout]);
    monoFilter->Update();
    vecInverseFFT->SetInput(monoFilter->GetOutput());
    vecInverseFFT->Update();
    phaseAnalyzer->SetInput(vecInverseFFT->GetOutput());
    phaseAnalyzer->Update();
    fftForwardPhaseFilter->SetInput(phaseAnalyzer->GetOutput(0));
    fftForwardPhaseFilter->Update();
    modifiedWavelets.push_back(fftForwardPhaseFilter->GetOutput());
    modifiedWavelets.back()->DisconnectPipeline();
    //TODO remove this visualize
#if ITK_VISUALIZE_TESTS != 0
    typedef itk::InverseFFTImageFilter<ComplexImageType> FFTInverseFilterType;
    typename FFTInverseFilterType::Pointer fftInv = FFTInverseFilterType::New();
    fftInv->SetInput(analysisWavelets[nout]);
    fftInv->Update();
    Testing::ViewImage(fftInv->GetOutput(),  "Wavelet coef " + n2s(nout) + "Original" );
    Testing::ViewImage(phaseAnalyzer->GetOutput(0), "Wavelet coef " + n2s(nout) + "PhaseAnalyzed" );
#endif
    }

  typedef itk::WaveletFrequencyInverse<ComplexImageType, ComplexImageType, WaveletFilterBankType> InverseWaveletType;
  typename InverseWaveletType::Pointer inverseWavelet = InverseWaveletType::New();
  inverseWavelet->SetHighPassSubBands(high_sub_bands);
  inverseWavelet->SetLevels(levels);
  inverseWavelet->SetInputs(modifiedWavelets);
  inverseWavelet->Update();

  typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType> InverseFFTFilterType;
  typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
  inverseFFT->SetInput(inverseWavelet->GetOutput());
  inverseFFT->Update();
#if ITK_VISUALIZE_TESTS != 0
    Testing::ViewImage(reader->GetOutput(), "Input Image");
    Testing::ViewImage(inverseFFT->GetOutput(), "Inverse Wavelet");
#endif

  typedef itk::ImageFileWriter<typename InverseFFTFilterType::OutputImageType>  WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImage);
  writer->SetInput(inverseFFT->GetOutput());

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error writing the image: " << std::endl;
    std::cerr << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

int itkRieszWaveletPhaseAnalysisTest(int argc, char *argv[])
{
  if( argc < 5 || argc > 6 )
    {
    std::cerr << "Usage: " << argv[0] <<
      " inputImage outputImage inputLevels inputBands [dimension]" << std::endl;
    return EXIT_FAILURE;
    }
  const string inputImage  = argv[1];
  const string outputImage = argv[2];
  const unsigned int inputLevels = atoi(argv[3]);
  const unsigned int inputBands  = atoi(argv[4]);

  unsigned int dimension = 3;
  if( argc == 6 )
  {
    dimension = atoi( argv[5] );
  }

  if( dimension == 2 )
  {
    return runRieszWaveletPhaseAnalysisTest<2>(inputImage, outputImage, inputLevels, inputBands);
  }
  else if( dimension == 3 )
  {
    return runRieszWaveletPhaseAnalysisTest<3>(inputImage, outputImage, inputLevels, inputBands);
  }
  else
  {
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
  }
}
