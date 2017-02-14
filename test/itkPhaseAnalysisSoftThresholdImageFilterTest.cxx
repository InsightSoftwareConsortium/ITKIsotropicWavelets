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
#include <cmath>
#include "itkPhaseAnalysisSoftThresholdImageFilter.h"
#include "itkMonogenicSignalFrequencyImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"

#include "itkVectorInverseFFTImageFilter.h"
//Visualize for dev/debug purposes. Set in cmake file. Require VTK
#ifdef ITK_VISUALIZE_TESTS
#include "itkViewImage.h"
#endif
using namespace std;
using namespace itk;

int itkPhaseAnalysisSoftThresholdImageFilterTest(int argc, char* argv[])
{
  if( argc != 3 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
    }
  const string inputImage  = argv[1];
  const string outputImage = argv[2];

  const unsigned int dimension = 3;
  typedef float                            PixelType;
  typedef itk::Image<PixelType, dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImage);
  reader->Update();
  reader->UpdateLargestPossibleRegion();

  // Perform FFT on input image.
  typedef itk::ForwardFFTImageFilter<ImageType> FFTForwardFilterType;
  FFTForwardFilterType::Pointer fftForwardFilter = FFTForwardFilterType::New();
  fftForwardFilter->SetInput(reader->GetOutput());
  fftForwardFilter->Update();
  typedef FFTForwardFilterType::OutputImageType ComplexImageType;

  // Get a Monogenic Vector. Other input to PhaseAnalysis could be derivatives.
  typedef itk::MonogenicSignalFrequencyImageFilter<ComplexImageType> MonogenicSignalFrequencyFilterType;
  MonogenicSignalFrequencyFilterType::Pointer monoFilter = MonogenicSignalFrequencyFilterType::New();
  monoFilter->SetInput(fftForwardFilter->GetOutput());
  monoFilter->Update();
  typedef MonogenicSignalFrequencyFilterType::OutputImageType VectorMonoOutputType;

  typedef itk::VectorInverseFFTImageFilter<VectorMonoOutputType> VectorInverseFFTType;
  VectorInverseFFTType::Pointer vecInverseFFT = VectorInverseFFTType::New();
  vecInverseFFT->SetInput(monoFilter->GetOutput());
  vecInverseFFT->Update();

  // Input to the PhaseAnalysisSoftThreshold
  typedef PhaseAnalysisSoftThresholdImageFilter<VectorInverseFFTType::OutputImageType> PhaseAnalysisSoftThresholdFilter;
  PhaseAnalysisSoftThresholdFilter::Pointer phaseAnalyzer = PhaseAnalysisSoftThresholdFilter::New();
  phaseAnalyzer->SetInput(vecInverseFFT->GetOutput());
  phaseAnalyzer->SetApplySoftThreshold(true);
  phaseAnalyzer->Update();
  PhaseAnalysisSoftThresholdFilter::OutputImageType::Pointer cosPhase = phaseAnalyzer->GetOutputCosPhase();
  PhaseAnalysisSoftThresholdFilter::OutputImageType::Pointer amp      = phaseAnalyzer->GetOutputAmplitude();
  PhaseAnalysisSoftThresholdFilter::OutputImageType::Pointer phase    = phaseAnalyzer->GetOutputPhase();

#ifdef ITK_VISUALIZE_TESTS
    Testing::ViewImage(cosPhase.GetPointer(), "PhaseAnalyzer(Soft) output:");
#endif

  return EXIT_SUCCESS;
}
