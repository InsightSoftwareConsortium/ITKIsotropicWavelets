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


#include "itkCastImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimoncelliIsotropicWavelet.h"
#include "itkWaveletFrequencyForward.h"
#include "itkWaveletFrequencyFilterBankGenerator.h"
#include "itkTestingMacros.h"


template< unsigned int VDimension >
int runSimoncelliIsotropicWaveletTest( char* inputImageFilename,
  char* outputImageFilename, unsigned int highPassSubBands )
{
  typedef float                                     InputPixelType;
  typedef unsigned char                             OutputPixelType;
  typedef itk::Point< InputPixelType, VDimension >  PointType;
  typedef itk::Image< InputPixelType, VDimension >  InputImageType;
  typedef itk::Image< OutputPixelType, VDimension > OutputImageType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFilename );

  TRY_EXPECT_NO_EXCEPTION( reader->Update() );


  typedef itk::SimoncelliIsotropicWavelet< InputPixelType, VDimension, PointType >
    SimoncelliIsotropicWaveletType;

  SimoncelliIsotropicWaveletType::Pointer simoncelliIsotropicWavelet =
    SimoncelliIsotropicWaveletType::New();

  TRY_EXPECT_EXCEPTION( simoncelliIsotropicWavelet->SetHighPassSubBands( 0 ) );

  simoncelliIsotropicWavelet->SetHighPassSubBands( highPassSubBands );
  TEST_SET_GET_VALUE( highPassSubBands, simoncelliIsotropicWavelet->GetHighPassSubBands() );


  // Perform FFT on input image
  typedef itk::ForwardFFTImageFilter< InputImageType > FFTFilterType;
  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();

  fftFilter->SetInput( reader->GetOutput() );


  typedef typename FFTFilterType::OutputImageType ComplexImageType;
  typedef itk::SimoncelliIsotropicWavelet<>       WaveletFunctionType;

  typedef itk::WaveletFrequencyFilterBankGenerator< ComplexImageType, WaveletFunctionType >
    WaveletFilterBankType;
  typedef itk::WaveletFrequencyForward< ComplexImageType, ComplexImageType, WaveletFilterBankType >
    ForwardWaveletType;

  typename ForwardWaveletType::Pointer forwardWavelet = ForwardWaveletType::New();

  forwardWavelet->SetInput( fftFilter->GetOutput() );

  TRY_EXPECT_NO_EXCEPTION( forwardWavelet->Update() );

  // Inverse FFT
  typedef itk::InverseFFTImageFilter< ComplexImageType, InputImageType > InverseFFTFilterType;
  typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();

  unsigned int subBand = 0;
  inverseFFT->SetInput( forwardWavelet->GetOutput( subBand ) );

  TRY_EXPECT_NO_EXCEPTION( inverseFFT->Update() );

  // Cast the output
  typedef itk::CastImageFilter< InputImageType, OutputImageType > CastImageFilterType;
  CastImageFilterType::Pointer caster = CastImageFilterType::New();
  caster->SetInput( inverseFFT->GetOutput() );

  // Write output
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( caster->GetOutput() );
  writer->SetFileName( outputImageFilename );

  TRY_EXPECT_NO_EXCEPTION( writer->Update() );


  std::cerr << "Test finished" << std::endl;
  return EXIT_SUCCESS;
}


int itkSimoncelliIsotropicWaveletTest( int argc, char *argv[] )
{
  if( argc != 5 )
    {
    std::cerr << "Usage: " << argv[0]
              << " inputImage outputImage dimension highPassSubBands" << std::endl;
    return EXIT_FAILURE;
    }
  const unsigned int ImageDimension = 3;
  typedef double                                          PixelType;
  typedef std::complex< PixelType >                       ComplexPixelType;
  typedef itk::Point< PixelType, ImageDimension >         PointType;
  typedef itk::Image< ComplexPixelType, ImageDimension >  ComplexImageType;

  // Exercise basic object methods
  // Done outside the helper function in the test because GCC is limited
  // when calling overloaded base class functions.
  typedef itk::SimoncelliIsotropicWavelet< PixelType, ImageDimension, PointType >
    SimoncelliIsotropicWaveletType;

  SimoncelliIsotropicWaveletType::Pointer simoncellidIsotropicWavelet =
    SimoncelliIsotropicWaveletType::New();
  EXERCISE_BASIC_OBJECT_METHODS( simoncellidIsotropicWavelet, SimoncelliIsotropicWavelet,
    IsotropicWaveletFrequencyFunction );


  char* inputImageFilename  = argv[1];
  char* outputImageFilename = argv[2];

  unsigned int highPassSubBands = atoi( argv[4] );

  unsigned int dimension = atoi( argv[3] );

  if( dimension == 2 )
    {
    return runSimoncelliIsotropicWaveletTest< 2 >( inputImageFilename,
      outputImageFilename, highPassSubBands );
    }
  else if( dimension == 3 )
    {
    return runSimoncelliIsotropicWaveletTest< 3 >( inputImageFilename,
      outputImageFilename, highPassSubBands );
    }
  else
    {
    std::cerr << "Test failed!" << std::endl;
    std::cerr << "Error: only 2 or 3 dimensions allowed, " << dimension << " selected." << std::endl;
    return EXIT_FAILURE;
    }
}
