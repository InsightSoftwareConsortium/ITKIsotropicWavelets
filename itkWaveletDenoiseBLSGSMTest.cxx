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

#include "itkWaveletDenoiseBLSGSM.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int itkWaveletDenoiseBLSGSMTest( int argc, char *argv[] )
{
  if ( argc != 7 )
    {
    std::cerr << "Usage: " << argv[0]
              << " inputImage outputImage inputLevels inputBands waveletFunction dimension" << std::endl;
    return EXIT_FAILURE;
    }
  const std::string inputImage  = argv[1];
  const std::string outputImage = argv[2];
  const unsigned int inputLevels = atoi( argv[3] );
  const unsigned int inputBands  = atoi( argv[4] );
  const std::string waveletFunction = argv[5];
  const unsigned int dimension = atoi( argv[6] );
}
