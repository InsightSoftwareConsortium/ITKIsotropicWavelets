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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkViewImage.h"

using namespace std;
using namespace itk;

int
runViewImage(int argc, char* argv[])
{
  if ( argc < 2 || argc == 4 || argc > 5 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImage [title] [win_size_x win_size_y] " << std::endl;
    return EXIT_FAILURE;
    }
  // Defaults
  std::string win_title = "itkViewImage";
  size_t win_x = 600;
  size_t win_y = 600;
  if ( argc >= 3 )
    {
    win_title = argv[2];
    }
  if ( argc == 5 )
    {
    win_x = atoi(argv[3]);
    win_y = atoi(argv[4]);
    }
  const string inputImage  = argv[1];

  const unsigned int dimension = 3;
  typedef float                              PixelType;
  typedef itk::Image< PixelType, dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImage);
  reader->Update();

  itk::Testing::ViewImage(reader->GetOutput(), win_title, win_x, win_y );

  return EXIT_SUCCESS;
}
