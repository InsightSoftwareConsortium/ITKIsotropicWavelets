cmake_minimum_required(VERSION 3.5)
project(IsotropicWaveletsExamples)

find_package(ITK REQUIRED COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  IsotropicWavelets
  ITKImageIO
  CONFIG
)
include(${ITK_USE_FILE})
include_directories(SYSTEM ${ITK_INCLUDE_DIRS})

add_executable(runWaveletCoeffsSpatialDomainImageFilter_example runWaveletCoeffsSpatialDomainImageFilter_example.cxx)
target_link_libraries(runWaveletCoeffsSpatialDomainImageFilter_example PUBLIC ${ITK_LIBRARIES})

add_executable(runWaveletCoeffsPhaseAnalyzisImageFilter_example runWaveletCoeffsPhaseAnalyzisImageFilter_example.cxx)
target_link_libraries(runWaveletCoeffsPhaseAnalyzisImageFilter_example PUBLIC ${ITK_LIBRARIES})

add_executable(runRieszWaveletAnalysis runRieszWaveletAnalysis.cpp)
target_link_libraries(runRieszWaveletAnalysis PUBLIC ${ITK_LIBRARIES})
