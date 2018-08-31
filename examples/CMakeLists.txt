cmake_minimum_required(VERSION 3.5)
project(IsotropicWavelets-riesz_phase_analysis)

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

add_executable(runRieszWaveletAnalysis runRieszWaveletAnalysis.cpp)
target_link_libraries(runRieszWaveletAnalysis PUBLIC ${ITK_LIBRARIES})
