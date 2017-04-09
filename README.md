[![DOI](https://zenodo.org/badge/67762635.svg)](https://zenodo.org/badge/latestdoi/67762635)

# ITKIsotropicWavelets
External Module for ITK, implementing Isotropic Wavelets and Riesz Filter for multiscale phase analysis.

Review in itk: http://review.source.kitware.com/#/c/21512/

# Commit message:
WIP: Add External Module IsotropicWavelets.

Module (external) that adds Isotropic Wavelet analysis.

##TODO:

- [ ] Add Generalized Riesz Filter Bank. Riesz of order N, including derivatives.
- [ ] Add Steering framework.
  - [ ] General case, U matrix from Chenouard, Unser.
  - [ ] Simoncelli Equiangular case
- [x] Add FrequencyBandImageFilter
- [x] Add Monogenic Signal Phase Analysis.
 - It now reproduces Held work as a brightness equalizator / local phase detector.
- [x] Add Simoncelli, Shannon, Held and Vow Isotropoic Wavelets.
- [x] Add Shrinker and Expander in spatial domain with no interpolation.
- [x] Add StructureTensor.
- [ ] Publish in InsightJournal about implementation and update handle.
- [ ] Add simple test to every wavelet (Vow,Held, Simoncelli, Shannon), instead of relying on the implicit testing with the WaveletBankGenerator.
- [ ] Report bug related to VNL with image size being multiple of 3 [21,21,9],
but generating exception about invalid size. More info in the cmake file for tests.

The work is inspired by the Monogenic Signal from literature, that uses
wavelets and riesz filter to provide a multiscale denoise and segmentation
mechanism.

The Riesz filter is a Hilbert transform for ND, that provides phase
information, ie feature detection, in every dimension.

Wavelets are really important in signal analysis, they are able to perform a
multiscale analysis of a signal.
Similar to a windowed FourierTransform, but with the advantage that the
spatial resolution (the window) can be modulated, retaining more
information from the original image.

In this implementation only IsotropicWavelet are considered. These are
non-separable wavelets that depend on the modulo of the frequency vector.
There are only 4 or 5 Mother Wavelets developed in the literature,
I implemented 2 of them here from respective papers (see specific docs for more info).
The mean advantage of IsotropicWavelets is that they are steerable, as
shown by Simoncelli, steering the wavelet at each location provides
adaptability to different signal, and can be used along PCA methods to
select the best matching 'steer' at each location and scale.

Input to filters in this module needs to be in the dual space (frequency).
For example, from the output of an forward FFT. The decision is made to
avoid performing multiple FFT.
Also a FrequencyShrinker and an Expander WITHOUT any interpolation, just
chopping and adding zeros has been added.

Because the layout of the frequencies after an FFT is implementation
dependent (FFTW and VNL should share the same layout, but python FFT
might be different, etc), I added an iterator to abstract this layout.
It has a function GetFrequencyIndex(), that facilitates implementation
of further frequency filters.
Right now this iterator has been tested with the option ITK_USES_FFTW,
but should work for the default VNL.


# Summary, files:

## Frequency Iterators:

```
itkFrequencyImageRegionConstIteratorWithIndex.h
itkFrequencyImageRegionIteratorWithIndex.h
```

## FrequencyFunctions
### Base and Derived Classes:

* `itkFrequencyFunction.h`

  * `itkIsotropicFrequencyFunction.h`

    * `itkIsotropicWaveletFrequencyFunction.h itkIsotropicWaveletFrequencyFunction.hxx`

### Wavelets Functions (IsotropicWaveletFrequencyFunction):

```
itkHeldIsotropicWavelet.h
itkHeldIsotropicWavelet.hxx

itkSimoncelliIsotropicWavelet.h
itkSimoncelliIsotropicWavelet.hxx

itkShannonIsotropicWavelet.h
itkShannonIsotropicWavelet.hxx

itkVowIsotropicWavelet.h
itkVowIsotropicWavelet.hxx
```

### Riesz Function (FrequencyFunction):

```
itkRieszFrequencyFunction.h
itkRieszFrequencyFunction.hxx
```

## Frequency Related Image Filters:

### Frequency Expand/Shrinkers

```
itkFrequencyExpandImageFilter.h
itkFrequencyExpandImageFilter.hxx
itkFrequencyShrinkImageFilter.h
itkFrequencyShrinkImageFilter.hxx

itkFrequencyExpandViaInverseFFTImageFilter.h
itkFrequencyExpandViaInverseFFTImageFilter.hxx
itkFrequencyShrinkViaInverseFFTImageFilter.h
itkFrequencyShrinkViaInverseFFTImageFilter.hxx
```

### MonogenicSignal Filter (Riesz Function in all dimensions)

```
itkMonogenicSignalFrequencyImageFilter.h
itkMonogenicSignalFrequencyImageFilter.hxx
```

### FrequencyBand Filter (pass or stop freq band)

```
itkFrequencyBandImageFilter.h
itkFrequencyBandImageFilter.hxx
```

## Generators (use functions to create ImageSources)

```
itkWaveletFrequencyFilterBankGenerator.h
itkWaveletFrequencyFilterBankGenerator.hxx

itkRieszFrequencyFilterBankGenerator.h
itkRieszFrequencyFilterBankGenerator.hxx
```

## Forward/Inverse Wavelet (ImageFilter, apply wavelet pyramid using generators)

```
itkWaveletFrequencyForward.h
itkWaveletFrequencyForward.hxx

itkWaveletFrequencyInverse.h
itkWaveletFrequencyInverse.hxx
```


## Wavelet independent:

Local estimator over a neighborhood. Get the linear combination of input that maximize the response at every pixel.

```
itkStructureTensor.h
itkStructureTensor.hxx
```

FFTPad but avoiding setting negative index, which is problematic working with neighborhoods.

```
itkFFTPadPositiveIndexImageFilter.h
itkFFTPadPositiveIndexImageFilter.hxx
```

### Regular shrinkers without interpolation

```
itkExpandWithZerosImageFilter.h
itkExpandWithZerosImageFilter.hxx
itkShrinkDecimateImageFilter.h
itkShrinkDecimateImageFilter.hxx
```

### Wrappers without new functionality:

```
itkVectorInverseFFTImageFilter.h
itkVectorInverseFFTImageFilter.hxx

itkZeroDCImageFilter.h
itkZeroDCImageFilter.hxx
```

### Helpers (Linear index to subindex array)

```
itkInd2Sub.h
```

### Phase Analysis:

```
itkPhaseAnalysisImageFilter.h
itkPhaseAnalysisImageFilter.hxx

itkPhaseAnalysisSoftThresholdImageFilter.h
itkPhaseAnalysisSoftThresholdImageFilter.hxx
```
