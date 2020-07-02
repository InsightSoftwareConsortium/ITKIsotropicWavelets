ITKIsotropicWavelets
=================================

.. image:: https://github.com/InsightSoftwareConsortium/ITKIsotropicWavelets/workflows/Build,%20test,%20package/badge.svg

.. image:: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_apis/build/status/InsightSoftwareConsortium.ITKIsotropicWavelets?branchName=master
    :target: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_build/latest?definitionId=5&branchName=master
    :alt:    Build Status

Overview
--------

.. |DOI| image:: https://zenodo.org/badge/67762635.svg
    :target: https://zenodo.org/badge/latestdoi/67762635

|DOI|

External Module for ITK, implementing Isotropic Wavelets and Riesz Filter for
multiscale phase analysis.

This document describes the implementation of the external module
`ITKIsotropicWavelets`, a multiresolution (MRA) analysis framework using
isotropic and steerable wavelets in the frequency domain. This framework
provides the backbone for state of the art filters for denoising, feature
detection or phase analysis in N-dimensions. It focuses on reusability, and
highly decoupled modules for easy extension and implementation of new filters,
and it contains a filter for multiresolution phase analysis.

The backbone of the multi-scale analysis is provided by an isotropic
band-limited wavelet pyramid, and the detection of directional features is
provided by coupling the pyramid with a generalized Riesz transform.

The generalized Riesz transform of order N behaves like a smoothed version of
the Nth order derivatives of the signal. Also, it is steerable: its components'
impulse responses can be rotated to any spatial orientation, reducing
computation time when detecting directional features.


For more information, see the `Insight Journal article <http://hdl.handle.net/10380/3558>`_::

  Hernandez-Cerdan P.
  Isotropic and Steerable Wavelets in N Dimensions. A multiresolution analysis framework
  The Insight Journal. January-December. 2016.
  http://hdl.handle.net/10380/3588
  http://www.insight-journal.org/browse/publication/986


Cite with::

  P. Hernandez-Cerdan, “Isotropic and Steerable Wavelets in N Dimensions. A multiresolution analysis framework for ITK,” arXiv:1710.01103 [cs], Oct. 2017.


Installation
------------

In ``C++``:

You need to `build ITK from source <https://itk.org/ITKSoftwareGuide/html/Book1/ITKSoftwareGuide-Book1ch2.html>`_
to use this module.

Since ITK version **4.13**, this module is available as a *Remote* module
in the ITK source code. 

Build it with the ``CMake`` option: ``Module_IsotropicWavelet``, this can be
switched on with a ``CMake`` graphical interface ``ccmake`` or directly from the
command line with: ``-DModule_IsotropicWavelet:BOOL=ON``

For **older** ITK versions (>4.10 required if ``BUILD_TEST=ON``), add it
manually as an *External* or *Remote* module to the ITK source code.

External::

  cd ${ITK_SOURCE_CODE}/Modules/External
  git clone https://github.com/phcerdan/ITKIsotropicWavelets


Remote:

Or create a file in ``${ITK_SOURCE_CODE}/Modules/Remote`` called
``IsotropicWavelets.remote.cmake`` (already there in ITK-4.13) with the
content::

  itk_fetch_module(IsotropicWavelets
    "IsotropicWavelets Extenal Module."
    GIT_REPOSITORY https://github.com/phcerdan/ITKIsotropicWavelets
    GIT_TAG master
    )

Python
^^^^^^

In ``Python``::

  pip install itk-isotropicwavelets



Components
----------

- [x] Steerable Pyramid in the frequency domain.
- [x] Undecimated Steerable Pyramid.
- [x] Generalized Riesz Filter Bank of order N (smoothed derivatives)
- [x] Steering framework (``RieszRotationMatrix``).
   - [NA] General case, U matrix from Chenouard, Unser.
   - [NA] Simoncelli Equiangular case.

- [x] ``FrequencyBandImageFilter``.
- [x] Monogenic Signal Phase Analysis.
   - Reproduces Held work as a brightness equalizator / local phase
     detector.

- [x] Simoncelli, Shannon, Held and Vow Isotropoic Wavelets.
- [x] Shrinker and Expander in spatial domain with no interpolation.
- [x] StructureTensor.
- [x] Simple test to every wavelet (Vow, Held, Simoncelli, Shannon),
  instead of relying on the implicit testing with the
  ``WaveletBankGenerator``.

The work is inspired by the monogenic signal from literature, that uses
wavelets and riesz filter to provide a multiscale denoise and segmentation
mechanism.

The Riesz filter is a Hilbert transform for ND, that provides phase
information, ie feature detection, in every dimension.

Wavelets are really important in signal analysis, they are able to perform a
multiscale analysis of a signal. Similar to a windowed Fourier Transform, but
with the advantage that the spatial resolution (the window) can be modulated,
retaining more information from the original image.

In this implementation only isotropic wavelets are considered. These are
wavelets that depend on the modulo of the frequency vector.
There are not many mother isotropic wavelets developed in the literature,
I implemented 4 of them here from respective papers (see specific docs for
more information). The main advantage of isotropic wavelets is that they are
steerable, as shown by Simoncelli, steering the wavelet at each location
provides adaptability to different signal, and can be used along PCA methods
to select the best matching "steer" at each location and scale.

Input to filters in this module needs to be in the dual space (frequency).
For example, from the output of an forward FFT. The decision is made to
avoid performing multiple FFT.

Also a `FrequencyShrinker` and an `Expander` WITHOUT any interpolation, just
chopping and adding zeros have been added.

Because the layout of the frequencies after an FFT is implementation
dependent (FFTW and VNL should share the same layout, but Python FFT
might be different, etc.), I added an iterator to abstract this layout.
It has a function GetFrequencyIndex(), that facilitates implementation
of further frequency filters.

Right now this iterator has been tested with the option ``ITK_USES_FFTW``,
but should work for the default VNL.

Summary of components and files:

Frequency Iterators
^^^^^^^^^^^^^^^^^^^

::

  itkFrequencyImageRegionConstIteratorWithIndex.h
  itkFrequencyImageRegionIteratorWithIndex.h
  itkFrequencyFFTLayoutImageRegionConstIteratorWithIndex.h
  itkFrequencyFFTLayoutImageRegionIteratorWithIndex.h
  itkFrequencyShiftedFFTLayoutImageRegionConstIteratorWithIndex.h
  itkFrequencyShiftedFFTLayoutImageRegionIteratorWithIndex.h


FrequencyFunctions
^^^^^^^^^^^^^^^^^^

Base and Derived Classes
''''''''''''''''''''''''

* ``itkFrequencyFunction.h``

  * ``itkIsotropicFrequencyFunction.h``

    * ``itkIsotropicWaveletFrequencyFunction.h itkIsotropicWaveletFrequencyFunction.hxx``

Wavelets Functions (``IsotropicWaveletFrequencyFunction``)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

::

  itkHeldIsotropicWavelet.h
  itkHeldIsotropicWavelet.hxx

  itkSimoncelliIsotropicWavelet.h
  itkSimoncelliIsotropicWavelet.hxx

  itkShannonIsotropicWavelet.h
  itkShannonIsotropicWavelet.hxx

  itkVowIsotropicWavelet.h
  itkVowIsotropicWavelet.hxx


Wavelets Generators (use functions to create ``ImageSources``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  itkWaveletFrequencyFilterBankGenerator.h
  itkWaveletFrequencyFilterBankGenerator.hxx


Riesz Function (``FrequencyFunction``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  itkRieszFrequencyFunction.h
  itkRieszFrequencyFunction.hxx


Riesz Generator (use functions to create ``ImageSources``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  itkRieszFrequencyFilterBankGenerator.h
  itkRieszFrequencyFilterBankGenerator.hxx


Frequency Related Image Filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Frequency Expand/Shrinkers
''''''''''''''''''''''''''

::

  itkFrequencyExpandImageFilter.h
  itkFrequencyExpandImageFilter.hxx
  itkFrequencyShrinkImageFilter.h
  itkFrequencyShrinkImageFilter.hxx

  itkFrequencyExpandViaInverseFFTImageFilter.h
  itkFrequencyExpandViaInverseFFTImageFilter.hxx
  itkFrequencyShrinkViaInverseFFTImageFilter.h
  itkFrequencyShrinkViaInverseFFTImageFilter.hxx


MonogenicSignal Filter (Riesz Function in all dimensions)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

::

  itkMonogenicSignalFrequencyImageFilter.h
  itkMonogenicSignalFrequencyImageFilter.hxx


FrequencyBand Filter (pass or stop freq band)
'''''''''''''''''''''''''''''''''''''''''''''

::

  itkFrequencyBandImageFilter.h
  itkFrequencyBandImageFilter.hxx


Forward/Inverse Wavelet (``ImageFilter``, apply wavelet pyramid using generators)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Decimated
'''''''''

::

  itkWaveletFrequencyForward.h
  itkWaveletFrequencyForward.hxx

  itkWaveletFrequencyInverse.h
  itkWaveletFrequencyInverse.hxx


Undecimated
'''''''''''

::

  itkWaveletFrequencyForwardUndecimated.h
  itkWaveletFrequencyForwardUndecimated.hxx

  itkWaveletFrequencyInverseUndecimated.h
  itkWaveletFrequencyInverseUndecimated.hxx


Wavelet independent
^^^^^^^^^^^^^^^^^^^

Local estimator over a neighborhood. Get the linear combination of inputs that
maximize the response at every pixel.

::

  itkStructureTensor.h
  itkStructureTensor.hxx


Regular shrinkers without interpolation
'''''''''''''''''''''''''''''''''''''''

::

  itkExpandWithZerosImageFilter.h
  itkExpandWithZerosImageFilter.hxx
  itkShrinkDecimateImageFilter.h
  itkShrinkDecimateImageFilter.hxx


Wrappers without new functionality
''''''''''''''''''''''''''''''''''

::

  itkVectorInverseFFTImageFilter.h
  itkVectorInverseFFTImageFilter.hxx

  itkZeroDCImageFilter.h
  itkZeroDCImageFilter.hxx


Helpers (Linear index to subindex array)
''''''''''''''''''''''''''''''''''''''''

::

  itkInd2Sub.h


Phase Analysis
''''''''''''''

::

  itkPhaseAnalysisImageFilter.h
  itkPhaseAnalysisImageFilter.hxx

  itkPhaseAnalysisSoftThresholdImageFilter.h
  itkPhaseAnalysisSoftThresholdImageFilter.hxx


Riesz Rotation Matrix (Steerable Matrix)
''''''''''''''''''''''''''''''''''''''''

::

  itkRieszRotationMatrix.h
  itkRieszRotationMatrix.hxx

  itkRieszUtilities.h
  itkRieszUtilities.cxx



License
-------

This software is distributed under the Apache 2.0 license. Please see the
*LICENSE* file for details.
