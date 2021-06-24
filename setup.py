# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sys

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

setup(
    name='itk-isotropicwavelets',
    version='0.6.1',
    author='Pablo Hernandez-Cerdan',
    author_email='pablo.hernandez.cerdan@outlook.com',
    packages=['itk'],
    package_dir={'itk': 'itk'},
    download_url=r'https://github.com/InsightSoftwareConsortium/ITKIsotropicWavelets',
    description=r'Steerable Wavelet Pyramid with Isotropic Wavelets and Riesz Functions',
    long_description='itk-isotropicwavelets provides a multiresolution analysis '
                     '(MRA) framework using isotropic and steerable wavelets in the frequency '
                     'domain.\n'
                     'This framework provides the backbone for state of the art filters '
                     'for denoising, feature detection or phase analysis in N-dimensions. '
                     'It focuses on reusability, and highly decoupled modules for easy extension '
                     'and implementation of new filters.\n'
                     'Please refer to:'
                     'P. Hernandez-Cerdan, “Isotropic and Steerable Wavelets in N Dimensions.\n'
                     'A multiresolution analysis framework”, Insight Journal, January-December 2016 '
                     'http://hdl.handle.net/10380/3558.',
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Android",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='Apache',
    keywords='ITK InsightToolkit Wavelet Pyramid Steerable Riesz',
    url=r'https://github.com/InsightSoftwareConsortium/ITKIsotropicWavelets',
    install_requires=[
        r'itk>=5.2.post3',
    ]
    )
