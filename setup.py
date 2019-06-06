"""
DICEseq - Dynamic Isoform spliCing Estimator via sequencing data
See: http://diceseq.sourceforge.net
"""

import sys
if sys.version_info >= (3,):
    sys.exit('diceseq requires Python 2')

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./diceseq/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy', 'pysam', 'matplotlib']

# from distutils.core import setup
# from Cython.Build import cythonize
# extensions = [Extension("*", ["diceseq/*/*.pyx"])]

setup(
    name='diceseq',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='DICEseq: Dynamic Isoform spliCing Estimator via sequencing data',
    long_description=long_description,

    # The project's main homepage.
    url='http://diceseq.sourceforge.net',

    # Author details
    author='Yuanhua Huang',
    author_email='yuanhua@ed.ac.uk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['splicing isoform quantification', 'time series RNA-seq',
              'Gaussian process', 'Markov chain Monte Carlo'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),#['diceseq', 'diceseq.utils', 'diceseq.models'],#

    entry_points={
          'console_scripts': [
              'diceseq = diceseq.diceseq:main',
              # 'dice-diff = diceseq.dice_diff:main',
              'dice-bias = diceseq.dice_bias:main',
              'dice-count = diceseq.dice_count:main',
              ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,

    py_modules = ['diceseq'],

    # ext_modules = cythonize(extensions)

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...

)