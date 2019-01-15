===============================================================
DICEseq: Dynamic Isoform spliCing Estimator via sequencing data
===============================================================

About DICEseq
=============

DICEseq (Dynamic Isoform spliCing Estimator via sequencing data) estimates the 
dynamics of isoform proportions jointly from time series RNA-seq experiments. 
DICEseq is a Bayesian method based on a mixture model whose mixing proportions 
represent isoform ratios. It incorporates the correlations from the temporal 
structure, by coupling the isoform proportions at different times through a 
latent Gaussian process (GP).

DICEseq provides following functions through command line:

1. ``diceseq``: estimate the isoform proportions and FPKM for time series data 
   jointly, or for a single time point. 

2. ``dice-count``: fetch reads counts for entile gene, or specific reads counts,
   e.g. junction reads, for genes with exact one intron. This is special design 
   mainly for yeast.

3. ``dice-bias``: estimate parameters for sequencing bias, including fragment 
   length distribution, reads sequence and position bias parameter. The output 
   file can be directly used for bias correction in ``diceseq``.

In addition, DICEseq package also provides interface of a set of functions and 
attributes as an object-oriented python module. Therefore, you could use some 
of the module e.g., ``SampleFile`` to visualize the samples in gzip file in a 
Gaussian process way, or ``BiasFile`` to visualize the bias parameters. Also, 
the ``gtf_utils`` provides a set of ways to load gtf file, choose the genes, or 
customize the coordinates of exons and introns, add and remove of specific 
transcripts.


Quick Start
===========

**Installation**: 

- ``pip install diceseq``
- or download this repository, and type ``python setup.py install``. 
- You may need to add ``--user`` if you don't have the root permission and you 
  don't use Anaconda_.

.. _Anaconda: https://www.continuum.io/anaconda-overview

**Arguments**

- Type command line ``diceseq -h``



Detailed Manual
===============

See the documentation_ on how to install, to use, to find the annotation data 
etc.

.. _documentation: http://diceseq.sourceforge.net


References
===========

Yuanhua Huang and Guido Sanguinetti. `Statistical modeling of isoform splicing 
dynamics from RNA-seq time series data 
<http://bioinformatics.oxfordjournals.org/content/32/19/2965.abstract>`_. \ **Bioinformatics**\, 2016, 32(19): 2965-2972.
