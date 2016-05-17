DICEseq Homepage
================

:Author: Yuanhua Huang
:Version: 0.2.0
:Last viewed: May 17, 2016

About DICEseq
-------------

DICEseq (Dynamic Isoform spliCing Estimator via sequencing data) estimates the dynamics of isoform proportions jointly from time series RNA-seq experiments. DICEseq is a Bayesian method based on a mixture model whose mixing proportions represent isoform ratios. It incorporates the correlations from the temporal structure, by coupling the isoform proportions at different times through a latent Gaussian process (GP).

DICEseq provides following functions through command line:

1. ``diceseq``: 
  * Estimate the isoform proportions and FPKM jointly or separately.

  * Calculate the Bayes factor to detect the differential dynamics of splicing profile.

2. ``dice-count``: 
  * Get the total reads counts of each gene. 

  * Get the specific reads counts, e.g., junction reads, for genes with exactly one intron. This is particularly design for yeast.

3. ``dice-bias``: 
  * Estimate fragment length distribution.

  * Estimate reads sequence and position bias parameter. The output file can be directly used for bias correction in ``diceseq``.

In addition, DICEseq package also provides interface of a set of functions and attributes as an object-oriented python module. Therefore, you could use some of the module e.g., ``SampleFile`` to visualize the samples in gzip file in a Gaussian process way, or ``BiasFile`` to visualize the bias parameters. Also, the ``gtf_utils`` provides a set of ways to load gtf file, choose the genes, or customize the coordinates of exons and introns, add and remove of specific transcripts.


Contents
--------

.. toctree::
   :maxdepth: 2

   install.rst
   usage.rst
   api.rst
   faq.rst
   release.rst


References
----------

Yuanhua Huang and Guido Sanguinetti. \ **Statistical modeling of isoform splicing dynamics from RNA-seq time series data**\. (`under review`_)

.. _under review: 


.. seealso::

   SourceForge package: including splicing annotation files
      http://sourceforge.net/projects/diceseq/files/annotation/

   MISO package: Mixture of Isoforms
      https://miso.readthedocs.io/en/fastmiso/

   RNA-seq aligner HISAT
      https://ccb.jhu.edu/software/hisat/index.shtml

   Information about pysam
      http://pysam.readthedocs.org

   The python language
      http://www.python.org
