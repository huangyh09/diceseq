DICEseq Homepage
================

:Author: Yuanhua Huang
:Version: 0.1.7
:Last viewed: December 21, 2015

About DICEseq
-------------

DICEseq (Dynamic Isoform spliCing Estimator via sequencing data) estimates the dynamics of isoform proportions jointly from time series RNA-seq experiments. DICEseq is a Bayesian method based on a mixture model whose mixing proportions represent isoform ratios. It incorporates the correlations from the temporal structure, by coupling the isoform proportions at different times through a latent Gaussian process (GP).

DICEseq provides following functions:

1. Estimate the isoform proportions jointly. The prior is GPs followed by a softmax functions transformation.
2. Estimate the isoform proportions separately. It is the same as MISO_, except prior.
3. Calculate the Bayes factor to detect the differential dynamics of splicing profile.
4. Get the total counts of each gene.
5. Get the specific counts, e.g., junction reads, for genes with exactly one intron.
6. Reads sequence and position bias correction and plot

.. _MISO: http://genes.mit.edu/burgelab/miso/

In addition to run the DICEseq functions from standard command line, DICEseq also provides interface of a set of functions and attributes as an object-oriented python module. Therefore, you could, for example customize the coordinates of exons and introns, add and remove of specific transcripts.


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

   Alternative splicing annotation (by DICEseq)
      http://sourceforge.net/projects/diceseq/files/annotation/

   RNA-seq aligner HISAT
      https://ccb.jhu.edu/software/hisat/index.shtml

   Information about pysam
      http://pysam.readthedocs.org

   The python language
      http://www.python.org
