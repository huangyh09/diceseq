DICEseq: Dynamic Isoform spliCing estimator via sequencing data
===============================================================

:Author: Yuanhua Huang
:Version: 0.1.2
:Last viewed: November 09, 2015

Different from most methods that quantifies the splicing isoforms statically, DICEseq estimates the dynamics of isoform proportions jointly from time series RNA-seq experiments. DICEseq is a Bayesian method based on a mixture model whose mixing proportions represent isoform ratios; however, DICEseq incorporates the correlations induced by the temporal structure by coupling the isoform proportions at different times through a latent Gaussian process (GP).

DICEseq offers a set of functions to run from the standard command line, and also provides interface of more functions and attributes as an object-oriented python module. Therefore, you could, for example customize the coordinates of exons and introns, add and remove of specific transcripts.

Contents
--------

.. toctree::
   :maxdepth: 2

   install.rst
   usage.rst
   api.rst
   models.rst
   faq.rst
   release.rst


References
----------

Yuanhua Huang and Guido Sanguinetti.\ **Statistical modeling of isoform splicing dynamics from RNA-seq time series data** \. [(under review)]()

.. seealso::

   Alternative splicing annotation for human and mouse (by MISO)
      http://miso.readthedocs.org/en/fastmiso/annotation.html

   RNA-seq aligner HISAT
      https://ccb.jhu.edu/software/hisat/index.shtml

   Information about pysam
      http://pysam.readthedocs.org

   The python language
      http://www.python.org
