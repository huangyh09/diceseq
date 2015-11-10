DICEseq: Dynamic Isoform spliCing Estimator via sequencing data
===============================================================

About DICEseq
-------------

Different from most methods that quantifies the splicing isoforms statically, DICEseq estimates the dynamics of isoform proportions jointly from time series RNA-seq experiments. DICEseq is a Bayesian method based on a mixture model whose mixing proportions represent isoform ratios; however, DICEseq incorporates the correlations induced by the temporal structure by coupling the isoform proportions at different times through a latent Gaussian process (GP).

DICEseq provides following functions:

1. Estimate the isoform proportions jointly. The prior is GPs followed by a softmax functions transform.
2. Estimate the isoform proportions separately. It is almost the same as MISO_, but with different prior distribution.
3. Calculate the Bayes factor to detect the differential dynamics of splicing profile.
4. Get the total counts of each gene.
5. Get the specific counts, e.g., junction reads, for genes with exactly one intron.
6. Generating simulation reads with given isoform proportions and sequencing coverage.
7. Reads sequence and position bias correction and plot (under study)

.. _MISO: http://genes.mit.edu/burgelab/miso/

In addition to run the DICEseq functions from standard command line, DICEseq also provides interface of a set of functions and attributes as an object-oriented python module. Therefore, you could, for example customize the coordinates of exons and introns, add and remove of specific transcripts.

More information
----------------

See the homepage_ of DICEseq for all information links, and the documentation_ on how to install, use, etc.

.. _homepage: http://diceseq.sourceforge.net
.. _documentation: http://diceseq.readthedocs.org


References
----------

Yuanhua Huang and Guido Sanguinetti. \ **Statistical modeling of isoform splicing dynamics from RNA-seq time series data**\. (`under review`_)

.. _under review: 
