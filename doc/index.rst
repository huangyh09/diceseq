DICEseq: Dynamic RNA splicing estimator via sequencing data
===========================================================

:Author: Yuanhua Huang
:Version: 0.0.1
:Last viewed: 01/09/2015

DICEseq is a Dynamic RNA splicing estimator via sequencing data. It could quantify alternative splicing isoform or RNA splicing ratio, i.e., proportion of mature mRNA from pre-mRNA. More generally, the program is designed object-oriented for transcripts. Thus, it supports arbitrary transcript, with customized the coordinates of exons and introns. It also supports addition and removal of specific transcripts.

This program also provides customized reads counting, especially boundary/junction reads, normalization, reads simulation, sequencing bias correction, etc.

Currently, the program for static quantification has been tested on the simulated and real data and showed good performance, but the dynamical model is still in preparation. 

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

Coming soon.

.. seealso::

   Alternative splicing annotation for human and mouse (by MISO)
      http://miso.readthedocs.org/en/fastmiso/annotation.html

   RNA-seq aligner HISAT
      https://ccb.jhu.edu/software/hisat/index.shtml

   Information about pysam
      http://pysam.readthedocs.org

   The python language
      http://www.python.org
