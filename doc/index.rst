DICEseq: Dynamic RNA splicing estimator via sequencing data
================================================

:Author: Yuanhua Huang
:Version: |version|
:Last viewed: June 20, 2015

DICEseq is a Dynamic RNA splicing estimator via sequencing data. It could quantify alternative splicing isoform or RNA splicing ratio, i.e., proportion of mature mRNA from pre-mRNA. More generally, the program is designed object-oriented for transcripts. Thus, it supports arbitrary transcript, with customised the coordinates of exons and introns. It also supports addition and removal of specific transcripts.

This program also provides customised reads counting, especially boundary/junction reads, normalisation, reads simulation, sequencing bias correction, etc.

Currently, the program for static quantification has been tested on the simulated and real data and showed good performance, but the dynamical model is still in preparation. 

Contents
--------

.. toctree::
   :maxdepth: 2

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
