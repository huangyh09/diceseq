============================
Installation and quick start
============================

Installation
============

* Required packages in python:
    * `h5py`, `scipy`, `numpy`, `matplotlib`, `pysam`
    * we suggest using Anaconda distribute (http://continuum.io/downloads), which includes most packages but `pysam`
    * you could install `pysam` by pypi in terminal, but you may need to re-start of the terminal for using it.

    ::
        
      pip install pysam

* You can install `DICEseq` simply via pypi in terminal (suggested):

  ::
        
    pip install diceseq

* Or you could download the source code via GitHub (latest version, suggested) or Sourceforge (any version) and run python setup in terminal:
    * GitHub: https://github.com/huangyh09/diceseq
    * Sourceforge: http://sourceforge.net/projects/diceseq/

  ::
        
    python setup.py install


Quick start
===========

* Run diceseq for isoform proportions:

  ::

    diceseq --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt  --bias_file=bias_file.hdf5 --ref_file=ref_file.fasta

* Run dice-count for reads counts:

  ::

    dice-count --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt

You could also see the file demo.sh in the source code package. And for more information, see the detailed usages page.

