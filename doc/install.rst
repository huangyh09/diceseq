============================
Installation and quick start
============================

Installation
============

* Required packages in python: `numpy`, `matplotlib`, `pysam`

  * we suggest using Anaconda_ distribute, which includes most packages but `pysam`
  * you could install `pysam` by pypi in terminal, or download_ and install as diceseq.

  .. _Anaconda: http://continuum.io/downloads
  .. _download: https://github.com/pysam-developers/pysam

* You can install `DICEseq` simply via pypi in terminal (suggested), or upgrade by add ``--upgrade`` as follows:

  ::

    pip install diceseq

    pip install --upgrade --no-deps diceseq


* Or you could download the source code via GitHub (latest version, suggested) or Sourceforge (any version) and run python setup in terminal:

  * GitHub: https://github.com/huangyh09/diceseq
  * Sourceforge: http://sourceforge.net/projects/diceseq/

  ::

    python setup.py install


Quick start with demo
=====================

The demo file `demo.sh <https://github.com/huangyh09/diceseq/blob/master/demo.sh>`_ and also the `data <https://github.com/huangyh09/diceseq/tree/master/data>`_ are included in `github <https://github.com/huangyh09/diceseq>`_

* Run diceseq for isoform proportions:

  * separated model

  ::

    anno_file=data/anno/yeast_RNA_splicing.gtf
    sam_dir=data/sam
    out_file=data/out/t1
    sam_list=$sam_dir/reads_t1.sorted.bam
    diceseq --anno_file=$anno_file --add_premRNA=True --sam_list=$sam_list --out_file=$out_file

  * joint model

  ::

    out_file=data/out/joint
    sam_list=$sam_dir/reads_t1.sorted.bam---$sam_dir/reads_t2.sorted.bam---$sam_dir/reads_t3.sorted.bam
    diceseq --anno_file=$anno_file --add_premRNA=True --sam_list=$sam_list --out_file=$out_file



* Run dice-count for reads counts:

  * total count

  ::

    anno_file=data/anno/yeast_RNA_splicing.gtf
    sam_file=data/sam/reads_t3.sorted.bam
    out_file=data/out/t3_cnt1.txt
    dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file

  * specific reads count (make sure it contains one intron only)

  ::
  
    out_file=data/out/t3_cnt2.txt
    dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file --total_only=False

For more options, see the detailed usages page.

