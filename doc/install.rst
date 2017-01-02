============
Installation
============

Easy install
============

* Required packages in python: `numpy`, `matplotlib`, `pysam`

  * we suggest using Anaconda_ distribute, which includes most packages (except 
    `pysam` here), and provides you a user specific environment, i.e., all your 
    python packages go to a single folder. Thus, you don't need the root to 
    install packages.

  * you could install `pysam` by pypi in terminal, or download_ and install as 
    diceseq.

  .. _Anaconda: http://continuum.io/downloads
  .. _download: https://github.com/pysam-developers/pysam

* You can install `DICEseq` simply via pypi in terminal (suggested), or upgrade 
  by add ``--upgrade`` as follows:

  ::

    pip install diceseq

    pip install --upgrade --no-deps diceseq


Source code
===========

* Alternatively, you also could download the source code via GitHub (latest 
  version, suggested) or Sourceforge (any version) and run python setup in 
  terminal:

  * GitHub: https://github.com/huangyh09/diceseq

  * Sourceforge: http://sourceforge.net/projects/diceseq/

  ::

    python setup.py install

* In any case, if had the permission error for installation as you are not 
  root, add ``--user``.
