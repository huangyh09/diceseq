=====================
Release notes (major)
=====================

Release 0.1.7 (21/12/2015)
==========================

* largely improved for the case with multiple (more than 10) isoforms, and introduced EM algorithm for buin-in
* extended for python 3
* improved bias estimates
* not rely on h5py anymore, by design diceseq format for bias parameter and samples
* introduced an out_utils.py for downstream analysis of diceseq results
* fixed some bugs on paired-end reads matching and processing
* removed reads simulation, as it can be done by spanki_

  .. _spanki: http://www.cbcb.umd.edu/software/spanki/


Release 0.1.6 (10/11/2015)
==========================

* slight change to be compatible with pysam 0.8.3

* change numpy.math.erf() to personally wroten erf(), as the former is not supported in numpy 1.8.2 and higher

* provide the option to not save the hdf5 copy

* remove dice-static for direct running


Release 0.1.0 (26/10/2015)
==========================

* probabilistic models
	* finished dynamic model: estimate jointly

* dirrect running file
	* diceseq: support separated and joint estimates
	* dice-static: not suggested anymore.

* bug fixed
	* reads simulation


Release 0.0.1 (01/09/2015)
==========================

* first version of DICEseq

* utils objects or functions
	* FastaFile: read and write fasta file
	* BiasFile: read and write bias parameters
	* ReadSet: get reads location index of genome region
	* sam_utils: load and preprocess RNA-seq reads
	* gtf_utils: load gtf file as list of Gene and Transcript objects
	* TranUnits: processing at isoform level (or transcript level)
	* TranSplice: processing on multiple isoforms (or gene level)

* probabilistic models
	* static model: estimate indivually
	* dynamic model: estimate jointly (under study)

* dirrect running file
	* diceseq
	* dice-static
	* dice-count
	* dice-simulate
	* dice-bias
	* dice-biasplot

