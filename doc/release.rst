=============
Release notes
=============


Release 0.2.2 (29/10/2016)
==========================

* fix some bugs
* simplify the arguments
* change default theta1 as 10
* change default save_sample (original sample_num) into 500
* keep transcript in diceseq results in the same order as in gtf file


Release 0.2.1 (16/06/2016)
==========================

* updated calculating TOTAL_COUNTS for pysam 0.90
* update documentary by adding more examples and tutorials


Release 0.2.0 (16/05/2016)
==========================

* fixed some bugs
* cleaned some codes


Release 0.1.9 (29/01/2016)
==========================

* fixed some bugs
* tested the estimate of bias parameters and fragment length distribution
* introduced the parallel computation, new you can run it on multiple cores
* changed the format of output file for diceseq, which is easier to read
* introducted the auto detection of using single-end or paired-end reads
* removed the EM algorithm for buin-in, as we could easily use multiple cores to reduce running time


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

