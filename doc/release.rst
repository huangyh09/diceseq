=============
Release notes
=============

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

