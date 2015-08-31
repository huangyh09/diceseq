==============
API of DICEseq
==============

Introduction
============

DICEseq is a python module that processes the RNA-seq data for basic purposes, e.g., getting reads counts and showing sequence bias, and modeles the mixture of transcripts (isoform) and the mRNA / pre-mRNA ratios. More uniquely, this module allows to statistically modelling the time series RNA-seq in studying the dynamics of isform regulations.

In addition, this module provides a set of object-oriented interface, with which users could develope their own codes for customized analysis. Also, a set of statistics functions are also easily available.


GTF file and isoforms
=====================

The objects of :class:`~diceseq.Transcript` and :class:`~diceseq.Gene` allow to load data from gtf annotation file. Usually, a gene could contain one or multiple transcripts, including a pre-mRNA. Each transcript could contains one or multiple exons.

.. autoclass:: diceseq.Transcript
   :members:

.. autoclass:: diceseq.Gene
   :members:

The objects of :class:`~diceseq.TranUnits` are based :class:`~diceseq.Transcript`, and allow setting sequecing data, bias parameters, etc. One could also add arbitrary transcript with customized the coordinates of exons and introns.

.. autoclass:: diceseq.TranUnits
   :members:

The objects of :class:`~diceseq.TranSplice` are based :class:`~diceseq.Gene`, and allow modelling the mixture of multiple transcripts (or isoforms) by using the observed sequecing reads. 

.. autoclass:: diceseq.TranSplice
   :members:



Reads set
=========

The objects of :class:`~diceseq.ReadSet` index the reads into specific location on genome, e.g., upstream exon, boundary, junction.

.. autoclass:: diceseq.ReadSet
   :members:

Bias files
==========

The objects of :class:`~diceseq.BiasFile` allow write and read sequence and position bias parameters in to a hdf5 file.

.. autoclass:: diceseq.BiasFile
   :members:

Fasta files
===========
The objects of :class:`~diceseq.FastaFile` allow load fasta file and get sequence from any part of the fasta file.

.. autoclass:: diceseq.FastaFile
   :members:

