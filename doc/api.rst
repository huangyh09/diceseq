===============================================
iSplice, intactive methods kit for RNA splicing
===============================================

Introduction
============

iSplice is a python module that take each isoform as an object and quantify its proportions.


API
===

Isoforms
--------

Objects of type :class:`~isplice.TranUnits` allow working on arbitrary isoform with customized the coordinates of exons and introns.

.. autoclass:: pysam.TranUnits
   :members:

An :class:`~pysam.TranSplice` represents a mixture of multiple isoforms.

.. autoclass:: pysam.TranSplice
   :members:

.. autoclass:: pysam.Transcript
   :members:

Reads set
---------

:class:`~pysam.ReadSet` index the reads into specific location on genome, e.g., upstream exon, boundary, junction.

.. autoclass:: pysam.ReadSet
   :members:

Bias files
----------

:class:`~pysam.BiasFile` allow write and read sequence and position bias parameters in to a hdf5 file.

.. autoclass:: pysam.BiasFile
   :members:

Fasta files
-----------

.. autoclass:: pysam.FastaFile
   :members:





