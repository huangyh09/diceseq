===
FAQ
===

Contacts
========
If you have any question about using DICEseq, you are welcome to discuss with 
the developer. Email: Yuanhua Huang <Y.Huang@ed.ac.uk>


Report bugs
===========
If you find any error or suspicious bug, I will appreciate your report.
Please write them in the github issues: 
https://github.com/huangyh09/diceseq/issues


DICEseq for single time point
=============================
Though DICEseq is developed for time series RNA-seq, it also performs very well
for single time point, like most RNA-seq experiment. With ``theta1=3.0``, you 
should get similar results as MISO.


input ``sam_list`` in `diceseq`
===============================

The format of ``sam_list`` is a string line, e.g.

::

  sam_list=$sam_dir/sam1_rep1.sorted.bam,$sam_dir/sam1_rep2.sorted.bam---$sam_dir/sam2.sorted.bam---$sam_dir/sam3.sorted.bam
  diceseq --anno_file=$anno_file --sam_list=$sam_list --out_file=$out_file

Please do not directly use a list file containing the sam files. Otherwise you 
need load it via bash command by yourself. Note: the delimiter ``,`` means 
merging the two samples, and ``---`` means different time points.


No output from DICEseq
======================
Currently, DICEseq requires gtf/gff3 file with specific order. 
gene-->transcript1-->exon1...exon_k-->transcript2...

If the annotation is not in such sorted or if there is no isoform annotation,
you probably won't get any output. 


Compatibilities
===============

compatibility with pysam
------------------------

If you are using `pysam` 0.8.2 or higher version, please use `diceseq` 0.1.2 or 
higher.


compatibility with numpy
------------------------

If you are using `numpy` 1.8.2 or higher version, please use `diceseq` 0.1.3 or 
higher.

