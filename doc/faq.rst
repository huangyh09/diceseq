===
FAQ
===

input ``sam_list`` in `diceseq`
===============================

The format of ``sam_list`` is a string line, e.g.

::

  sam_list=$sam_dir/sam1_rep1.sorted.bam,$sam_dir/sam1_rep2.sorted.bam---$sam_dir/sam2.sorted.bam---$sam_dir/sam3.sorted.bam
  diceseq --anno_file=$anno_file --sam_list=$sam_list --out_file=$out_file

Please do not directly use a list file containing the sam files. Otherwise you need load it via bash command by yourself. Note: the delimiter ``,`` means merging the two samples, and ``---`` means different time points.


compatibility with pysam
------------------------

If you are using `pysam` 0.8.2 or higher version, please use `diceseq` 0.1.2 or higher.


compatibility with numpy
------------------------

If you are using `numpy` 1.8.2 or higher version, please use `diceseq` 0.1.3 or higher.

