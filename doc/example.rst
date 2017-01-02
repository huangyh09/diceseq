========
Examples
========

The demo file `demo.sh <https://github.com/huangyh09/diceseq/blob/master/demo.sh>`_ and also the `data <https://github.com/huangyh09/diceseq/tree/master/data>`_ are included in `github <https://github.com/huangyh09/diceseq>`_

diceseq demo
============

* Run diceseq for isoform proportions (add ``--add_premRNA`` if you want pre-mRNA as an extra transcript):

  * separated model

  ::

    anno_file=data/anno/yeast_RNA_splicing.gtf
    sam_dir=data/sam
    out_file=data/out/t1
    sam_list=$sam_dir/reads_t1.sorted.bam
    diceseq -a $anno_file -s $sam_list -o $out_file --add_premRNA

  * joint model

  ::

    out_file=data/out/joint
    sam_list=$sam_dir/reads_t1.sorted.bam---$sam_dir/reads_t2.sorted.bam---$sam_dir/reads_t3.sorted.bam
    diceseq -a $anno_file -s $sam_list -o $out_file --add_premRNA


dice-count demo
===============

* Run dice-count for reads counts:

  * total count

  ::

    anno_file=data/anno/yeast_RNA_splicing.gtf
    sam_file=data/sam/reads_t3.sorted.bam
    out_file=data/out/t3_cnt1.txt
    dice-count -a $anno_file -s $sam_list -o $out_file

  * specific reads count, e.g., junctions (make sure it contains one intron only)

  ::
  
    out_file=data/out/t3_cnt2.txt
    dice-count  -a $anno_file -s $sam_list -o $out_file --junction

For more options, see the detailed usages page.

