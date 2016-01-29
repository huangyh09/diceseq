=================
Usages of DICEseq
=================

This page contains the details on how to use the functions that DICEseq provides. Before using diceseq and dice-count, you need the annotation file, which you could download from sepcific database, e.g., Ensembl_, but you may need to change it to the format_ that diceseq supports.

.. _Ensembl: http://www.ensembl.org/info/data/ftp/index.html 



preprocess
==========

.. _format:

annotation format
-----------------

The basic format is gtf_. In order to simply process the gtf file, we used a few formats of key words. The default is similar as Ensembl format. Namely, the order is ``gene`` --> ``transcript 1`` --> ``exon ...`` --> ``transcript 2`` --> ``exon ...`` , and the attributes includes ``gene_id``, ``gene_name``, ect. as follows

::

  XV  ensembl gene  93395 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
  XV  ensembl transcript  93395 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
  XV  ensembl exon  94291 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
  XV  ensembl exon  93395 93843 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
  XI  ensembl gene  431906  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
  XI  ensembl transcript  431906  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
  XI  ensembl exon  431906  432034  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
  XI  ensembl exon  432433  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";

.. _gtf: http://www.ensembl.org/info/website/upload/gff.html

alignment
---------

There are quite a fewer aligner that allows mapping reads to genome reference with big gaps, mainly caused by splicing. You could use STAR_, TOPHAT_, but I would suggest HISAT_ here, which is fast and returns reads with good aligment quality.

You could run it like this (based on HISAT 0.1.5), which including alignment, sort and index:

::

  ($hisatDir/hisat -x $hisatRef -1 $fq_dir/"$file"_1.fq.gz -2 $fq_dir/"$file"_2.fq.gz --no-unal | samtools view -bS -> $out_dir/$file.bam) 2> $out_dir/$file.err
  samtools sort $out_dir/$file.bam $out_dir/$file.sorted
  samtools index $out_dir/$file.sorted.bam

.. _STAR: https://code.google.com/p/rna-star/
.. _TOPHAT: https://ccb.jhu.edu/software/tophat/index.shtml
.. _HISAT: https://ccb.jhu.edu/software/hisat/index.shtml


diceseq
=======

This command allows you to estimate isoform proportions jointly (or separately if you only input one time point each time). It also allows to merging multiple replicates. You could run it like this:

::

  my_sam_list=t1_rep1.sorted.bam,t1_rep2.sorted.bam---t1_rep1.sorted.bam---t3_rep1.sorted.bam

  diceseq --anno_file=anno_file.gtf --sam_list=$my_sam_list --out_file=out_file

There are more parameters for setting (``diceseq -h`` always give the version you are using):

* ``--anno_file`` or ``-a`` (default=None): The annotation file in gtf format.
* ``--anno_source`` (default=Ensembl): The annotation source of the gtf file.
* ``--sam_list`` or ``-s`` (default=None): The indexed alignement file in bam/sam format, use ``,`` for replicates and ``---`` for time points, e.g., my_sam1_rep1.sorted.bam,my_sam1_rep2.sorted.bam---my_sam2.sorted.bam.
* ``--ref_file`` or ``-r`` (default=None): The genome reference file in fasta format. This is necessary for bias correction, otherwise uniform mode will be used.
* ``--out_file`` or ``-o`` (default=diceseq_out): The prefix of the output file. There will be two output file, one in plain text format, the other in gzip format.
* ``--bias_file`` or ``-b`` (default=""): The files for bias parameter. You could use one bise file for all time points, or use T bias files for each time point, by ``---``, e.g., file1.bias---file2.bias
* ``--bias_mode`` (default=unif): The bias mode: unif, end5, end3 or both. Without ``ref_file`` or ``--bias_file``, it will be changed into unif.

* ``--time_seq`` or ``-t`` (default=None): The time for the input samples, e.g., 0,1,2,3, the default values will be the index of all time points, i.e., 0,1,...
* ``--sample_num`` (default=0): The number of MCMC samples to save, 0 for no such file. Advice: lower than 3/4 of `min_run`, e.g, 500.

* ``--nproc`` (default=4): The number of subprocesses.
* ``--add_premRNA``: Add the pre-mRNA as a transcript.
* ``--no_twice``: No quick estimate of the variance, but use fixed.
* ``--print_detail``: Print the detail of the sampling.

* ``--mate_mode`` (default=pair): The mode for using paired-end reads: auto, pair, single.
* ``--auto_min`` (default=200): The minimum pairs of read mates if mate_mode=auto.
* ``--fl_mean`` (default=None, i.e., auto dected): The mean of fragment length.
* ``--fl_std`` (default=None, i.e., auto dected): The standard deviation of fragment length.

* ``--max_run`` (default=5000): The maximum iterations for the MCMC sampler.
* ``--min_run`` (default=1000): The minimum iterations for the MCMC sampler.
* ``--gap_run`` (default=100): The increase gap of iterations for the MCMC sampler.
* ``--theta1`` (default=3.0): The fixed hyperparameter theta1 for the GP model.
* ``--theta2`` (default=None): The fixed hyperparameter theta2 for the GP model. The default will cover 1/3 of the duration.


Suggestions on setting hyperparameter :math:`\theta_2`: if you want :math:`\theta_2` cover :math:`\eta \in (0,1)` of duration, then you should set :math:`\theta_2=(\eta(t_{max}-t_{min}))^2`. The default is :math:`\eta = 1/3`. Generally, we suggest using a small :math:`\theta_2`, e.g., covering less than 1/3 length, while it really depends on the time scale.


dice-count
==========

This command allows you to calculate the reads counts in an aligned + sorted + indexed sam (or bam) file with an annotation file in gtf format. It allows calculating the total counts for each gene, but also specific counts of different segments (e.g., junction, exon, and intron) if a gene has exactly one intron. You could run it like this:

::

  dice-count --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt

There are more parameters for setting (``dice-count -h`` always give the version you are using):

* ``--anno_file`` or ``-a`` (default=None): The annotation file in gtf format.
* ``--anno_source`` (default=Ensembl): The annotation source of the gtf file.
* ``--sam_file`` or ``-s`` (default=None): The indexed alignement file in bam/sam format;
* ``--out_file`` or ``-o`` (default=dice_count.txt): the counts in plain text file;

* ``--duplicate``: Keep duplicate reads.
* ``--partial``: Keep reads partial in the region.
* ``--single_end``: Use the reads as single-end.
* ``--junction``: return junction and boundary reads, only for gene with one exon-intron-exon structure; other wise return total counts for the whole gene.

* ``--nproc`` (default=4): The number of subprocesses.
* ``--mapq_min`` (default=10): The minimum mapq for reads.
* ``--mismatch_max`` (default=5): The maximum mismatch for reads.
* ``--rlen_min`` (default=1): The mimimum length of reads.

