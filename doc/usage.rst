=================
Usages of DICEseq
=================

This page contains the details on how to use the functions that DICEseq provides. Before using diceseq and dice-count, you need the annotation file, which you could download from sepcific database, e.g., Ensembl_, but you may need to change it to the format_ that diceseq supports.

.. _Ensembl: http://www.ensembl.org/info/data/ftp/index.html 



Preprocess
==========

.. _format:

Annotation format
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

Alignment
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

  diceseq --anno_file=anno_file.gtf --sam_list=my_sam_list --out_file=out_file

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format.
* ``--sam_list`` or ``-s`` (default=None): the indexed alignement file in bam/sam format, use "," for replicates and "---" for time points, e.g., my_sam1_rep1.sorted.bam,my_sam1_rep2.sorted.bam---my_sam2.sorted.bam.
* ``--time`` or ``-t`` (default=None): The time for the input samples, e.g., 0,1,2,3, the default values will be the index of all time points, i.e., 0,1,...
* ``--ref_file`` or ``-r`` (default=None): the genome reference file in fasta format. This is necessary for bias correction, otherwise uniform mode will be used.
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use. It is the gene id in the gtf annotation. Default is all genes in annotations.
* ``--out_file`` or ``-o`` (default=diceseq_out): The prefix of the output file. There will be two output file, one in plain text format, the other in hdf5 format.
* ``--bias_file`` or ``-b`` (default=None): the parameter file for bias in hdf5 format.
* ``--bias_mode`` (default=unif): The bias mode: unif (uniform), end3, end5, both.
* ``--sample_num`` (default=500): The number of MCMC samples to save.
* ``--add_premRNA`` (default=False): Whether adding pre-mRNA or not.
* ``--theta1_fix`` (default=3.0): The fixed hyperparameter theta1 for the GP model.
* ``--theta2_fix`` (default=None): The fixed hyperparameter theta2 for the GP model. The default will cover 1/3 of the duration.
* ``--is_twice`` (default=True): Whether estimate the rates twice with a quick check first. It is useful for ensuring the 30-50% acceptances in MH sampler.

Suggestions on setting hyperparameter :math:`\theta_2`: if you want :math:`\theta_2` cover :math:`\eta \in (0,1)` of duration, then you should set :math:`\theta_2=(\eta(t_{max}-t_{min}))^2`. The default is :math:`\eta = 1/3`.

theta_2 setting (with covering duration as \eta, default: \eta=1/3): 

  .. math::

    \theta_2=(\eta(t_{max}-t_{min}))^2



dice-count
==========

This command allows you to calculate the reads counts in an aligned + sorted + indexed sam (or bam) file with an annotation file in gtf format. It allows calculating the total counts for each gene, but also specific counts of different segments (e.g., junction, exon, and intron) if a gene has exactly one intron. You could run it like this:

::

  dice-count --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format;
* ``--anno_source`` (default="Ensembl"): the annotation source of the gtf file, current supporting Ensemble, SGD, MISO, and Sander (a personalized way);
* ``--sam_file`` or ``-s`` (default=None): the indexed alignement file in bam/sam format;
* ``--total_reads`` (default=None): the total aligned reads for calculating RPKM;
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use. It is the gene id in the gtf annotation. Default is all genes in annotations;
* ``--out_file`` or ``-o`` (default=dice_count.txt): the counts in plain text file;
* ``--rm_duplicate`` (default=True): remove duplicate reads or not;
* ``--inner_only`` (default=True): only include the reads inside or not;
* ``--mapq_min`` (default=10): the minimum mapq for reads;
* ``--mismatch_max`` (default=5): the maximum mismatch for reads;
* ``--rlen_min`` (default=1): the mimimum length of reads;
* ``--is_mated`` (default=True): process reads as paired-end or not;
* ``--total_only`` (default=True): provide total reads count only (for a whole gene); if False, then the specific reads for the exon-intron-exon structure will be provide;
* ``--biotype_rm`` (default=None): the exclusive biotype(s); e.g., snRNA or snRNA---tRNA.
* ``--biotype_only`` (default=None): the only used biotype(s); e.g., snRNA or snRNA---tRNA.


dice-simulate
=============

This command allows generating simulated reads. You could run it like this:

::

  dice-simulate --anno_file=anno_file.gtf --out_file=out_file --ref_file=ref_file.fasta

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format;
* ``--sam_file`` or ``-s`` (default=None): the indexed alignement file in bam/sam format, e.g., my_sam1.sorted.bam---my_sam2.sorted.bam;
* ``--ref_file`` or ``-r`` (default=None): the genome reference file in fasta format;
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use. It is the gene id in the gtf annotation. Default is all genes in annotations.
* ``--out_file`` or ``-o`` (default=out_file): The prefix of output file with simulated reads;
* ``--bias_file`` or ``-b``(default=None): the parameter file for bias in hdf5 format;
* ``--bias_mode`` (default=unif): The bias mode: unif (uniform), end3, end5, both;
* ``--add_premRNA`` (default=False): Whether adding pre-mRNA or not.
* ``--RPK`` (default=1000): The all used sequence depths, e.g., 100,200,400 and 100;
* ``--ratio`` (default=0.5): The all ratios of the first C-1 isoform, e.g., 0.3,0.6,0.8 and 0.3;
* ``--noise`` (default=0.001): The noise in the reads number for each isoform;
* ``--rlen`` (default=100): The length of reads;
* ``--fl_mean`` (default=200): The mean length of fragment;
* ``--fl_sigma`` (default=20): The stand variance of fragment length;

