=================
Usages of DICEseq
=================

diceseq
=======

This command allows you to estimate isoform proportions jointly (or separately if you only input one time point each time). It also allows to merging multiple replicates. You could run it like this:

  ::

    diceseq --anno_file=anno_file.gtf --sam_list=my_sam_list --out_file=out_file

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format.
* ``--sam_list`` or ``-s`` (default=None): the indexed alignement file in bam/sam format, use "," for replicates and "---" for time points, e.g., my_sam1_rep1.sorted.bam,my_sam1_rep2.sorted.bam---my_sam2.sorted.bam.
* ``--time`` or ``-t`` (default=None): The time for the input samples, e.g., 0,1,2,3, the default values will be the index of all time points, i.e., 0,1,...
* ``--ref_file`` or ``-r`` (default=None): the genome reference file in fasta format. This is necessary for bias correction, otherwise uniform mode will be used.
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use. It is the gene id in the gtf annotation. Default is all genes in annotations.
* ``--out_file`` or ``-o`` (default=diceseq_out): The prefix of the output file. There will be two output file, one in plain text format, the other in hdf5 format.
* ``--bias_file`` or ``-b``(default=None): the parameter file for bias in hdf5 format.
* ``--bias_mode`` (default=unif): The bias mode: unif (uniform), end3, end5, both.
* ``--sample_num`` (default=500): The number of MCMC samples to save.
* ``--add_premRNA`` (default=False): Whether adding pre-mRNA or not.
* ``--theta1_fix`` (default=3.0): The fixed hyperparameter theta1 for the GP model.
* ``--theta2_fix`` (default=None): The fixed hyperparameter theta2 for the GP model. The default will cover 1/3 of the duration.
* ``--is_twice`` (default=True): Whether estimate the rates twice with a quick check first. It is useful for ensuring the 30-50% acceptances in MH sampler.


Suggestions on setting hyperparameter $\theta_2$: given a coverage of duration $p$, then $\theta_2=(p(t_{max}-t_{min}))^2$.


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

