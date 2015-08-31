=================
Usages of DICEseq
=================

dice-count
==========

This command allows you to calculate the reads counts in an aligned + sorted + indexed sam (or bam) file with an annotation file in gtf format. It allows getting the counts for each gene, but also seperating the counts into different segments (e.g., junction, exon, and intron) if the gene just contains two exons. You could run it like this:

::

  dice-count --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format;
* ``--anno_source`` (default="Ensembl"): the annotation source of the gtf file, current supporting Ensemble, SGD, MISO, and Sander (a personalized way);
* ``--sam_file`` or ``-s`` (default=None): the indexed alignement file in bam/sam format;
* ``--total_reads`` (default=None): the total aligned reads for calculating RPKM;
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use;
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


dice-static
===========

This command allows you to calculate the proportions of each isoform or the ratio of mRNA and pre-mRNA. This is static model, which considers the data for only one time point, but still allows to merging multiple replicates. You could run it like this:

  ::

    dice-static --anno_file=anno_file.gtf --sam_file=sam_file.bam --out_file=out_file.txt  --bias_file=bias_file.hdf5 --ref_file=ref_file.fasta

There are more parameters for setting:

* ``--anno_file`` or ``-a`` (default=None): the annotation file in gtf format;
* ``--sam_file`` or ``-s`` (default=None): the indexed alignement file in bam/sam format, e.g., my_sam1.sorted.bam---my_sam2.sorted.bam;
* ``--ref_file`` or ``-r`` (default=None): the genome reference file in fasta format;
* ``--gene_file`` or ``-g`` (default=None): the list of genes in use;
* ``--out_file`` or ``-o`` (default=estimated_rate.psi): The output file in txt format;
* ``--bias_file`` or ``-b``(default=None): the parameter file for bias in hdf5 format;
* ``--bias_mode`` (default=unif): The bias mode: unif (uniform), end3, end5, both;


More ...
===========

Comming soon.