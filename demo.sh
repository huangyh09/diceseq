#!/bin/bash
#a demo for running diceseq

anno_file=data/yeast/Saccharomyces_cerevisiae.R64-1-1.77_filtered.gtf
gene_file=data/yeast/gene_name_use.lst
bias_file=data/yeast/merged_bias.tophat.hdf5
sam_file=data/yeast/readp_100_5.sorted.bam
out_file=data/yeast/out.txt

# estimate isoform proportions with static model
diceseq --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file  --bias_file=$bias_file --gene_file=$gene_file #--ref_file=$ref_file

