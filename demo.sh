#!/bin/bash
#a demo for running diceseq

anno_file=data/yeast/Saccharomyces_cerevisiae.R64-1-1.77_filtered.gtf

# demo 1: running the main program of diceseq to estimate the isoform proportions.
sam_file=data/yeast/readp_100_3.sorted.bam---data/yeast/readp_100_5.sorted.bam
out_file=data/yeast/ratios.txt
dice-static --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file

# demo 2: get the counts
sam_file=data/yeast/readp_100_3.sorted.bam
out_file=data/yeast/count.txt

dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file #--total_only=False