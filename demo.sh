#!/bin/bash
#a demo for running diceseq

anno_file=data/anno/yeast_RNA_splicing.gtf

# demo 1: dice-count
sam_file=data/sam/reads_t3.sorted.bam

## total count
out_file=data/out/t3_cnt1.txt
dice-count -a $anno_file -s $sam_file -o $out_file -p 4

## specific reads count
out_file=data/out/t3_cnt2.txt
dice-count -a $anno_file -s $sam_file -o $out_file -p 4 --junction --overhang=3


# demo 2: diceseq
sam_dir=data/sam

## separated model
out_file=data/out/t1
sam_list=$sam_dir/reads_t1.sorted.bam
diceseq -a $anno_file -s $sam_list -o $out_file -p 4 --add_premRNA

## joint model
out_file=data/out/joint
sam_list=$sam_dir/reads_t1.sorted.bam---$sam_dir/reads_t2.sorted.bam---$sam_dir/reads_t3.sorted.bam
diceseq -a $anno_file -s $sam_list -o $out_file -p 4 --add_premRNA --time_seq=1.5,2.5,5 

