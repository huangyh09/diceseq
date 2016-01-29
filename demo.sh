#!/bin/bash
#a demo for running diceseq

anno_file=data/anno/yeast_RNA_splicing.gtf

# demo 1: dice-count
sam_file=data/sam/reads_t3.sorted.bam

## total count
out_file=data/out/t3_cnt1.txt
dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file --nproc=4

## specific reads count
out_file=data/out/t3_cnt2.txt
dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file --junction --nproc=4


# demo 2: diceseq
sam_dir=data/sam

## separated model
out_file=data/out/t1
sam_list=$sam_dir/reads_t1.sorted.bam
diceseq --anno_file=$anno_file --add_premRNA --sam_list=$sam_list --out_file=$out_file

## joint model
out_file=data/out/joint
sam_list=$sam_dir/reads_t1.sorted.bam---$sam_dir/reads_t2.sorted.bam---$sam_dir/reads_t3.sorted.bam
diceseq --anno_file=$anno_file --add_premRNA --sam_list=$sam_list --out_file=$out_file


