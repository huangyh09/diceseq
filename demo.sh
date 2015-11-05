#!/bin/bash
#a demo for running diceseq

anno_file=data/anno/yeast_RNA_splicing.gtf
gene_file=data/anno/gene_id_use.lst
bias_file=data/bias/bias_weights.hdf5
out_file=data/out/out

# demo 1: dice-count
sam_file=data/sam/reads_t3.sorted.bam

# total count
out_file=data/out/t3_cnt1.txt
dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file

# specific reads count
out_file=data/out/t3_cnt2.txt
dice-count --anno_file=$anno_file --sam_file=$sam_file --out_file=$out_file --total_only=False


# demo 2: diceseq
sam_dir=data/sam

# separated model
out_file=data/out/t1
sam_list=$sam_dir/reads_t1.sorted.bam
diceseq --anno_file=$anno_file --add_premRNA=True --sam_list=$sam_list --out_file=$out_file

# joint model
out_file=data/out/joint
sam_list=$sam_dir/reads_t1.sorted.bam---$sam_dir/reads_t2.sorted.bam---$sam_dir/reads_t3.sorted.bam
diceseq --anno_file=$anno_file --add_premRNA=True --sam_list=$sam_list --out_file=$out_file



# demo 1: reads simulation
# ref_file=/media/Data/research/splicing/data/Annotation/Saccharomyces_cerevisiae.R64-1-1.75.fa
# out_dir=/media/Data/research/splicing/data/dice/demo/fastq
# dice-simulate --anno_file=$anno_file --anno_source=miso --ref_file=$ref_file --RPK=100 --ratio=0.10,0.15,0.45,0.85 --rlen=75 --out_file=$out_dir/reads_unif

# # reads alignment
# fasta_file=/media/Data/research/splicing/data/Annotation/yeast/Saccharomyces_cerevisiae.R64-1-1.75.fa
# hisatRef=/media/Data/research/splicing/data/Annotation/yeast/hisatRef/Saccharomyces_cerevisiae.R64-1-1.75
# hisatDir=/home/yuanhua/Tool/HISAT/hisat-0.1.5-beta

# ratios=( 0.10 0.15 0.45 0.85 )
# fq_dir=/media/Data/research/splicing/data/dice/demo/fastq
# out_dir=/media/Data/research/splicing/data/dice/demo/hisatSam

# for j in "${ratios[@]}";
# do
# 	file=reads_unif_100_$j
	
# 	($hisatDir/hisat -x $hisatRef -1 $fq_dir/$file"_1.fq.gz" -2 $fq_dir/$file"_2.fq.gz" --no-unal | samtools view -bS -> $out_dir/$file.bam) 2> $out_dir/$file.err
# 	samtools sort $out_dir/$file.bam $out_dir/$file.sorted
# 	samtools index $out_dir/$file.sorted.bam
# 	rm $out_dir/$file.bam
# done

