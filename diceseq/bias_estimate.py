# This file is to esitmate the sequecne and positon bias parameters.

import sys
import h5py
import numpy as np
from optparse import OptionParser
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads

def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="The indexed alignement file in bam/sam format")
    parser.add_option("--refseq_file", "-r", dest="refseq_file", default=None,
        help="The fasta file of genome reference sequences")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="parameters.hdf5", help="The output file in hdf5 format")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    if options.sam_file == None:
        print "Error: need --sam_file for reads indexed and aliged reads."
        sys.exit(1)
    
    sam_file    = options.sam_file
    out_file    = options.out_file
    anno_file   = options.anno_file
    refseq_file = options.refseq_file

    #part 0.1: load the data
    biasFile  = BiasFile()
    fastaFile = FastaFile(refseq_file)
    samFile   = load_samfile(sam_file)
    anno      = load_annotation(anno_file)

    #part 1. filtering the annotations
    idx     = (anno["exon_num"] == 1) * (anno["biotype"] == "protein_coding")
    exons   = anno["exons"][idx]
    gene_id = anno["gene_id"][idx]
    chrom   = anno["chrom"][idx]
    strand  = anno["strand"][idx]
    tran_start = anno["tran_start"][idx]
    tran_stop  = anno["tran_stop"][idx]
    tran_len   = np.abs(tran_stop - tran_start) + 1

    #part 2.0: initialize a bias file
    biasFile.set_percentile(tran_len, 5)
    print biasFile.percentile
    print sum(idx)

    #part 2.1: get the bias and uniform weights
    for i in range(len(gene_id)):
        # i=154#2297
        print i, gene_id[i], chrom[i], tran_start[i], tran_stop[i]
        
        # 2.1 get the reads and seg_len
        reads = fetch_reads(samFile, chrom[i], tran_start[i], tran_stop[i], 
                            rm_duplicate=True, inner_only=True, mapq_min=10,
                            mismatch_max=10, rlen_min=1, is_mated=False)
        if len(reads["reads1u"])+len(reads["reads2u"])==0:
            print "No reads for %s" %gene_id[i]
            continue
        else:
            print (len(reads["reads1u"]), len(reads["reads2u"]))

        tlen = tran_len[i]
        for j in range(2):
            if   j == 0:
                end = 5
                reads_use = reads["reads1u"];
                if strand[i] == "+" or strand[1] == "1":
                    forward = True
                else:
                    forward = False
            elif j == 1:
                end = 3
                reads_use = reads["reads2u"]
                if strand[i] == "+" or strand[1] == "1":
                    forward = False
                else:
                    forward = True
            if len(reads_use) == 0: continue

            rleng  = 0
            weight = 1.0 / len(reads_use)
            for r in reads_use:
                rleng += r.rlen
                # tlen = tran_len[i] - r.rlen + 1 # comment this for whole length
                if forward == True:
                    pos = r.pos
                    seq = fastaFile.get_seq(chrom[i], pos-8, pos+12)
                    pos = pos - tran_start[i]
                else:
                    pos = r.aend
                    seq = fastaFile.get_seq(chrom[i], pos-12, pos+8)[::-1]
                    pos = tran_stop[i] - pos
                biasFile.set_both_bias(seq, pos, tlen, weight, end, "bias")
            rleng = int(round(rleng / (len(reads_use)+0.0)))

            if forward == True:
                sites = np.arange(tran_start[i], tran_stop[i]-rleng+1+1)
            else:
                sites = np.arange(tran_start[i]+rleng, tran_stop[i]+1)
            weight = 1.0 / sites.shape[0]
            for k in range(sites.shape[0]):
                # tlen = tran_len[i] - rleng + 1 # comment this for whole length
                pos = sites[k]
                if forward == True:
                    seq = fastaFile.get_seq(chrom[i], pos-8, pos+12)
                    pos = pos - tran_start[i]
                else:
                    seq = fastaFile.get_seq(chrom[i], pos-12, pos+8)[::-1]
                    pos = tran_stop[i] - pos
                biasFile.set_both_bias(seq, pos, tlen, weight, end, "unif")

    biasFile.save_file(out_file)

if __name__ == "__main__":
    main()
