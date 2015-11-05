# This is a direct running file to calculate the specific counts of RNA
# splicing, and currently it is mainly for intronic splicing, rather than
# alternative splicing. There are 8 types for reads for an exon1-intron-
# exon2 structure: (0) exon1, (1) exon1-intron boundary, (2) intron, (3) 
# intron-exon2 boundary, (4) exon2, (5) exon1-exon2 junction, (6) exon1-
# intron-exon2, (7) exon1-exon2 unsure. In addition, it also provides the
# normalized counts RPK (reads per kilo-base).

import sys
import h5py
import pysam
import numpy as np
from optparse import OptionParser
from utils.reads_utils import ReadSet
from utils.gtf_utils import load_annotation
from utils.sam_utils import load_samfile, fetch_reads

def main():
    print "Welcome to dice-count!"

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--anno_source", dest="anno_source", default="Ensembl",
        help="The annotation source of the gtf file")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="The indexed alignement file in bam/sam format")
    parser.add_option("--total_reads", dest="total_reads", default=None,
        help="The total reads for calculate RPKM.")
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list of genes in use.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="dice_count.txt", help="The counts in plain text file")
    parser.add_option("--rm_duplicate", dest="rm_duplicate",
        default="True", help="remove duplicate reads or not.")
    parser.add_option("--inner_only", dest="inner_only",
        default="True", help="only unclude the reads inside or not.")
    parser.add_option("--mapq_min", dest="mapq_min",
        default="10", help="the minimum mapq for reads.")
    parser.add_option("--mismatch_max", dest="mismatch_max",
        default="5", help="the maximum mismatch for reads.")
    parser.add_option("--rlen_min", dest="rlen_min",
        default="1", help="the mimimum length of reads.")
    parser.add_option("--is_mated", dest="is_mated",
        default="True", help="process reads as paired-end or not.")
    parser.add_option("--total_only", dest="total_only", default="True",
        help="provide total reads count only; if False the specific reads" +
        " for the exon-intron-exon structure.")
    parser.add_option("--biotype_rm", dest="biotype_rm", default=None,
        help="The exclusive biotype.")
    parser.add_option("--biotype_only", dest="biotype_only", default=None,
        help="The only used biotype.")

    (options, args) = parser.parse_args()
    anno_source = options.anno_source
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file, anno_source)
    if options.sam_file == None:
        print "Error: need --sam_file for reads indexed and aliged reads."
        sys.exit(1)
    else:
        samFile = load_samfile(options.sam_file)

    if options.gene_file == None:
        gene_list = anno["gene_id"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")
    
    out_file     = options.out_file
    is_mated     = bool(options.is_mated == "True")
    mapq_min     = int(options.mapq_min)
    rlen_min     = int(options.rlen_min)
    inner_only   = bool(options.inner_only == "True")
    total_only   = bool(options.total_only == "True")
    biotype_rm   = options.biotype_rm
    total_reads  = options.total_reads
    biotype_only = options.biotype_only
    rm_duplicate = bool(options.rm_duplicate == "True")
    mismatch_max = int(options.mismatch_max)

    cnt_name = ["ex1", "ex1_int", "int", "int_ex2", "ex2", "ex1_ex2_junc",
                "ex1_int_ex2", "ex1_ex2_vague"]

    if total_reads is None: RPKM_symbol = "_RPK"
    else: RPKM_symbol = "_RPKM"

    fid = open(out_file, "w")
    if total_only != True:
        head_line = "gene_id\tgene_name\tbiotype"
        for i in range(len(cnt_name)):
            head_line += "\t" + cnt_name[i] + "_NUM"
        for i in range(len(cnt_name)):
            head_line += "\t" + cnt_name[i] + RPKM_symbol
    else:
        head_line = "gene_id\tgene_name\tbiotype\tcount\t" + RPKM_symbol[1:]
    fid.writelines(head_line + "\n")
    
    g_cnt = 0
    for g in range(gene_list.shape[0]):
        i = np.where(np.array(anno["gene_id"]) == gene_list[g])[0][0]

        if (biotype_rm is not None and 
            biotype_rm.split("---").count(anno["biotype"][i]) == 1): continue
        if (biotype_only is not None and 
            biotype_only.split("---").count(anno["biotype"][i]) == 0): continue

        _gene   = anno["gene_id"][i]
        _strand = anno["strand"][i]
        _chrom  = str(anno["chrom"][i])
        _start  = anno["gene_start"][i]
        _stop   = anno["gene_stop"][i]

        if total_only == True:
            reads = fetch_reads(samFile, _chrom, _start, _stop, rm_duplicate,
                            inner_only, mapq_min, mismatch_max,
                            rlen_min, is_mated)

            count = len(reads["reads1u"]) + len(reads["reads2u"])
            count += len(reads["reads2"])
            if total_reads is None:
                RPK   = count * 1000.0 / abs(_stop - _start + 1)
            else:
                RPK   = count * 10**9 / (abs(_stop - _start + 1) * float(total_reads))
            count = [count]
            RPK   = [RPK]
            g_cnt += 1
        elif (anno["genes"][i].tranNum > 0 and 
              anno["genes"][i].trans[0].exonNum == 2):
            reads = fetch_reads(samFile, _chrom, _start, _stop, rm_duplicate,
                            inner_only, mapq_min, mismatch_max,
                            rlen_min, is_mated)

            rdSet = ReadSet(reads["reads1u"])
            rdSet.get_loc_idx(anno["genes"][i].trans[0].exons, _strand)
            count = rdSet.loc_idx.sum(axis=0)
            RPK   = rdSet.RPK_use.sum(axis=0)

            rdSet = ReadSet(reads["reads2u"])
            rdSet.get_loc_idx(anno["genes"][i].trans[0].exons, _strand)
            count += rdSet.loc_idx.sum(axis=0)
            RPK   += rdSet.RPK_use.sum(axis=0)

            rdSet = ReadSet(reads["reads1"], reads["reads2"])
            rdSet.get_loc_idx(anno["genes"][i].trans[0].exons, _strand)
            count += rdSet.loc_idx.sum(axis=0)
            RPK   += rdSet.RPK_use.sum(axis=0)

            if total_reads is not None:
                RPK = RPK * 10**6 / float(total_reads)
            g_cnt += 1
        else: continue

        a_line = anno["gene_id"][i] + "\t" + anno["gene_name"][i] + "\t"
        a_line += anno["biotype"][i] + "\t"
        a_line += "\t".join(["%d" %num for num in list(count)]) + "\t" 
        a_line += "\t".join(["%.2f" %num for num in list(RPK)]) + "\n"
        fid.writelines(a_line)

        if g_cnt % 100 == 0 and g_cnt != 0:
            print "%d genes have been processed." %g_cnt
    fid.close()
    print "%d genes have been processed. Done!" %g_cnt

if __name__ == "__main__":
    main()
