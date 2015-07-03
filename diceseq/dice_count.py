# This is a direct running file to calculate the specific counts of RNA
# splicing, and currently it is mainly for intronic splicing, rather than
# alternative splicing. There are 7 types for reads for an exon1-intron-
# exon2 structure: (0) exon1, (1) exon1-intron boundary, (2) intron, (3) 
# intron-exon2 boundary, (4) exon2, (5) exon1-exon2 junction, (6) exon1-
# exon2 paired reads. In addition, it also provides the normalized counts
# RPK (reads per kilo-base).

import h5py
import pysam
import numpy as np
from optparse import OptionParser
from utils.reads_utils import ReadSet
from utils.gtf_utils import load_annotation
from utils.sam_utils import load_samfile, fetch_reads

if __name__ == "__main__":
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="The indexed alignement file in bam/sam format")
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


    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    if options.sam_file == None:
        print "Error: need --sam_file for reads indexed and aliged reads."
        sys.exit(1)
    if options.gene_file == None:
        gene_list = anno["gene_name"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")
    
    sam_file  = options.sam_file
    out_file  = options.out_file
    anno_file = options.anno_file
    rm_duplicate = bool(options.rm_duplicate)
    inner_only   = bool(options.inner_only)
    is_mated     = bool(options.is_mated)
    mapq_min     = int(options.mapq_min)
    rlen_min     = int(options.rlen_min)
    mismatch_max = int(options.mismatch_max)
    
    samFile   = load_samfile(sam_file)
    anno      = load_annotation(anno_file)

    print "Welcome to dice-count! %d genes are under counting." %len(gene_list)

    cnt_name = ["exon1", "exon1_intron", "intron", "intron_exon2", "exon2",
                "exon1_exon2_junc", "exon1_exon2_vague"]

    fid = open(out_file, "w")
    head_line = "gene"
    for i in range(len(cnt_name)):
        head_line += "\t" + cnt_name[i] + "_cnt"
    for i in range(len(cnt_name)):
        head_line += "\t" + cnt_name[i] + "_RPK"
    fid.writelines(head_line + "\n")
    
    for g in range(gene_list.shape[0]):
        i = np.where(np.array(anno["gene_name"]) == gene_list[g])[0][0]

        _gene   = anno["gene_id"][i]
        _exons  = anno["exons"][i]
        _strand = anno["strand"][i]
        _chrom  = str(anno["chrom"][i])
        _start  = anno["tran_start"][i]
        _stop   = anno["tran_stop"][i]

        reads = fetch_reads(samFile, _chrom, _start, _stop, rm_duplicate,
                            inner_only, mapq_min, mismatch_max,
                            rlen_min, is_mated)
        
        rdSet = ReadSet(reads["reads1u"])
        rdSet.get_loc_idx(_exons, _strand)
        count = rdSet.loc_idx.sum(axis=0)
        RPK   = rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads2u"])
        rdSet.get_loc_idx(_exons, _strand)
        count += rdSet.loc_idx.sum(axis=0)
        RPK   += rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads1"], reads["reads2"])
        rdSet.get_loc_idx(_exons, _strand)
        count += rdSet.loc_idx.sum(axis=0)
        RPK   += rdSet.RPK_use.sum(axis=0)

        a_line = gene_list[g] + 
        a_line += "\t".join(["%d" %num for num in list(count[g,:])]) + "\t" 
        a_line += "\t".join(["%.3f" %num for num in list(RPK[g,:])]) + "\n"
        fid.writelines(a_line)
    fid.close()
