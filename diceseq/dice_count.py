# This is a direct running file to calculate the specific counts of RNA
# splicing, and currently it is mainly for intronic splicing, rather than
# alternative splicing. There are 8 types for reads for an exon1-intron-
# exon2 structure: (0) exon1, (1) exon1-intron boundary, (2) intron, (3) 
# intron-exon2 boundary, (4) exon2, (5) exon1-exon2 junction, (6) exon1-
# intron-exon2, (7) exon1-exon2 unsure. In addition, it also provides the
# normalized counts RPKM (reads per kilo-base per million reads).

import sys
import time
import pysam
import subprocess
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from .utils.reads_utils import ReadSet
from .utils.gtf_utils import loadgene #load_annotation
from .utils.sam_utils import load_samfile, fetch_reads

FID = None 
PROCESSED = 0
TOTAL_GENE = 0
START_TIME = time.time()

def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME, FID
    if RV is None:
        return RV
    else:
        FID.writelines(RV)
    
    PROCESSED += 1
    bar_len = 20
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[dice-count] [%s] %.1f%% done in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()
    return RV

def get_count(gene, sam_file, total_reads, rm_duplicate, inner_only, mapq_min, 
    mismatch_max, rlen_min, is_mated, total_only, overhang):
    samFile = load_samfile(sam_file)
    reads = fetch_reads(samFile, gene.chrom, gene.start, gene.stop, 
        rm_duplicate, inner_only, mapq_min, mismatch_max, rlen_min, is_mated)

    if total_only:
        count = [len(reads["reads1u"])+len(reads["reads2u"])+len(reads["reads1"])]
        RPKM  = [count[0] * 10**9 / abs(gene.start - gene.stop) / total_reads]
    elif gene.tranNum == 1 and gene.trans[0].exonNum == 2:
        rdSet = ReadSet(reads["reads1u"])
        rdSet.get_loc_idx(gene.trans[0].exons, gene.strand, overhang)
        count = rdSet.loc_idx.sum(axis=0)
        RPKM  = rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads2u"])
        rdSet.get_loc_idx(gene.trans[0].exons, gene.strand, overhang)
        count += rdSet.loc_idx.sum(axis=0)
        RPKM  += rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads1"], reads["reads2"])
        rdSet.get_loc_idx(gene.trans[0].exons, gene.strand, overhang)
        count += rdSet.loc_idx.sum(axis=0)
        RPKM  += rdSet.RPK_use.sum(axis=0)

        RPKM = RPKM * 10**6 / total_reads
    else: return None

    gLen = str(abs(gene.stop - gene.start) + 1)
    a_line =  "\t".join([gene.geneID, gene.geneName, gene.biotype, gLen]) + "\t"
    a_line += "\t".join(["%d" %num for num in list(count)]) + "\t" 
    a_line += "\t".join(["%.2e" %num for num in list(RPKM)]) + "\n"
    return a_line

def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="Sorted and indexed bam/sam files")
    parser.add_option("--out_file", "-o", dest="out_file",  
        default="dice_count.tsv", help="The counts in tsv file")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    # group.add_option("--anno_type", dest="anno_type", default="GTF",
    #     help="Type of annotation file: GTF, GFF3, UCSC_table "
    #     "[default: %default]")

    group.add_option("--mapq_min", dest="mapq_min", default="10", 
        help="Minimum mapq for reads. [default: %default]")
    group.add_option("--mismatch_max", dest="mismatch_max", default="5", 
        help="Maximum mismatch for reads. [default: %default]")
    group.add_option("--rlen_min", dest="rlen_min", default="1", 
        help="Minimum length for reads. [default: %default]")
    group.add_option("--overhang", dest="overhang", default="1", 
        help="Minimum overhang on junctions. [default: %default]")

    group.add_option("--duplicate", action="store_true", 
        dest="duplicate_use", default=False, help="keep duplicate reads; "
        "otherwise not")
    group.add_option("--partial", action="store_true", dest="partial_use",
        default=False, help="keep reads partial in the region; otherwise not")
    group.add_option("--single_end", action="store_true", dest="single_end",
        default=False, help="use reads as single-end; otherwise paired-end")
    group.add_option("--junction", action="store_true", dest="junction_reads",
        default=False, help=("return junction and boundary reads, only for "
        "gene with one exon-intron-exon structure; otherwise no junction."))
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to dice-count!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[dice-count] Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[dice-count] loading annotation file...")
        sys.stdout.flush()    
        # anno = load_annotation(options.anno_file, options.anno_type)
        # genes = anno["genes"]
        genes = loadgene(options.anno_file)
        sys.stdout.write("\r[dice-count] loading annotation file... Done.\n")
        sys.stdout.flush()
    if options.sam_file == None:
        print("[dice-count] Error: need --sam_file for aliged & indexed reads.")
        sys.exit(1)
    else:
        sam_file = options.sam_file

    total_reads = 0
    pysam_stats = pysam.idxstats(sam_file)
    if type(pysam_stats) is not list:
        pysam_stats = pysam_stats.split("\n")
    for tp in pysam_stats: 
        tmp = tp.strip().split("\t")
        if len(tmp) >= 3:
            total_reads += float(tp.strip().split("\t")[2])

    nproc = int(options.nproc)
    overhang = int(options.overhang)
    mapq_min = int(options.mapq_min)
    rlen_min = int(options.rlen_min)
    mismatch_max = int(options.mismatch_max)

    is_mated = (options.single_end == False)
    inner_only = (options.partial_use == False)
    total_only = (options.junction_reads == False)
    rm_duplicate = (options.duplicate_use == False)
    
    global FID, TOTAL_GENE
    FID = open(options.out_file, "w")
    if total_only == False:
        head_line = "gene_id\tgene_name\tbiotype\tgene_length"
        cnt_name = ["ex1", "ex1_int", "int", "int_ex2", "ex2", "ex1_ex2_junc",
            "ex1_int_ex2", "ex1_ex2_vague"]
        for i in range(len(cnt_name)):
            head_line += "\t" + cnt_name[i] + "_NUM"
        for i in range(len(cnt_name)):
            head_line += "\t" + cnt_name[i] + "_FPKM"
        cnt = 0
        for g in genes:
            if g.tranNum == 1 and g.trans[0].exonNum == 2: cnt += 1
        TOTAL_GENE = cnt
    else:
        TOTAL_GENE = len(genes)
        head_line = "gene_id\tgene_name\tbiotype\tgene_length\tcount\tFPKM"
    FID.writelines(head_line + "\n")

    print("[dice-count] running dice-count for %d genes with %d cores..." %(TOTAL_GENE, nproc))

    if nproc <= 1:
        for g in genes:
            if total_only == False and (g.tranNum > 1 or g.trans[0].exonNum != 2):
                continue
            RV = get_count(g, sam_file, total_reads, rm_duplicate, inner_only, 
                mapq_min, mismatch_max, rlen_min, is_mated, total_only, overhang)
            show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        for g in genes:
            if total_only == False and (g.tranNum > 1 or g.trans[0].exonNum != 2):
                continue
            pool.apply_async(get_count, (g, sam_file, total_reads, rm_duplicate, 
                inner_only, mapq_min, mismatch_max, rlen_min, is_mated, total_only,
                overhang),  callback=show_progress)
        pool.close()
        pool.join()
    FID.close()
    print("")


if __name__ == "__main__":
    main()
