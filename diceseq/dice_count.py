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
from .utils.run_utils import sort_dice_file
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

def check_reads_end(reads, start, stop, filter=False, TSS_method="mean"):
    """
    Check whether reads partially mapped and get TSS and TTS

    Parameters
    ----------
    reads: dictionary
        "reads1u", "reads2u", "reads1", "reads2"
    start: int
        gene start site
    stop: int
        gene stop site
    filter: bool
        whether to filter reads too far from gene (for partially mapped only)
    TSS_method: string
        method to calculate the TSS/TTS: mean, max, 95
    """
    RV = {}
    cnt5, cnt3 = 0, 0
    dis5, dis3 = 0, 0
    dis5_all, dis3_all = [], []
    ReadsAll = [[reads["reads1u"], []],
                [reads["reads2u"], []],
                [reads["reads1" ], reads["reads2"]]]

    for s in range(3):
        read1_out = []
        read2_out = []
        reads_set = ReadsAll[s]
        for i in range(len(reads_set[0])):
            r1 = reads_set[0][i]
            r1len = r1.qlen
            _gap1 = start - r1.pos
            _gap2 = r1.aend - stop

            if s < 2:
                r2len =  0
                _gap3 = -1
                _gap4 = -1
            else:
                r2 = reads_set[1][i]
                r2len = r2.qlen
                _gap3 = start - r2.pos
                _gap4 = r2.aend - stop              

            # flitering
            if filter is True:
                if max(_gap1, _gap2) > r1len or max(_gap3, _gap4) > r2len:
                    continue
                read1_out.append(r1)
                if s == 2: 
                    read2_out.append(r2)

            # TSS gap and TTS gap
            if max(_gap1, _gap3) > 0:
                cnt5 += 1
                dis5 += max(_gap1, _gap3)
                dis5_all.append(max(_gap1, _gap3))
            if max(_gap2, _gap4) > 0:
                cnt3 += 1
                dis3 += max(_gap2, _gap4)
                dis3_all.append(max(_gap2, _gap4))

        # flitering
        if filter is True:
            if s == 0:
                RV["reads1u"] = read1_out
            elif s == 1: 
                RV["reads2u"] = read1_out
            elif s == 2: 
                RV["reads1"] = read1_out
                RV["reads2"] = read2_out

    # if cnt5 > 0: dis5 = dis5 / (cnt5+0.0)
    # if cnt3 > 0: dis3 = dis3 / (cnt3+0.0)
    if len(dis5_all) > 3:
        cnt = len(dis5_all)
        (values,counts) = np.unique(dis5_all, return_counts=True)
        idx = np.argmax(counts)
        dis5_max = values[idx]
        dis5_avg = np.mean(dis5_all)
        dis5_C95 = np.sort(dis5_all)[int(cnt*0.95)]
        idx = np.argsort(dis5_all)[int(cnt*0.025):int(cnt*0.975)]
        dis5_m95 = np.mean(np.array(dis5_all)[idx])
    else:
        dis5_max = 0.0
        dis5_avg = 0.0
        dis5_C95 = 0.0
        dis5_m95 = 0.0
    if len(dis3_all) > 3:
        cnt = len(dis3_all)
        (values,counts) = np.unique(dis3_all, return_counts=True)
        idx = np.argmax(counts)
        dis3_max = values[idx]
        dis3_avg = np.mean(dis3_all)
        dis3_C95 = np.sort(dis3_all)[int(len(dis3_all)*0.95)]
        idx = np.argsort(dis3_all)[int(cnt*0.025):int(cnt*0.975)]
        dis3_m95 = np.mean(np.array(dis3_all)[idx])
    else:
        dis3_max = 0.0
        dis3_avg = 0.0
        dis3_C95 = 0.0
        dis3_m95 = 0.0
        
    if filter is False: 
        RV = reads
    if TSS_method == "peak":
        dis5, dis3 = dis5_max, dis3_max
    elif TSS_method == "edge":
        dis5, dis3 = dis5_C95, dis3_C95
    else:
        dis5, dis3 = dis5_avg, dis3_avg
    
    return RV, dis5, dis3


def get_count(gene, sam_file, total_reads, rm_duplicate, inner_only, mapq_min, 
    mismatch_max, rlen_min, is_mated, total_only, overhang, gap5, gap3, 
    TSS_method, show_gap):
    """Get reads counts

    For inner only: 
        all reads from start to stop, considering gap5 and gap3
    For partial: 
        same as inner only, excluding reads whole outside genebody
    """
    # adjust TSS and TTS
    if gene.strand == "-":
        gap5, gap3 = gap3, gap5

    # fetch reads
    samFile = load_samfile(sam_file)
    if inner_only == True:
        reads = fetch_reads(samFile, gene.chrom, gene.start-gap5, gene.stop+gap3,
            rm_duplicate, inner_only, mapq_min, mismatch_max, rlen_min, 
            is_mated)
        if gap5 == 0 and gap3 == 0:
            dis5, dis3 = 0, 0
        else:
            reads, dis5, dis3 = check_reads_end(reads, gene.start, gene.stop, 
                                                False, TSS_method)
    else:
        reads = fetch_reads(samFile, gene.chrom, gene.start, gene.stop, 
            rm_duplicate, inner_only, mapq_min, mismatch_max, rlen_min, 
            is_mated)
        reads, dis5, dis3 = check_reads_end(reads, gene.start, gene.stop, True, 
                                            TSS_method)

    if total_only:
        count = [len(reads["reads1u"])+len(reads["reads2u"])+len(reads["reads1"])]
        RPKM  = [count[0] * 10**9 / abs(gene.start - gene.stop) / total_reads]
    elif gene.tranNum == 1 and gene.trans[0].exonNum == 2:
        exons = gene.trans[0].exons + 0
        exons[0,  0] -= gap5
        exons[-1,-1] += gap3

        rdSet = ReadSet(reads["reads1u"])
        rdSet.get_loc_idx(exons, gene.strand, overhang)
        count = rdSet.loc_idx.sum(axis=0)
        RPKM  = rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads2u"])
        rdSet.get_loc_idx(exons, gene.strand, overhang)
        count += rdSet.loc_idx.sum(axis=0)
        RPKM  += rdSet.RPK_use.sum(axis=0)

        rdSet = ReadSet(reads["reads1"], reads["reads2"])
        rdSet.get_loc_idx(exons, gene.strand, overhang)
        count += rdSet.loc_idx.sum(axis=0)
        RPKM  += rdSet.RPK_use.sum(axis=0)

        RPKM = RPKM * 10**6 / total_reads
    else: return None

    gLen = str(abs(gene.stop - gene.start) + 1)
    a_line =  "\t".join([gene.geneID, gene.geneName, gene.biotype, gLen]) + "\t"
    a_line += "\t".join(["%d" %num for num in list(count)]) + "\t" 
    a_line += "\t".join(["%.2e" %num for num in list(RPKM)])
    # if gap5 != 0 or gap3 != 0 or inner_only == False:
    if show_gap:
        if gene.strand == "-": 
            dis5, dis3 = dis3, dis5
        a_line += "\t%.2f\t%.2f" %(dis5, dis3)

    return a_line + "\n"

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
    group.add_option("--junction", action="store_true", dest="junction_reads",
        default=False, help=("return junction and boundary reads, only for "
        "gene with one exon-intron-exon structure; otherwise no junction."))
    parser.add_option_group(group)

    group = OptionGroup(parser, "Arguments for reads quality")
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
    group.add_option("--single_end", action="store_true", dest="single_end",
        default=False, help="use reads as single-end; otherwise paired-end")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Arguments for TSS and TTS ")
    group.add_option("--partial", action="store_true", dest="partial_use",
        default=False, help="keep reads partial in the region; otherwise not")
    group.add_option("--gap5", "-5", dest="gap5", type="int", default=0,
        help="Extended length on 5-end of the gene, --partial mode free. "
        "[default: %default]")
    group.add_option("--gap3", "-3", dest="gap3", type="int", default=0,
        help="Extended length on 3-end of the gene, --partial mode free. "
        "[default: %default]")
    group.add_option("--gapFile", dest="gap_file", default=None,
        help="File for gene-id, gap5 and gap3. Replace --gap5 and --gap3.")
    group.add_option("--TSSmethod", dest="TSS_method", default="mean",
        help="Method to estimate TSS and TTS from outside reads: 1. mean (mean "
        "reads density), 2. peak (maximum of reads density) or 3. edge (edge "
        "of 95%% reads) [default: %default]")
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

    gap3 = options.gap3
    gap5 = options.gap5
    nproc = int(options.nproc)
    overhang = int(options.overhang)
    mapq_min = int(options.mapq_min)
    rlen_min = int(options.rlen_min)
    TSS_method = options.TSS_method
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

    gap_file = options.gap_file
    if gap_file is not None:
        gap_info = np.genfromtxt(gap_file, dtype="str", skip_header=1)


    if gap5 != 0 or gap3 != 0 or inner_only == False or gap_file is not None:
        show_gap = True
        head_line += "\tTSS_gap\tTTS_gap"
    else:
        show_gap = False

    FID.writelines(head_line + "\n")

    print("[dice-count] running dice-count for %d genes with %d cores..." 
        %(TOTAL_GENE, nproc))

    gene_ids = []
    if nproc <= 1:
        for g in genes:
            if total_only == False and (g.tranNum > 1 or g.trans[0].exonNum!=2):
                continue
            if gap_file is not None:
                ii = np.where(g.geneID == gap_info[:,0])[0]
                if len(ii) >= 1:
                    gap5 = int(float(gap_info[ii[0], 1]))#
                    gap3 = int(float(gap_info[ii[0], 2]))

            RV = get_count(g, sam_file, total_reads, rm_duplicate, inner_only, 
                mapq_min, mismatch_max, rlen_min, is_mated, total_only, 
                overhang, gap5, gap3, TSS_method, show_gap)
            show_progress(RV)
            gene_ids.append(g.geneID)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        for g in genes:
            if total_only == False and (g.tranNum > 1 or g.trans[0].exonNum!=2):
                continue
            if gap_file is not None:
                ii = np.where(g.geneID == gap_info[:,0])[0]
                if len(ii) >= 1:
                    gap5 = int(float(gap_info[ii[0], 1]))#*0.5
                    gap3 = int(float(gap_info[ii[0], 2]))

            pool.apply_async(get_count, (g, sam_file, total_reads, rm_duplicate, 
                inner_only, mapq_min, mismatch_max, rlen_min, is_mated, 
                total_only, overhang, gap5, gap3, TSS_method, show_gap),  
                callback=show_progress)
            gene_ids.append(g.geneID)
        pool.close()
        pool.join()
        
    FID.close()
    
    sort_dice_file(options.out_file, gene_ids)
    print("")


if __name__ == "__main__":
    main()
