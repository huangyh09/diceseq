# This file is to esitmate the sequecne and positon bias parameters.

import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from .utils.tran_utils import TranUnits
from .utils.gtf_utils import loadgene #load_annotation
from .utils.bias_utils import BiasFile, FastaFile
from .utils.sam_utils import load_samfile, fetch_reads

PROCESSED = 0
USED_GENE = 0
TOTAL_GENE = 0
TOTAL_READ = 0
BIAS_FILE = BiasFile()
START_TIME = time.time()

def show_progress(RV):
    global PROCESSED, BIAS_FILE, TOTAL_READ, USED_GENE

    PROCESSED += 1
    bar_len = 20
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[dice-bias] [%s] %.1f%% done in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()

    if RV is None: return None
    USED_GENE += 1
    TOTAL_READ += RV["BF"].read_num
    BIAS_FILE.add_bias_file(RV["BF"])
    return RV

def get_bias(BF, g, ref_seq, samFile, threshold=0.01):
    """to get the bias parameters from a transcript `t`."""
    samFile = load_samfile(samFile)
    ref_seq = FastaFile(ref_seq)

    RV = {}
    reads = fetch_reads(samFile, g.chrom, g.start, g.stop, 
                        rm_duplicate=True, inner_only=True, mapq_min=10,
                        mismatch_max=10, rlen_min=1, is_mated=True)
    rcount = len(reads["reads1"])+len(reads["reads1u"])+len(reads["reads2u"])
    if rcount < threshold * g.trans[0].tranL:
        # print("Coverage is only RPK=%.1f on %s, skipped." 
        #     %(rcount * 1000.0 / g.trans[0].tranL, g.trans[0].tranID))
        #RV["BF"] = BF
        #RV["flen"] = np.array([])
        return None
    # print(g.geneID, g.chrom, g.strand, g.start, g.stop)
    # print(len(reads["reads1"]), len(reads["reads1u"]), len(reads["reads2u"]))

    t = TranUnits(g.trans[0])
    t.set_sequence(ref_seq)
    seqs = t.seq
    tLen = t.ulen
    flen = np.array([], "float")
    idx5 = np.array([], "float")
    idx3 = np.array([], "float")

    if len(reads["reads1"]) > 0:
        t.set_reads(reads["reads1"], reads["reads2"], "unif")
        idx5 = np.append(idx5, t.idx5)
        idx3 = np.append(idx3, t.idx3)
        flen = np.append(flen, t.flen[t.Rmat])
        #print(len(flen))

    if len(reads["reads1u"]) > 0:
        t.set_reads(reads["reads1u"], [], "unif")
        idx5 = np.append(idx5, t.idx5)
        # idx3 = np.append(idx3, t.idx3)
        # flen = np.append(flen, t.flen[t.Rmat])

    if len(reads["reads2u"]) >0:
        t.set_reads([], reads["reads2u"], "unif")
        idx3 = np.append(idx3, t.idx3)
        # idx5 = np.append(idx5, t.idx5)
        # flen = np.append(flen, t.flen[t.Rmat])

    _i5 = idx5 == idx5 #> -1000 #
    _i3 = idx3 == idx3 #> -1000 #
    idx5 = idx5[_i5]
    idx3 = idx3[_i3]

    for i in range(tLen):
        ipos = i + 20
        _seq5 = t.seq[ipos-8  : ipos+13]
        _seq3 = t.seq[ipos-12 : ipos+9 ][::-1]
        _pos5 = i
        _pos3 = i #tLen - i -1
        if t.strand != "+" and t.strand != "1":
            _seq5, _seq3 = _seq3, _seq5
            # _pos5, _pos3 = _pos3, _pos5
            _pos5, _pos3 = tLen-1-_pos3, tLen-1-_pos5
        if len(idx5) > 0 and np.mean(idx5==i) > 0:
            BF.set_both_bias(_seq5, _pos5, tLen, np.mean(idx5==i), 5, "bias")
        if len(idx3) > 0 and np.mean(idx3==i) > 0: 
            BF.set_both_bias(_seq3, _pos3, tLen, np.mean(idx3==i), 3, "bias")
        BF.set_both_bias(_seq5, _pos5, tLen, 1.0/tLen, 5, "unif")
        BF.set_both_bias(_seq3, _pos3, tLen, 1.0/tLen, 3, "unif")
    BF.read_num = len(flen)
    BF.flen_sum1 = np.sum(flen)
    BF.flen_sum2 = np.sum(flen**2)

    RV["BF"] = BF
    #RV["flen"] = flen
    return RV


def main():
    # import warnings
    # warnings.filterwarnings('error')

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file with single-transcript genes")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="Sorted and indexed bam/sam files")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="Genome reference sequences in FASTA format")
    parser.add_option("--out_file", "-o", dest="out_file",  
        default="output.bias", help="Output file in BIAS format")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    # group.add_option("--anno_type", dest="anno_type", default="GTF",
    #     help="Type of annotation file: GTF, GFF3, UCSC_table "
    #     "[default: %default]")
    group.add_option("--num_max", dest="num_max", default=None,
        help="The maximum number of genes for bias estimate [default: inf].")
    group.add_option("--percentile", type="int", nargs=4, dest="percentile",
        default=[500,800,1300,2200], help=("N-1 arguments for N percentiles "
        "of transcript length. [default: 500 800 1300 2200]"))
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to dice-bias!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[dice-bias] Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[dice-bias] loading annotation file...")
        sys.stdout.flush()    
        # anno = load_annotation(options.anno_file, options.anno_type)
        # genes = anno["genes"]
        genes = loadgene(options.anno_file)
        sys.stdout.write("\r[dice-bias] loading annotation file... Done.\n")
        sys.stdout.flush()
    if options.sam_file == None:
        print("[dice-bias] Error: need --sam_file for aliged & indexed reads.")
        sys.exit(1)
    else:
        sam_file = options.sam_file
    
    nproc = int(options.nproc)
    out_file = options.out_file
    ref_file = options.ref_file

    #part 1. transcript length
    global TOTAL_GENE
    tran_len_all = []
    for g in genes:
        for t in g.trans:
            tran_len_all.append(t.tranL)
    BF = BiasFile()
    BF.set_percentile(np.array(tran_len_all))
    BIAS_FILE.set_percentile(np.array(tran_len_all))
    TOTAL_GENE = np.sum(anno["tran_num"]==1)

    if options.num_max is None:
        num_max = TOTAL_GENE
    else:
        num_max = min(int(options.num_max), TOTAL_GENE)

    print("[dice-bias] isoform length: 0--%d--%d--%d--%d--inf" 
        %tuple(BF.percentile[1:,0]))
    print("[dice-bias] running dice-bias for %d one-isoform genes with %d cores..." 
        %(TOTAL_GENE, nproc))

    #part 2. estimate
    cnt = 0
    if nproc <= 1:
        for g in genes:
            if g.tranNum > 1: continue
            cnt += 1
            if cnt >= num_max: break
            RV = get_bias(BF, g, ref_file, sam_file, 0.01)
            show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        for g in genes:
            if g.tranNum > 1: continue
            cnt += 1
            if cnt >= num_max: break
            pool.apply_async(get_bias, (BF, g, ref_file, sam_file, 0.01), 
                callback=show_progress)
        pool.close()
        pool.join()

    BIAS_FILE.save_file(out_file)
    print("[dice-bias] \n%d reads from %d genes were used in bias estimate." 
        %(TOTAL_READ, USED_GENE))


if __name__ == "__main__":
    main()
