# This file is to esitmate the sequecne and positon bias parameters.

import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser
from .utils.tran_utils import TranUnits
from .utils.gtf_utils import load_annotation
from .utils.bias_utils import BiasFile, FastaFile
from .utils.sam_utils import load_samfile, fetch_reads


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
        print("Coverage is only RPK=%.1f on %s, skipped." 
            %(rcount * 1000.0 / g.trans[0].tranL, g.trans[0].tranID))
        RV["BF"] = BF
        RV["flen"] = np.array([])
        return RV
    else:
        print(g.geneID, g.chrom, g.strand, g.start, g.stop)
        print(len(reads["reads1"]), len(reads["reads1u"]), len(reads["reads2u"]))

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

    if len(reads["reads1u"]) > 0:
        t.set_reads(reads["reads1u"], [], "unif")
        idx5 = np.append(idx5, t.idx5)
        # idx3 = np.append(idx3, t.idx3)
        # flen = np.append(flen, t.flen[t.Rmat])

    if len(reads["reads2u"]) >0:
        t.set_reads([], reads["reads2u"], "unif")
        # idx5 = np.append(idx5, t.idx5)
        idx3 = np.append(idx3, t.idx3)
        # flen = np.append(flen, t.flen[t.Rmat])
    
    _i5 = idx5 == idx5
    _i3 = idx3 == idx3
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

    RV["BF"] = BF
    RV["flen"] = flen
    return RV

def bias_add(BF, BFnew):
    BF.pos5_bias += BFnew.pos5_bias
    BF.pos3_bias += BFnew.pos3_bias
    BF.pos5_unif += BFnew.pos5_unif
    BF.pos3_unif += BFnew.pos3_unif
    for i in range(len(BF.chain_len)):
        BF.seq5_bias[str(i)] += BFnew.seq5_bias[str(i)]
        BF.seq3_bias[str(i)] += BFnew.seq3_bias[str(i)]
        BF.seq5_unif[str(i)] += BFnew.seq5_unif[str(i)]
        BF.seq3_unif[str(i)] += BFnew.seq3_unif[str(i)]
    return BF

def main():
    print("Welcome to dice-bias!")
    import warnings
    warnings.filterwarnings('error')

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format.")
    parser.add_option("--anno_source", dest="anno_source", default="Ensembl",
        help="The annotation source of the gtf file [default: %default]")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="The sourted & indexed bam/sam file.")
    parser.add_option("--refseq_file", "-r", dest="refseq_file", default=None,
        help="The fasta file of genome reference sequences.")
    parser.add_option("--nproc", dest="nproc", default="4",
        help="The number of subprocesses [default: %default].")
    parser.add_option("--out_file", "-o", dest="out_file", default="out.bias",
        help="The output file in BIAS format.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file, options.anno_source)
        genes = anno["genes"]
    if options.sam_file == None:
        print("Error: need --sam_file for reads indexed and aliged reads.")
        sys.exit(1)
    else:
        # samFile = load_samfile(options.sam_file)
        sam_file = options.sam_file
    
    nproc = int(options.nproc)
    biasFile = BiasFile()
    out_file = options.out_file
    refseq_file = options.refseq_file
    # fastaFile = FastaFile(options.refseq_file)

    #part 1.1. all transcript length
    
    tran_len_all = []
    for g in genes:
        for t in g.trans:
            tran_len_all.append(t.tranL)
    biasFile.set_percentile(np.array(tran_len_all))
    print(biasFile.percentile)
    print("%.0f out of %.0f genes have one isoform for bias parameter estimate."
          %(np.sum(anno["tran_num"]==1), len(anno["tran_num"])))

    start_time = time.time()
    flen_all = np.array([],"float")
    pool = multiprocessing.Pool(processes=nproc)
    result = []
    cnt = 0
    for g in genes:
        if g.tranNum > 1: continue
        # if cnt == 100: break
        cnt += 1

        result.append(pool.apply_async(get_bias, (biasFile, g, refseq_file, sam_file, 0.01, )))
    pool.close()
    pool.join()
    cnt = 0
    for res in result:
        RV = res.get()
        if len(RV["flen"]) == 0: continue
        cnt += 1
        flen_all = np.append(flen_all, RV["flen"])
        biasFile = bias_add(biasFile, RV["BF"])
    
    # result = []
    # for g in genes:
    #     if g.tranNum > 1: continue
    #     if cnt == 100: break
    #     cnt += 1
    #     result.append((get_bias(biasFile, g, refseq_file, sam_file, 0.01)))
    # for RV in result:
    #     flen_all = np.append(flen_all, RV["flen"])
    #     biasFile = bias_add(biasFile, RV["BF"])

    if len(flen_all) == 0:
        biasFile.flen_mean = 0
        biasFile.flen_std = 0
    else:
        biasFile.flen_mean = np.mean(flen_all)
        biasFile.flen_std = np.std(flen_all)
    biasFile.save_file(out_file)

    print("%d effect genes" %cnt)
    print("--- %.3f seconds ---" % (time.time() - start_time))
    

if __name__ == "__main__":
    main()
