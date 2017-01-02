# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import os
import sys
import gzip
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from .utils.gtf_utils import loadgene
from .utils.run_utils import get_psi, sort_dice_file

FID1 = None 
FID2 = None
PROCESSED = 0
TOTAL_GENE = 0
TOTAL_READ = []
START_TIME = time.time()
    
def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME, FID1, FID2
    if RV is None: 
        return RV
    else:
        FID1.writelines(RV["dice_line"])
        if FID2 is not None: FID2.writelines(RV["sample_line"])
    
    PROCESSED += 1
    bar_len = 20
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[DICEseq] [%s] %.1f%% done in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()
    return RV


def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts in GTF or GFF3")
    parser.add_option("--sam_list", "-s", dest="sam_list", default=None,
        help=("Sorted and indexed bam/sam files, use ',' for replicates "
        "and '---' for time points, e.g., T1_rep1.bam,T1_rep2.bam---T2.bam"))
    parser.add_option("--time_seq", "-t", dest="time_seq", default=None,
        help="The time for the input samples [Default: 0,1,2,...]")
    parser.add_option("--out_file", "-o", dest="out_file", default="output", 
        help="Prefix of the output files with full path")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")    
    group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
        default=False, help="Add the pre-mRNA as a transcript")
    
    group.add_option("--fLen", type="float", nargs=2, dest="frag_leng",
        default=[None,None], help=("Two arguments for fragment length: "
        "mean and standard diveation, default: auto-detected"))
    group.add_option("--bias", nargs=3, dest="bias_args",
        default=["unif","None","None"], help=("Three argments for bias "
        "correction: BIAS_MODE,REF_FILE,BIAS_FILE(s). BIAS_MODE: unif, end5, "
        "end3, both. REF_FILE: the genome reference file in fasta format. "
        "BIAS_FILE(s): bias files from dice-bias, use '---' for time specific "
        "files, [default: unif None None]"))

    group.add_option("--thetas", nargs=2, dest="thetas", default=[3,"None"], 
        help=("Two arguments for hyperparameters in GP model: theta1,theta2. "
        "default: [3 None], where theta2 covers 1/3 duration."))
    group.add_option("--mcmc", type="int", nargs=4, dest="mcmc_run",
        default=[0,20000,1000,100], help=("Four arguments for in MCMC "
        "iterations: save_sample,max_run,min_run,gap_run. Required: "
        "save_sample =< 3/4*mim_run. [default: 0 20000 1000 100]"))
    # SAVE_NUM: the number of samples for saving out; 
    # MAX_NUM,MIN_NUM: the maximum and the minmum samples;
    # GAP_NUM: after min_num, the gap_run added till convergency.

    # group.add_option("--anno_type", dest="anno_type", default="GTF",
    #     help="Type of annotation file: GTF, GFF3, UCSC_table "
    #     "[default: %default]")

    parser.add_option_group(group)


    ##### FOR DEVELOPMENT #####
    # parser.add_option("--mate_mode", dest="mate_mode", default="pair",
    #     help=("The mode for using paired-end reads: auto, pair, single "
    #     "[default: %default]."))
    # parser.add_option("--auto_min", dest="auto_min", default="200",
    #     help=("The minimum pairs of read mates in auto mode. "
    #     "[default: %default]."))

    # parser.add_option("--print_detail", action="store_true", dest="print_detail", 
    #     default=False, help="print the detail of the sampling.")
    # parser.add_option("--no_twice", action="store_true", dest="no_twice", 
    #     default=False, help="No quick estimate of the variance, but use fixed.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to diceseq!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[DICEseq] Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[DICEseq] loading annotation file...")
        sys.stdout.flush()    
        # anno = load_annotation(options.anno_file, options.anno_type)
        # genes = anno["genes"]
        genes = loadgene(options.anno_file)
        sys.stdout.write("\r[DICEseq] loading annotation file... Done.\n")
        sys.stdout.flush()
        
        global TOTAL_GENE
        TOTAL_GENE = len(genes)

    if options.sam_list == None:
        print("[DICEseq] Error: need --sam_list for aliged & indexed reads.")
        sys.exit(1)
    else:
        sam_list = options.sam_list.split("---")
        global TOTAL_READ
        for i in range(len(sam_list)):
            sam_list[i] = sam_list[i].split(",")
            _cnt = 0
            for ss in sam_list[i]:
                if not os.path.isfile(ss):
                    print("Error: No such file\n    -- %s" %ss)
                    sys.exit(1)
                pysam_stats = pysam.idxstats(ss)
                if type(pysam_stats) is not list:
                    pysam_stats = pysam_stats.split("\n")
                for tp in pysam_stats: 
                    tmp = tp.strip().split("\t")
                    if len(tmp) >= 3:
                        _cnt += float(tmp[2])
            TOTAL_READ.append(_cnt)

    no_twice = False
    auto_min = 200
    mate_mode = "pair"
    print_detail = False

    nproc = options.nproc
    out_file = options.out_file
    add_premRNA = options.add_premRNA
    FLmean, FLstd = options.frag_leng
    sample_num, Mmax, Mmin, Mgap = options.mcmc_run
    
    if options.time_seq is None: 
        X = np.arange(len(sam_list))
    else: 
        X = np.array(options.time_seq.split(","), "float")

    theta1, theta2 = options.thetas
    theta1 = float(theta1)
    if theta2 is None or ["None", "Auto", "auto"].count(theta2) == 1:
        theta2 = ((max(X) - min(X) + 0.1) / 3.0)**2
    elif theta2 == 'learn':
        theta2 = None
    else:
        theta2 = max(0.00001, float(theta2))

    bias_mode, ref_file, bias_file = options.bias_args
    if bias_mode == "unif":
        ref_file = None
        bias_file = None
    elif ref_file is "None": 
        ref_file = None
        bias_file = None
        bias_mode = "unif"
        print("[DICEseq] No reference sequence, change to uniform mode.")
    elif bias_file is "None":
        ref_file = None
        bias_file = None
        bias_mode = "unif"
        print("[DICEseq] No bias parameter file, change to uniform mode.")
    else:
        bias_file = bias_file.split("---")


    global FID1, FID2
    # if not os.path.exists(os.path.dirname(out_file)):
    #     try:
    #         os.makedirs(os.path.dirname(out_file))
    #     except OSError as exc: # Guard against race condition
    #         if exc.errno != errno.EEXIST:
    #             raise
    FID1 = open(out_file + ".dice", "w")
    headline = "tran_id\tgene_id\tlogLik\ttransLen"
    for i in range(len(X)):
        _t = str(X[i])
        headline += "\tFPKM_T%s\tratio_T%s\tratio_lo_T%s\tratio_hi_T%s" %(_t,
            _t, _t, _t)
    FID1.writelines(headline + "\n")
    if sample_num > 0:
        FID2 = gzip.open(out_file + ".sample.gz", "w")
        FID2.writelines("# MCMC samples for latent Y\n")
        FID2.writelines("# @gene|transcripts|theta2\n")
        FID2.writelines("# y_c1t1,y_c2t1;y_c1t2,y_c2t2;...\n")

    print("[DICEseq] running diceseq for %d genes with %d cores..." %(
        TOTAL_GENE, nproc))

    tran_ids = []
    if nproc <= 1:
        for g in genes:
            if add_premRNA: g.add_premRNA()
            for t in g.trans: tran_ids.append(t.tranID)
            RV = get_psi(g, sam_list, ref_file,  bias_file, bias_mode, X, Mmax, 
                 Mmin, Mgap, theta1, theta2, no_twice, sample_num, print_detail,
                 FLmean, FLstd, mate_mode, auto_min, TOTAL_READ)
            show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        for g in genes:
            if add_premRNA: g.add_premRNA()
            for t in g.trans: tran_ids.append(t.tranID)
            pool.apply_async(get_psi, (g, sam_list, ref_file,  bias_file, 
                bias_mode, X, Mmax, Mmin, Mgap, theta1, theta2, no_twice, 
                sample_num, print_detail, FLmean, FLstd, mate_mode, 
                auto_min, TOTAL_READ), callback=show_progress)
        pool.close()
        pool.join()
    FID1.close()
    if FID2 is not None: FID2.close()
    
    sort_dice_file(out_file+".dice", tran_ids)
    print("")


if __name__ == "__main__":
    main()
    