# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import sys
import gzip
import time
import subprocess
import numpy as np
import multiprocessing
from optparse import OptionParser
from .models.model_GP import Psi_GP_MH
from .utils.sam_utils import load_samfile
from .utils.gtf_utils import load_annotation
from .utils.bias_utils import BiasFile, FastaFile
from .utils.tran_utils import TranUnits, TranSplice

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
    bar_len = 30
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[%s] %.2f%% processed in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()
    return RV

def get_CI(data, percent=0.95):
    if len(data.shape) == 0:
        data = data.reshape(-1,1)
    RV = np.zeros((data.shape[1],2))
    CI_idx = int(data.shape[0] * (1-percent)/2)
    for k in range(data.shape[1]):
        temp = np.sort(data[:,k])
        RV[k,:] = [temp[-CI_idx], temp[CI_idx]]
    return RV

def harmonic_mean_log(data, log_in=True, log_out=True):
    if log_in is not True:
        data = np.log(data)
    RV = np.log(len(data)) + np.min(data, axis=0)
    RV = RV - np.log(np.sum(np.exp(np.min(data, axis=0) - data), axis=0))
    if log_out is not True: RV = np.exp(RV)
    return RV

def Psi2Y(Psi):
    Psi = Psi / np.sum(Psi)
    Y = np.zeros(len(Psi))
    Y[:-1] = np.log(Psi[:-1] / Psi[-1])
    return Y

def get_psi(gene, sam_list, ref_file,  biasFile, bias_mode, X, M, initial,
            gap, theta1, theta2, no_twice, sample_num, pout):
    samFiles = []
    for s in sam_list:
        _sam = []
        ss_list = s.split(",")
        for ss in ss_list:
            _sam.append(load_samfile(ss))
        samFiles.append(_sam)
    if ref_file is not None: fastaFile = FastaFile(ref_file)

    R_all, len_iso_all, prob_iso_all, _count = [], [], [], []
    for j in range(len(samFiles)):
        t = TranSplice(gene)
        for _sam in samFiles[j]:
            t.set_reads(_sam)
        if bias_mode != "unif":
            t.set_sequence(fastaFile)
            t.set_bias(biasFile)

        t.get_ready(bias_mode, 200, 20)
        R_all.append(t.Rmat)
        if bias_mode == "unif": 
            len_iso_all.append(t.efflen_unif)
            prob_iso_all.append(t.proU)
        else: 
            len_iso_all.append(t.efflen_bias)
            prob_iso_all.append(t.proB)
        # print(len(t.read1p), len(t.read1u), len(t.read2u))
        # _count.append(len(t.read1p) + len(t.read1u) + len(t.read2u))
        _count.append(np.sum(t.Rmat.sum(axis=1)>0))
    count = np.array(_count)
    
    print_info = ("%s: %d transcript, %d reads." %
        (gene.geneID, gene.tranNum, np.sum(count)))
    # if pout: print("\n" + print_info)#

    if R_all[0].shape[1] == 1:
        psi_mean = np.ones((1,len(X)))
        psi_95 = [np.ones((len(X), 2))]
        lik_marg = 0.0
        theta_mean = None
        sample_all = []

    elif R_all[0].shape[1] > 1:
        var = 0.1 * np.ones(R_all[0].shape[1] - 1)
        Ymean = np.zeros((R_all[0].shape[1], len(X)))
        if no_twice == False:
            _var = np.zeros(R_all[0].shape[1]-1)
            for j in range(len(R_all)):
                _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all[j:j+1],
                    len_iso_all[j:j+1], prob_iso_all[j:j+1], X[j:j+1], 
                    Ymean[:,j:j+1], var, theta1, theta2, 100, 100, 50)
                _var += np.var(_Y[int(_m/4):, :-1, 0], axis=0)
            var = (_var + 0.000001) / len(R_all)
        # print(var)

        _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all, len_iso_all, 
            prob_iso_all, X, Ymean, var, theta1, theta2, M, initial, gap)
        print_info += (" %d acceptances in %d iterations." %(_cnt, _m))

        psi_95 = []
        for c in range(_Psi.shape[1]):
            psi_95.append(get_CI(_Psi[int(_m/4):, c, :], 0.95))
        psi_mean = _Psi[int(_m/4):, :, :].mean(axis=0)
        lik_marg = harmonic_mean_log(_Lik[int(_m/4):])
        theta_mean = _theta[int(_m/4):,:].mean(axis=0)
        sample_all = _Y[-sample_num:,:-1,:]
    if pout: print("\n" + print_info)

    # for output
    gene_info = gene.get_gene_info()
    transLen = []
    for tr in gene.trans:
        transLen.append(tr.tranL)
    transLen = np.array(transLen)

    diceL = ""
    for c in range(len(gene.trans)):
        _line = "%s\t%s\t%.1e\t%d" %(gene.trans[c].tranID, gene.geneID,
            lik_marg, gene.trans[c].tranL)
        for t in range(len(X)):
            fsi = (transLen * psi_mean[:,t]) / np.sum(transLen * psi_mean[:,t])
            FPKM = count[t] * fsi[c] * 10**9 / transLen[c] / TOTAL_READ[t]
            _line += "\t%.1f\t%.3f\t%.3f\t%.3f" %(FPKM, psi_mean[c,t],
                psi_95[c][t,1], psi_95[c][t,0])
        diceL += _line + "\n"

    _theta2 = []
    for tt in range(len(theta_mean)):
        _theta2.append("%.2e" %_theta[int(_m/4):,tt].mean())
    sampL = "@%s|%s|%s\n" %(gene_info[0], gene_info[-1], ",".join(_theta2))
    for ss in range(sample_num):
        for tt in range(_Y.shape[2]):
            _line_time = []
            for cc in range(_Y.shape[1]):
                _line_time.append("%.2e" %_Y[int(_m/2)+ss, cc, tt])
            sampL += ",".join(_line_time)
            if tt < _Y.shape[2] - 1: 
                sampL += ";"
            else:
                sampL += "\n"
    # show_progress()
    
    RV = {}
    RV["dice_line"] = diceL
    RV["sample_line"] = sampL
    return RV


def main():
    # import warnings
    # warnings.filterwarnings('error')
    print("Welcome to diceseq!")

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format.")
    parser.add_option("--anno_source", dest="anno_source", default="Ensembl",
        help="The annotation source of the gtf file [default: %default].")
    parser.add_option("--sam_list", "-s", dest="sam_list", default=None,
        help="the indexed alignement file in bam/sam format, use ',' for\
        replicates and '---' for time points, e.g.,\
        my_sam1_rep1.sorted.bam,my_sam1_rep2.sorted.bam---my_sam2.sorted.bam.")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format. This is necessary\
        for bias correction, otherwise uniform mode will be used.")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The file for bias parameter.")
    parser.add_option("--bias_mode", dest="bias_mode", default="unif", 
        help="The bias mode: unif, end5, end3 or both. without ref_file, it \
        will be changed into unif [default: %default].")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="diceseq_out", help="The prefix of the output file. There will\
        be two files: one in plain text format, the other in gzip format.")
    parser.add_option("--sample_num", dest="sample_num", default="0",
        help="The number of MCMC samples to save, 0 for no such file. Advice: \
        lower than 3/4 of `min_run`, e.g, 500 [default: %default].")
    parser.add_option("--time_seq", "-t", dest="time_seq", default=None,
        help="The time for the input samples, e.g., 0,1,2,3, the default\
        values will be the index of all time points, i.e., 0,1,...")

    parser.add_option("--nproc", dest="nproc", default="4",
        help="The number of subprocesses [default: %default].")
    parser.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
        default=False, help="add the pre-mRNA as a transcript.")
    parser.add_option("--no_twice", action="store_true", dest="no_twice", 
        default=False, help="No quick estimate of the variance, but use fixed.")
    parser.add_option("--print_detail", action="store_true", dest="print_detail", 
        default=False, help="print the detail of the sampling.")

    parser.add_option("--theta1", dest="theta1", default="3.0",
        help="The fixed hyperparameter theta1 for the GP model \
        [default: %default].")
    parser.add_option("--theta2", dest="theta2", default=None,
        help="The fixed hyperparameter theta2 for the GP model \
        The default will cover 1/3 of the duration. [default: %default].")
    parser.add_option("--max_run", dest="max_run", default="5000",
        help="The maximum iterations for the MCMC sampler \
        [default: %default].")
    parser.add_option("--min_run", dest="min_run", default="1000",
        help="The minimum iterations for the MCMC sampler \
        [default: %default].")
    parser.add_option("--gap_run", dest="gap_run", default="100",
        help="The increase gap of iterations for the MCMC sampler \
        [default: %default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\rloading annotation file...")
        sys.stdout.flush()    
        anno = load_annotation(options.anno_file, options.anno_source)
        sys.stdout.write("\rloading annotation file... Done.\n")
        sys.stdout.flush()
        genes = anno["genes"]
        global TOTAL_GENE
        TOTAL_GENE = len(genes)
    if options.sam_list == None:
        print("Error: need --sam_list for reads indexed and aliged reads.")
        sys.exit(1)
    else:
        sam_list = options.sam_list.split("---")
        global TOTAL_READ
        for s in sam_list:
            _cnt = 0
            ss_list = s.split(",")
            for ss in ss_list:
                bashCommand = "samtools view -c " + ss
                process = subprocess.Popen(bashCommand.split(), 
                    stdout=subprocess.PIPE)
                _cnt += float(process.communicate()[0].split()[0])
            TOTAL_READ.append(_cnt)

    no_twice = options.no_twice
    add_premRNA = options.add_premRNA
    print_detail = options.print_detail

    M = int(options.max_run)
    gap = int(options.gap_run)
    initial = int(options.min_run)
    if options.time_seq is None: 
        X = np.arange(len(sam_list))
    else: 
        X = np.array(options.time_seq.split(","), "float")
    theta1 = float(options.theta1)
    if options.theta2 is None:
        theta2 = ((max(X) - min(X) + 0.1) / 3.0)**2
    else:
        theta2 = float(options.theta2)
    
    nproc = int(options.nproc)
    ref_file  = options.ref_file
    out_file  = options.out_file
    bias_file = options.bias_file
    bias_mode = options.bias_mode
    sample_num = int(options.sample_num)

    if ref_file is None: 
        if bias_mode != "unif":
            bias_mode = "unif"
            print("No reference sequence, so we change to uniform mode.")

    if bias_file is None: 
        biasFile = None
        if bias_mode != "unif":
            bias_mode = "unif"
            print("No bias parameter file, so we change to uniform mode.")
    else: biasFile = BiasFile(bias_file)

    global FID1, FID2
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

    print("running diceseq for %d genes with %d cores..." %(TOTAL_GENE, nproc))

    if nproc <= 1:
        for g in genes:
            if add_premRNA: g.add_premRNA()
            RV = get_psi(g, sam_list, ref_file,  biasFile, bias_mode, X, M, 
                 initial, gap, theta1, theta2, no_twice, sample_num, print_detail)
            show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        for g in genes:
            if add_premRNA: g.add_premRNA()
            pool.apply_async(get_psi, (g, sam_list, ref_file,  biasFile, bias_mode,
                 X, M, initial, gap, theta1, theta2, no_twice, sample_num, 
                 print_detail), callback=show_progress)
        pool.close()
        pool.join()
    FID1.close()
    if FID2 is not None: FID2.close()
    print("")


if __name__ == "__main__":
    main()
    