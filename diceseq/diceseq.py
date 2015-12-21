# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import sys
import gzip
import numpy as np
from optparse import OptionParser
from .models.model_GP import Psi_GP_MH, Psi_EM, EM_filter, EM_bootstrap
from .models.model_static import Psi_MCMC_MH
from .utils.gtf_utils import load_annotation
from .utils.bias_utils import BiasFile, FastaFile
from .utils.sam_utils import load_samfile, fetch_reads
from .utils.tran_utils import TranUnits, TranSplice

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

def main():
    import warnings
    warnings.filterwarnings('error')
    print("Welcome to diceseq GP model!")

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--anno_source", dest="anno_source", default="Ensembl",
        help="The annotation source of the gtf file")
    parser.add_option("--sam_list", "-s", dest="sam_list", default=None,
        help="the indexed alignement file in bam/sam format, use ',' for\
        replicates and '---' for time points, e.g.,\
        my_sam1_rep1.sorted.bam,my_sam1_rep2.sorted.bam---my_sam2.sorted.bam")
    parser.add_option("--time", "-t", dest="time", default=None,
        help="The time for the input samples, e.g., 0,1,2,3, the default\
        values will be the index of all time points, i.e., 0,1,...")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format. This is necessary\
        for bias correction, otherwise uniform mode will be used.")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The file for bias parameter")
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list genes in use. It is the gene id in the gtf annotation.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="diceseq_out", help="The prefix of the output file. There will\
        be two files: one in plain text format, the other in gzip format")
    # parser.add_option("--save_sample", dest="save_sample", default="True",
    #     help="whether save the sample file as well.")
    parser.add_option("--sample_num", dest="sample_num", default="500",
        help="The number of MCMC samples to save.")
    parser.add_option("--bias_mode", dest="bias_mode", default="unif", 
        help="The bias mode")
    parser.add_option("--add_premRNA", dest="add_premRNA", default="False",
        help="Whether adding pre-mRNA or not.")
    parser.add_option("--theta1_fix", dest="theta1_fix", default="3.0",
        help="The fixed hyperparameter theta1 for the GP model.")
    parser.add_option("--theta2_fix", dest="theta2_fix", default="None",
        help="The fixed hyperparameter theta2 for the GP model.")
    parser.add_option("--is_twice", dest="is_twice", default="True",
        help="Whether estimate the rates twice with a quick check first.")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file, options.anno_source)
    if options.sam_list == None:
        print("Error: need --sam_list for reads indexed and aliged reads.")
        sys.exit(1)
    else:
        sam_list = options.sam_list.split("---")
        samFiles = []
        for s in sam_list:
            _sam = []
            ss_list = s.split(",")
            for ss in ss_list:
                _sam.append(load_samfile(ss))
            samFiles.append(_sam)

    if options.gene_file == None:
        gene_list = anno["gene_id"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")
    
    time      = options.time
    ref_file  = options.ref_file
    out_file  = options.out_file
    bias_file = options.bias_file
    bias_mode = options.bias_mode
    is_twice  = options.is_twice == "True"
    sample_num = int(options.sample_num)
    add_premRNA = options.add_premRNA == "True"
    # save_sample = options.save_sample == "True"
    theta1 = float(options.theta1_fix)
    try: theta2 = float(options.theta2_fix)
    except ValueError: theta2 = None

    if ref_file is None: 
        if bias_mode != "unif":
            bias_mode = "unif"
            print("No reference sequence, so we change to uniform mode.")
    else: fastaFile = FastaFile(ref_file)

    if bias_file is None: 
        if bias_mode != "unif":
            bias_mode = "unif"
            print("No bias parameter file, so we change to uniform mode.")
    else: biasFile = BiasFile(bias_file)

    # 1. run the model
    if time is None: X = np.arange(len(sam_list))
    else: X = np.array(time.split(","), "float")
    if theta2 is None: theta2 = ((max(X) - min(X) + 0.1) / 3.0)**2

    fid1 = open(out_file + ".dice", "w")
    headline = "gene_id\ttranscripts\ttransLen\tlogLik"
    for i in range(len(X)):
        _t = str(X[i])
        headline += "\tcount_T%s\tratio_T%s\t95CI_T%s" %(_t, _t, _t)
    fid1.writelines(headline + "\n")

    fid2 = gzip.open(out_file + ".sample.gz", "w")
    fid2.writelines("# MCMC samples for latent Y\n")
    fid2.writelines("# @gene|transcripts|theta2\n")
    fid2.writelines("# y_c1t1,y_c2t1;y_c1t2,y_c2t2\n")

    g_cnt = 0
    for g in range(gene_list.shape[0]):
        g_cnt += 1
        i = np.where(np.array(anno["gene_id"]) == gene_list[g])[0][0]

        if add_premRNA == True: 
            anno["genes"][i].add_premRNA()

        gene_info = anno["genes"][i].get_gene_info()
        _transLen = []
        for t in anno["genes"][i].trans:
            _transLen.append("%d" %t.tranL)
        transLen = ",".join(_transLen)

        R_all, len_iso_all, prob_iso_all, _count = [], [], [], []
        for j in range(len(samFiles)):
            t = TranSplice(anno["genes"][i])
            for _sam in samFiles[j]:
                t.set_reads(_sam)
            if bias_mode != "unif":
                t.set_sequence(fastaFile)
                t.set_bias(biasFile)

            t.get_ready(bias_mode)
            R_all.append(t.Rmat)
            if bias_mode == "unif": 
                len_iso_all.append(t.efflen_unif)
                prob_iso_all.append(t.proU)
            else: 
                len_iso_all.append(t.efflen_bias)
                prob_iso_all.append(t.proB)
            # _count.append(len(t.read1p) + len(t.read1u) + len(t.read2u))
            _count.append(np.sum(t.Rmat.sum(axis=1)>0))
        count = np.array(_count)
        
        print_info = ("%s: %d transcript, %d reads." %
            (anno["gene_id"][i], anno["genes"][i].tranNum, np.sum(count)))

        if R_all[0].shape[1] == 1:
            psi_mean = np.ones((1,len(X)))
            psi_95 = [np.ones((len(X), 2))]
            lik_marg = 0.0
            theta_mean = None
            sample_all = []

        elif R_all[0].shape[1] > 1:
            KP_FLAG, Psi_EM = EM_filter(R_all, len_iso_all, prob_iso_all, 2)

            for t in range(len(R_all)):
                R_all[t] = R_all[t][:,KP_FLAG]
                len_iso_all[t] = len_iso_all[t][KP_FLAG]
                prob_iso_all[t] = prob_iso_all[t][:,KP_FLAG]
            print_info += " KP_NUM: %d." %(sum(KP_FLAG))
            Ymean = np.zeros((R_all[0].shape[1], len(X)))
            for tt in range(len(X)):
                Ymean[:,tt] = Psi2Y(Psi_EM[tt, KP_FLAG] + 
                                    np.random.rand(np.sum(KP_FLAG)) * 0.0001)

            M = 20000
            gap = 100
            var = 0.1 * np.ones(np.sum(KP_FLAG)-1)
            initial = 1000 + gap * len(X)

            if is_twice == True:
                _var = np.zeros(R_all[0].shape[1]-1)
                for j in range(len(R_all)):
                    _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all[j:j+1],
                        len_iso_all[j:j+1], prob_iso_all[j:j+1], X[j:j+1], 
                        Ymean[:,j:j+1], var, theta1, theta2, 100, 100, 50)
                    _var += np.var(_Y[int(_m/4):, :-1, 0], axis=0)
                var = (_var + 0.000000001) / len(R_all)
                #print(var)

            _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all, len_iso_all, 
                prob_iso_all, X, Ymean, var, theta1, theta2, M, initial, gap)
            print_info += (" %d acceptances in %d iterations." %(_cnt, _m))

            psi_mean = np.zeros((len(KP_FLAG),len(X)))
            psi_95 = []
            _idx = -1
            for c in range(len(KP_FLAG)):
                if KP_FLAG[c]:
                    _idx += 1
                    psi_95.append(get_CI(_Psi[int(_m/4):,_idx,:], 0.95))
                else: 
                    psi_95.append(np.zeros((len(X), 2)))

            psi_mean[KP_FLAG,:] = _Psi[int(_m/4):,:,:].mean(axis=0)
            lik_marg = harmonic_mean_log(_Lik[int(_m/4):])
            theta_mean = _theta[int(_m/4):,:].mean(axis=0)
            sample_all = _Y[-sample_num:,:-1,:]

            _tran_name = []
            for tt in range(len(gene_info[-1].split(","))):
                if KP_FLAG[tt]: _tran_name.append(gene_info[-1].split(",")[tt])
            _theta2 = []
            for tt in range(len(theta_mean)):
                _theta2.append("%.2e" %_theta[int(_m/4):,tt].mean())
            _head = "@%s|%s|%s\n" %(gene_info[0], ",".join(_tran_name), ",".join(_theta2))
            _samples = ""
            for ss in range(sample_num):
                for tt in range(_Y.shape[2]):
                    _line_time = []
                    for cc in range(_Y.shape[1]):
                        _line_time.append("%.2e" %_Y[int(_m/2)+ss, cc, tt])
                    _samples += ",".join(_line_time)
                    if tt < _Y.shape[2] - 1: 
                        _samples += ";"
                    else:
                        _samples += "\n"
            fid2.writelines(_head)
            fid2.writelines(_samples)

        #print(print_info)

        aline = (gene_info[0] + "\t" + gene_info[-1] + "\t" + transLen + "\t%.1e" %lik_marg)
        for i in range(len(X)):
            _ratio,_ci95 = [], []
            for c in range(len(psi_mean[:,i])):
                _ratio.append("%.3f" %psi_mean[c,i])
                _ci95.append("%.3f:%.3f" %(psi_95[c][i,1], psi_95[c][i,0]))
            _ratio,_ci95 = ",".join(_ratio), ",".join(_ci95)
            aline += "\t%d\t%s\t%s" %(count[i], _ratio, _ci95)
        fid1.writelines(aline + "\n")

        if g_cnt % 10 == 0 and g_cnt != 0:
            print("%d genes have been processed." %g_cnt)
    print("%d genes have been processed. Done!" %g_cnt)
    fid1.close()
    fid2.close()


if __name__ == "__main__":
    main()
    