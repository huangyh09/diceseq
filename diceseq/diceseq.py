# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import sys
import h5py
import pysam
import numpy as np
from optparse import OptionParser
from models.model_GP import Psi_GP_MH
from models.model_static import Psi_MCMC_MH
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads
from utils.tran_utils import TranUnits, Transcript, TranSplice

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

def main():
    print "Welcome to diceseq GP model!"

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
        help="The parameter file for bias in hdf5 format")
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list genes in use. It is the gene id in the gtf annotation.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="diceseq_out", help="The prefix of the output file. There will\
        be two files: one in plain text format, the other in hdf5 format")
    parser.add_option("--out_h5", dest="out_h5", default="True",
        help="whether save a hdf5 file as well.")
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
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file, options.anno_source)
    if options.sam_list == None:
        print "Error: need --sam_list for reads indexed and aliged reads."
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
    out_h5    = options.out_h5 == "True"
    is_twice  = options.is_twice == "True"
    sample_num = int(options.sample_num)
    add_premRNA = options.add_premRNA == "True"
    theta1 = float(options.theta1_fix)
    try: theta2 = float(options.theta2_fix)
    except ValueError: theta2 = None

    if ref_file is None: 
        if bias_mode != "unif":
            bias_mode = "unif"
            print "No reference sequence, so we change to uniform mode."
    else: fastaFile = FastaFile(ref_file)

    if bias_file is None: 
        if bias_mode != "unif":
            bias_mode = "unif"
            print "No bias parameter file, so we change to uniform mode."
    else: biasFile = BiasFile(bias_file)

    # 1. run the model
    y_mean, y_var = [], []
    psi_75, psi_95 = [], []
    gene_info, count = [], []
    psi_mean, psi_var = [], []
    theta_mean, pro_mean = [], []
    lik_marg, sample_all = [], []

    if time is None: X = np.arange(len(sam_list))
    else: X = np.array(time.split(","), "float")
    if theta2 is None: theta2 = ((max(X) - min(X) + 0.1) / 3.0)**2

    fid = open(out_file + ".dice", "w")
    headline = "gene_id\ttranscripts\tlogLik"
    for i in range(len(X)):
        _t = str(X[i])
        headline += "\tcount_T%s\tratio_T%s\t95CI_T%s" %(_t, _t, _t)
    fid.writelines(headline + "\n")

    g_cnt = 0
    for g in range(gene_list.shape[0]):
        g_cnt += 1
        i = np.where(np.array(anno["gene_id"]) == gene_list[g])[0][0]

        if add_premRNA == True: 
            anno["genes"][i].add_premRNA()

        gene_info.append(anno["genes"][i].get_gene_info())

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
            _count.append(len(t.read1p) + len(t.read1u) + len(t.read2u))
        count.append(np.array(_count))

        M = 50000       
        var = 0.5
        gap = 500
        Ymean = np.zeros((R_all[0].shape[1], len(X)))
        initial = 1000 + gap * len(X)

        if is_twice == True:
            _var = 0
            for j in range(len(R_all)):
                _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all[j:j+1],
                    len_iso_all[j:j+1], prob_iso_all[j:j+1], X[j:j+1], 
                    Ymean[:,j:j+1], var, theta1, theta2, 100, 100, 50)
                # Ymean[:,j] = _Y[_m/4:,:,0].mean(axis=0)
                _var += np.var(_Y[_m/4:, :-1, 0], axis=0).mean()
            if _var / len(R_all) is not None and _var / len(R_all) > 0:
                var = _var / len(R_all)

        _Psi, _Y, _theta, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all, len_iso_all, 
            prob_iso_all, X, Ymean, var, theta1, theta2, M, initial, gap)
        #print "%s done. %d acceptances in %d iterations." %(anno["gene_id"][i],_cnt,_m)

        temp75 = []
        temp95 = []
        for c in range(_Psi.shape[1]):
            temp75.append(get_CI(_Psi[_m/4:,c,:], 0.75))
            temp95.append(get_CI(_Psi[_m/4:,c,:], 0.95))

        psi_75.append(temp75)
        psi_95.append(temp95)
        y_var.append(_Y[_m/4:,:,:].var(axis=0))
        y_mean.append(_Y[_m/4:,:,:].mean(axis=0))
        psi_var.append(_Psi[_m/4:,:,:].var(axis=0))
        psi_mean.append(_Psi[_m/4:,:,:].mean(axis=0))
        pro_mean.append(_Pro[_m/4:].mean())
        lik_marg.append(harmonic_mean_log(_Lik[_m/4:]))
        theta_mean.append(_theta[_m/4:,:].mean(axis=0))
        sample_all.append(_Y[-sample_num:,:-1,:])

        aline = gene_info[g][0] + "\t" + gene_info[g][-1] + "\t%.1e" %lik_marg[g]
        for i in range(len(X)):
            _ratio,_ci95 = [], []
            for c in range(len(psi_mean[g][:,i])-1):
                _ratio.append("%.3f" %psi_mean[g][c,i])
                _ci95.append("%.3f:%.3f" %(psi_95[g][c][i,1], psi_95[g][c][i,0]))
            _ratio,_ci95 = ",".join(_ratio), ",".join(_ci95)
            aline += "\t%d\t%s\t%s" %(count[g][i], _ratio, _ci95)
        fid.writelines(aline + "\n")

        if g_cnt % 10 == 0 and g_cnt != 0:
            print "%d genes have been processed." %g_cnt
    print "%d genes have been processed. Done!" %g_cnt
    fid.close()

    if out_h5:
        fout = h5py.File(out_file + ".hdf5", "w")
        fout["X"] = X
        fout["count"] = count
        fout["y_var"] = y_var
        fout["y_mean"] = y_mean
        fout["psi_75"] = psi_75
        fout["psi_95"] = psi_95
        fout["psi_var"] = psi_var
        fout["psi_mean"] = psi_mean
        fout["pro_mean"] = pro_mean
        fout["lik_marg"] = lik_marg
        fout["gene_info"] = gene_info
        fout["sample_all"] = sample_all
        fout["theta_mean"] = theta_mean
        fout.close()

if __name__ == "__main__":
    main()
    