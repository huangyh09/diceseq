import os
import numpy as np

import pyximport; pyximport.install()
from ..models.model_GP import Psi_GP_MH
from .out_utils import id_mapping
from .sam_utils import load_samfile
from .gtf_utils import load_annotation
from .bias_utils import BiasFile, FastaFile
from .tran_utils import TranUnits, TranSplice


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
    RV = RV - np.log(np.sum(np.exp(np.min(data, axis=0) - data), axis=0)) #make sure avoiding log(0)
    if log_out is not True: RV = np.exp(RV)
    return RV

def Psi2Y(Psi):
    Psi = Psi / np.sum(Psi)
    Y = np.zeros(len(Psi))
    Y[:-1] = np.log(Psi[:-1] / Psi[-1])
    return Y


def sort_dice_file(fname, tran_ids):
    data = np.loadtxt(fname, dtype=str, delimiter="\t")
    idx = id_mapping(tran_ids, data[1:, 0])
    fid = open(fname, "w")
    fid.writelines("\t".join(list(data[0,:])) + "\n")
    for i in idx:
        fid.writelines("\t".join(list(data[i+1,:])) + "\n")
    fid.close()
    

def get_psi(gene, sam_list, ref_file,  bias_file, bias_mode, X, M, initial, gap,
            theta1, theta2, no_twice, sample_num, pout, FLmean, FLstd, 
            mate_mode, auto_min, TOTAL_READ):
    samFiles = []
    for sam_time in sam_list:
        _sam = []
        for sam_rep in sam_time:
            _sam.append(load_samfile(sam_rep))
        samFiles.append(_sam)
    if ref_file is not None: 
        fastaFile = FastaFile(ref_file)

    R_all, len_iso_all, prob_iso_all, _count = [], [], [], []
    for i in range(len(samFiles)):
        t = TranSplice(gene)
        for j in range(len(samFiles[i])):
            t.set_reads(samFiles[i][j])
            if bias_mode != "unif":
                if len(bias_file) == len(samFiles):
                    _bias = BiasFile(bias_file[i])
                elif len(bias_file) == 1:
                    _bias = BiasFile(bias_file[0])
                else: continue
                t.set_sequence(fastaFile)
                t.set_bias(_bias, "seq") #under development
                if FLmean is None and _bias.flen_mean != 0: 
                    FLmean = _bias.flen_mean
                if FLstd is None and _bias.flen_std != 0: 
                    FLstd = _bias.flen_std

        t.get_ready(bias_mode, FLmean, FLstd, mate_mode, auto_min)
        R_all.append(t.Rmat)
        if bias_mode == "unif": 
            len_iso_all.append(t.efflen_unif)
            prob_iso_all.append(t.proU)
        else: 
            #len_iso_all.append(t.efflen_bias)
            len_iso_all.append(t.efflen_unif) #under development
            prob_iso_all.append(t.proB)

        # for count only
        if mate_mode != "pair":
            t.get_ready("unif", FLmean, FLstd, "pair", auto_min)
        _count.append(np.sum(t.Rmat.sum(axis=1)>0))
        # print(len(t.read1p), len(t.read1u), len(t.read2u))
    count = np.array(_count)
    
    print_info = ("%s: %d transcript, %d reads." %
        (gene.geneID, gene.tranNum, np.sum(count)))

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
                _Psi, _Y, _theta2, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all[j:j+1],
                    len_iso_all[j:j+1], prob_iso_all[j:j+1], X[j:j+1], 
                    Ymean[:,j:j+1], var, theta1, theta2, 100, 100, 50)
                _var += np.var(_Y[int(_m/4):, :-1, 0], axis=0)
            var = (_var + 0.000001) / len(R_all)
            theta2_var = np.var(_theta2[int(_m/4):, 0])/10.0 + 0.001
        # print(var)

        _Psi, _Y, _theta2, _Pro, _Lik, _cnt, _m = Psi_GP_MH(R_all, len_iso_all, 
            prob_iso_all, X, Ymean, var, theta1, theta2, M, initial, gap, 
            theta2_std=theta2_var)
        print_info += (" %d acceptances in %d iterations." %(_cnt, _m))

        psi_95 = []
        for c in range(_Psi.shape[1]):
            psi_95.append(get_CI(_Psi[int(_m/4):, c, :], 0.95))
        psi_mean = _Psi[int(_m/4):, :, :].mean(axis=0)
        lik_marg = harmonic_mean_log(_Lik[int(_m/4):])
        theta2_avg = _theta2[int(_m/4):,:].mean(axis=0)
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
            _line += "\t%.2e\t%.2e\t%.2e\t%.2e" %(FPKM, psi_mean[c,t],
                psi_95[c][t,1], psi_95[c][t,0])
        diceL += _line + "\n"

    Theta2 = []
    for tt in range(_theta2.shape[1]):
        Theta2.append("%.2e" %_theta2[int(_m/4):,tt].mean())
    sampL = "@%s|%s|%s\n" %(gene_info[0], gene_info[-1], ",".join(Theta2))
    for ss in range(min(int(0.5*_m), sample_num)):
        for tt in range(_Y.shape[2]):
            _line_time = []
            for cc in range(_Y.shape[1]):
                _line_time.append("%.2e" %_Y[int(_m/2)+ss, cc, tt])
            sampL += ",".join(_line_time)
            if tt < _Y.shape[2] - 1: 
                sampL += ";"
            else:
                sampL += "\n"
    
    RV = {}
    RV["dice_line"] = diceL
    RV["sample_line"] = sampL
    return RV
