# This file is for generating reads for simulation

import sys
import gzip
import pysam
import numpy as np
from optparse import OptionParser
from diceseq.utils.tran_utils import TranUnits
from diceseq.utils.gtf_utils import load_annotation
from diceseq.utils.bias_utils import BiasFile, FastaFile

def GP_K(X, theta, x=None):
    N = len(X)
    if x is None:
        K = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                K[i,j] = theta[0] * np.exp(-0.5/theta[1]*(X[i]-X[j])**2)
    else:
        K = np.zeros(N)
        for i in range(N):
            K[i] = theta[0] * np.exp(-0.5/theta[1]*(X[i]-x)**2)
    return K

def YtoPsi(y):
    return np.exp(y) / (np.exp(y) + 1)

def YtoPsi(y, axis=0):
    """softmax function: transfer y to psi
    """
    size = np.array(y.shape)
    size[axis] = -1
    return np.exp(y) / (np.exp(y).sum(axis=axis).reshape(size))


def generate_reads(tran, N, rlen=100, flen=200, sigma=20, bias_mode="unif"):
    """generate the simulation reads. Only implimented the single-end
    reads for 5'end or 3'end. bias_mode could be 'unif', 'end5', 'end3'
    and 'both'. """
    reads1, reads2 = [], []
    if tran.strand == "1" or tran.strand == "+":
        idx5 = np.arange(0, tran.ulen-rlen+1)
        idx3 = np.arange(rlen-1, tran.ulen)
    else:
        idx5 = np.arange(rlen-1, tran.ulen)
        idx3 = np.arange(0, tran.ulen-rlen+1)

    cnt = -1
    probs = np.ones((len(idx5)+1)*len(idx5)/2)
    idx5s = np.ones((len(idx5)+1)*len(idx5)/2, "int")
    idx3s = np.ones((len(idx5)+1)*len(idx5)/2, "int")
    for i in range(len(idx5)):
        for j in range(i,len(idx3)):
            cnt += 1
            if tran.strand == "1" or tran.strand == "+":
                x = idx3[j] - idx5[i] + 1
                idx5s[cnt] = idx5[i]
                idx3s[cnt] = idx3[j]
            else:
                x = idx5[j] - idx3[i] + 1
                idx5s[cnt] = idx5[j]
                idx3s[cnt] = idx3[i]
            if x < 0: print("Errors! flen is negative!")
            probs[cnt] = normal_pdf(x, flen, sigma)
    if bias_mode == "end5" or bias_mode == "both":
        probs *= tran.bias5[idx5s]
    if bias_mode == "end3" or bias_mode == "both":
        probs *= tran.bias3[idx3s]
    probs = probs / np.sum(probs)
    count = np.random.multinomial(N, probs)

    for i in range(count.shape[0]):
        if count[i] == 0: continue
        if tran.strand == "1" or tran.strand == "+":
            fwd = tran.seq[idx5s[i]+20 : idx3s[i]+20+1]
        else:
            fwd = tran.seq[idx3s[i]+20 : idx5s[i]+20+1]
        rev = []
        rev[:] = fwd[::-1]
        for c in range(len(rev)):
            if   rev[c] == "A": rev[c] = "T"
            elif rev[c] == "T": rev[c] = "A"
            elif rev[c] == "G": rev[c] = "C"
            elif rev[c] == "C": rev[c] = "G"
        rev = "".join(rev)
        for j in range(count[i]):
            if tran.strand == "1" or tran.strand == "+":
                reads1.append(fwd[:rlen])
                reads2.append(rev[:rlen])
            else:
                reads1.append(rev[:rlen])
                reads2.append(fwd[:rlen])
    RV = {}
    RV["reads1"] = reads1
    RV["reads2"] = reads2
    return RV


def main():
    print "Welcome to dice-simulate. Reads are under generating."
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--anno_source", dest="anno_source", default="Ensembl",
        help="The annotation source of the gtf file: miso, ensemble")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The parameter file for bias in hdf5 format")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="simulate", help="The output file with simulated reads")
    parser.add_option("--add_premRNA", dest="add_premRNA", default="False",
        help="Whether adding pre-mRNA or not.")
    parser.add_option("--RPK", dest="RPK", default="1000", 
        help="The sequence depth, default=1000")
    parser.add_option("--ratio", dest="ratio", default="0.5", 
        help="The ratio of the first isoform, default=0.5.")
    parser.add_option("--times", dest="times", default=None, 
        help="The time points for generating reads with random ratios." + 
        "When it is not None, --ratio will be no useful, default=None.")
    parser.add_option("--noise", dest="noise", default="0",
        help="The noise in the reads number for each isoform, default=0.")
    parser.add_option("--rlen", dest="rlen", default="100", 
        help="The length of reads, default=100")
    parser.add_option("--fl_mean", dest="fl_mean", default="200", 
        help="The mean length of fragment, default=200")
    parser.add_option("--fl_sigma", dest="fl_sigma", default="20", 
        help="The stand variance of fragment length, default=20")
    parser.add_option("--bias_mode", dest="bias_mode", default="both", 
        help="The bias mode; unif, end5, end3, both")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file, options.anno_source)

    if options.ref_file == None:
        print "Error: need --ref_file for reference sequence."
        sys.exit(1)
    else:
        fastaFile = FastaFile(options.ref_file)

    bias_mode = options.bias_mode
    if options.bias_file is None:
        biasFile = None
        if bias_mode != "unif":
            bias_mode = "unif"
            print "No bias parameter file, so we change to uniform mode."
    else: biasFile = BiasFile(options.bias_file)

    RPK = np.array(options.RPK.split(","), "int")
    if options.times is not None:
        T = np.arange(int(options.times))
    else: T = np.arange(0)
    Ymean = np.zeros(8)
    theta = np.array([3, 5.0])
    cov = GP_K(T, theta)
    #     ratio = YtoPsi(np.random.multivariate_normal(Ymean, cov))
    # else:
    #     ratio = np.array(options.ratio.split(","), "float")
    #     T = np.arange(len(ratio))
    # Np = np.zeros((len(RPK), len(ratio)), "int")
    # Nm = np.zeros((len(RPK), len(ratio)), "int")
    # cnt = np.zeros((len(RPK), len(ratio)), "int")

    rlen = int(options.rlen)
    noise = float(options.noise)
    fl_mean = float(options.fl_mean)
    fl_sigma = float(options.fl_sigma)
    out_file = options.out_file
    add_premRNA = options.add_premRNA == "True"

    fid0 = open(out_file + "_ratios.txt", "w")
    headline = "gene_id\ttrans_id"
    for t in range(len(T)):
        headline += "\tT%d" %T[t]
    fid0.writelines(headline + "\n")
    fid0.close()

    qua = []
    qua[:] = ">" * rlen
    qua = "".join(qua)
    for i in range(len(anno["gene_id"])):
        if i == 0: ftype = "w"
        else: ftype = "a+"
        _gene = anno["genes"][i]
        if add_premRNA == True: _gene.add_premRNA()

        unitp = TranUnits(anno["genes"][i].trans[1])
        unitm = TranUnits(anno["genes"][i].trans[0])
        unitp.set_sequence(fastaFile)
        unitm.set_sequence(fastaFile)
        unitp.set_bias(biasFile)
        unitm.set_bias(biasFile)

        if options.times is not None:
            ratio = YtoPsi(np.random.multivariate_normal(Ymean, cov))
        fid0 = open(out_file + "_ratios.txt", "a+")
        aline = _gene.geneID
        for t in range(len(T)):
            aline += "\t%.5f" %ratio[t]
        fid0.writelines(aline + "\n")
        fid0.close()

        for tr in range(len(_gene.trans)):

        for m in range(len(RPK)):
            for n in range(len(ratio)):
                _Np = max(0, (unitp.ulen) / 1000.0 * RPK[m] * (1-ratio[n]))
                _Nm = max(0, (unitm.ulen) / 1000.0 * RPK[m] * ratio[n])
                _aa = np.random.rand()
                _bb = np.random.rand()
                if _aa <= _Np - int(_Np): _Np = int(_Np) + 1.0
                else: _Np = int(_Np) + 0.0
                if _bb <= _Nm - int(_Nm): _Nm = int(_Nm) + 1.0
                else: _Nm = int(_Nm) + 0.0
                Np[m,n] = max(0, np.round(np.random.normal(_Np, max(0.01,noise*_Np))))
                Nm[m,n] = max(0, np.round(np.random.normal(_Nm, max(0.01,noise*_Nm))))
        readp = unitp.generate_reads(Np.max(),rlen,fl_mean,fl_sigma,bias_mode)
        readm = unitm.generate_reads(Nm.max(),rlen,fl_mean,fl_sigma,bias_mode)

        for m in range(len(RPK)):
            for n in range(len(ratio)):
                idxp = np.random.permutation(Np.max())[:Np[m,n]]
                idxm = np.random.permutation(Nm.max())[:Nm[m,n]]
                
                fid1 = gzip.open(out_file + "_RPK%d_T%d_1.fq.gz" %(RPK[m], n), ftype)
                fid2 = gzip.open(out_file + "_RPK%d_T%d_2.fq.gz" %(RPK[m], n), ftype)
                for r in idxp:
                    cnt[m,n] += 1
                    fid1.writelines("@simulated_reads_%d_1\n" %cnt[m,n])
                    fid2.writelines("@simulated_reads_%d_2\n" %cnt[m,n])
                    fid1.writelines(readp["reads1"][r] + "\n+\n" + qua + "\n")
                    fid2.writelines(readp["reads2"][r] + "\n+\n" + qua + "\n")
                for r in idxm:
                    cnt[m,n] += 1
                    fid1.writelines("@simulated_reads_%d_1\n" %cnt[m,n])
                    fid2.writelines("@simulated_reads_%d_2\n" %cnt[m,n])
                    fid1.writelines(readm["reads1"][r] + "\n+\n" + qua + "\n")
                    fid2.writelines(readm["reads2"][r] + "\n+\n" + qua + "\n")
                fid1.close()
                fid2.close()
        if i == 0 or (i+1) % 10 == 0 or (i+1) == len(anno["gene_id"]):
            print "%d out %d genes have been finished." % (i+1, len(anno["gene_id"]))

if __name__ == "__main__":
    main()

