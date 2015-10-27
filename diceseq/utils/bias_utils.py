# This module is to estimate the sequence or position biases. It provides
# ways to access and save the bias parameters.

import h5py
import numpy as np

class FastaFile:
    """docstring for FastaFile"""
    def __init__(self, fasta_file):
        fid = open(fasta_file,"r")
        all_lines = fid.readlines()
        fid.close()
        seq, self.ref, self.seq = "", [], []
        for line in all_lines:
            line = line.split("\n")[0]
            if len(self.ref) == 0 and line[0] != ">": continue
            if line.startswith(">"):
                self.ref.append(line[1:].split(" ")[0])
                if seq == "": continue
                self.seq.append(seq)
                seq = ""
            else:
                seq = seq + line
        self.seq.append(seq)

    def get_seq(self, qref, start, stop):
        """get the sequence in a given region, the start is from 1"""
        try:
            idx = self.ref.index(qref.split("chr")[-1])
        except ValueError:
            try:
                idx = self.ref.index("chr" + qref.split("chr")[-1])
            except ValueError:
                print "No reference id as the query: %s" %qref
                return None
        try:
            RV = self.seq[idx][start-1 : stop]
        except ValueError:
            print "Wrong start or stop position: %d, %d" %(start, stop)
            return None
        return RV

class BiasFile:
    """docstring for BiasFile"""
    def __init__(self, bias_file=None):
        """get the bias parameters from the hdf5 file"""
        if bias_file is None:
            self.initial_bias()
            return
        f = h5py.File(bias_file, "r")
        self.chain_len  = list(f["chain_len"])
        self.percentile = np.array(f["percentile"])
        # self.pos5_bias  = np.array(f["pos5_bias"]) / np.array(f["pos5_unif"])
        # self.pos3_bias  = np.array(f["pos3_bias"]) / np.array(f["pos3_unif"])
        self.pos5_bias  = np.array(f["pos5_bias"])
        self.pos3_bias  = np.array(f["pos3_bias"])
        self.pos5_unif  = np.array(f["pos5_unif"])
        self.pos3_unif  = np.array(f["pos3_unif"])
        self.pos5_prob  = self.pos5_bias / self.pos5_unif
        self.pos3_prob  = self.pos3_bias / self.pos3_unif
        
        idx5 = np.where(self.pos5_prob < 0)
        idx3 = np.where(self.pos3_prob < 0)
        self.pos5_prob[idx5] = 0
        self.pos3_prob[idx3] = 0

        self.base_chain = {}
        self.seq5_bias, self.seq3_bias = {}, {}
        self.seq5_unif, self.seq3_unif = {}, {}
        self.seq5_prob, self.seq3_prob = {}, {}
        for i in range(21):
            self.seq5_bias[str(i)] = np.array(f["seq5_bias/" + str(i)])
            self.seq3_bias[str(i)] = np.array(f["seq3_bias/" + str(i)])
            self.seq5_unif[str(i)] = np.array(f["seq5_unif/" + str(i)])
            self.seq3_unif[str(i)] = np.array(f["seq3_unif/" + str(i)])
            self.seq5_prob[str(i)] = (self.seq5_bias[str(i)] / 
                                      self.seq5_unif[str(i)])
            self.seq3_prob[str(i)] = (self.seq3_bias[str(i)] / 
                                      self.seq3_unif[str(i)])
            self.base_chain[str(i)]= list(f["base_chain/" + str(i)])
        f.close()

    def get_both_bias(seq, loc, ulen, end_num=5):
        """get the bias from the bias parameters"""
        prob = (self.get_seq_bias(seq, end_num) * 
                self.get_pos_bias(loc, ulen, end_num))
        return prob

    def get_seq_bias(self, seq, end_num):
        """get the sequence bias score"""
        if end_num == 5:
            parameters = self.seq5_prob
        elif end_num == 3:
            parameters = self.seq3_prob
        else:
            print "wrong end_num: %s" %str(end_num)
            return None

        prob = 1.0
        for j in range(len(seq)):
            _len = self.chain_len[j]
            _bas = seq[j-_len+1 : j+1]
            _idx = self.base_chain[str(j)].index(_bas)
            prob = prob * parameters[str(j)][_idx]
        return prob

    def get_pos_bias(self, loc, ulen, end_num):
        """get the position bias score, the loc is base pair distance
        from the 5'end of the units"""
        if end_num == 5:
            parameters = self.pos5_prob
        elif end_num == 3:
            parameters = self.pos3_prob
        else:
            print "wrong end_num: %s" %str(end_num)
            return None

        bin1 = (ulen >= self.percentile[:,0]) * (ulen <= self.percentile[:,1])
        bin2 = 20.0 * loc / ulen
        prob = parameters[bin1, bin2]
        return prob

    def set_percentile(self, ulen, K=5):
        """set the percentiles by input the lengths of unitsets, i.e., ulen, 
        and number of percentiles, K."""
        perc_gap = np.linspace(0, 100, K+1)
        _percent = np.percentile(ulen, list(perc_gap))
        self.percentile = np.zeros((K, 2))
        for i in range(K):
            self.percentile[i, 0] = int(_percent[i])+1
            self.percentile[i, 1] = int(_percent[i+1])
            if i == 0:
                self.percentile[i,0] = 0
            elif i==4:
                self.percentile[i,1] = float("inf")

    def set_base_chain(self):
        """set the sub-base chain for the variable-length Markov model (VLMM),
        which was proposed by Reberts et al, Genome Biology, 2011: 
        Figure2 in supp 3. http://genomebiology.com/2011/12/3/r22/"""
        b1 = ["A","T","G","C"]
        b2, b3 = [], []
        for i in b1:
            for j in b1:
                b2.append(j+i)
                for k in b1:
                    b3.append(k+j+i)
        base_comb = [b1, b2, b3]

        self.chain_len = [1]*4 + [2]*3 + [3]*10 + [2]*2 + [1]*2
        self.base_chain = {}
        for i in range(21):
            self.base_chain[str(i)] = base_comb[self.chain_len[i]-1]

    def initial_bias(self):
        self.set_base_chain()
        self.pos5_bias = np.zeros((5, 20))
        self.pos3_bias = np.zeros((5, 20))
        self.pos5_unif = np.zeros((5, 20))
        self.pos3_unif = np.zeros((5, 20))

        self.seq5_bias, self.seq3_bias = {}, {}
        self.seq5_unif, self.seq3_unif = {}, {}
        for i in range(len(self.chain_len)):
            self.seq5_bias[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq3_bias[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq5_unif[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq3_unif[str(i)] = np.zeros(4**self.chain_len[i])

    def set_both_bias(self, seq, loc, ulen, weight, end_num=5, mode="bias"):
        """get the bias from the bias parameters"""
        self.set_seq_bias(seq, weight, end_num, mode)
        self.set_pos_bias(loc, ulen, weight, end_num, mode)

    def set_seq_bias(self, seq, weight, end_num=5, mode="bias"):
        """get the sequence bias score"""
        for j in range(len(seq)):
            _len = self.chain_len[j]
            _bas = seq[j-_len+1 : j+1]
            _idx = self.base_chain[str(j)].index(_bas)
            if end_num == 5:
                if mode == "bias":
                    self.seq5_bias[str(j)][_idx] += weight
                elif mode == "unif":
                    self.seq5_unif[str(j)][_idx] += weight
            else:
                if mode == "bias":
                    self.seq3_bias[str(j)][_idx] += weight
                elif mode == "unif":
                    self.seq3_unif[str(j)][_idx] += weight

    def set_pos_bias(self, loc, ulen, weight, end_num=5, mode="bias"):
        """get the position bias score, the loc is base pair distance
        from the 5'end of the units"""
        bin1 = (ulen >= self.percentile[:,0]) * (ulen <= self.percentile[:,1])
        bin2 = int(20.0 * loc / (ulen + 0.0001))
        if end_num == 5:
            if mode == "bias":
                self.pos5_bias[bin1, bin2] += weight
            elif mode == "unif":
                self.pos5_unif[bin1, bin2] += weight
        else:
            if mode == "bias":
                self.pos3_bias[bin1, bin2] += weight
            elif mode == "unif":
                self.pos3_unif[bin1, bin2] += weight

    def save_file(self, out_file):
        f = h5py.File(out_file, "w")
        f["percentile"] = self.percentile
        f["pos5_bias"] = self.pos5_bias
        f["pos3_bias"] = self.pos3_bias
        f["pos5_unif"] = self.pos5_unif
        f["pos3_unif"] = self.pos3_unif
        f["chain_len"] = self.chain_len
        for i in range(21):
            f["seq5_bias/"+str(i)]  = self.seq5_bias[str(i)]
            f["seq3_bias/"+str(i)]  = self.seq3_bias[str(i)]
            f["seq5_unif/"+str(i)]  = self.seq5_unif[str(i)]
            f["seq3_unif/"+str(i)]  = self.seq3_unif[str(i)]
            f["base_chain/"+str(i)] = self.base_chain[str(i)]
        f.close()


def load_refseq(refseq_file):
    """
    Load the fasta file
    Parameters
    ----------
    refseq_file: a fasta file name
    Returns
    -------
    ref_chr: list, string
        the name of chromesomes
    ref_seq: list, string
        the sequence of chromesomes
    """
    #1. open file and read lines
    fid = open(refseq_file,"r")
    all_lines = fid.readlines()
    fid.close()
    #2. process all lines
    ref_chr, ref_seq = [], []
    seq = ""
    for line in all_lines:
        line = line.split()[0]
        if line.startswith( ">" ):
            ref_chr.append(line[1:])
            if seq == "":
                continue
            ref_seq.append(seq)
            seq = ""
        else:
            seq = seq + line
    ref_seq.append(seq)
    #3.reture chrom and sequence
    return ref_chr, ref_seq

    
def get_percentile(tran_len, K=5):
    """
    Calculate the percentiles with K parts
    Parameters
    ----------
    tran_len : numpy array, (N,)
        the length of N transcripts, >0
    K : int
        the number of parts for percentiles
    Returns
    -------
    percentile: numpy array, (K, 2)
        percentile of the transcript length
    """
    perc_gap = np.linspace(0, 100, K+1)
    _percent = np.percentile(tran_len, list(perc_gap))
    percentile = np.zeros((K,2))
    for i in range(K):
        percentile[i,0] = int(_percent[i])+1
        percentile[i,1] = int(_percent[i+1])
        if i == 0:
            percentile[i,0] = 0
        elif i==4:
            percentile[i,1] = float("inf")
    return percentile

def position_bias(end_pos, percentile, tran_len):
    """
    Calculate the weights of proportion bias for 5'/3' fragment ends
    Parameters
    ----------
    end_pos : numpy array, (N,)
        5'/3'-end positions; N: total number of reads
    percentile : numpy array, (K, 2)
        the given percentiles, K is set as 5
    tran_len : numpy array, (N,)
        the length of N transcripts, >0
    Returns
    -------
    parameters: numpy array, (5, 20)
        percentages of end postions in each of the 20 bins
    """
    parameters = np.zeros((5,20), "float")
    try:
        if len(tran_len) == 1:
            tran_len = np.ones(len(end_pos), "int") * tran_len[0]
    except:
        tran_len = np.ones(len(end_pos), "int") * tran_len
    for i in range(len(end_pos)):
        if tran_len[i] <= 0:
            break
        idx_bin1 = (tran_len[i] >= percentile[:,0]) * (tran_len[i] <= percentile[:,1])
        idx_bin2 = 20.0 * end_pos[i] / tran_len[i]
        # print end_pos[i], tran_len[i]
        if idx_bin2 < 20:
            idx_bin2 = int(idx_bin2)
        elif idx_bin2 == 20:
            idx_bin2 = 19
        elif idx_bin2 > 20:
            continue
        parameters[idx_bin1, idx_bin2] += 1
    return parameters

def get_base_chain():
    """
    discribe the VLMM
    Set the sub-base chain for the variable-length Markov model
    Returns
    -------
    chain_len: list of int, len = 21
        the length of each sub-base chain
    base_chain: list of string, len = 21
        the sub-base chain for each 21 positions with different length
    """
    b1 = ["A","T","G","C"]
    b2, b3 = [], []
    for i in b1:
        for j in b1:
            b2.append(j+i)
            for k in b1:
                b3.append(k+j+i)
    base_comb = [b1, b2, b3]

    chain_len = [1]*4 + [2]*3 + [3]*10 + [2]*2 + [1]*2
    base_chain = {}
    for i in range(21):
        base_chain[str(i)] = base_comb[chain_len[i]-1]
    return chain_len, base_chain

def sequence_bias(ref_seq, end_pos, chain_len, base_chain, strand, end):
    """
    Calculate the weights of sequence bias for 5'/3' fragment ends
    Parameters
    ----------
    ref_seq: string
        the sequence of a chromesome
    end_pos : numpy array, (N,)
        N 5'/3'-end positions on the chromesome
    chain_len: list of int, len = 21
        the length of each sub-base chain
    base_chain: list of string, len = 21
        the sub-base chain for each 21 positions with different length
    strand: "+/-" or "1/-1", string
        the strand of the transcripts
    end: 5 or 3, int
        the end of the fragment
    Returns
    -------
    parameters: structure, 21 items with variable-length list as element
        weights of each sub-base chain on specific position
    """
    parameters = {}
    for i in range(21):
        parameters[str(i)] = np.zeros(4**chain_len[i])

    for i in range(len(end_pos)):
        if end == 5 and (strand == "+" or strand == "1"):
            if end_pos[i]<8 or end_pos[i]>len(ref_seq)-13-1:
                continue
            _seq = ref_seq[end_pos[i]-8 : end_pos[i]+13]
        else:
            if end_pos[i]<12 or end_pos[i]>len(ref_seq)-9-1:
                continue
            _seq = ref_seq[end_pos[i]-12 : end_pos[i]+9][::-1]
        for j in range(len(_seq)):
            _bases = _seq[j-chain_len[j]+1 : j+1]
            _base_idx = base_chain[str(j)].index(_bases)
            parameters[str(j)][_base_idx] += 1
    return parameters