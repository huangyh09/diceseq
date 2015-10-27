# This module is to process different transcript state, and set information
# for the transcript states, and finally get ready to run MCMC sampling.

import h5py
import numpy as np
from sam_utils import load_samfile, fetch_reads
from bias_utils import FastaFile, BiasFile
from gtf_utils import Gene, Transcript

def normal_pdf(x, mu, sigma):
    RV = 1 / (sigma*np.sqrt(2*np.pi)) * np.exp(-1.0/2*((x-mu)/sigma)**2)
    return RV

class TranUnits:
    """docstring for TranUnits"""
    def __init__(self, transcript):
        self.chrom  = transcript.chrom
        self.strand = transcript.strand
        self.units  = transcript.exons
        self.loci   = np.array([],"int")
        for i in range(self.units.shape[0]):
            _loci   = np.arange(self.units[i,0], self.units[i,1]+1)
            self.loci = np.append(self.loci, _loci)
        self.ulen   = len(self.loci)

    def set_sequence(self, fastaFile=None):
        """set the sequence of the transcript units, with 20 bases longer in 
        each end"""
        if fastaFile is not None:
            self.seq = fastaFile.get_seq(self.chrom, self.units[0,0] - 20, 
                                                     self.units[0,0] - 1)
            for i in range(self.units.shape[0]):
                self.seq += fastaFile.get_seq(self.chrom, self.units[i,0], 
                                                          self.units[i,1])
            self.seq += fastaFile.get_seq(self.chrom, self.units[i,1] + 1, 
                                                      self.units[i,1] + 20)
        else:
            print "This is a Null sequence file."

    def set_bias(self, biasFile, mode="both"):
        """set the bias for all loci, with the bias modes of sequence / 
        position / both; make sure setting quence before using sequence
        or both modes"""
        if biasFile is None:
            # print "This is a Null bias file."
            return
            
        self.bias_method = mode
        self.bias5 = np.ones(self.loci.shape[0])
        self.bias3 = np.ones(self.loci.shape[0])

        if ["seq", "sequence", "both"].count(mode) > 0:
            for i in range(len(self.bias5)):
                ipos = i + 20
                if self.strand == "+" or self.strand == "1":
                    _seq5 = self.seq[ipos-8  : ipos+13]
                    _seq3 = self.seq[ipos-12 : ipos+9 ][::-1]
                else:
                    _seq5 = self.seq[ipos-12 : ipos+9 ][::-1]
                    _seq3 = self.seq[ipos-8  : ipos+13]

                self.bias5[i] *= biasFile.get_seq_bias(_seq5, 5)
                self.bias3[i] *= biasFile.get_seq_bias(_seq3, 3)

        if ["pos", "position", "both"].count(mode) > 0:
            for i in range(len(self.bias5)):
                if self.strand == "+" or self.strand == "1":
                    _pos = i
                else:
                    _pos = len(self.loci) - 1 - i

                self.bias5[i] *= biasFile.get_pos_bias(_pos, self.ulen, 5)
                self.bias3[i] *= biasFile.get_pos_bias(_pos, self.ulen, 3)

    # def generate_reads(self, N, rlen=50, bias_mode="unif"):
    #     """generate the simulation reads. Only implimented the single-end
    #     reads for 5'end or 3'end. bias_mode could be 'unif', 'end5', 'end3'
    #     and 'both'. """
    #     reads = []
    #     if self.strand == "1" or self.strand == "+":
    #         idx5 = np.arange(0, self.ulen-rlen+1)
    #         idx3 = np.arange(rlen-1, self.ulen)
    #     else:
    #         idx5 = np.arange(rlen-1, self.ulen)
    #         idx3 = np.arange(0, self.ulen-rlen+1)

    #     probs = np.ones(self.ulen-rlen+1)
    #     if bias_mode == "end5" or bias_mode == "both":
    #         probs *= self.bias5[idx5]
    #     if bias_mode == "end3" or bias_mode == "both":
    #         probs *= self.bias3[idx3]
    #     probs = probs / np.sum(probs)
    #     count = np.random.multinomial(N, probs)
    #     #print count

    #     for i in range(count.shape[0]):
    #         r = self.seq[i+20 : i+20+rlen]
    #         if self.strand == "-1" or self.strand == "-":
    #             rev = []
    #             rev[:] = r[::-1]
    #             for c in range(len(rev)):
    #                 if   rev[c] == "A": rev[c] = "T"
    #                 elif rev[c] == "T": rev[c] = "A"
    #                 elif rev[c] == "G": rev[c] = "C"
    #                 elif rev[c] == "C": rev[c] = "G"
    #             r = "".join(rev)
    #         for j in range(count[i]):
    #             reads.append(r)
    #     return reads

    def generate_reads(self, N, rlen=100, flen=200, sigma=20, bias_mode="both"):
        """generate the simulation reads. Only implimented the single-end
        reads for 5'end or 3'end. bias_mode could be 'unif', 'end5', 'end3'
        and 'both'. """
        reads1, reads2 = [], []
        if self.strand == "1" or self.strand == "+":
            idx5 = np.arange(0, self.ulen-rlen+1)
            idx3 = np.arange(rlen-1, self.ulen)
        else:
            idx5 = np.arange(rlen-1, self.ulen)
            idx3 = np.arange(0, self.ulen-rlen+1)

        cnt = -1
        probs = np.ones((len(idx5)+1)*len(idx5)/2)
        idx5s = np.ones((len(idx5)+1)*len(idx5)/2, "int")
        idx3s = np.ones((len(idx5)+1)*len(idx5)/2, "int")
        for i in range(len(idx5)):
            for j in range(i,len(idx3)):
                cnt += 1
                if self.strand == "1" or self.strand == "+":
                    x = idx3[j] - idx5[i] + 1
                    idx5s[cnt] = idx5[i]
                    idx3s[cnt] = idx3[j]
                else:
                    x = idx5[j] - idx3[i] + 1
                    idx5s[cnt] = idx5[j]
                    idx3s[cnt] = idx3[i]
                if x < 0: print "Errors! flen is negative!"
                probs[cnt] = (1 / (sigma*np.sqrt(2*np.pi)) * 
                              np.exp(-0.5*((x-flen)/sigma)**2))
        if bias_mode == "end5" or bias_mode == "both":
            probs *= self.bias5[idx5s]
        if bias_mode == "end3" or bias_mode == "both":
            probs *= self.bias3[idx3s]
        probs = probs / np.sum(probs)
        count = np.random.multinomial(N, probs)
        #print count

        for i in range(count.shape[0]):
            if count[i] == 0: continue
            if self.strand == "1" or self.strand == "+":
                fwd = self.seq[idx5s[i]+20 : idx3s[i]+20+1]
            else:
                fwd = self.seq[idx3s[i]+20 : idx5s[i]+20+1]
            rev = []
            rev[:] = fwd[::-1]
            for c in range(len(rev)):
                if   rev[c] == "A": rev[c] = "T"
                elif rev[c] == "T": rev[c] = "A"
                elif rev[c] == "G": rev[c] = "C"
                elif rev[c] == "C": rev[c] = "G"
            rev = "".join(rev)
            for j in range(count[i]):
                if self.strand == "1" or self.strand == "+":
                    reads1.append(fwd[:rlen])
                    reads2.append(rev[:rlen])
                else:
                    reads1.append(rev[:rlen])
                    reads2.append(fwd[:rlen])
        RV = {}
        RV["reads1"] = reads1
        RV["reads2"] = reads2
        return RV

    def get_index(self, loc):
        RV = [-1, -1]
        if   loc < self.units[0,0]:   RV[:] = [0, -2]
        elif loc > self.units[-1,-1]: RV[:] = [self.ulen-1, -3]
        else:
            cnt = 0
            for i in range(self.units.shape[0]):
                if loc >= self.units[i,0] and loc <= self.units[i,1]:
                    RV = [cnt + loc - self.units[i,0], i]
                    break
                else:
                    cnt += self.units[i,1] - self.units[i,0] + 1
        return RV

    def get_read_info(self, r1, r2, part_in=True):
        if r1 is None and r2 is None:
            return None
        idx5, idx3 = None, None
        mapq1, idx51, idx31 = 0.0, None, None
        mapq2, idx52, idx32 = 0.0, None, None

        if r1 is not None:
            mapq1 = 1.0 - 10 ** (0 - r1.mapq / 10.0)
            if self.strand == "+" or self.strand == "1":
                idx51 = self.get_index(r1.pos)
                idx31 = self.get_index(r1.aend - 1)
            else:
                idx31 = self.get_index(r1.pos)
                idx51 = self.get_index(r1.aend - 1)
            if  idx51[1] == -1 or idx31[1] == -1: return None
            elif idx51[1] >= 0 and idx31[1] >= 0:
                if (abs(idx51[0] - idx31[0]) + 1 > r1.qlen + 3 or 
                    abs(idx51[0] - idx31[0]) + 1 < r1.qlen - 3): return None

        if r2 is not None:
            mapq2 = 1.0 - 10 ** (0 - r2.mapq / 10.0)
            if self.strand == "+" or self.strand == "1":
                idx52 = self.get_index(r2.pos)
                idx32 = self.get_index(r2.aend - 1)
            else:
                idx32 = self.get_index(r2.pos)
                idx52 = self.get_index(r2.aend - 1)
            if  idx52[1] == -1 or idx32[1] == -1: return None
            elif idx52[1] >= 0 and idx32[1] >= 0:
                if (abs(idx52[0] - idx32[0]) + 1 > r2.qlen + 3 or 
                    abs(idx52[0] - idx32[0]) + 1 < r2.qlen - 3): return None

        if r1 is None: 
            flen = abs(idx52[0] - idx32[0]) + 1
            if idx32[1] >= 0: idx3 = idx32[0]
        elif r2 is None:
            flen = abs(idx51[0] - idx31[0]) + 1
            if idx51[1] >= 0: idx5 = idx51[0]
        else:
            flen = abs(idx32[0] - idx51[0]) + 1
            if idx51[1] >= 0: idx5 = idx51[0]
            if idx32[1] >= 0: idx3 = idx32[0]
            
        RV = {}
        RV["idx5"] = idx5
        RV["idx3"] = idx3
        RV["flen"] = flen
        RV["prob"] = max(mapq1, mapq2)
        return RV

    def set_reads(self, reads1=[], reads2=[], bias_mode="both"):
        """identify whether a read (pair) or is in this units, and return the 
        identity of the reads in this units, the units specific fragment  
        length bias scores. The bias score can be based on both (mode=both) 
        ends or single end5 (mode=end5) or single end3 (mode=end5), uniform 
        (mode=unif). Make sure the loc of read1 is smaller than read2."""
        if len(reads1) == 0 and len(reads2) == 0:
            print "Please input paired-end reads or singled end reads!"
            sys.exit(1)
        elif (len(reads1) * len(reads2)) != 0 and len(reads2) != len(reads1):
            print "Please input the same number of both mates of the reads!"
            sys.exit(1)

        self.rcnt = max(len(reads1), len(reads2))
        self.Rmat = np.ones(self.rcnt, "bool")
        self.flen = np.ones(self.rcnt, "float")
        self.proB = np.ones(self.rcnt, "float")
        self.proU = np.ones(self.rcnt, "float")
        self.bias_mode = bias_mode
        self.efflen_bias = 0
        self.efflen_unif = 0

        for i in range(self.rcnt):
            if len(reads1) == 0: r1 = None
            else: r1 = reads1[i]
            if len(reads2) == 0: r2 = None
            else: r2 = reads2[i]

            rinfo = self.get_read_info(r1, r2)
            if rinfo is None: 
                self.Rmat[i] = False
                self.flen[i] = None
                self.proB[i] = None
                self.proU[i] = None
            else:
                self.Rmat[i] = True
                self.flen[i] = rinfo["flen"]
                self.proB[i] = rinfo["prob"]
                self.proU[i] = rinfo["prob"]
                if self.bias_mode == "unif": continue
                elif self.bias_mode != "end3" and rinfo["idx5"] is not None: 
                    self.proB[i] *= self.bias5[rinfo["idx5"]]
                elif self.bias_mode != "end5" and rinfo["idx3"] is not None: 
                    self.proB[i] *= self.bias3[rinfo["idx3"]]

        # fragement distribution
        flen  = self.flen[self.Rmat]
        mu    = np.mean(flen)
        sigma = np.std(flen)
        x     = np.arange(1, self.ulen+1)
        uniqL = np.unique(flen)
        Total = sum(self.Rmat) + 0.0
        self.probs = np.zeros(self.ulen)
        if   Total == 0: self.probs[0] = 1.0
        elif uniqL.shape[0] >= 10:
            self.probs[:] = (1 / (sigma*np.sqrt(2*np.pi)) * 
                             np.exp(-0.5*((x-mu)/sigma)**2))
            if sum(self.probs) != 0: self.probs = self.probs / sum(self.probs)
        else:
            for i in uniqL: 
                self.probs[i-1] = sum(flen==i) / Total
        # effective length
        self.biasLen = np.zeros(self.ulen)
        for i in range(1, self.ulen+1):
            self.efflen_unif += self.probs[i-1] * (self.ulen-i+1)
            if self.bias_mode == "unif": continue
            if self.probs[i-1] == 0:     continue

            for j in range(self.ulen - i + 1):
                if self.strand == "+" or self.strand == "1":
                    pos5, pos3 = j, j+i-1
                else:
                    pos3, pos5 = j, j+i-1
                if   self.bias_mode == "end5": _bias = self.bias5[pos5]
                elif self.bias_mode == "end3": _bias = self.bias3[pos3]
                else : _bias = self.bias5[pos5] * self.bias3[pos3]
                self.biasLen[i-1] += _bias
            self.efflen_bias += self.probs[i-1] * self.biasLen[i-1]
        # print self.flen
        # reads probability
        for i in range(self.rcnt):
            if self.Rmat[i] == False: continue
            fL = self.flen[i]
            self.proU[i] *= self.probs[fL-1] / (self.ulen - fL + 1)                
            if self.bias_mode != "unif":
                self.proB[i] *= (self.probs[fL-1] / self.biasLen[fL-1])


# class Transcript:
#     def __init__(self,tran_id,chrom,strand,start,stop,exons):
#         """a general purpose transcript object with the basic information and
#         the sequencing reads.
#         """
#         self.tranID = tran_id
#         self.chrom  = chrom
#         self.strand = strand
#         self.start  = start
#         self.stop   = stop
#         self.exons  = np.sort(np.array(exons).reshape(-1,2), axis=0)
#         self.seglen = np.array([self.exons[0,1]-self.exons[0,0] + 1,
#                                 self.exons[1,0]-self.exons[0,1] - 1,
#                                 self.exons[1,1]-self.exons[1,0] + 1])
#         if ["-","-1","0",0,-1].count(self.strand) > 0:
#             self.seglen = self.seglen[::-1]


class TranSplice:
    def __init__(self, Gene):
        self.gene = Gene
        self.chrom = Gene.chrom
        self.start = Gene.start
        self.stop = Gene.stop
        self.initial_unitSet()

        self.read1p = []
        self.read2p = []
        self.read1u = []
        self.read2u = []
    def initial_unitSet(self, add_premRNA=False):
        """set the initial units: pre-mRNA & mature mRNA"""
        if add_premRNA == True:
            self.gene.add_premRNA()
        self.unitSet = []
        for i in range(len(self.gene.trans)):
            self.unitSet.append(TranUnits(self.gene.trans[i]))
        # self.unitSet["u1"] = TranUnits(self.chrom, self.strand, 
        #                               [self.start, self.stop])
        # self.unitSet["u2"] = TranUnits(self.chrom, self.strand, self.exons)

    def add_units(self, tranunits):
        """"add a new state, which is a TranUnits object."""
        self.unitSet.append(tranunits)

    def set_sequence(self, fastafile):
        """get the sequence from the genome sequence in FastaFile object"""
        for i in len(self.unitSet):
            self.unitSet[i].set_sequence(fastafile)

    def set_bias(self, biasfile):
        """get the bias parameters from the BiasFile object"""
        for i in len(self.unitSet):
            self.unitSet[i].set_bias(biasfile)
        self.bias_in = True

    def set_reads(self, samfile, rm_duplicate=True, inner_only=True,
                  mapq_min=0, mismatch_max=5, rlen_min=1, is_mated=True):
        """fetch the reads on this transcript from sam file"""
        reads = fetch_reads(samfile, self.chrom, self.start, self.stop,  
                            rm_duplicate, inner_only, mapq_min, 
                            mismatch_max, rlen_min,   is_mated)

        self.read1p += reads["reads1"]
        self.read2p += reads["reads2"]
        self.read1u += reads["reads1u"]
        self.read2u += reads["reads2u"]

    def get_ready(self, bias_mode="both"):
        """get the location index of the transcript, need implimentation
        in future. Then, we could remove the Rmat, and flen from ReadSet
        object, and set the Rmat and flen by the info of states."""
        rcnt  = [len(self.read1p), len(self.read1u), len(self.read2u)]
        unit_cnt  = len(self.unitSet)
        self.Rmat = np.ones((sum(rcnt), unit_cnt), "bool")
        self.flen = np.ones((sum(rcnt), unit_cnt), "float")
        self.proB = np.ones((sum(rcnt), unit_cnt), "float")
        self.proU = np.ones((sum(rcnt), unit_cnt), "float")
        self.efflen_bias = np.zeros(unit_cnt, "float")
        self.efflen_unif = np.zeros(unit_cnt, "float")

        for i in range(unit_cnt):
            _units = self.unitSet[i]
            for j in range(3):
                _idx = np.arange(sum(rcnt[:j]), sum(rcnt[:j+1]))
                _reads1, _reads2 = [], []
                if   j == 0: _reads1, _reads2 = self.read1p, self.read2p
                elif j == 1: _reads1 = self.read1u
                elif j == 2: _reads2 = self.read2u
                if (len(_reads1) + len(_reads2)) == 0: continue

                _units.set_reads(_reads1, _reads2, bias_mode)
                self.Rmat[_idx, i] = _units.Rmat
                self.flen[_idx, i] = _units.flen
                self.proB[_idx, i] = _units.proB
                self.proU[_idx, i] = _units.proU

                self.efflen_unif[i] += sum(_units.Rmat) * _units.efflen_unif
                self.efflen_bias[i] += sum(_units.Rmat) * _units.efflen_bias

            self.efflen_unif[i] /= sum(self.Rmat[:, i])
            self.efflen_bias[i] /= sum(self.Rmat[:, i])

    def quick_start(self, samfile=None, fastafile=None, biasfile=None):
        """update the state_read after change the transcript information or 
        adding or removing states"""
        self.initial_unitSet()
        if samfile   is not None: self.set_reads(samfile)
        if fastafile is not None: self.set_sequence(fastafile)
        if biasfile  is not None: self.set_bias(biasfile)
        print len(self.read1p), len(self.read1u), len(self.read2u)

        self.get_ready()
        print "ready"
        # try:
        #     self.get_ready()
        # except NameError:
        #     print "Cannot get ready now; please input necessary info first!"

