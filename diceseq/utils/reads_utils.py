# This module is to process the reads, and the location index of reads 
# in each state, i.e.,  exons, juntion between exons.

import numpy as np

class ReadSet(object):
    """docstring for ReadSet"""
    def __init__(self, reads1, reads2=[]):
        self.length = len(reads1)
        if len(reads1) != len(reads2): self.size = 1
        if len(reads1) == len(reads2): self.size = 2

        self.uloc = np.ones((self.length,self.size)) #up-stream location(small)
        self.dloc = np.ones((self.length,self.size)) #down-stream location(big)
        self.qlen = np.ones((self.length,self.size)) #aligned reads length
        self.rlen = np.ones((self.length,self.size)) #orignal reads length
        self.prob = np.ones((self.length,self.size)) #reads report probability
        self.mapp = np.ones((self.length,self.size)) #reads mapping probability

        for s in range(self.size):
            if s == 0: reads = reads1
            if s == 1: reads = reads2
            for i in range(self.length):
                r = reads[i]
                self.uloc[i,s] = r.pos
                self.dloc[i,s] = r.aend - 1 # qlen=dloc-uloc-1
                self.qlen[i,s] = r.qlen
                self.rlen[i,s] = r.rlen
                self.mapp[i,s] = 1.0 - 10 ** (0 - r.mapq / 10.0)
                qual_bases = []
                qual_bases[:] = r.qqual
                #Note: still need improvement on the mismatched bases
                for _qbase in qual_bases:
                    self.prob[i,s] = self.prob[i,s]*(1-10**(-ord(_qbase)/10.0))
                # if r.rlen < (r.aend - r.pos): print(r.rlen, r.aend, r.pos)

    def get_loc_idx(self, exons, strand, overhang=1):
        """get the index on exon or intron"""
        self.strand = strand
        self.exons  = np.sort(np.array(exons).reshape(-1,2), axis=0)
        self.seglen = np.array([self.exons[0,1]-self.exons[0,0] + 1,
                                self.exons[1,0]-self.exons[0,1] - 1,
                                self.exons[1,1]-self.exons[1,0] + 1])
        
        loc_idx = np.zeros((self.length, 8), bool)
        len_use = np.zeros((self.length, 8))

        flen  = self.dloc.max(axis=1) - self.uloc.min(axis=1) + 1
        rlen1 = self.dloc.min(axis=1) - self.uloc.min(axis=1) + 1
        rlen2 = self.dloc.max(axis=1) - self.uloc.max(axis=1) + 1

        # exon1
        loc_idx[:,0] = self.dloc.max(axis=1) <= self.exons[0,1]
        len_use[:,0] = self.seglen[0] - flen + 1

        # exon1-intron1 boundary, SSite5
        loc_idx[:,1] = ((self.uloc.min(axis=1) <= self.exons[0,1]-overhang+1) *
                        (self.dloc.max(axis=1) <= self.exons[1,0]-1) *
                        (self.dloc.max(axis=1) >= self.exons[0,1]+overhang))
        len_use[:,1] = np.minimum(flen, min(self.seglen[0], self.seglen[1]))

        # intron1
        loc_idx[:,2] = ((self.uloc.min(axis=1) >= self.exons[0,1]+1) *
                        (self.dloc.max(axis=1) <= self.exons[1,0]-1))
        len_use[:,2] = self.seglen[1] - flen + 1

        # intron1-exon2 boundary, SSite3
        loc_idx[:,3] = ((self.uloc.min(axis=1) >= self.exons[0,1]+1) *
                        (self.uloc.min(axis=1) <= self.exons[1,0]-overhang) *
                        (self.dloc.max(axis=1) >= self.exons[1,0]+overhang-1))
        len_use[:,3] = np.minimum(flen, min(self.seglen[1], self.seglen[2]))

        # exon2
        loc_idx[:,4] = self.uloc.min(axis=1) >= self.exons[1,0]
        len_use[:,4] = self.seglen[2] - flen + 1

        # exon1-exon2 junction, junction
        loc_idx[:,5] = ((self.uloc <= self.exons[0,1]-overhang+1) *
                        (self.dloc >= self.exons[1,0]+overhang-1) *
                        (self.dloc-self.uloc-self.qlen > self.seglen[1]-5)
                        ).sum(axis=1)>0
        len_use[:,5] = np.minimum(np.maximum(rlen1, rlen2), 
                                  min(self.seglen[0], self.seglen[2]))

        # exon1-intron-exon2
        loc_idx[:,6] = ((self.uloc.min(axis=1) <= self.exons[0,1]) *
                        (self.dloc.max(axis=1) >= self.exons[1,0]) *
                        ((self.uloc.max(axis=1) <= self.exons[1,0]-overhang) +
                         (self.dloc.min(axis=1) >= self.exons[0,1]+overhang))*
                        (True - loc_idx[:,5]))
        len_use[:,6] = min(self.seglen[0], self.seglen[2])
        
        # exon1-exon2 unsure, exon_both
        loc_idx[:,7] = ((self.uloc.max(axis=1) >= self.exons[1,0]) *
                        (self.dloc.min(axis=1) <= self.exons[0,1]))
        len_use[:,7] = np.minimum(self.seglen[0]-rlen1,self.seglen[2]-rlen2)+1

        idx = np.where(len_use <= 0)
        len_use[idx] = 1.0

        # revers strand
        if ["-","-1","0",0,-1].count(self.strand) > 0: 
            self.seglen = self.seglen[::-1]
            for i in range(loc_idx.shape[0]):
                loc_idx[i,:5] = loc_idx[i,:5][::-1]
                len_use[i,:5] = len_use[i,:5][::-1]

        self.loc_idx  = loc_idx
        self.len_use  = len_use
        self.RPK_use  = loc_idx / len_use
        self.loc_name = ["exon1", "SSite5", "intron", "SSite3", "exon2",
                         "jucntion", "ex1_intn_ex2", "ex1_ex2_unsure"]

        # Rmat and flen are both for pre-mRNA (p) and mature mRNA (m)
        self.Rmat = np.zeros((self.length, 2), dtype="bool")
        self.Rmat[:,0] = (self.loc_idx[:,:5].sum(axis=1) + 
                          self.loc_idx[:,6:].sum(axis=1)) #p
        self.Rmat[:,1] = (self.loc_idx[:,0] + self.loc_idx[:,4] +
                          self.loc_idx[:,5] + self.loc_idx[:,7]) #m

        idx = self.loc_idx[:,5] + self.loc_idx[:,7] #junctions
        self.flen = np.zeros((self.length, 2))
        self.flen[:, 0]   = self.dloc.max(axis=1) - self.uloc.min(axis=1)
        self.flen[:, 1]   = self.dloc.max(axis=1) - self.uloc.min(axis=1)
        self.flen[idx, 1] = self.flen[idx, 1] - self.seglen[1]

        