# This module is to parse the gtf file, whoes format can be found at:
# www.ensembl.org/info/website/upload/gff.html

# Note: here we need the order of feature lines like this: gene --> transcript
# --> exon 1 --> exon 2 ..., which is not required for the gtf format. But this
# assumption is usually right for the gtf file downloaded from Ensembl database.

import sys
import numpy as np

class Transcript:
    def __init__(self,chrom,strand,start,stop,tran_id,tran_name="*",biotype="*"):
        """a general purpose transcript object with the basic information.
        """
        self.chrom  = chrom
        self.strand = strand
        self.start  = int(start)
        self.stop   = int(stop)
        self.tranID = tran_id
        self.exons  = np.zeros((0,2), "int")
        self.seglen = None
        self.tranL  = 0
        self.exonNum = 0
        self.biotype = biotype
        self.tranName = tran_name
        

    def add_exon(self,chrom,strand,start,stop):
        if strand != self.strand or chrom != self.chrom:
            print("The exon has different chrom or strand to the transcript.")
            return
        _exon = np.array([start, stop], "int").reshape(1,2)
        self.exons = np.append(self.exons, _exon, axis=0)
        self.exons = np.sort(self.exons, axis=0)
        self.tranL += abs(int(stop) - int(start) + 1)
        self.exonNum += 1


        self.seglen = np.zeros(self.exons.shape[0] * 2 - 1, "int")
        self.seglen[0] = self.exons[0,1]-self.exons[0,0] + 1
        for i in range(1, self.exons.shape[0]):
            self.seglen[i*2-1] = self.exons[i,0]-self.exons[i-1,1] - 1
            self.seglen[i*2] = self.exons[i,1]-self.exons[i,0] + 1

        if ["-","-1","0",0,-1].count(self.strand) > 0:
            self.seglen = self.seglen[::-1]

class Gene:
    def __init__(self,chrom,strand,start,stop,gene_id,gene_name="*",biotype="*"):
        self.chrom  = chrom
        self.strand = strand
        self.start  = int(start)
        self.stop   = int(stop)
        self.geneID = gene_id
        self.trans  = []
        self.tranNum = 0
        self.biotype = biotype
        self.geneName = gene_name
        
    def add_transcipt(self, transcript):
        self.trans.append(transcript)
        self.tranNum += 1

    def get_gene_info(self):
        RV = [self.geneID, self.geneName, self.chrom, self.strand, self.start,
              self.stop, self.biotype]
        _trans = []
        for t in self.trans:
            _trans.append(t.tranID)
        RV.append(",".join(_trans))
        return RV

    def add_premRNA(self):
        _tran = Transcript(self.chrom, self.strand, self.start, self.stop, 
                           self.geneID+".p", self.geneName, self.biotype)
        _tran.add_exon(self.chrom, self.strand, self.start, self.stop)
        self.trans.append(_tran)
        self.tranNum += 1
        
    def get_exon_max_num(self):
        exonMax = 0
        for _tran in self.trans:
            exonMax = max(exonMax, _tran.exonNum)
        return exonMax

    def gene_ends_update(self):
        for t in self.trans:
            self.start = min(self.start, np.min(t.exons))
            self.stop  = max(self.stop,  np.max(t.exons))

def ensembl_gtf(anno_in):
    genes = []
    _gene = None
    for _line in anno_in:
        # comment lines
        if _line[0] == "#" or _line[0] == ">" : continue
            
        a_line = _line.split("\t")
        if len(a_line) < 9: continue
        elif a_line[2] == "gene":
            if _gene is not None:
                genes.append(_gene)
                _gene = None

            _gene_name, _gene_id, _biotype = "*", "*", "*"
            idx = a_line[8].find("gene_name")
            if idx > -1:
                _gene_name = a_line[8][idx:].split('"')[1].split("\n")[0]
            idx  = a_line[8].find("gene_id")
            if idx > -1:
                _gene_id = a_line[8][idx:].split('"')[1].split("\n")[0]
            idx = a_line[8].find("gene_biotype")
            if idx == -1: idx = a_line[8].find("gene_type")
            if idx > -1:
                _biotype = a_line[8][idx:].split('"')[1].split("\n")[0]

            _gene = Gene(a_line[0], a_line[6], a_line[3], a_line[4],
                         _gene_id, _gene_name, _biotype)

        elif a_line[2] == "transcript":
            _tran_name, _tran_id, _biotype = "*", "*", "*"

            idx = a_line[8].find("transcript_name")
            if idx > -1:
                _tran_name = a_line[8][idx:].split('"')[1].split("\n")[0]
            idx = a_line[8].find("transcript_id")
            if idx > -1:
                _tran_id = a_line[8][idx:].split('"')[1].split("\n")[0]
            idx = a_line[8].find("gene_biotype")
            if idx > -1:
                _biotype = a_line[8][idx:].split('"')[1]
            _tran = Transcript(a_line[0], a_line[6], a_line[3], a_line[4],
                               _tran_id, _tran_name, _biotype)

            if _gene is not None:
                _gene.add_transcipt(_tran)
            else:
                print("Gene is not ready before transcript.")

        elif a_line[2] == "exon":
            if a_line[0] != _gene.trans[-1].chrom:
                print("The exon is on the different chrom from the transcript.")
                continue
            if a_line[6] != _gene.trans[-1].strand:
                print("The exon is on the different strand from the transcript.")
                continue
            if _gene is not None and len(_gene.trans) > 0:
                _gene.trans[-1].add_exon(a_line[0],a_line[6],a_line[3],a_line[4])
                # _gene.gene_ends_update()
            else:
                print("Gene or transcript is not ready before exon.")

    if _gene is not None: genes.append(_gene)
    return genes


def sander_gtf(anno_in):
    genes = []
    _gene = None
    for _line in anno_in:
        # comment lines
        if _line[0] == "#" or _line[0] == ">" : continue
            
        a_line = _line.split("\t")
        if len(a_line) < 9: continue
        elif a_line[2] == "exon":
            _gene_id = "#"
            idx  = a_line[8].find("gene_id")
            if idx > -1:
                _gene_id = a_line[8][idx:].split('"')[1]
            if len(genes) > 0 and _gene_id == genes[-1].trans[-1].tranID :
                genes[-1].trans[-1].add_exon(a_line[0], a_line[6],
                                             a_line[3], a_line[4])
                genes[-1].gene_ends_update()
            else:
                _biotype = a_line[1]
                _gene_name, _gene_id = "*", "*"
                idx = a_line[8].find("gene_name")
                if idx > -1:
                    _gene_name = a_line[8][idx:].split('"')[1].split("\n")[0]
                idx  = a_line[8].find("gene_id")
                if idx > -1:
                    _gene_id = a_line[8][idx:].split('"')[1].split("\n")[0]
                _gene = Gene(a_line[0], a_line[6], a_line[3], a_line[4],
                             _gene_id, _gene_name, _biotype)

                _tran_name, _tran_id = "*", "*"
                idx = a_line[8].find("transcript_name")
                if idx > -1:
                    _tran_name = a_line[8][idx:].split('"')[1].split("\n")[0]
                idx = a_line[8].find("transcript_id")
                if idx > -1:
                    _tran_id = a_line[8][idx:].split('"')[1].split("\n")[0]
                idx = a_line[8].find("gene_biotype")
                _tran = Transcript(a_line[0], a_line[6], a_line[3], a_line[4],
                                   _tran_id, _tran_name, _biotype)

                _tran.add_exon(a_line[0], a_line[6], a_line[3], a_line[4])

                _gene.add_transcipt(_tran)
                genes.append(_gene)
    return genes


def sgd_gtf(anno_in):
    genes = []
    _gene = None
    cnt =+ 0
    for _line in anno_in:
        # comment lines
        if _line[0] == "#":continue
        elif _line[0] == ">" : break
            
        a_line = _line.split("\t")
        if len(a_line) < 9: continue
        elif a_line[2].find("gene") > -1:
            if _gene is not None:
                genes.append(_gene)
                _gene = None
            
            if a_line[2].find("gene") == 0:
                _biotype = "protein_coding"
            else:
                _biotype = a_line[2].split("_")[0]
            _gene_name, _gene_id = "*", "*"
            idx = a_line[8].find("Name")
            if idx > -1:
                _gene_name = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]
            idx  = a_line[8].find("ID")
            if idx > -1:
                _gene_id = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]
            _gene = Gene(a_line[0], a_line[6], a_line[3], a_line[4],
                         _gene_id, _gene_name, _biotype)

        # elif a_line[2] == "mRNA":
        #     _tran_name, _tran_id, _biotype = "*", "*", "*"
        #     idx = a_line[8].find("Name")
        #     if idx > -1:
        #         _tran_name = a_line[8][idx:].split(";")[0].split("=")[1]
        #     idx = a_line[8].find("ID")
        #     if idx > -1:
        #         _tran_id = a_line[8][idx:].split(";")[0].split("=")[1]
            # _tran = Transcript(a_line[0], a_line[6], a_line[3], a_line[4],
            #                    _tran_id, _tran_name, _biotype)

            _tran_name, _tran_id = _gene_name, _gene_id
            _tran = Transcript(a_line[0], a_line[6], a_line[3], a_line[4],
                               _tran_id, _tran_name, _biotype)

            if _gene is not None:
                _gene.add_transcipt(_tran)
            else:
                print("Gene is not ready before transcript.")

        elif a_line[2] == "CDS" or a_line[2] == "noncoding_exon" :
            if _gene is not None and len(_gene.trans) > 0:
                _gene.trans[-1].add_exon(a_line[0],a_line[6],a_line[3],a_line[4])
                # _gene.gene_ends_update()
            else:
                print("Gene or transcript is not ready before exon.")

    if _gene is not None: genes.append(_gene)
    return genes


def miso_gtf(anno_in):
    genes = []
    _gene = None
    for _line in anno_in:
        # comment lines
        if _line[0] == "#" or _line[0] == ">" : continue
            
        a_line = _line.split("\t")
        if len(a_line) < 9: continue
        elif a_line[2] == "gene":
            if _gene is not None:
                genes.append(_gene)
                _gene = None

            _gene_name, _gene_id, _biotype = "*", "*", "*"
            idx = a_line[8].find("Name")
            if idx > -1:
                _gene_name = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]
            idx  = a_line[8].find("ID=")
            if idx > -1:
                _gene_id = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]
            _gene = Gene(a_line[0], a_line[6], a_line[3], a_line[4],
                         _gene_id, _gene_name, _biotype)

        elif a_line[2] == "mRNA":
            _tran_name, _tran_id, _biotype = "*", "*", "*"
            idx = a_line[8].find("ID=")
            if idx > -1:
                _tran_name = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]       
            idx = a_line[8].find("ID=")
            if idx > -1:
                # print(a_line[8][idx:])
                _tran_id = a_line[8][idx:].split(";")[0].split("=")[1].split("\n")[0]
            _tran = Transcript(a_line[0], a_line[6], a_line[3], a_line[4],
                               _tran_id, _tran_name, _biotype)

            if _gene is not None:
                _gene.add_transcipt(_tran)
            else:
                print("Gene is not ready before transcript.")

        elif a_line[2] == "exon":
            if _gene is not None and len(_gene.trans) > 0:
                _gene.trans[-1].add_exon(a_line[0],a_line[6],a_line[3],a_line[4])
                # _gene.gene_ends_update()
            else:
                print("Gene or transcript is not ready before exon.")
    
    if _gene is not None: genes.append(_gene)
    return genes

def load_annotation(anno_gtf, source="Ensembl"):
    """load gene annotation in convenient structures;
    we support several source of gtf files: Ensembl, SGD, MISO, Sander.
    """
    # Note: we need the start and stop site of genes but only transcripts,
    # especially there are multiple transcripts on a gene.

    fid = open(anno_gtf, "r")
    anno_in = fid.readlines()
    fid.close()
    
    gene_name, gene_id, genes = [], [], []
    strand, exon_num, biotype = [], [], []
    chrom, gene_start, gene_stop = [], [], []
    tran_num = []

    if source == "SGD" or source == "sgd":
        genes = sgd_gtf(anno_in)
    elif source == "MISO" or source == "miso":
        genes = miso_gtf(anno_in)
    elif source == "Sander" or source == "sander":
        genes = sander_gtf(anno_in)
    else:
        genes = ensembl_gtf(anno_in)

    for g in genes:
        exon_num.append(g.get_exon_max_num())
        chrom.append(g.chrom)
        strand.append(g.strand)
        gene_start.append(g.start)
        gene_stop.append(g.stop)
        gene_id.append(g.geneID)
        gene_name.append(g.geneName)
        biotype.append(g.biotype)
        tran_num.append(g.tranNum)

    RV = {}
    RV['genes'] = genes
    RV['chrom'] = np.array(chrom)
    RV['strand'] = np.array(strand)
    RV['gene_id'] = np.array(gene_id)
    RV['biotype'] = np.array(biotype)
    RV['exon_num'] = np.array(exon_num)
    RV['tran_num'] = np.array(tran_num)
    RV['gene_name'] = np.array(gene_name)
    RV['gene_stop'] = np.array(gene_stop)
    RV['gene_start'] = np.array(gene_start)
    return RV
    
