
# This module is to parse the gtf file, whoes format can be found at:
# www.ensembl.org/info/website/upload/gff.html

# Note: here we assume that the lines of a gene or transcript are put together
# which is not required for the gtf format. This assumption is usually right 
# for the gtf file downloaded from Ensembl database.

import sys
import numpy as np

def load_annotation(anno_gtf):
    """load gene annotation in convenient structures"""
    # Note: we need the start and stop site of genes but only transcripts,
    # especially there are multiple transcripts on a gene.

    gene_name, gene_id, tran_id = [], [], []
    chrom, tran_start, tran_stop = [], [], []
    strand, exons, exon_num, biotype = [], [], [], []
    
    fid = open(anno_gtf, "r")
    anno_in = fid.readlines()
    fid.close()

    _exons = []
    for i in range(len(anno_in))[::-1]: # from down to top
        # comment lines
        if len(anno_in[i]) != 0 and anno_in[i][0] == "#":
            break

        a_line = anno_in[i].split("\t")
        if a_line[2] == "exon":
            _exons.append(a_line[3:5])

        if a_line[2] == "transcript":
            _gene_name, _gene_id, _tran_id, _biotype = "*", "*", "*", "*"

            idx = a_line[8].find("gene_name")
            if idx > -1:
                _gene_name = a_line[8][idx:].split('"')[1]
            
            idx  = a_line[8].find("gene_id") 
            if idx > -1:
                _gene_id = a_line[8][idx:].split('"')[1]
                
            idx = a_line[8].find("transcript_id")
            if idx > -1:
                _tran_id = a_line[8][idx:].split('"')[1]

            idx = a_line[8].find("gene_biotype")
            if idx > -1:
                _biotype = a_line[8][idx:].split('"')[1]


            if _tran_id != "*":
                exon_num.append(len(_exons))
                exons.append(np.array(_exons[::-1], "int"))
                chrom.append(a_line[0])
                strand.append(a_line[6])
                tran_start.append(a_line[3])
                tran_stop.append(a_line[4])
                tran_id.append(_tran_id)
                gene_id.append(_gene_id)
                biotype.append(_biotype)
                gene_name.append(_gene_name)
            _exons = []

    RV = {}
    RV['chrom'] = np.array(chrom[::-1])
    RV['exons'] = np.array(exons[::-1])
    RV['strand'] = np.array(strand[::-1])
    RV['tran_id'] = np.array(tran_id[::-1])
    RV['gene_id'] = np.array(gene_id[::-1])
    RV['biotype'] = np.array(biotype[::-1])
    RV['exon_num'] = np.array(exon_num[::-1])
    RV['gene_name'] = np.array(gene_name[::-1])
    RV['tran_stop'] = np.array(tran_stop[::-1], "int")
    RV['tran_start'] = np.array(tran_start[::-1], "int")

    return RV