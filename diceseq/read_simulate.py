# This file is for generating reads for simulation

import sys
import gzip
import h5py
import pysam
import numpy as np
from optparse import OptionParser
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads
from utils.tran_utils import TranUnits, Transcript, TranSplice


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
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list genes in use. Default is all genes in annotation.")
    parser.add_option("--add_premRNA", dest="add_premRNA", default="False",
        help="Whether adding pre-mRNA or not.")
    parser.add_option("--RPK", dest="RPK", default="1000", 
        help="The all used sequence depths, e.g., 100,200,400 and 100")
    parser.add_option("--ratio", dest="ratio", default="0.5", 
        help="The ratio of the first isoform, default=0.5.")
    parser.add_option("--noise", dest="noise", default="0.001",
        help="The noise in the reads number for each isoform, default=0.001.")
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

    if options.gene_file == None:
        gene_list = anno["gene_id"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")

    bias_mode = options.bias_mode
    if options.bias_file is None:
        biasFile = None
        if bias_mode != "unif":
            bias_mode = "unif"
            print "No bias parameter file, so we change to uniform mode."
    else: biasFile = BiasFile(options.bias_file)

    RPK = options.RPK.split(",")
    ratio = options.ratio.split(",")
    Np = np.zeros((len(RPK), len(ratio)), "int")
    Nm = np.zeros((len(RPK), len(ratio)), "int")
    cnt = np.zeros((len(RPK), len(ratio)), "int")

    rlen = int(options.rlen)
    noise = float(options.noise)
    fl_mean = float(options.fl_mean)
    fl_sigma = float(options.fl_sigma)
    add_premRNA = options.add_premRNA == "True"

    qua = []
    qua[:] = "~" * rlen
    qua = "".join(qua)
    for g in range(len(gene_list)):
        if g == 0: ftype = "w"
        else: ftype = "a+"
        i = np.where(np.array(anno["gene_id"]) == gene_list[g])[0][0]
        if add_premRNA == True: anno["genes"][i].add_premRNA()

        unitp = TranUnits(anno["genes"][i].trans[1])
        unitm = TranUnits(anno["genes"][i].trans[0])
        unitp.set_sequence(fastaFile)
        unitm.set_sequence(fastaFile)
        unitp.set_bias(biasFile)
        unitm.set_bias(biasFile)

        for m in range(len(RPK)):
            for n in range(len(ratio)):
                _Np = (unitp.ulen) / 1000.0 * float(RPK[m]) * (1-float(ratio[n]))
                _Nm = (unitm.ulen) / 1000.0 * float(RPK[m]) * float(ratio[n])
                _aa = np.random.rand()
                _bb = np.random.rand()
                if _aa <= _Np - int(_Np): _Np = int(_Np) + 1.0
                else: _Np = int(_Np) + 0.0
                if _bb <= _Nm - int(_Nm): _Nm = int(_Nm) + 1.0
                else: _Nm = int(_Nm) + 0.0
                
                Np[m,n] = max(0, np.round(np.random.normal(_Np, noise*_Np)))
                Nm[m,n] = max(0, np.round(np.random.normal(_Nm, noise*_Nm)))
        readp = unitp.generate_reads(Np.max(),rlen,fl_mean,fl_sigma,bias_mode)
        readm = unitm.generate_reads(Nm.max(),rlen,fl_mean,fl_sigma,bias_mode)

        for m in range(len(RPK)):
            for n in range(len(ratio)):
                idxp = np.random.permutation(Np.max())[:Np[m,n]]
                idxm = np.random.permutation(Nm.max())[:Nm[m,n]]
                fid1 = gzip.open(options.out_file + "_%s_%s_1.fq.gz" %(RPK[m], ratio[n]), ftype)
                fid2 = gzip.open(options.out_file + "_%s_%s_2.fq.gz" %(RPK[m], ratio[n]), ftype)
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
        if g == 0 or (g+1) % 10 == 0 or (g+1) == len(gene_list):
            print "%d out %d genes have been finished." % (g+1, len(gene_list))

if __name__ == "__main__":
    main()
