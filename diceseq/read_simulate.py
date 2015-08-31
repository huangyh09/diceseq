# This file is for generating reads for simulation


import h5py
import pysam
import numpy as np
from optparse import OptionParser
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads
from utils.tran_utils import TranUnits, Transcript, TranSplice


def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The parameter file for bias in hdf5 format")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="simulate.fastq", help="The output file with simulated reads")
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list genes in use.")
    parser.add_option("--RPK", dest="RPK", default="1000", 
        help="The sequence depth")
    parser.add_option("--ratio", dest="ratio", default="0.5", 
        help="The ratio of mRNA")
    parser.add_option("--rlen", dest="rlen", default="100", 
        help="The length of reads")
    parser.add_option("--bias_mode", dest="bias_mode", default="both", 
        help="The bias mode; unif, end5, end3, both")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file)

    if options.gene_file == None:
        gene_list = anno["gene_id"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")

    ref_file  = options.ref_file
    out_file  = options.out_file
    bias_file = options.bias_file
    RPK       = int(options.RPK)
    rlen      = int(options.rlen)
    ratio     = float(options.ratio)
    bias_mode = options.bias_mode
    
    biasFile  = BiasFile(bias_file)
    fastaFile = FastaFile(ref_file)
    fid1      = open(out_file + "_1.fq", "w")
    fid2      = open(out_file + "_2.fq", "w")

    cnt = 0
    qua = []
    qua[:] = "~" * rlen
    qua = "".join(qua)
    for g in gene_list:
        i = np.where(np.array(anno["gene_id"]) == g)[0][0]
        _gene   = anno["gene_id"][i]
        _exons  = anno["exons"][i]
        _strand = anno["strand"][i]
        _chrom  = str(anno["chrom"][i]) 
        _start  = anno["tran_start"][i]
        _stop   = anno["tran_stop"][i]

        unitp = TranUnits(_chrom, _strand, [_start, _stop])
        
        unitp.set_sequence(fastaFile)
        unitp.set_bias(biasFile)
        Np    = (unitp.ulen - rlen + 1) / 1000.0 * RPK * (1-ratio)
        readp = unitp.generate_reads(Np, rlen, 200, 20, bias_mode)

        unitm = TranUnits(_chrom, _strand, _exons)
        unitm.set_sequence(fastaFile)
        unitm.set_bias(biasFile)
        Nm    = (unitm.ulen - rlen + 1) / 1000.0 * RPK * ratio
        readm = unitm.generate_reads(Nm, rlen, 200, 20, bias_mode)

        print g, Np, Nm
        for r in range(len(readp["reads1"])):
            cnt += 1
            fid1.writelines("@simulated_reads_%d_1\n" %cnt)
            fid2.writelines("@simulated_reads_%d_2\n" %cnt)
            fid1.writelines(readp["reads1"][r] + "\n+\n" + qua + "\n")
            fid2.writelines(readp["reads2"][r] + "\n+\n" + qua + "\n")
        for r in range(len(readm["reads1"])):
            cnt += 1
            fid1.writelines("@simulated_reads_%d_1\n" %cnt)
            fid2.writelines("@simulated_reads_%d_2\n" %cnt)
            fid1.writelines(readm["reads1"][r] + "\n+\n" + qua + "\n")
            fid2.writelines(readm["reads2"][r] + "\n+\n" + qua + "\n")
    fid1.close()
    fid2.close()

if __name__ == "__main__":
    main()
