# This is a main file to run the isplice software, which will return the 
# splicing ratio for each gene.

import sys
import h5py
import pysam
import numpy as np
from optparse import OptionParser
from models.model_static import Psi_MCMC_MH
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads
from utils.tran_utils import TranUnits, Transcript, TranSplice

def get_CI(data, percent=0.95):
    if len(data.shape) == 0:
        data = data.reshape(-1,1)
    for k in range(data.shape[1]):
        data[:,k] = np.sort(data[:,k])
        CI_idx = int(data.shape[0] * (1-percent)/2)
        RV = data[-CI_idx,:] - data[CI_idx,:]
    return RV

def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help="The indexed alignement file in bam/sam format")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The parameter file for bias in hdf5 format")
    parser.add_option("--bias_mode", dest="bias_mode", default="unif", 
        help="The bias mode")
    parser.add_option("--gene_file", "-g", dest="gene_file", default=None,
        help="The list genes in use.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="estimated_rate.psi", help="The output file in txt format")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    else:
        anno = load_annotation(options.anno_file)
    if options.sam_file == None:
        print "Error: need --sam_file for reads indexed and aliged reads."
        sys.exit(1)
    else:
        sam_list = options.sam_file.split("---")
        samFiles = []
        for s in sam_list:
            samFiles.append(load_samfile(s))

    if options.gene_file == None:
        gene_list = anno["gene_id"]
    else:
        gene_list = np.loadtxt(options.gene_file, dtype="str")
    
    ref_file  = options.ref_file
    out_file  = options.out_file
    bias_file = options.bias_file
    bias_mode = options.bias_mode

    if ref_file is None: fastaFile = None
    else: fastaFile = FastaFile(ref_file)

    if bias_file is None: biasFile = None
    else: biasFile = BiasFile(bias_file)

    print "Welcome to diceseq static model! %d genes are under estimating." %len(gene_list)

    fid = open(out_file, "w")
    head_line = "gene_id\tgene_name\t"
    head_line = head_line + "count\tpsi_iso1\tpsi_iso2\tCI95_iso1\tCI95_iso2\n"
    fid.writelines(head_line)
    
    # 1 set the parameters for Psi inference MCMC method
    M = 100000
    cov_fixed = np.array([[0.15]])
    alpha_diri = np.ones(2) / 2

    # 2. run the model
    Psi_all = np.zeros((gene_list.shape[0],2))
    CIs_all = np.zeros((gene_list.shape[0],2))
    for g in range(gene_list.shape[0]):
        i = np.where(np.array(anno["gene_id"]) == gene_list[g])[0][0]
        Psi_all[g,:] = None

        _gene   = anno["gene_id"][i]
        # _exons  = anno["exons"][i]
        _strand = anno["strand"][i]
        _chrom  = str(anno["chrom"][i])
        _start  = anno["gene_start"][i]
        _stop   = anno["gene_stop"][i]
        # print g, _gene, _chrom, _start, _stop, _strand

        # t = TranSplice(_gene, _chrom, _strand, _start, _stop, _exons)
        # t.initial_unitSet()
        t = TranSplice(anno["genes"][i])
        t.initial_unitSet(add_premRNA=True)
        for samFile in samFiles:
            t.set_reads(samFile)

        if bias_mode != "unif":
            if fastaFile is None:
                print ("we re-set bias mode as uniform as reference " +
                "sequence are not available.")
                bias_mode = "unif"
            elif biasFile is None:
                print ("we re-set bias mode as uniform as bias " +
                "parameters are not available.")
                bias_mode = "unif"
            else:
                t.set_sequence(fastaFile)
                t.set_bias(biasFile)

        t.get_ready(bias_mode)
        count = len(t.read1p) + len(t.read1u) + len(t.read2u)
        if count == 0:
            print "No reads for %s" %anno["gene_id"][i]
            tmp = []; tmp[:] = ["nan"] * 5; tmp[0] = "0"
            a_line = anno["gene_id"][i] + "\t" + anno["gene_name"][i]
            a_line = a_line + "\t" + "\t".join(tmp) + "\n"
            fid.writelines(a_line)
            continue

        elif bias_mode == "unif":
            Psi_tmp = Psi_MCMC_MH(t.Rmat, t.efflen_unif, t.proU, cov_fixed,
                alpha_diri, M)
        else:
            Psi_tmp = Psi_MCMC_MH(t.Rmat, t.efflen_bias, t.proB, cov_fixed,
                alpha_diri, M)

        Psi_all[g,:] = np.mean(Psi_tmp[Psi_tmp.shape[0]/2:,:], axis=0)
        CIs_all[g,:] = get_CI(Psi_tmp[Psi_tmp.shape[0]/2:,:])

        # 3 save output
        rInfo = np.append(Psi_all[g,:], CIs_all[g,:])
        rInfo = list(np.append(count, rInfo))
        aLine = anno["gene_id"][i] + "\t" + anno["gene_name"][i]
        aLine = aLine + "\t" + "\t".join(["%.3f" %x for x in rInfo]) + "\n"
        fid.writelines(aLine)
    fid.close()

if __name__ == "__main__":
    main()
    
