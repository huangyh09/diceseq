# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import h5py
import pysam
import numpy as np
from optparse import OptionParser
from models.model_GP import Psi_GP_MH
from models.model_static import Psi_MCMC_MH
from utils.gtf_utils import load_annotation
from utils.bias_utils import BiasFile, FastaFile
from utils.sam_utils import load_samfile, fetch_reads
from utils.tran_utils import TranUnits, Transcript, TranSplice

def get_CI(data, percent=0.95):
    if len(data.shape) == 0:
        data = data.reshape(-1,1)
    RV = np.zeros(data.shape[1])
    CI_idx = int(data.shape[0] * (1-percent)/2)
    for k in range(data.shape[1]):
        temp = np.sort(data[:,k])
        RV[k] = temp[-CI_idx] - temp[CI_idx]
    return RV

def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format")
    parser.add_option("--sam_list", "-s", dest="sam_list", default=None,
        help="The indexed alignement file in bam/sam format")
    parser.add_option("--ref_file", "-r", dest="ref_file", default=None,
        help="The genome reference file in fasta format")
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The parameter file for bias in hdf5 format")
    parser.add_option("--gene_file", "-g", dest="gene_file",
        help="The list genes in use.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="estimated_rate.psi", help="The output file in txt format")
    parser.add_option("--bias_mode", dest="bias_mode", default="unif", 
        help="The bias mode")
    parser.add_option("--merge_mode", dest="merge_mode", default="single", 
        help="The merge mode: single or multiple")

    (options, args) = parser.parse_args()
    if options.anno_file == None:
        print "Error: need --anno_file for annotation."
        sys.exit(1)
    if options.sam_list == None:
        print "Error: need --sam_list for reads indexed and aliged reads."
        sys.exit(1)
    
    ref_file  = options.ref_file
    sam_list  = options.sam_list
    out_file  = options.out_file
    bias_file = options.bias_file
    anno_file = options.anno_file
    gene_file = options.gene_file
    bias_mode = options.bias_mode
    merge_mode = options.merge_mode

    if ref_file is None: fastaFile = None
    else: fastaFile = FastaFile(ref_file)

    if bias_file is None: biasFile = None
    else: biasFile = BiasFile(bias_file)

    sam_list  = sam_list.split("---")
    anno      = load_annotation(anno_file)
    gene_list = np.loadtxt(gene_file, dtype="str")

    samFiles  = []
    for i in range(len(sam_list)):
        samFiles.append(load_samfile(sam_list[i]))
    print len(anno["gene_name"])
    T = len(sam_list)
    K = 2
    
    # 1. set the parameters
    M = 100000
    cov_Y = np.ones(2) * 0.1
    sca_theta = np.ones(2) * 0.1
    
    # 2. run the model
    Psi_all = np.zeros((gene_list.shape[0], K, T))
    CIs_all = np.zeros((gene_list.shape[0], K, T))
    for g in range(gene_list.shape[0]):
        #g = 90 # "RPL28"#"APS3" # 27 RPL28; 43 RPL39
        i = np.where(np.array(anno["gene_name"]) == gene_list[g])[0][0]
        Psi_all[g, :, :] = None

        _gene   = anno["gene_id"][i]
        # _exons  = anno["exons"][i]
        _strand = anno["strand"][i]
        _chrom  = str(anno["chrom"][i])
        _start  = anno["gene_start"][i]
        _stop   = anno["gene_stop"][i]

        # print gene_list[g], _chrom, _exons
        R_all, len_iso_all, prob_iso_all = [], [], []
        for j in range(len(sam_list)):
            # t = TranSplice(_gene, _chrom, _strand, _start, _stop, _exons)
            t = TranSplice(anno["genes"][i])
            # t.quick_start(samFile, fastaFile, biasFile)
            t.initial_unitSet(add_premRNA=True)
            t.set_reads(samFiles[j])
            if bias_mode != "unif": 
                t.set_sequence(fastaFile)
                t.set_bias(biasFile)
            t.get_ready(bias_mode)
            R_all.append(t.Rmat)
            if bias_mode == "unif": 
                len_iso_all.append(t.efflen_unif)
                prob_iso_all.append(t.proU)
            else: 
                len_iso_all.append(t.efflen_bias)
                prob_iso_all.append(t.proB)

        if merge_mode == "multiple":
            X = np.arange(len(R_all))
            Psi_tmp = Psi_GP_MH(R_all, len_iso_all, prob_iso_all, X, cov_Y, sca_theta, M)
            Psi_all[g,:,:] = np.mean(Psi_tmp[Psi_tmp.shape[0]/2:,:,:], axis=0)
            for s in range(T):
                CIs_all[g,:,s] = get_CI(Psi_tmp[Psi_tmp.shape[0]/2:,:,s])
        else:
            for s in range(T):
                X = [0]
                Psi_tmp = Psi_GP_MH(R_all[s:s+1], len_iso_all[s:s+1], prob_iso_all[s:s+1], X, cov_Y, sca_theta, M)
                Psi_all[g,:,s] = np.mean(Psi_tmp[Psi_tmp.shape[0]/2:,:,0], axis=0)
                CIs_all[g,:,s] = get_CI(Psi_tmp[Psi_tmp.shape[0]/2:,:,0])
        
        print Psi_all[g, :, :]
        print CIs_all[g, :, :]

    fout = h5py.File(out_file, "w")
    fout["Psi_all"] = Psi_all
    fout["CIs_all"] = CIs_all
    fout.close()


if __name__ == "__main__":
    main()
    