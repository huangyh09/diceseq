# This is a direct running file to plot the parameters of sequencing bias.

import h5py
import pysam
import pylab as pl
import numpy as np
from optparse import OptionParser
from utils.reads_utils import ReadSet
from utils.gtf_utils import load_annotation
from utils.sam_utils import load_samfile, fetch_reads

def lognorm_pdf(x, mu, sigma):
    return 1 / (sigma*x*np.sqrt(2*np.pi)) * np.exp(-1/2*((np.log(x)-mu)/sigma)**2)

def norm_pdf(x, mu, sigma):
    return 1 / (sigma*np.sqrt(2*np.pi)) * np.exp(-1/2*((x-mu)/sigma)**2)

def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--bias_file", "-b", dest="bias_file", default=None,
        help="The bias parameters file in hdf5 format")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="./sequence_position_bias.pdf", help="The name of the figure.")

    (options, args) = parser.parse_args()
    if options.bias_file == None:
        print "Error: need --bias_file for bias parameters."
        sys.exit(1)
    
    bias_file = options.bias_file
    out_file  = options.out_file

    # loading data
    f = h5py.File(bias_file, "r")
    percentile = np.array(f["percentile"])
    pos5_bias = np.array(f["pos5_bias"]) / np.array(f["pos5_unif"])
    pos3_bias = np.array(f["pos3_bias"]) / np.array(f["pos3_unif"])
    chain_len = list(f["chain_len"])

    idx5 = np.where(pos5_bias != pos5_bias)
    idx3 = np.where(pos3_bias != pos3_bias)
    pos5_bias[idx5] = 0
    pos3_bias[idx3] = 0

    base_chain = {}
    seq5_bias, seq3_bias = {}, {}
    seq5_unif, seq3_unif = {}, {}
    for i in range(21):
        seq5_bias[str(i)] = np.array(f["seq5_bias/"+str(i)]) / np.array(f["seq5_unif/"+str(i)])
        seq3_bias[str(i)] = np.array(f["seq3_bias/"+str(i)]) / np.array(f["seq3_unif/"+str(i)])
        seq5_unif[str(i)] = np.array(f["seq5_unif/"+str(i)])
        seq3_unif[str(i)] = np.array(f["seq3_unif/"+str(i)])
        base_chain[str(i)] = np.array(f["base_chain/"+str(i)])

        print i, np.array(f["seq5_unif/"+str(i)]).sum(), np.array(f["seq5_bias/"+str(i)]).sum()
    f.close()

    pl.rcParams['pdf.fonttype'] = 42
    pl.rcParams['font.sans-serif'] = "Arial"

    # #fragment distribution
    # fig = pl.figure()
    # xx = np.arange(80, 400)
    # # pl.subplot(1,2,1)
    # # for i in range(5):
    # #     mu, sigma = flen_bias["normal"][i,:]
    # #     yy = norm_pdf(xx, mu, sigma)
    # #     pl.plot(xx, yy, linewidth=2.0, label="bin %d" %i)
    # # pl.legend(loc="best")
    # # pl.title("normal distribution")

    # # pl.subplot(1,2,2)
    # for i in range(5):
    #     mu, sigma = lognormal[i,:]
    #     yy = lognorm_pdf(xx, mu, sigma)
    #     pl.plot(xx, yy, linewidth=2.0, label="bin %d" %i)
    # pl.legend(loc="best")
    # pl.title("lognormal distribution")

    #position bias
    fig = pl.figure()
    pl.subplot(2,2,1)
    pl.plot(np.arange(20), np.ones(20), '--k')
    for i in range(5):
        if i == 4:
            _label="bin%d: " %(i+1) + ">=%d bp" %int(percentile[i,0])
        else :
            _label="bin%d: " %(i+1) + "%d-%d bp" %(int(percentile[i,0]), int(percentile[i,1]))
        pl.plot(np.arange(20), pos5_bias[i,:], linewidth=2.0, label=_label)
    pl.legend(loc="best")
    pl.xlabel("fractional transcription position")
    pl.ylabel("position bias weight")
    pl.ylim(0,2)

    pl.subplot(2,2,2)
    pl.plot(np.arange(20), np.ones(20), '--k')
    for i in range(5):
        if i == 4:
            _label="bin%d: " %(i+1) + ">=%d bp" %int(percentile[i,0])
        else :
            _label="bin%d: " %(i+1) + "%d-%d bp" %(int(percentile[i,0]), int(percentile[i,1]))
        pl.plot(np.arange(20), pos3_bias[i,:], linewidth=2.0, label= _label)
    pl.legend(loc="best")
    pl.xlabel("fractional transcription position")
    pl.ylabel("position bias weight")
    pl.ylim(0,2)
    # fig.set_size_inches(12,4.5)
    # fig.savefig(out_dir + "/position_bias.pdf",dpi=300, bbox_inches='tight')

    #sequence bias
    # fig = pl.figure()
    base = ["A", "T", "G", "C"]
    pl.subplot(2,2,3)
    pl.plot(np.arange(21)-8, np.ones(21), '--k')
    pl.plot(np.zeros(2), np.array([0, 2.0]), '--k', linewidth=2.0)
    percent = np.zeros((4,21))
    for i in range(4):
        for j in range(21):
            _seq_bias = np.array(seq5_bias[str(j)])
            percent[i,j] = np.sum(_seq_bias[i*4**(chain_len[j]-1): (i+1)*4**(chain_len[j]-1)]) / 4**(chain_len[j]-1) 
        pl.plot(np.arange(21)-8, percent[i,:], ":o", linewidth=1.0, label="%s" %base[i])

    pl.legend(loc="best")
    pl.xlabel("offset from 5' fragment end")
    pl.ylabel("sequence bias weight")
    pl.xlim(-8,12)
    pl.ylim(0.5,2)

    pl.subplot(2,2,4)
    pl.plot(np.arange(21)-12, np.ones(21), '--k')
    pl.plot(np.zeros(2), np.array([0, 2.0]), '--k', linewidth=2.0)
    percent = np.zeros((4,21))
    for i in range(4):
        for k in range(21):
            j = 20 - k
            percent[i,j] = np.sum(np.array(seq3_bias[str(j)])[i*4**(chain_len[j]-1):(i+1)*4**(chain_len[j]-1)]) / 4**(chain_len[j]-1)
        pl.plot(np.arange(21)-12, percent[i,:], ":o", linewidth=1.0, label="%s" %base[i])
    pl.legend(loc="best")
    pl.xlabel("offset from 3' fragment end")
    pl.ylabel("sequence bias weight")
    pl.xlim(-12,8)
    pl.ylim(0.5,2)
    fig.set_size_inches(12,9.5)
    fig.savefig(out_file, dpi=300, bbox_inches='tight')
    #pl.show()

if __name__ == "__main__":
    main()
    