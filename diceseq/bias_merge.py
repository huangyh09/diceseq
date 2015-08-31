# This file is to merge bias parameters from all bias files.

import h5py
import numpy as np
from optparse import OptionParser
from utils.bias_utils import BiasFile

def main():
   #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--bias_list", dest="bias_list", default=None,
        help="The list of bias file path, separated by ---")
    parser.add_option("--out_file", dest="out_file", default=None,
        help="The merged bias file for output.")

    #part 0.1. load the data
    (options, args) = parser.parse_args()
    file_list = list_file.split("---")
    out_file  = options.out_file

    #part 1.1. copy a file as output file
    N = len(file_list)
    biasFile  = BiasFile(file_list[0])
    biasFile.pos5_bias /= N
    biasFile.pos3_bias /= N
    biasFile.pos5_unif /= N
    biasFile.pos3_unif /= N
    for i in range(len(biasFile.chain_len)):
        biasFile.seq5_bias[str(i)] /= N
        biasFile.seq3_bias[str(i)] /= N
        biasFile.seq5_unif[str(i)] /= N
        biasFile.seq3_unif[str(i)] /= N

    #part 1.2. merge other files
    for n in np.arange(1, N):
        _biasFile = BiasFile(file_list[n])

        biasFile.pos5_bias += _biasFile.pos5_bias / N
        biasFile.pos3_bias += _biasFile.pos3_bias / N
        biasFile.pos5_unif += _biasFile.pos5_unif / N
        biasFile.pos3_unif += _biasFile.pos3_unif / N
        for i in range(len(biasFile.chain_len)):
            biasFile.seq5_bias[str(i)] += _biasFile.seq5_bias[str(i)] / N
            biasFile.seq3_bias[str(i)] += _biasFile.seq3_bias[str(i)] / N
            biasFile.seq5_unif[str(i)] += _biasFile.seq5_unif[str(i)] / N
            biasFile.seq3_unif[str(i)] += _biasFile.seq3_unif[str(i)] / N

    #part 2. save file
    biasFile.save_file(out_file) 

if __name__ == "__main__":
    main()
