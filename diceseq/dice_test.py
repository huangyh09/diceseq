# This file is to test the performance of the HM sampler,
# for both single and multiple time points.

import time
import numpy as np
from optparse import OptionParser
from models.model_GP import Psi_GP_MH

if __name__ == "__main__":
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--mode", dest="mode", default="single",
        help="The mode of the model, single or multiple.")
    parser.add_option("--model", "-m", dest="model", default="GP",
        help="The model in use: GP, static, or dynamic.")

    (options, args) = parser.parse_args()
    mode  = options.mode
    model = options.model 

    print "Welcome to diceseq sampler test!"
    print "Current model is %s with mode %s." %(model, mode)


    R_all = []
    Psi_all = []
    len_iso_all = []
    prob_iso_all = []

    xx = np.arange(0, 10, 1)
    y1 = 1 - np.exp(-0.25*(xx+0.5))

    for i in range(0, xx.shape[0]):
        Psi = np.array([y1[i], 1-y1[i]])

        N  = 50
        perc = np.array([[2.0/8, 2.0/7, 2.0/5], [3.0/8, 3.0/7, 0], [1.0/8, 0, 1.0/5], [2.0/8, 2.0/7, 2.0/5]]) * np.array([8,7,5])
        perc = np.array([[1.0/6, 1.0/9], [0, 3.0/9], [5.0/6, 5.0/9]]) * np.array([6,9])
        Num = N * (Psi * perc).sum(axis=1)
        N = Num.sum()

        R_mat = np.zeros((N, 2), "bool")
        R_mat[Num[:0].sum():Num[:1].sum(),:] = [True, True]
        R_mat[Num[:1].sum():Num[:2].sum(),:] = [False, True]
        R_mat[Num[:2].sum():Num[:3].sum(),:] = [True, True]

        len_isos  = np.array([600, 900])
        prob_isos = np.zeros((N, 2))
        prob_isos[Num[:0].sum():Num[:1].sum(),:] = 1.0 / len_isos
        prob_isos[Num[:1].sum():Num[:2].sum(),:] = 1.0 / len_isos * np.array([0,1])
        prob_isos[Num[:2].sum():Num[:3].sum(),:] = 1.0 / len_isos

        Psi_all.append(Psi)
        R_all.append(R_mat)
        len_iso_all.append(len_isos)
        prob_iso_all.append(prob_isos)

    print "%d total simulated reads are used here." %Num.sum()

    if model == "GP":
        start_time = time.time()
        M = 100000
        cov_Y = np.ones(2) * 0.1
        sca_theta = np.ones(2) * 0.1
        if mode == "multiple":
            X = np.arange(len(R_all[:10]))
            Psi_all = Psi_GP_MH(R_all[:10], len_iso_all[:10], prob_iso_all[:10], X, cov_Y, sca_theta, M)
            temp = np.mean(Psi_all[Psi_all.shape[0]/2:,0,:], axis=0)
        else:
            temp = np.array([])
            for i in range(10):
                X = [0]
                Psi_all = Psi_GP_MH(R_all[i:i+1], len_iso_all[i:i+1], prob_iso_all[i:i+1], X, cov_Y, sca_theta, M)
                temp = np.append(temp, np.mean(Psi_all[Psi_all.shape[0]/2:,0,:], axis=0))

    print "The true isoform proportions are:"
    print y1
    print "The estimated isoform proportions:"
    print temp
    print("--- running %.2f seconds ---" % (time.time() - start_time))

