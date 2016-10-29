# This module is an interface of the diceseq output,
# with which you could easily access to the results,
# also visualize the data.

import sys
import gzip
import numpy as np
import pylab as pl

class SampleFile:
    """docstring for SampleFile"""
    def __init__(self, sample_file):
        fid = gzip.open(sample_file, "r")
        all_lines = fid.readlines()
        fid.close()
        sample = []
        self.gene, self.tran = [], []
        self.theta, self.sample = [], []
        for line in all_lines:
            line = line.split("\n")[0]
            if len(self.gene) == 0 and line[0] != "@": continue
            if line.startswith("@"):
                self.gene.append(line[1:].split("|")[0])
                self.tran.append(line[1:].split("|")[1].split(","))
                self.theta.append(np.array(line[1:].split("|")[2].split(","),"float"))
                if sample == []: continue
                self.sample.append(np.array(sample, "float"))
                sample = []
            else:
                _sample = line.split(";")
                T = len(line.split(";"))
                for i in range(len(_sample)):
                    _sample[i] = _sample[i].split(",")
                sample.append(_sample)
        self.sample.append(np.array(sample, "float"))

    def GP_visual(self, gquery, T, xx, thetas=None):
        """visualization function for Gaussian process"""
        if self.gene.count(gquery) == 0:
            print("No %s in the sample file." %gquery)
            sys.exit(1)
        else: qidx = self.gene.index(gquery)
        _tran = self.tran[qidx]
        _sample = self.sample[qidx]
        
        yy = np.zeros((_sample.shape[0], xx.shape[0], len(_tran)))
        for c in range(len(_tran)):
            if c < len(_tran)-1:
                if thetas is None: 
                    theta1 = 3.0
                    theta2 = self.theta[qidx][c]
                    Theta = [theta1, theta2]
                elif len(np.array(thetas).shape) == 1:
                    Theta = thetas
                else:
                    Theta = thetas[c]    
            yy[:,:,c] = GP_sample_predict(Theta, T, _sample[:,:,c], xx)

        psi_mean = YtoPsi(_sample, axis=2).mean(axis=0)
        psi_pred = YtoPsi(yy, axis=2).mean(axis=0)
        psi_CI95 = np.zeros((len(T), len(_tran), 2))
        psi_pred95 = np.zeros((xx.shape[0], len(_tran), 2))
        for t in range(len(T)):
            psi_CI95[t,:,:] = get_CI(YtoPsi(_sample[:,t,:], axis=1))
        for t in range(xx.shape[0]):
            psi_pred95[t,:,:] = get_CI(YtoPsi(yy[:,t,:], axis=1))

        np.random.seed(10)
        for c in range(len(_tran)):
            _color = np.random.rand(3,1)
            psi_err = [psi_mean[:,c]-psi_CI95[:,c,1], psi_CI95[:,c,0]-psi_mean[:,c]]
            pl.errorbar(T, psi_mean[:,c], psi_err, fmt='.', color=_color, markersize=7)
            pl.plot(xx, psi_pred[:,c], '-', color=_color, label=_tran[c])
            pl.fill(np.concatenate([xx, xx[::-1]]), np.concatenate([psi_pred95[:,c,0],psi_pred95[:,c,1][::-1]]), 
                    alpha=.15, fc=_color, ec='None')

        pl.xlabel('$t$')
        pl.ylabel('$\psi(t)$')
        pl.ylim(-0.02, 1.02)
        pl.xticks(T, T)
        pl.title(self.gene[qidx])
        pl.legend(loc="best")
        
class DiceFile:
    """docstring for DiceFile"""
    def __init__(self, file):
        dice_data = np.loadtxt(file, skiprows=1, delimiter="\t", dtype="str")
        idx = np.arange(int((dice_data.shape[1]-4) / 4))
        self.tran = dice_data[:,0]
        self.gene = dice_data[:,1]
        self.FPKM = dice_data[:,idx*4+4].astype("float")
        self.LogLike = dice_data[:,2].astype("float")
        self.tranLen = dice_data[:,3].astype("float")
        self.IsoRatio = dice_data[:,idx*4+5].astype("float")

        self.IsoCI95s = np.zeros((dice_data.shape[0], len(idx), 2))
        self.IsoCI95s[:,:,0] = dice_data[:,idx*4+6].astype("float")
        self.IsoCI95s[:,:,1] = dice_data[:,idx*4+7].astype("float")

def id_mapping(IDs1, IDs2):
    """
    Mapping IDs2 to IDs1, both of which should only contain unique ids.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = np.argsort(IDs1)
    idx2 = np.argsort(IDs2)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RV_idx1)
    RV_idx = np.array(RV_idx2)[origin_idx]
    return RV_idx

def get_CI(data, percent=0.95):
    """calculate the confidence intervals
    """
    if len(data.shape) == 0:
        data = data.reshape(-1,1)
    RV = np.zeros((data.shape[1],2))
    CI_idx = int(data.shape[0] * (1-percent)/2)
    for k in range(data.shape[1]):
        temp = np.sort(data[:,k])
        RV[k,:] = [temp[-CI_idx], temp[CI_idx]]
    return RV

def GP_K(X, theta, x=None):
    """calculate the covariance matrix
    """
    N = len(X)
    if x is None:
        K = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                K[i,j] = theta[0] * np.exp(-0.5/theta[1]*(X[i]-X[j])**2)
    else:
        K = np.zeros(N)
        for i in range(N):
            K[i] = theta[0] * np.exp(-0.5/theta[1]*(X[i]-x)**2)
    return K

def GP_predict(Theta, X, Y, Yvar, xx):
    """predict by GP with its parameters
    """
    if len(np.array(xx).shape) == 0:
        xx = np.array([xx])
    yy = np.zeros(xx.shape[0])
    ss = np.zeros(xx.shape[0])

    K0 = GP_K(X, Theta)
    for i in range(K0.shape[0]):
        K0[i,i] += Yvar[i]
    K2 = Theta[0]

    for i in range(xx.shape[0]):
        K1 = GP_K(X, Theta, xx[i])
        yy[i] = np.dot(np.dot(K1, np.linalg.inv(K0)), Y)
        ss[i] = K2 - np.dot(np.dot(K1, np.linalg.inv(K0)), np.transpose(K1))
    return yy, ss

def GP_sample_predict(Theta, X, Yall, xx):
    """predict by GP with its all samples
    """
    if len(np.array(xx).shape) == 0:
        xx = np.array([xx])
    yy = np.zeros((Yall.shape[0], xx.shape[0]))

    K0 = GP_K(X, Theta)
    K2 = Theta[0]

    for i in range(xx.shape[0]):
        idx = np.where(abs(xx[i] - X)<10**(-10))[0]
        if len(idx) == 1:
            yy[:,i] = Yall[:,idx[0]]
            continue
        K1 = GP_K(X, Theta, xx[i])
        for j in range(Yall.shape[0]):
            _y = np.dot(np.dot(K1, np.linalg.inv(K0)), Yall[j,:])
            _v = K2 - np.dot(np.dot(K1, np.linalg.inv(K0)), np.transpose(K1))
            yy[j,i] = np.random.normal(_y, np.sqrt(_v))
    return yy

def YtoPsi(y, axis=0):
    """softmax function: transfer y to psi
    """
    size = np.array(y.shape)
    size[axis] = -1
    return np.exp(y) / (np.exp(y).sum(axis=axis).reshape(size))