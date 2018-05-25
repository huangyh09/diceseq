# This file contains some random but maybe useful stuff

import numpy as np
from .out_utils import id_mapping

def get_fraction(gene_ids, FPKM, ignoreNan=False):
    """Get the fraction from FPKM"""
    idx0 = 0
    frac = np.zeros(len(FPKM))
    for i in range(len(FPKM)+1):
        if i >= len(FPKM) or gene_ids[idx0] != gene_ids[i]:
            FPKM_sum = float(np.sum(FPKM[idx0:i]))
            if FPKM_sum == 0:
                if ignoreNan == True:
                    frac[idx0:i] = None
                else:
                    frac[idx0:i] = 1.0/(i-idx0)
            else:
                frac[idx0:i] = FPKM[idx0:i] / FPKM_sum
            idx0 = i
    return frac

def loadresult(file_name, tran_ids, gene_ids=None, method="dice", 
    ignoreNan=False):
    """
    Load results for isoforms from different methods.
    ------------------------------------------------

    Methods supported: DICEseq, BRIE, MISO, Cufflinks, Kallisto, flux, spanki
    """
    frac = np.zeros(len(tran_ids))
    FPKM = np.zeros(len(tran_ids))
    CI95 = np.zeros(len(tran_ids)) + 1.0
    if ["dice", "brie"].count(method) == 1:
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        idxT = id_mapping(tran_ids, data[:,0])
        frac = data[idxT, 5].astype(float)
        FPKM = data[idxT, 4].astype(float)
        
    elif method == "miso":
        tran_miso = []
        frac_miso = np.array([])
        CI95_miso = np.array([])
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        for gn in range(data.shape[0]):
            _tran_ids = data[gn,4].split(",")
            for tr in range(len(_tran_ids)):
                #_tran_ids[tr] = ".".join(_tran_ids[tr].split(".")[:2])[1:]
                _tran_ids[tr] = _tran_ids[tr].split("_")[0][1:-2]

            _frac = np.array(data[gn,1].split(","),"float")
            _CI95 = np.array(data[gn,3].split(","),"float")
            _CI95 = _CI95 - np.array(data[gn,2].split(","),"float")
            if len(_tran_ids) == 2:
                _frac = np.append(_frac, 1-_frac[-1])
                _CI95 = np.append(_CI95, _CI95[-1])
            tran_miso += _tran_ids
            frac_miso = np.append(frac_miso, _frac)
            CI95_miso = np.append(CI95_miso, _CI95)
            # FPKM_miso = [float(x.split(":")[1]) for x in data[gn, 6].split(",")]
        
        tran_miso = np.array(tran_miso, dtype="str")
        idxT = id_mapping(tran_ids, tran_miso)        
        for j in range(len(FPKM)):
            if idxT[j] >= 0: 
                FPKM[j] = frac_miso[idxT[j]]
                CI95[j] = CI95_miso[idxT[j]]
            else:
                FPKM[j] = 0.0
                CI95[j] = 1.0
        frac = get_fraction(gene_ids, FPKM, ignoreNan)
                
    elif method == "kallisto":
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        idxT = id_mapping(tran_ids, [x.split("|")[1] for x in data[:,0]])
        FPKM = data[idxT, 4].astype(float)
        frac = get_fraction(gene_ids, FPKM, ignoreNan)

    elif method == "rsem":
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        idxT = id_mapping(tran_ids, [x.split("|")[1] for x in data[:,0]])
        FPKM = data[idxT, 6].astype(float)
        frac = get_fraction(gene_ids, FPKM, ignoreNan)
        
    elif method == "cufflinks":
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        idxT = id_mapping(tran_ids, data[:,0])
        for j in range(len(FPKM)):
            if idxT[j] >= 0: 
                FPKM[j] = data[idxT[j], 9].astype(float)
            else:
                if ignoreNan == False:
                    FPKM[j] = 0.0
                else:
                    FPKM[j] = None
        frac = get_fraction(gene_ids, FPKM, ignoreNan)
                
    elif method == "flux":
        data = np.genfromtxt(file_name, skip_header=0, dtype="str")
        idxT = id_mapping(tran_ids, data[:,1])
        FPKM = data[idxT, 5].astype(float)
        frac = get_fraction(gene_ids, FPKM, ignoreNan)
        
    elif method == "spanki":
        data = np.genfromtxt(file_name, skip_header=1, dtype="str")
        idxT = id_mapping(tran_ids, data[:,1])
        FPKM = data[idxT, 3].astype(float)
        frac = get_fraction(gene_ids, FPKM, ignoreNan)
        
    else:
        print("Unsupported method: %")
        frac[:], FPKM[:] = None, None
    
    return frac, FPKM