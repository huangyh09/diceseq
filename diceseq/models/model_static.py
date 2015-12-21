import numpy as np

def logistic_normal_pdf(D, mu, cov, log=False):
    """
    Calculate the probability density of logistic-normal distribution, see
    Atchison J, Shen S M. Logistic-normal distributions: Some properties and 
    uses[J]. Biometrika, 1980, 67(2): 261-272.
    http://www.jstor.org/stable/2335470

    Parameters
    ----------
    D   : 1-D array_like, (K, )
        The variable for calculating the probability density. Range: [0,1].
    mu  : 1-D array_like, (K-1, )
        Mean of the logistic multinormal distribution. Range: [-inf,inf].
    cov : 1-D array_like, (K-1, K-1)
        Covariance of for the logistic multinormal distribution.
    log : bool
        If true, the return value is at log scale.

    Returns
    -------
    RV  : numpy float
        The probability density of D.
    """
    if len(D.shape) > 1:
        print("Only one array is calculated each time.")
        return None
    D_logit = np.log(D[:-1]/D[-1]) - mu
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    if cov_det < 0:
        print("The det of covariance is negative, please check!")
        return None
    RV = (-0.5*np.log(2*np.pi*cov_det) - np.sum(np.log(D)) - 
           0.5*np.dot(np.dot(D_logit, cov_inv), D_logit))
    if log == False: RV = np.exp(RV)
    return RV


def Dirichlet_pdf(D, alpha, log=False):
    """
    Calculate the probability density of Dirichlet distribution.

    Parameters
    ----------
    D : 1-D array_like, (K, )
        The variable for calculating the probability density. Range: [0,1].
    alpha : 1-D array_like, (K, )
        The parameter array for the Dirichlet distribution.
    log : bool
        If true, the return value is at log scale.

    Returns
    -------
    RV : numpy float
        The probability density of D.
    """
    if len(D.shape) > 1:
        print("Only one array is calculated each time.")
        return None
    part2 = 0
    for i in range(len(alpha)):
        part2 += np.math.lgamma(alpha[i])
    RV = np.math.lgamma(np.sum(alpha)) - part2 + np.sum((alpha-1)*np.log(D))
    if log == False: RV = np.exp(RV)
    return RV

def Geweke_Z(X, first=0.1, last=0.5):
    """
    Geweke diagnostics for MCMC chain convergence.
    See Geweke J. Evaluating the accuracy of sampling-based approaches to the 
    calculation of posterior moments[M]. Minneapolis, MN, USA: Federal Reserve 
    Bank of Minneapolis, Research Department, 1991.
    and https://pymc-devs.github.io/pymc/modelchecking.html#formal-methods

    Parameters
    ----------
    X : 1-D array_like, (N, )
        The uni-variate MCMC sampled chain for convergence diagnostic.
    first : float
        The proportion of first part in Geweke diagnostics.
    last : float
        The proportion of last part in Geweke diagnostics.

    Returns
    -------
    Z : float
        The Z score of Geweke diagnostics.
    """
    N = X.shape[0]
    A = X[:first*N]
    B = X[last*N:]
    Z = abs(A.mean() - B.mean()) / np.sqrt(np.var(A) + np.var(B))
    return Z

def Psi_MCMC_MH(R_mat, len_isos, prob_isos, cov_fixed, alpha_Dir, M=100000, start=3000, gap=1000):
    """
    Calculate the proportion of each isoforms with all reads via MCMC samplling
    Parameters
    ----------
    R_mat : 2-D array_like, (N, K)
        N reads identities of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.
    cov_fixed : 2-D array_like, (K-1, K-1)
        The fixed covariance of Logistic-Normal distribution.
    alpha_Dir : 1-D array_like, (K,)
        The parameter of Dirichlet distribution.
    M : int
        The iteration times of MH sampler, default=100.

    Returns
    -------
    RV : 2-D array_like, (M, K)
        M sampled Psi for abundance of each isoforms
    """
    # 1. check input data
    idx = (len_isos != len_isos)
    len_isos[idx] = 0.0
    prob_isos[:,idx] = 0.0
    R_mat[:,idx] = False

    idx = np.where(R_mat != R_mat)
    R_mat[idx] = False

    idx = np.where(prob_isos != prob_isos)
    prob_isos[idx] = 0.0

    idx = (R_mat.sum(axis=1) > 0) * (prob_isos.sum(axis=1) > 0)
    R_mat = R_mat[idx,:]
    prob_isos = prob_isos[idx,:]

    N, K = R_mat.shape

    # MCMC random initializations
    mu_now  = np.random.rand(K-1) - 0.5
    mu_now  = np.random.multivariate_normal(mu_now, cov_fixed)
    psi_now = np.zeros(K)
    psi_now[:-1] = np.exp(mu_now) / (1 + np.sum(np.exp(mu_now)))
    psi_now[-1]  = 1-np.sum(psi_now[:-1])
    fsi_now = len_isos*psi_now / np.sum(len_isos*psi_now)
    P_now   = (np.sum(np.log(np.dot(R_mat*prob_isos, fsi_now))) + 
               Dirichlet_pdf(psi_now, alpha_Dir, log=True))
    
    # MCMC running
    psi_try = np.zeros(K)
    Psi_all = np.zeros((M, K))
    for i in range(M):
        # step 1: propose a value
        mu_try  = np.random.multivariate_normal(mu_now, cov_fixed)
        psi_try[:-1] = np.exp(mu_try) / (1.0 + np.sum(np.exp(mu_try)))
        psi_try[-1]  = 1 - np.sum(psi_try[:-1])
        fsi_try = len_isos*psi_try / np.sum(len_isos*psi_try)
        P_try   = (np.sum(np.log(np.dot(R_mat*prob_isos, fsi_try))) + 
                   Dirichlet_pdf(psi_try, alpha_Dir, log=True))
        
        Q_now = logistic_normal_pdf(psi_now, mu_try, cov_fixed, log=True)
        Q_try = logistic_normal_pdf(psi_try, mu_now, cov_fixed, log=True)

        # step 2: calculate the MH ratio, accept or reject the proposal
        alpha = min(np.exp(P_try+Q_now-P_now-Q_try), 1)
        if alpha is None: print("alpha is none!")
        elif np.random.rand(1) < alpha:
            mu_now  = mu_try  + 0
            psi_now = psi_try + 0
            fsi_now = fsi_try + 0
            P_now   = P_try   + 0
        Psi_all[i,:] = psi_now

        #step 3. convergence diagnostics
        if i >= start and i % gap == 0:
            conv = 1
            for k in range(K):
                if Geweke_Z(Psi_all[:i, k]) > 2: 
                    conv = 0
                    break
            if conv == 1:
                Psi_all = Psi_all[:i,:]
                break

    return Psi_all


def Psi_junction(R_mat, len_juncts):
    """
    Calculate the proportion of each isoforms with junction reads only
    Parameters
    ----------
    R_mat : 2D matrix, (N, K)
        Parameters; N and K are the number of reads and isoforms, respectively.
    len_juncts : 1D array, (K,)
        Parameters; the used junction length of each isoforms
    Returns
    -------
    RV : 1D array, (K,)
    The the proportion of each isoforms
    """
    N, K = R_mat.shape
    if R_mat.shape[1] != len_juncts.shape[0]:
        print("The number of isoforms in the reads and junction lengths are different!")
        return [None] * K
    if np.sum(R_mat) == 0:
        print("No junction reads!")
        return [None] * K
    Dens = np.sum(R_mat, axis=0) / len_juncts.astype("float")
    RV = Dens / np.sum(Dens)
    return RV

def Psi_analytic(R_mat, len_isos):
    """
    Calculate the proportion of each isoforms with all reads
    Parameters
    ----------
    R_mat : 2D matrix, (N, K)
        Parameters; N and K are the number of reads and isoforms, respectively.
    len_isos : 1D array, (K,)
        Parameters; the used length of each isoforms
    Returns
    -------
    RV : 1D array, (K,)
    The the proportion of each isoforms
    """
    # transfer the parameters
    p_p, p_m = 1.0/len_isos

    # counting the reads of each type 
    Nc = np.sum(R_mat[:,0] * R_mat[:,1])
    Np = np.sum(R_mat[:,0] * (R_mat[:,1]==False))
    Nm = np.sum(R_mat[:,1] * (R_mat[:,0]==False))

    # The coefficients of a quadratic equation
    A = (-Np-Nm-Nc) * (p_p-p_m)
    B = Np*(p_p - 2*p_m) - Nm*p_m + Nc*(p_p - p_m)
    C = Np * p_m

    # # The analytical solution the quadratic equation on psi_freq
    Psi_freq = np.zeros(2)
    Psi_freq[0] = (-B - np.sqrt(B**2-4*A*C)) / (2.0*A)
    Psi_freq[1] = 1 - Psi_freq[0]
    RV = Psi_freq*len_isos[::-1] / np.sum(Psi_freq*len_isos[::-1])
    return RV

