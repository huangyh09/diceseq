import numpy as np

def normal_pdf(x, mu, cov, log=True):
    """
    Calculate the probability density of Gaussian (Normal) distribution.
    Parameters
    ----------
    x : float, 1-D array_like (K, ), or 2-D array_like (K, N)
        The variable for calculating the probability density.
    mu : float or 1-D array_like, (K, )
        The mean of the Gaussian distribution.
    cov : float or 2-D array_like, (K, K)
        The variance or the covariance matrix of the Gaussian distribution.
    log : bool
        If true, the return value is at log scale.
    Returns
    -------
    pdf : numpy float
        The probability density of x. 
        if N==1, return a float
        elif N>1, return an array
    """
    if len(np.array(mu).shape) == 0:
        x = np.array(x).reshape(-1,1)
    elif len(np.array(x).shape) <= 1:
        x = np.array(x).reshape(1, -1)
    x = x - np.array(mu)
    N, K = x.shape
    if len(np.array(cov).shape) < 2:
        cov = np.array(cov).reshape(-1,1)
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    if cov_det <= 0:
        print("Warning: the det of covariance is not positive!")
        return None
    pdf_all = np.zeros(N)
    pdf_part1 = -(K*np.log(2*np.pi) + np.log(cov_det)) / 2.0
    for i in range(N):
        pdf_all[i] = pdf_part1 - np.dot(np.dot(x[i,:], cov_inv), x[i,:]) / 2.0
    if log == False: pdf_all = np.exp(pdf_all)
    if N == 1: pdf_all = pdf_all[0]
    return pdf_all
    

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




def Iso_read_check(R_mat, len_isos, prob_isos):
    """
    Check the input data for isoform quantification.

    Parameters
    ----------
    R_mat : 2-D array_like, (N, K)
        N reads identities of belonging to K isoforms.
    prob_isos : 2-D array_like, (N, K)
        N reads probablity of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.

    Returns
    -------
    prob_isos : 2-D array_like, (N, K)
        N reads probablity of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.
    """
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

    return R_mat, prob_isos, len_isos


### adaptive MCMC
def Psi_MH(R_mat, len_isos, prob_isos, cov=None, alpha_Dir=None, M=10000,
    initial=500, gap=100):
    """
    Calculate the proportion of each isoforms with all reads via MCMC samplling
    Parameters
    ----------
    R_mat : 2-D array_like, (N, K)
        N reads identities of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.
    cov : 2-D array_like, (K-1, K-1)
        The fixed covariance of multivariate normal distribution of Y.
    alpha_Dir : 1-D array_like, (K,)
        The hyper-parameter of Dirichlet distribution, the prior.
    M : int
        The maximum iteration times of MH sampler, default=10000.
    initial : int
        The minmum iteration times of MH sampler, default=1000.
    gap : 
        The gap iteration times of MH sampler, default=100.

    Returns
    -------
    Psi_all : 2-D array_like, (M, K)
        M sampled Psi for abundance of each isoforms
    Y_all : 2-D array_like, (M, K)
        M sampled Psi for abundance of each isoforms
    cnt : float
        The accepted times of iterations
    """
    # 1. check input data
    R_mat, prob_isos, len_isos = Iso_read_check(R_mat, len_isos, prob_isos)
    prob_isos = R_mat * prob_isos

    N, K = prob_isos.shape
    if cov is None:
        cov = 1.5 * np.identity(K-1) / (K-1)
    if alpha_Dir is None:
        alpha_Dir = np.ones(K)

    # MCMC random initializations
    Y_now = np.zeros(K)
    psi_now = np.exp(Y_now) / np.sum(np.exp(Y_now))
    fsi_now = len_isos*psi_now / np.sum(len_isos*psi_now)
    P_now = Dirichlet_pdf(psi_now, alpha_Dir, log=True)
    P_now += np.sum(np.log(np.dot(prob_isos, fsi_now)))
    
    # MCMC running
    cnt = 0
    Y_avg = 0.0
    Y_var = 0.0
    Y_try = np.zeros(K)
    Y_all = np.zeros((M, K))
    psi_try = np.zeros(K)
    Psi_all = np.zeros((M, K))
    for m in range(M):
        P_try, Q_now, Q_try = 0, 0, 0
        # step 1: propose a value
        Y_try[:K-1] = np.random.multivariate_normal(Y_now[:K-1], cov)
        Y_try[Y_try < -700] = -700
        Y_try[Y_try > 700 ] = 700

        Q_now += normal_pdf(Y_now[:K-1], Y_try[:K-1], cov)
        Q_try += normal_pdf(Y_try[:K-1], Y_now[:K-1], cov)

        psi_try = np.exp(Y_try) / np.sum(np.exp(Y_try))
        fsi_try = len_isos*psi_try / np.sum(len_isos*psi_try)
        P_try += Dirichlet_pdf(psi_try, alpha_Dir, log=True)
        P_try += np.sum(np.log(np.dot(prob_isos, fsi_try)))
        
        # step 2: calculate the MH ratio, accept or reject the proposal
        alpha = np.exp(min(P_try+Q_now-P_now-Q_try, 0))
        # print alpha, np.exp(P_try), np.exp(Q_now), np.exp(P_now), np.exp(Q_try)
        if alpha is None: 
            print("alpha is none!")
        elif np.random.rand(1) < alpha:
            cnt += 1
            Y_now  = Y_try  + 0
            psi_now = psi_try + 0
            fsi_now = fsi_try + 0
            P_now   = P_try   + 0
        Y_all[m,:] = Y_now
        Psi_all[m,:] = psi_now

        #step 3. convergence diagnostics
        if m >= initial and m % gap == 0:
            conv = True
            for k in range(K):
                Z = Geweke_Z(Psi_all[:m, k])
                if Z is None or Z > 2: 
                    conv = False
                    break
            if conv:
                Y_all = Y_all[:m,:]
                Psi_all = Psi_all[:m,:]
                break

        _avg1 = (Y_avg*m + Y_now) / (m+1.0)
        Y_var = (Y_var*m + (Y_now-_avg1)*(Y_now-Y_avg)) / (m+1.0)
        Y_avg = _avg1 + 0.0

        # adaptive MCMC
        if m >= 11: #and m % gap == 0:
            cov = np.cov(Y_all[:m,:-1].T) + np.diag(np.ones(K-1))*0.001
            cov = cov * 5.0 / (K-1) / (1 + prob_isos.shape[0]/500.0)

            # cov1 = np.diag(Y_var[:-1]+0.001)
            # cov2 = np.cov(Y_all[:m,:-1].T)
            # beta = 0.05 #min((K-1)/20.0, 1.0)
            # cov = 5.5 / (K-1) * (beta*cov1 + (1-beta)*cov2) / (1 + prob_isos.shape[0]/500.0)
            
            # print("accept ratio: %.3f" %(cnt/(m+1.0)))

    print("Total accept ratio: %.3f for %d isoforms with %d reads." 
        %(cnt/(m+1.0), prob_isos.shape[1], prob_isos.shape[0]))
    return Psi_all, Y_all, cnt





##### EM
def Psi_EM(R_mat, len_isos, prob_isos, M=100):
    # check input data
    R_mat, prob_isos, len_isos = Iso_read_check(R_mat, len_isos, prob_isos)
    prob_isos = R_mat * prob_isos

    # random initialization
    Psi = np.random.dirichlet(np.ones(prob_isos.shape[1]))
    Fsi = Psi * len_isos / np.sum(Psi * len_isos)
    Lik = np.log(np.dot(prob_isos, Fsi)).sum()
    
    cnt = 0
    Lik_all = np.zeros(M)
    for i in range(M):
        cnt += 1
        # E step
        Num = (prob_isos*Fsi)
        Num = Num / Num.sum(axis=1).reshape(-1,1)
        Num = Num.sum(axis=0)
        
        # M step
        Pro = Num / len_isos
        Psi = Pro / sum(Pro)
        Fsi = Psi * len_isos / sum(Psi * len_isos)
        Lik = np.log(np.dot(prob_isos, Fsi)).sum()
        Lik_all[i] = Lik

        # check convergence
        if i > 0 and (Lik - Lik_all[i-1]) < np.log(1+10**(-5)):
            break
    return Psi, Lik_all