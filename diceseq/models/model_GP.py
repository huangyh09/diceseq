"""
Inference model with Gaussian process for sampling isoform proportions.
"""

import numpy as np

def erf(x):
    """error function approximation, with maximum error 1.5e-7.
    """
    # save the sign of x
    if x >= 0: sign = 1
    else: sign = -1
    x = abs(x)
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)


def gamma_pdf(x, k, theta, log=True):
    """
    Calculate the probability density of Gamma distribution.

    Parameters
    ----------
    x : float
        The variable for calculating the probability density.
    k : float
        The shape of the Gamma distribution.
    theta : float
        The scale of the Gamma distribution.
    log : bool
        If true, the return value is at log scale.

    Returns
    -------
    pdf : numpy float
        The probability density of x.
    """
    if k == 0: 
        print("The shape of Gamma distribution cannot be 0, please check!")
        return None
    pdf = -np.math.lgamma(k) - k*np.log(theta)
    pdf += (k-1)*np.log(x) - x/theta
    if log == False: pdf = np.exp(pdf)
    return pdf

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

def trun_normal_pdf(x, mu, sigma, a, b, log=True):
    """
    Calculate the probability density of Truncated Normal distribution.

    Parameters
    ----------
    x : float
        The variable for calculating the probability density.
    mu : float
        The mean of the Gaussian distribution.
    sigma : float
        The standard variance of the Gaussian distribution.
    a : float
        The lower bounder of the Truncated Normal distribution
    b : float
        The upper bounder of the Truncated Normal distribution
    log : bool
        If true, the return value is at log scale.

    Returns
    -------
    pdf : float
        The probability density of x.
    """
    x = x - mu
    a = a - mu
    b = b - mu
    
    pdf = np.exp(-0.5 * (x/sigma)**2) / (sigma * np.sqrt(2 * np.pi))
    cdf_a = (1 + erf(a / sigma / np.sqrt(2))) / 2.0
    cdf_b = (1 + erf(b / sigma / np.sqrt(2))) / 2.0
    
    pdf = pdf / abs(cdf_b - cdf_a)
    if log == True: pdf = np.log(pdf)
    return pdf

def GP_K(X, theta):
    """
    Covariance of Gaussian process generator.
    It is based on a common squared-exponential kernel, with two parameters. 

    Parameters
    ----------
    X : 1-D array_like, (N, )
        The x-axis of the Gaussian process, e.g., time points.
    theta : 1-D array_like, (2,)
        The array of the two parameters of the squared-exponential kernel.

    Returns
    -------
    K : 2-D array_like, (N, N)
        The covariance matrix of the N points at x-axis.
    """
    N = len(X)
    K = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            K[i,j] = theta[0] * np.exp(-0.5 * (X[i]-X[j])**2 / theta[1])
    return K

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
    A = X[:int(first*N)]
    B = X[int(last*N):]
    if np.sqrt(np.var(A) + np.var(B)) == 0:
        Z = None
    else:
        Z = abs(A.mean() - B.mean()) / np.sqrt(np.var(A) + np.var(B))
    return Z


def Psi_GP_MH(R_mat, len_isos, prob_isos, X=None, Ymean=None, var=None, 
              theta1=3.0, theta2=None, M=20000, initial=1000, gap=500, 
              randomS=None, theta2_std=1.0, theta2_low=0.00001, theta2_up=100):
    """
    Estimate the proportion of C isoforms at T time points with all reads 
    by MCMC samplling (MH algorithm) combined with a GP prior.

    Parameters
    ----------
    R_mat : list of 2-D array_like, of length T
        A set of reads identities of belonging to C isoforms
    len_isos : list of 2-D array_like, of length T
        A set of effective length of C isoforms
    prob_isos : list of 2-D array_like, of length T
        A set of probablity for isoform specific reads
    X : 1-D array_like, (T, )
        An array of time points.
    Ymean : 2-D array_like, (C, T)
        The means for Y.
    var : 1-D array_like, (C-1, )
        An array of variation of each y.
    theta1 : float
        The fixed hyper-parameter theta1
    theta2 : float
        The fixed hyper-parameter theta2. If it is None, then sample it.
    theta2_std : float
        The jump std of hyper-parameter theta2 for each dimension, default=1.0
    theta2_low : float
        The lower bound of Truncated Normal distribution for sampling theta2.
    theta2_up : float
        The upper bound of Truncated Normal distribution for sampling theta2.
    randomS : float
        The fixed seeds for random number generator. None means ignoring it.
    M : int
        the maximum iterations of in MCMC sampler, default=100000
    initial : int
        the minmum iterations of in MCMC sampler, default=3000
    gap : int
        the gap iterations of in MCMC sampler, default=1000

    Returns
    -------
    Psi_all : 3-D array_like, (m, C, T)
        The the proposed m proportion of C isoforms of T time points
    Y_all : 3-D array_like, (m, C, T)
        The the proposed m latent y for C isoforms of T time points
    theta2_all : 2-D array_like, (m, C-1)
        The the proposed m hyper-parameter theta2 for C-1 isoforms
    Pro_all : 1-D array_like, (m,)
        The the probability for accepted proposals
    Lik_all : 1-D array_like, (m,)
        The the probability for accepted proposals
    cnt : int
        The number of acceptances
    m : int
        The number of iterations
    """
    T = len(len_isos)
    C = len(len_isos[0])
    if X is None: X = np.arange(T)
    if Ymean is None: Ymean = np.zeros((C,T))
    if randomS is not None: np.random.seed(randomS)
    for t in range(T):
        idx = (len_isos[t] != len_isos[t])
        len_isos[t][idx] = 0.0
        prob_isos[t][:,idx] = 0.0
        R_mat[t][:,idx] = False

        idx = np.where(R_mat[t] != R_mat[t])
        R_mat[t][idx] = False

        idx = np.where(prob_isos[t] != prob_isos[t])
        prob_isos[t][idx] = 0.0

        idx = (R_mat[t].sum(axis=1) > 0) * (prob_isos[t].sum(axis=1) > 0)
        R_mat[t] = R_mat[t][idx,:]
        prob_isos[t] = prob_isos[t][idx,:]

    # step 0: MCMC fixed initializations
    if var is None: 
        var = 0.05 * np.ones(C-1)
    theta_now = np.zeros((C-1, 2))
    theta_now[:,0] = theta1
    if theta2 is not None: 
        theta_now[:,1] = theta2
    else: 
        theta_now[:,1] = 0.1 * (np.max(T)-np.min(T)+0.001)**2 #0.75
    Y_now = Ymean + 0.0
    Ymean = np.zeros((C,T))

    psi_now = np.zeros((C, T))
    fsi_now = np.zeros((C, T))
    for t in range(T):
        psi_now[:,t] = np.exp(Y_now[:,t]) / np.sum(np.exp(Y_now[:,t]))
        fsi_now[:,t] = len_isos[t]*psi_now[:,t]/np.sum(len_isos[t]*psi_now[:,t])
    
    P_now, L_now = 0, 0
    cov_now = np.zeros((T, T, C-1))
    for c in range(C-1):
        cov_now[:,:,c] = GP_K(X, theta_now[c,:])
        P_now += normal_pdf(Y_now[c,:], Ymean[c,:], cov_now[:,:,c])
        
    for t in range(T):
        P_now += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_now[:, t])).sum()
        L_now += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_now[:, t])).sum()
        
    # MCMC running
    Y_try = np.zeros((C, T))
    Y_all = np.zeros((M, C, T))
    psi_try = np.zeros((C, T))
    fsi_try = np.zeros((C, T))
    Psi_all = np.zeros((M, C, T))
    cov_try = np.zeros((T, T, C-1))
    theta_try = np.zeros((C-1, 2))
    theta2_all = np.zeros((M, C-1))
    theta_try[:, 0] = theta1
    if theta2 is not None:
        theta_try[:,1] = theta2
        cov_try[:,:,:] = GP_K(X, theta_try[0,:]).reshape(T,T,1)
    
    cnt = 0
    Pro_all = np.zeros(M)
    Lik_all = np.zeros(M)
    for m in range(M):
        P_try, L_try, Q_now, Q_try = 0, 0, 0, 0
        
        # step 1: propose a value
        for c in range(C-1):
            # sample single theta2 for all isoforms
            if theta2 is None and c==0:
                theta_try[:,1] = np.random.normal(theta_now[c,1], theta2_std)
                while theta_try[c,1]<theta2_low or theta_try[c,1]>theta2_up:
                    theta_try[:,1] = np.random.normal(theta_now[c,1],theta2_std)
                cov_try[:,:,c] = GP_K(X, theta_try[c,:])

                Q_now += trun_normal_pdf(theta_now[c,1], theta_try[c,1], 
                    theta2_std, theta2_low, theta2_up)
                Q_try += trun_normal_pdf(theta_try[c,1], theta_now[c,1], 
                    theta2_std, theta2_low, theta2_up)

            cov_jmp = cov_try[:,:,c] * var[c] * 5 / (T * C * theta1)
            Y_try[c,:] = np.random.multivariate_normal(Y_now[c,:], cov_jmp)
            Q_now += normal_pdf(Y_now[c,:], Y_try[c,:], cov_jmp)
            Q_try += normal_pdf(Y_try[c,:], Y_now[c,:], cov_jmp)
            P_try += normal_pdf(Y_try[c,:], Ymean[c,:], cov_try[:,:,c])
            
        for t in range(T):
            psi_try[:,t] = np.exp(Y_try[:,t]) / np.sum(np.exp(Y_try[:,t]))
            fsi_try[:,t] = (len_isos[t]*psi_try[:,t] / 
                np.sum(len_isos[t]*psi_try[:,t]))
            _lik_list = np.dot(R_mat[t]*prob_isos[t], fsi_try[:,t])
            # if min(_lik_list) <= 0:
            #     P_try, L_try = -np.inf, -np.inf
            # else:
            #     P_try += np.log(_lik_list).sum()
            #     L_try += np.log(_lik_list).sum()
            P_try += np.log(_lik_list).sum()
            L_try += np.log(_lik_list).sum()
        
        # step 2: calculate the MH ratio; accept or reject the proposal
        alpha = np.exp(min(P_try+Q_now-P_now-Q_try, 0))
        if alpha is None:
            print("alpha is none!")
        elif np.random.rand(1) < alpha:
            #print alpha
            cnt += 1
            P_now = P_try + 0.0
            L_now = L_try + 0.0
            Y_now = Y_try + 0.0
            cov_now = cov_try + 0.0
            psi_now = psi_try + 0.0
            fsi_now = fsi_try + 0.0
            theta_now = theta_try + 0.0            

        Pro_all[m]     = P_now
        Lik_all[m]     = L_now
        Y_all[m,:,:]   = Y_now
        Psi_all[m,:,:] = psi_now
        theta2_all[m,:] = theta_now[:,1]
        
        #step 3. convergence diagnostics
        if m >= initial and m % gap == 0:
            conv = 1
            for c in range(C-1):
                for t in range(T):
                    # Z = Geweke_Z(Y_all[:m,c,t])
                    Z = Geweke_Z(Psi_all[:m,c,t])
                    if Z is None or Z > 2: 
                        conv = 0
                        break
                    #print("psi converged!")
                if theta2 is None:
                    Z = Geweke_Z(theta2_all[:m, c])
                    if Z is None or Z > 2: 
                        conv = 0
                        break
                if conv == 0: break
            if conv == 1:
                Pro_all = Pro_all[:m,]
                Lik_all = Lik_all[:m,]
                Y_all   = Y_all[:m,:,:]
                Psi_all = Psi_all[:m,:,:]
                theta2_all = theta2_all[:m,:]
                break

    # if m >= initial and conv == 0:
    #     print("Warning: Not converged. Need a longer MCMC chain.")
    return Psi_all, Y_all, theta2_all, Pro_all, Lik_all, cnt, m
