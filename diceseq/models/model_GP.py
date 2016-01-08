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
    x : float or 1-D array_like, (K, )
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
    """
    x = np.array(x) - np.array(mu)
    if len(np.array(cov).shape) < 2:
        cov = np.array(cov).reshape(-1,1)
    if len(np.array(x).shape) < 1:
        x = np.array(x).reshape(-1)
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    if cov_det < 0:
        print("The det of covariance is negative, please check!")
        return None
    pdf = (-0.5*np.log(2*np.pi*cov_det) - 0.5*np.dot(np.dot(x, cov_inv), x))
    if log == False: pdf = np.exp(pdf)
    return pdf

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
              randomS=None, theta2_std=1.0, theta2_low=0.01, theta2_up=100):
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
    var : float
        The average variation of all y.
    theta : 1-D array_like or list, (2, )
        The fixed hyper-parameter theta1, and initial theta2, default=[3.0,5.0]
    theta2_std : float
        The jump std of hyper-parameter theta2 for each dimension, default=1.0
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
     :
    Y_all : 3-D array_like, (m, C, T)
        The the proposed m latent y for C isoforms of T time points
    theta_all : 3-D array_like, (m, C-1, 2)
        The the proposed m hyper-parameters for C-1 isoforms
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
    theta_now = np.zeros((C-1, 2))
    theta_now[:,0] = theta1
    theta_now[:,1] = theta2
    if var is None: var = 0.05 * np.ones(C-1)
    if theta2 is not None: theta_now[:,1] = theta2
    else: theta_now[:,1] = 0.75
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
    theta_try = np.zeros((C-1, 2))
    theta_all = np.zeros((M, C-1, 2))
    Y_try   = np.zeros((C, T))
    Y_all   = np.zeros((M, C, T))
    psi_try = np.zeros((C, T))
    fsi_try = np.zeros((C, T))
    Psi_all = np.zeros((M, C, T))
    cov_try = np.zeros((T, T, C-1))
    
    cnt = 0
    Pro_all = np.zeros(M)
    Lik_all = np.zeros(M)
    for m in range(M):
        P_try, L_try, Q_now, Q_try = 0, 0, 0, 0
        
        # step 1: propose a value
        for c in range(C-1):
            theta_try[c,0] = theta1
            if theta2 is not None: theta_try[c,1] = theta2
            else:
                theta_try[c,1] = np.random.normal(theta_now[c,1], theta2_std)
                while theta_try[c,1]<theta2_low or theta_try[c,1]>theta2_up:
                    theta_try[c,1] = np.random.normal(theta_now[c,1],theta2_std)
            cov_try[:,:,c] = GP_K(X, theta_try[c,:])
            Y_try[c,:] = np.random.multivariate_normal(Y_now[c,:], 
                cov_try[:,:,c]/theta1*var[c]*5/T/C)
                        
            P_try += normal_pdf(Y_try[c,:], Ymean[c,:], cov_try[:,:,c])
            Q_now += normal_pdf(Y_now[c,:], Y_try[c,:], 
                cov_try[:,:,c]/theta1*var[c]*5/T/C)
            Q_try += normal_pdf(Y_try[c,:], Y_now[c,:], 
                cov_now[:,:,c]/theta1*var[c]*5/T/C)
            Q_now += trun_normal_pdf(theta_now[c,1], theta_try[c,1], 
                theta2_std, theta2_low, theta2_up)
            Q_try += trun_normal_pdf(theta_try[c,1], theta_now[c,1], 
                theta2_std, theta2_low, theta2_up)
            
        for t in range(T):
            psi_try[:,t] = np.exp(Y_try[:,t]) / np.sum(np.exp(Y_try[:,t]))
            fsi_try[:,t] = (len_isos[t]*psi_try[:,t] / 
                np.sum(len_isos[t]*psi_try[:,t]))
            P_try += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_try[:,t])).sum()
            L_try += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_try[:,t])).sum()
        
        # step 2: calculate the MH ratio and accept or reject the proposal
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
        theta_all[m,:] = theta_now
        
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
                Z = Geweke_Z(theta_all[:m,c,1])
                if Z is None or Z > 2: conv == 0
                if conv == 0: break
            if conv == 1:
                Pro_all = Pro_all[:m,]
                Lik_all = Lik_all[:m,]
                Y_all   = Y_all[:m,:,:]
                Psi_all = Psi_all[:m,:,:]
                theta_all = theta_all[:m,:]
                break
    return Psi_all, Y_all, theta_all, Pro_all, Lik_all, cnt, m

def EM_filter(R_mat, len_isos, prob_isos, min_num=2):
    T = len(len_isos)
    C = len(len_isos[0])
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

    KP_FLAG = np.ones(C, 'bool')
    count_tmp = np.zeros(T)
    psi_tmp = np.zeros((T, C))
    for t in range(T):
        if np.sum(R_mat[t]) == 0: continue
        count_tmp = np.sum(R_mat[t][:,:].sum(axis=1) > 0)
        psi_tmp[t,:], _Lik_all = Psi_EM(R_mat[t], len_isos[t], prob_isos[t])

    idx = np.argsort(psi_tmp.mean(axis=0))
    for i in range(len(idx)-min_num):
        if np.max(psi_tmp[:,idx[i]]) < 0.001 and np.mean(count_tmp) > 100:
            KP_FLAG[idx[i]] = False
    return KP_FLAG, psi_tmp


def Psi2Y(Psi):
    Psi = Psi / np.sum(Psi)
    Y = np.zeros(len(Psi))
    Y[:-1] = np.log(Psi[:-1] / Psi[-1])
    return Y

def EM_bootstrap(R_mat, len_isos, prob_isos, bootN=10):
    T = len(len_isos)
    C = len(len_isos[0])
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

    Y_boot = np.zeros((bootN, T, C))
    psi_boot = np.ones((bootN, T, C))/(C+0.0)
    for i in range(bootN):
        for t in range(T):
            _N = R_mat[t].shape[0]
            if _N == 0: continue
            idx = np.random.randint(_N, size=_N)
            if np.sum(R_mat[t][idx,:]) == 0: continue
            psi_boot[i,t,:], _Lik_all = Psi_EM(R_mat[t][idx,:], len_isos[t], prob_isos[t][idx,:])
            Y_boot[i,t:] = Psi2Y(psi_boot[i,t,:] + 0.000001)
    var_boot = Y_boot.var(axis=0).mean(axis=0)
    # var_boot = psi_boot.var(axis=0).mean(axis=0) * 36
    
    return Y_boot.mean(axis=0), var_boot
    # return psi_boot.mean(axis=0), var_boot

def Psi_EM(R_mat, len_isos, prob_isos, M=100):
    # check input data
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

    # random initialization
    Psi = np.random.dirichlet(np.ones(R_mat.shape[1]))
    Fsi = Psi * len_isos / np.sum(Psi * len_isos)
    Lik = np.log(np.dot(R_mat*prob_isos, Fsi)).sum()
    
    cnt = 0
    Lik_all = np.zeros(M)
    for i in range(M):
        cnt += 1
        # E step
        Num = (R_mat*prob_isos*Fsi)
        Num = Num / Num.sum(axis=1).reshape(-1,1)
        Num = Num.sum(axis=0)
        
        # M step
        Pro = Num / len_isos
        Psi = Pro / sum(Pro)
        Fsi = Psi * len_isos / sum(Psi * len_isos)
        Lik = np.log(np.dot(R_mat*prob_isos, Fsi)).sum()
        Lik_all[i] = Lik

        # check convergence
        if i > 0 and (Lik - Lik_all[i-1]) < np.log(1+10**(-5)):
            break
    return Psi, Lik_all