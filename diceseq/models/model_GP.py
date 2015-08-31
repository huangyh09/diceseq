import numpy as np

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
        print "The shape of Gamma distribution cannot be 0, please check!"
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
        print "The det of covariance is negative, please check!"
        return None
    pdf = (-0.5*np.log(2*np.pi*cov_det) - 0.5*np.dot(np.dot(x, cov_inv), x))
    if log == False: pdf = np.exp(pdf)
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
            K[i,j] = theta[0] * np.exp(-0.5/theta[1]*(X[i]-X[j])**2)
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
    A = X[:first*N]
    B = X[last*N:]
    Z = abs(A.mean() - B.mean()) / np.sqrt(np.var(A) + np.var(B))
    return Z

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
        print "The number of isoforms in the reads and junction lengths are different!"
        return [None] * K
    if np.sum(R_mat) == 0:
        print "No junction reads!"
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

def Psi_GP_MH(R_mat, len_isos, prob_isos, X, cov_Y, sca_theta, M=100000, start=3000, gap=1000):
    """
    Estimate the proportion of K isoforms at T time points with all reads 
    by MCMC samplling (MH algorithm) combined with a GP prior.

    Parameters
    ----------
    R_mat : list of 2-D array_like, of length T
        A set of reads identities of belonging to K isoforms
    len_isos : list of 2-D array_like, of length T
        A set of effective length of K isoforms
    prob_isos : list of 2-D array_like, of length T
        A set of probablity for isoform specific reads
    X : 1-D array_like, (T, )
        An array of time points.
    cov_Y : 1-D array_like, (K, )
        An array of jump variance of each dimension of Y.
    sca_theta : 1-D array_like, (2, )
        An array of jump scale of each parameters of GP theta
    M : int
        the maximum iterations of in MCMC sampler, default=100000
    start : int
        the minmum iterations of in MCMC sampler, default=3000
    gap : int
        the gap iterations of in MCMC sampler, default=1000

    Returns
    -------
    RV : 3-D array_like, (M, K, T)
        The the proposed M proportion of K isoforms of T time points
    """
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

    # MCMC random initializations
    theta_now = np.zeros((C, 2))
    Y_now = np.zeros((C, T))
    for c in range(C):
        theta_now[c,0] = np.random.gamma(1.0/sca_theta[0], sca_theta[0])
        theta_now[c,1] = np.random.gamma(1.0/sca_theta[1], sca_theta[1])
        Y_now[c,:] = np.random.normal(0, cov_Y[c], T)
    psi_now = np.zeros((C, T))
    fsi_now = np.zeros((C, T))
    for t in range(T):
        psi_now[:,t] = np.exp(Y_now[:,t]) / np.sum(np.exp(Y_now[:,t]))
        fsi_now[:,t] = len_isos[t]*psi_now[:,t] /  np.sum(len_isos[t]*psi_now[:,t])
    
    P_now = 0
    cov_now = np.zeros((T, T, C))
    for c in range(C):
        cov_now[:,:,c] = GP_K(X, theta_now[c,:])
        P_now += normal_pdf(Y_now[c,:], np.zeros(T), cov_now[:,:,c]) #this may be a problem, try to set priors
    for t in range(T):
        P_now += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_now[:, t])).sum()
        
    # MCMC running
    theta_try = np.zeros((C, 2))
    theta_all = np.zeros((M, C, 2))
    Y_try   = np.zeros((C, T))
    Y_all   = np.zeros((M, C, T))
    psi_try = np.zeros((C, T))
    fsi_try = np.zeros((C, T))
    Psi_all = np.zeros((M, C, T))
    cov_try = np.zeros((T, T, C))
    
    cnt = 0
    Pro_all = np.zeros(M)
    for m in range(M):
        P_try, Q_now, Q_try = 0, 0, 0
        # step 1: propose a value
        for c in range(C):
            for j in range(2):
                theta_try[c,j] = np.random.gamma(theta_now[c,j]/sca_theta[j], sca_theta[j]) # make sure this is good
                while theta_try[c,j] == 0:
                    theta_try[c,j] = np.random.gamma(theta_now[c,j]/sca_theta[j], sca_theta[j]) # make sure this is good
                # the two important paramters jump in the same way. Is it good?
                Q_now += gamma_pdf(theta_now[c,j], theta_try[c,j]/sca_theta[j], sca_theta[j])
                Q_try += gamma_pdf(theta_try[c,j], theta_now[c,j]/sca_theta[j], sca_theta[j])
        
        for t in range(T):
            for c in range(C):
                Y_try[c,t] = np.random.normal(Y_now[c,t], cov_Y[c]) # the jump at different time point is the same (good?)
                Q_now += normal_pdf(Y_now[c,t], Y_try[c,t], cov_Y[c])
                Q_try += normal_pdf(Y_try[c,t], Y_now[c,t], cov_Y[c])
            psi_try[:,t] = np.exp(Y_try[:,t]) / np.sum(np.exp(Y_try[:,t]))
            fsi_try[:,t] = len_isos[t]*psi_try[:,t] /  np.sum(len_isos[t]*psi_try[:,t])

        for c in range(C):
            cov_try[:,:,c] = GP_K(X, theta_try[c,:])
            P_try += normal_pdf(Y_try[c,:], np.zeros(T), cov_try[:,:,c])
        for t in range(T):
            P_try += np.log(np.dot(R_mat[t]*prob_isos[t], fsi_try[:,t])).sum()
        
        # step 2: calculate the MH ratio and decide to accapt or reject
        alpha = min(np.exp(P_try+Q_now-P_now-Q_try), 1)
        if alpha is None: print "alpha is none!"        
        elif np.random.rand(1) < alpha:
            cnt += 1
            Y_now = Y_try + 0.0
            P_now = P_try + 0.0
            cov_now = cov_try + 0.0
            theta_now = theta_try + 0.0
            psi_now = psi_try + 0.0
            fsi_now = fsi_try + 0.0

        Pro_all[m] = P_now
        Y_all[m,:,:] = Y_now
        Psi_all[m,:,:] = psi_now
        theta_all[m,:,:] = theta_now
        
        #step 3. convergence diagnostics
        if m >= start and m % gap == 0:
            conv = 1
            for c in range(C):
                for t in range(T):
                    if Geweke_Z(Psi_all[:m,c,t]) > 2: 
                        conv = 0
                        break
                if conv == 0: break
            if conv == 1:
                Psi_all, Pro_all = Psi_all[:m,:,:], Pro_all[:m,]
                break   

        # if m >= start and m % gap == 0:
        #     conv = 1
        #     for c in range(C):
        #         if Geweke_Z(theta_all[:m,c,0]) > 2: conv = 0
        #         if Geweke_Z(theta_all[:m,c,1]) > 2: conv = 0
        #         if conv == 0: break
        #         for t in range(T):
        #             if Geweke_Z(Y_all[:m,c,t]) > 2: 
        #                 conv = 0
        #                 break
        #     if conv == 1:
        #         Psi_all, Pro_all = Psi_all[:m,:,:], Pro_all[:m,]
        #         break    
            
    print "The numbers of iterations and acceptances are %d and %d." %(m, cnt)
    return Psi_all

