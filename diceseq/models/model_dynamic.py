import numpy as np

def normal_pdf(x, mu, cov, log=True):
    """Only for one sample.
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
    RV = (-0.5*np.log(2*np.pi*cov_det) - 0.5*np.dot(np.dot(x, cov_inv), x))
    if log == False: RV = np.exp(RV)
    return RV


def logistic_normal_pdf(D, mu, cov, log=True):
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
        print "Only one array is calculated each time."
        return None
    D_logit = np.log(D[:-1]/D[-1]) - mu
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    if cov_det < 0:
        print "The det of covariance is negative, please check!"
        return None
    RV = (-0.5*np.log(2*np.pi*cov_det) - np.sum(np.log(D)) - 
           0.5*np.dot(np.dot(D_logit, cov_inv), D_logit))
    if log == False: RV = np.exp(RV)
    return RV


def Dirichlet_pdf(D, alpha, log=True):
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
        print "Only one array is calculated each time."
        return None
    part2 = 0
    for i in range(len(alpha)):
        part2 += np.math.lgamma(alpha[i])
    RV = np.math.lgamma(np.sum(alpha)) - part2 + np.sum((alpha-1)*np.log(D))
    if log == False: RV = np.exp(RV)
    return RV


def Geweke_diag():
    """The Geweke diagnostic of the convergency of MCMC samplling.
    """
    pass


def Psi_dynamic(R_mat, len_isos, prob_isos, cov_try_A, cov_rej_A, cov_try_psi,
                cov_rej_psi, A_prior, alpha_Dir, M=100):
    """
    Estimate the proportion of K isoforms at T time points with all reads 
    by MCMC samplling (MH algorithm)

    Parameters
    ----------
    R_mat : list of 2-D array_like, of length T
        A set of reads identities of belonging to K isoforms
    len_isos : list of 2-D array_like, of length T
        A set of effective length of K isoforms
    prob_isos : list of 2-D array_like, of length T
        A set of probablity for isoform specific reads
    cov_A : 2-D array_like, (K-1, K-1)
        the fixed hyper parameters; std variance of Normal distribution
        for transform matrix
    cov_psi : 2-D array_like, (K-1, K-1)
        the fixed hyper parameters; covariance of Logistic-Normal distribution
    A_prior : 2-D array_like, (K-1, K-1)
        the fixed hyper parameters; mean of Normal distribution for matrix A
    alpha_Dir : 1-D array_like, (K,)
        the fixed hyper parameter for Dirichlet distribution
    M : int
        the iteration times of in MCMC samplling, default=100

    Returns
    -------
    RV : 3-D array_like, (M, K, T)
        The the proposed M proportion of K isoforms of T time points
    """
    T = len(len_isos)
    K = len(len_isos[0])
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
    A_now   = np.identity(K-1)
    mu_now  = np.ones((K-1, T+1)) / K
    psi_now = np.zeros((K, T+1))
    fsi_now = np.zeros((K, T+1))
    for t in range(T+1):
        psi_now[:-1,t] = np.exp(mu_now[:,t]) / (1+np.sum(np.exp(mu_now[:,t])))
        psi_now[ -1,t] = 1-np.sum(psi_now[:-1, t])
        if t > 0: fsi_now[:, t]  = (len_isos[t-1]*psi_now[:,t] / 
                                    np.sum(len_isos[t-1]*psi_now[:,t]))

    # MCMC running
    A_try   = np.zeros((K-1, K-1))
    A_all   = np.zeros((M, K-1, K-1))
    mu_try  = np.zeros((K-1, T+1))
    psi_try = np.zeros((K, T+1))
    fsi_try = np.zeros((K, T+1))
    Psi_all = np.zeros((M, K, T))
    
    cnt = 0
    Pro_all = np.zeros(M)
    for m in range(M):
        Q_now, Q_try = 0, 0
        P_now, P_try = 0, 0
        # step 1: propose a value
        for i in range(K-1):
            for j in range(K-1):
                A_try[i,j] = np.random.normal(A_now[i,j], cov_try_A[i,j])
                Q_now += normal_pdf(A_now[i,j], A_try[i,j], cov_rej_A[i,j])
                Q_try += normal_pdf(A_try[i,j], A_now[i,j], cov_rej_A[i,j])
                P_now += normal_pdf(A_now[i,j], A_prior[i,j], cov_rej_A[i,j])
                P_try += normal_pdf(A_try[i,j], A_prior[i,j], cov_rej_A[i,j])

        for t in range(T+1):
            mu_try[:, t] = np.random.multivariate_normal(mu_now[:,t], cov_try_psi)
            psi_try[:-1,t] = np.exp(mu_try[:,t])/(1+np.exp(mu_try[:,t]).sum())
            psi_try[-1, t] = 1 - np.sum(psi_try[:-1, t])
            if t > 0: fsi_try[:, t]  = (len_isos[t-1] * psi_try[:,t] / 
                                        np.sum(len_isos[t-1] * psi_try[:,t]))

            Q_now += logistic_normal_pdf(psi_now[:, t], mu_try[:, t], cov_rej_psi)
            Q_try += logistic_normal_pdf(psi_try[:, t], mu_now[:, t], cov_rej_psi)
            if t == 0:
                P_now += Dirichlet_pdf(psi_now[:, t], alpha_Dir)
                P_try += Dirichlet_pdf(psi_try[:, t], alpha_Dir)
            else:
                P_now += normal_pdf(mu_now[:, t], np.dot(mu_now[:,t-1], A_now), cov_rej_psi)
                P_try += normal_pdf(mu_try[:, t], np.dot(mu_try[:,t-1], A_try), cov_rej_psi)
                P_now += np.log(np.dot(R_mat[t-1]*prob_isos[t-1], fsi_now[:, t])).sum()
                P_try += np.log(np.dot(R_mat[t-1]*prob_isos[t-1], fsi_try[:, t])).sum()

        # step 2: calculate the MH ratio at log scale
        alpha = min(np.exp(P_try+Q_now-P_now-Q_try), 1)

        if alpha is None:
            print "alpha is none!"
            Psi_all[m,:,:] = psi_now[:,1:]
            A_all[m,:,:] = A_now
            continue
        
        # step 3. accept or reject the proposal
        mu = np.random.rand(1)
        if mu < alpha:
            cnt += 1
            #print m, P_now, P_try
            A_now   = A_try + 0.0
            mu_now  = mu_try + 0.0
            psi_now = psi_try + 0.0
            fsi_now = fsi_try + 0.0
        Psi_all[m,:,:] = psi_now[:,1:]
        A_all[m,:,:] = A_now
        Pro_all[m] = P_now
    print cnt
    return Psi_all, A_all, Pro_all




def get_best_transition(X):
    K, T = X.shape
    A = np.zeros((K, K))
    B = np.zeros((K, K))
    C = np.zeros((K, K))
    for t in range(T-1):
        for i in range(K):
            for j in range(K):
                A[i,j] += X[i,t+1] * X[j,t]
            B[:,i]  += X[i,t]**2
    A = A / B
    for t in range(T-1):
        for i in range(K):
            C[i,i] += (X[i,t+1]**2 - 2*X[i,t+1]*np.dot(A, X[:,t])[i] + np.dot(A, X[:,t])[i]**2)
    C = C / (T-1)
    return A, C

def Psi_dynamic_collapse(R_mat, len_isos, prob_isos, cov_try_psi, alpha_Dir, M=100):
    T = len(len_isos)
    K = len(len_isos[0])
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
    # mu_now  = np.ones((K-1, T+1)) / K
    mu_now  = np.random.rand(K-1, T+1)
    mu_now  = mu_now / (mu_now.sum(axis=0) + np.random.rand(T+1))
    psi_now = np.zeros((K, T+1))
    fsi_now = np.zeros((K, T+1))
    for t in range(T+1):
        psi_now[:-1,t] = np.exp(mu_now[:,t]) / (1+np.sum(np.exp(mu_now[:,t])))
        psi_now[ -1,t] = 1-np.sum(psi_now[:-1, t])
        if t > 0: fsi_now[:, t]  = (len_isos[t-1]*psi_now[:,t] / 
                                    np.sum(len_isos[t-1]*psi_now[:,t]))
    A_now, C_now = get_best_transition(mu_now)

    # MCMC running
    A_try   = np.zeros((K-1, K-1))
    A_all   = np.zeros((M, K-1, K-1))
    C_try   = np.zeros((K-1, K-1))
    mu_try  = np.zeros((K-1, T+1))
    psi_try = np.zeros((K, T+1))
    fsi_try = np.zeros((K, T+1))
    Psi_all = np.zeros((M, K, T))
    
    cnt = 0
    Pro_all = np.zeros(M)
    for m in range(M):
        Q_now, Q_try = 0, 0
        P_now, P_try = 0, 0
        
        # step 1: propose a value
        for t in range(T+1):
            mu_try[:, t] = np.random.multivariate_normal(mu_now[:,t], cov_try_psi)
            psi_try[:-1,t] = np.exp(mu_try[:,t])/(1+np.exp(mu_try[:,t]).sum())
            psi_try[-1, t] = 1 - np.sum(psi_try[:-1, t])
            if t > 0: fsi_try[:, t]  = (len_isos[t-1] * psi_try[:,t] / 
                                        np.sum(len_isos[t-1] * psi_try[:,t]))
        A_try, C_try = get_best_transition(mu_try)
        
        for t in range(T+1):
            Q_now += logistic_normal_pdf(psi_now[:, t], mu_try[:, t], cov_try_psi)
            Q_try += logistic_normal_pdf(psi_try[:, t], mu_now[:, t], cov_try_psi)
            if t == 0:
                P_now += Dirichlet_pdf(psi_now[:, t], alpha_Dir)
                P_try += Dirichlet_pdf(psi_try[:, t], alpha_Dir)
            else:
                P_now += normal_pdf(mu_now[:, t], np.dot(mu_now[:,t-1], A_now), C_now)
                P_try += normal_pdf(mu_try[:, t], np.dot(mu_try[:,t-1], A_try), C_try)
                P_now += np.log(np.dot(R_mat[t-1]*prob_isos[t-1], fsi_now[:, t])).sum()
                P_try += np.log(np.dot(R_mat[t-1]*prob_isos[t-1], fsi_try[:, t])).sum()
        
        # step 2: calculate the MH ratio at log scale
        alpha = min(np.exp(P_try+Q_now-P_now-Q_try), 1)
        if alpha is None:
            print "alpha is none!"
            Psi_all[m,:,:] = psi_now[:,1:]
            A_all[m,:,:] = A_now
            continue
        
        # step 3. accept or reject the proposal
        mu = np.random.rand(1)
        if mu < alpha:
            cnt += 1
            #print m, P_now, P_try
            A_now   = A_try + 0.0
            C_now   = C_try + 0.0
            mu_now  = mu_try + 0.0
            psi_now = psi_try + 0.0
            fsi_now = fsi_try + 0.0
        Psi_all[m,:,:] = psi_now[:,1:]
        A_all[m,:,:] = A_now
        Pro_all[m] = P_now
    print cnt
    return Psi_all, A_all, Pro_all