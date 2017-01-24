import numpy as np
from scipy.stats import gaussian_kde, mvn
from .model_GP import normal_pdf, GP_K, Geweke_Z

def logistic(y, C=2):
    """logistic function: transfer y to psi
    """
    return np.exp(y) / (C - 1.0 + np.exp(y))

def logit(psi, C=2):
    """logistic function: transfer psi to y
    """
    return np.log((C-1.0) * psi / (1-psi))
    
def miso_BF(sample1, sample2, max_diff=0.0, bootstrap_size=None, 
    max_bf=1e12, min_unique=5, log=True):
    """
    Bayes factor calculator in MISO.
    
    Assumptions
    -----------
    prior : 
        uniform distribution for psi in both conditions.
        analytically calculating prior distribution of 
        difference, with shape like a triangle.
    posterior :
        smooth the posterior of difference by a Gaussian
        density estimate.
    
    Parameters
    ----------
    sample1 : array_like or list
        samples of psi in condition 1.
    sample2 : array_like or list
        samples of psi in condition 2.
    max_diff : float or int
        null hypothesis: abs(condtion1 - conditon2) <= max_diff.
        max_diff ranges from 0 to 1.
    bootstrap_size : int or None
        if None, don't do bootstrapping.
        if int, resampling with size of bootstrap_size.
    max_bf : float or int
        upper bound of Bayes factor.
    min_unique : int
        minimum number of unique samples in both conditions
        if any unique samples is smaller than min_unique , 
        supporting Null hypothesis.
        
    Returns
    -------
    bayes_factor : float
        The Bayes factor supporting alternative hypothesis
        abs(condtion1 - conditon2) > max_diff
    """
    if min(len(np.unique(sample1)), len(np.unique(sample1))) < min_unique:
        print("Improper sampling of psi! Supporting Null hypothesis.")
        return 0 #None
    
    if bootstrap_size is None:
        idx1 = np.arange(len(sample1))
        idx2 = np.arange(len(sample2))
    else:
        idx1 = np.random.randint(len(sample1), size=bootstrap_size)
        idx2 = np.random.randint(len(sample2), size=bootstrap_size)
    samp_diff = sample1[idx1] - sample2[idx2]
    
    prior_density = lambda x: 1 + x if x <= 0 else 1 - x
    posterior_density = gaussian_kde(samp_diff, bw_method='scott')
    
    #TODO: the decimals make a huge difference!
    #Maybe the inverval can be less sensitive!
    if max_diff == 0:
        diff_prior = np.log(prior_density(max_diff))
        diff_posterior = np.log(posterior_density.evaluate(max_diff)[0])
    else:
        diff_prior = np.log((prior_density(0.0) + prior_density(max_diff)) * (max_diff))
        diff_posterior = np.log(posterior_density.integrate_box_1d(-max_diff, max_diff))

    if diff_posterior == -np.inf:
        bayes_factor = max_bf
    elif diff_posterior == np.inf:
        bayes_factor = -max_bf
    else:
        bayes_factor = diff_prior - diff_posterior
    if log is False:
        bayes_factor = np.exp(bayes_factor)

    return bayes_factor


def dice_BF(sample1, sample2, bootstrap_size=None, max_bf=1e12, 
    min_unique=5, log=True):
    """
    Bayes factor calculator in DICE-diff.
    
    Assumptions
    -----------
    prior : 
        uniform distribution for psi in both conditions.
        analytically calculating prior distribution of 
        difference, with shape like a triangle.
    posterior :
        smooth the posterior of difference by a Gaussian
        density estimate.
    
    Parameters
    ----------
    sample1 : array_like or list
        samples of psi in condition 1.
    sample2 : array_like or list
        samples of psi in condition 2.
    max_diff : float or int
        null hypothesis: abs(condtion1 - conditon2) <= max_diff.
        max_diff ranges from 0 to 1.
    bootstrap_size : int or None
        if None, don't do bootstrapping.
        if int, resampling with size of bootstrap_size.
    max_bf : float or int
        upper bound of Bayes factor.
    min_unique : int
        minimum number of unique samples in both conditions
        if any unique samples is smaller than min_unique , 
        supporting Null hypothesis.
        
    Returns
    -------
    bayes_factor : float
        The Bayes factor supporting alternative hypothesis
        abs(condtion1 - conditon2) > max_diff
    """
    if len(np.array(sample1).shape) == 1:
        sample1 = sample1.reshape(-1,1)
        sample2 = sample2.reshape(-1,1)
    k = sample1.shape[1]
    # if min(len(np.unique(sample1)), len(np.unique(sample1))) < min_unique:
    #     print("Improper sampling of psi! Supporting Null hypothesis.")
    #     return 0 #None
    
    if bootstrap_size is None:
        idx1 = np.arange(len(sample1))
        idx2 = np.arange(len(sample2))
    else:
        idx1 = np.random.randint(len(sample1), size=bootstrap_size)
        idx2 = np.random.randint(len(sample2), size=bootstrap_size)
    samp_diff = sample1[idx1,:] - sample2[idx2,:]
    posterior_density = gaussian_kde(samp_diff.T, bw_method='scott')

    #TODO: the decimals make a huge difference!
    #Maybe the inverval can be less sensitive!
    mu = np.zeros(k)
    cov = GP_K(np.arange(k), [3.0, 5.0])
    diff_prior_log = normal_pdf(mu, mu, 2*cov, log=True)
    diff_posterior_log = np.log(posterior_density.evaluate(mu)[0])
    
    if diff_posterior_log == -np.inf:
        bayes_factor = max_bf
    elif diff_posterior_log == np.inf:
        bayes_factor = -max_bf
    else:
        bayes_factor = diff_prior_log - diff_posterior_log
    if log is False:
        bayes_factor = np.exp(bayes_factor)

    return bayes_factor


def dicediff_BF(samples1, samples2, bio_cov=None, bootstrap_size=None, 
                max_diff=0.0, max_bf=1e12, min_unique=5, log=True,
                post_mode="sample", prior_mode="GP", prior_param=[3, 5],
                time_corr=True, is_logit=True):
    """
    Bayes factor calculator in DICE-diff.
    
    Assumptions
    -----------
    prior : uniform or GP
        uniform: uniform distribution for psi in both conditions.
        analytically calculating prior distribution of 
        difference, with shape like a triangle.
    posterior : mean, kde, or sample
        smooth the posterior of difference by a Gaussian
        density estimate.
    
    Parameters
    ----------
    sample1 : array_like or list
        samples of psi in condition 1.
    sample2 : array_like or list
        samples of psi in condition 2.
    max_diff : float or int
        null hypothesis: abs(condtion1 - conditon2) <= max_diff.
        max_diff ranges from 0 to 1.
    bootstrap_size : int or None
        if None, don't do bootstrapping.
        if int, resampling with size of bootstrap_size.
    max_bf : float or int
        upper bound of Bayes factor.
    min_unique : int
        minimum number of unique samples in both conditions
        if any unique samples is smaller than min_unique , 
        supporting Null hypothesis.
        
    Returns
    -------
    bayes_factor : float
        The Bayes factor supporting alternative hypothesis
        abs(condtion1 - conditon2) > max_diff

    #TODO: the decimals make a huge difference!
    #Maybe the inverval can be less sensitive!
    """
    if type(samples1) is not list:
        samples1 = [samples1]
    if type(samples2) is not list:
        samples2 = [samples2]
    for i in range(len(samples1)):
        if len(np.array(samples1[i]).shape) == 1:
            samples1[i] = samples1[i].reshape(-1,1)
    for i in range(len(samples2)):
        if len(np.array(samples2[i]).shape) == 1:
            samples2[i] = samples2[i].reshape(-1,1)
    k = samples1[0].shape[1]

    for samp in samples1 + samples2:
        for j in range(samp.shape[1]):
            if len(np.unique(samp[:,j])) < min_unique:
                print("Improper sampling of psi! Supporting Null hypothesis.")
                return 0 #None
            
    if is_logit is False:
        for s in range(len(samples1)):
            samples1[s] = logistic(samples1[s])
        for s in range(len(samples2)):
            samples2[s] = logistic(samples2[s])
    
    samp_mu = np.zeros((len(samples1)*len(samples2), k))
    samp_cov = np.zeros((k, k))
    samp_diff = np.zeros((0, k))
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            if bootstrap_size is None:
                # idx1 = np.random.permutation(len(samples1[i]))
                # idx2 = np.random.permutation(len(samples2[j]))
                idx1 = np.arange(len(samples1[i]))
                idx2 = np.arange(len(samples2[j]))
            else:
                idx1 = np.random.randint(len(samples1[i]), size=bootstrap_size)
                idx2 = np.random.randint(len(samples2[j]), size=bootstrap_size)
            samp_temp = samples1[i][idx1,:] - samples2[j][idx2,:]
            samp_cov += np.cov(samp_temp.T)
            samp_mu[i*len(samples1)+j,:] = np.mean(samp_temp, axis=0)
            samp_diff = np.append(samp_diff, samp_temp, axis=0)
        
    # define biological variance
    if bio_cov is None:
        if samp_mu.shape[0] == 1:
            bio_cov = np.diag(np.zeros(k))
        else:
            bio_cov = np.diag(np.var(samp_mu, axis=0))
    sum_cov = samp_cov / samp_mu.shape[0] + bio_cov
    
    mu = np.zeros(k)
    cov = GP_K(np.arange(k), prior_param)
    if time_corr is False:
        cov = np.diag(np.diag(cov))
        sum_cov = np.diag(np.diag(sum_cov))
        
    if max_diff == 0:
        if prior_mode.lower() == "uniform":
            diff_prior_log = 0.0 # check it!!!
        else: #GP
            diff_prior_log = normal_pdf(mu, mu, 2*cov, log=True)
        
        if post_mode.lower() == "kde":
            posterior_density = gaussian_kde(samp_diff.T, bw_method='scott')
            diff_posterior_log = np.log(posterior_density.evaluate(mu)[0])
        elif post_mode.lower() == "mean":
            diff_posterior_log = np.log(np.mean(normal_pdf(
                samp_diff.mean(axis=0), mu, sum_cov, log=False)))
        else: #sample
            diff_posterior_log = np.log(np.mean(normal_pdf(
                samp_diff, mu, sum_cov, log=False)))
    else: #integral
        upp = np.ones(k) * max_diff
        diff_prior_log = np.log(mvn.mvnun(-upp, upp, mu, 2*cov)[0])
        
        temp = np.zeros(samp_diff.shape[0])
        for i in range(temp.shape[0]):
            temp[i] = np.log(mvn.mvnun(-upp, upp, samp_diff[i,:], sum_cov)[0])
        diff_posterior_log = np.sum(temp)
    
    if diff_posterior_log == -np.inf:
        bayes_factor = max_bf
    elif diff_posterior_log == np.inf:
        bayes_factor = -max_bf
    else:
        bayes_factor = diff_prior_log - diff_posterior_log
    if log is False:
        bayes_factor = np.exp(bayes_factor)

    return bayes_factor



def get_BioVar(samp_mean1, samp_mean2, share=True, diagonal=True):
    """get the biological variance.
    """
    N = samp_mean1.shape[0]
    T = samp_mean1.shape[1]
    R1 = samp_mean1.shape[2]
    R2 = samp_mean2.shape[2]
    mean_n_diff = np.zeros((N, T*2+1))
    mean_n_diff[:,:T] = samp_mean1.mean(axis=2)
    mean_n_diff[:,T:T*2] = samp_mean2.mean(axis=2)

    _mean = mean_n_diff[:,:T] + mean_n_diff[:,T:T*2]
    _diff = mean_n_diff[:,:T] - mean_n_diff[:,T:T*2]
    mean_n_diff[:,T*2] = np.sqrt(np.sum(_diff**2, axis=1))

    samp_diff_all = np.zeros((N,T,R1*R2))
    for i in range(R1):
        for j in range(R2):
            samp_diff_all[:,:,i*R1+j] = samp_mean1[:,:,i] - samp_mean2[:,:,j]
    bio_var = np.zeros((N, T, T))
    if R1*R2 > 1:
        for i in range(T):
            if diagonal is True:
                bio_var[i,:,:] = np.diag(np.var(samp_diff_all[i,:,:], axis=1))
            else:
                bio_var[i,:,:] = np.cov(samp_diff_all[i,:,:])

    if share is True and R1*R2 > 1:
        sort_idx = np.argsort(np.sum(_mean**2, axis=1))
        for i in range(N):
            idx1 = max(0, i-250)
            idx2 = min(N, i+250)
            idx_use = sort_idx[idx1:idx2]
            bio_var[sort_idx[i],:,:] = np.mean(bio_var[idx_use,:,:], axis=0)

    return bio_var, mean_n_diff
