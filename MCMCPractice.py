#Monte Carlo Practise
#Chris Conway 10/11/2016
import numpy as np
import matplotlib.pyplot as plt
import emcee

filepath = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/'

x = np.array([201, 244, 47 ,287, 203, 58, 210, 202, 198, 158, 165, 201, 157, 131, 166, 160, 186, 125, 218, 146])
y = np.array([592, 401, 583, 402, 495, 173, 479, 504, 510, 416, 393, 442, 317, 311, 400, 337, 423, 334, 533, 344])
yerr = np.array([61, 25, 38, 15, 21, 15, 27, 14, 30, 16, 14, 25, 52, 16, 34, 31, 42, 26, 16, 22])
xerr = np.array([9, 4, 11, 7, 5, 9, 4, 4, 11, 7, 5, 5, 5, 6, 6, 5, 9, 8, 6, 5])
rho = np.array([-0.84, 0.31, 0.64, -0.27, -0.33, 0.67, -0.02, -0.05, -0.84, -0.69, 0.3, -0.46, -0.03, 0.5, 0.73, -0.52, 0.9, 0.4, -0.78, -0.56])

#Without outliers
x_1= np.delete(x, range(4), None) #array as above with outliers deleted
y_1= np.delete(y, range(4), None)#array as above with outliers deleted
yerr_1= np.delete(yerr, range(4), None)#array as above with outliers deleted
xerr_1= np.delete(xerr, range(4), None)#array as above with outliers deleted

def ls_fit(x, y, yerr, exercise, plot, quad):  # least squares straight line plot
    if quad == True:
        A = np.vstack((np.ones_like(x), x, x**2)).T
    else:
        A = np.vstack((np.ones_like(x), x)).T #matrix for linear equations
    C = np.diag(yerr * yerr)# diagonal covariance across y values
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A))) # covariance matrix of b and m
    b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))

    #print( exercise + " standard m uncertainty ", cov[1,1])
    if plot == True:
        #plot data with best fit overlayed
        plt.errorbar(x, y, yerr, fmt='.k')
        xlin = np.linspace(np.min(x), np.max(x), 1000) # x values for plotting
        if quad == True:
            plt.plot(xlin, q_ls*xlin**2 + xlin*m_ls + b_ls)
        else:
            plt.plot(xlin, xlin*m_ls + b_ls)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("HBL " + exercise)
        plt.show()
        plt.clf()
    return m_ls, b_ls

#Exercise 1
#ls_fit(x_1, y_1, yerr_1, 'Ex1', plot=True, quad=False)

#Excersise 2
#repeat process of exercise 1 with data that includes outliers
#ls_fit(x,y,yerr,'Ex2',plot = True, quad=False)
#Excerise 3
#plot with outliers assuming a quadractic rather than linear relationship
#ls_fit(x_1,y_1,yerr_1,'Ex3',plot = True, quad=True)

#Multigaussian Monte Carlo Funcitons
#function for log of likelihood
def lnlike(theta, x, y, yerr):
    m, b, P_b, Y_b, lnV_b = theta
    model = m * x + b
    fg = np.log((1-P_b)/(yerr*(2*np.pi)**0.5)) -0.5*((y-model)/yerr)**2
    bg = np.log(P_b/((2*np.pi*(yerr**2+np.exp(lnV_b)))**0.5)) -0.5*(1/(yerr**2+np.exp(lnV_b)))*(y-Y_b)**2
    return np.sum(np.logaddexp(fg, bg))

#function for log of priors
def lnprior(theta):
    m, b, P_b, Y_b, lnV_b = theta
    if (0.0 < m <5.0 and -100 < b  <500 and 0.0 < P_b < 1.0 and 0<Y_b< 600 and 0 < lnV_b < 10):
        return 0.00
    return -np.inf
# function for log of posterior
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

m_ls, b_ls =  ls_fit(x,y,yerr,'Ex2',plot = False, quad=False)
MCinit = np.array([m_ls, b_ls, 0.1, 400, 5]) # initial estimate for MCMC

import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
result = op.minimize(nll, MCinit, args=(x, y, yerr))
m_ml, b_ml, Pb_ml, Yb_ml, lnVb_ml = result["x"]

#print(m_ml, b_ml, Pb_ml, Yb_ml, lnVb_ml)

ndim, nwalkers = 5, 300
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

def Ex6and7(x, y, yerr, ndim, nwalkers, pos):
    mls, bls = ls_fit(x,y,yerr,'Ex2',plot = False, quad=False)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)
    """
    samplesa = sampler.chain
    for j in range(ndim):
        plt.figure(j)
        for i in range(nwalkers):
            plt.plot(range(500), samplesa[i, : , j])
        plt.show()
    """

    samples = sampler.chain[:, 100:, :].reshape((-1, ndim)) # remove first 100 (burn in) #ask will to explain this line
    plt.figure(1)
    plt.subplot(221)
    plt.hist2d(samples[:,1],samples[:,0] , bins=50)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")

    plt.subplot(222)
    for m, b, P_b, Y_b, lnV_b in samples[np.random.randint(len(samples), size=10)]:
        plt.plot(xlin, m*xlin +b, color="k", alpha=0.1)
    plt.errorbar(x, y, yerr=yerr, fmt=".k")
    plt.xlabel("x")
    plt.ylabel("y")

    plt.subplot(223)
    plt.hist(samples[:,2])
    plt.xlabel("Pb")

    yerr_2 = yerr/2
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr_2))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    plt.subplot(224)
    plt.hist(samples[:,2])
    plt.xlabel("Pb")

    plt.savefig(filepath + 'HBL_ex6', bbox_inches='tight')
    plt.clf()
#Ex6and7(x,y,yerr,ndim,nwalkers,pos)
#-----------------------------------------------------------------------------------
#Ex8
def Ex8(x,y,yerr):
    mjk_trials = []
    #delete one data point at a time and recalculate m
    for i in range(np.size(y)):
        y_jk = np.delete(y,i) # y jacknife
        x_jk = np.delete(x,i)
        yerr_jk = np.delete(yerr, i)
        dummy_m, dummy_b = ls_fit(x_jk,y_jk,yerr_jk,'Ex8',plot = False, quad=False)
        mjk_trials.append(dummy_m)
    mjk_trials = np.array(mjk_trials)
    mjk = np.sum(mjk_trials)/np.size(mjk_trials)
    sigma_mjk = ((np.size(mjk_trials)-1)/np.size(mjk_trials))*np.sum( (mjk_trials-mjk)**2 )
    print(mjk, sigma_mjk)
#Ex8(x,y,yerr)
def Ex9(x,y,yerr):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_1, y_1, yerr_1))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    plt.figure(1)
    plt.subplot(211)
    plt.hist2d(samples[:,1],samples[:,0] , bins=100)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_1, y_1, yerr_1/2))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    plt.subplot(212)
    plt.hist2d(samples[:,1],samples[:,0] , bins=100)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")
    plt.show()
#Ex9(x_1,y_1,yerr_1)
def Chi_sq(x,y,yerr):
    m_ls, b_ls = ls_fit(x,y,yerr, 'Ex10', plot = False, quad = False)
    X= np.array([b_ls,m_ls]).reshape(2,1)
    A = np.vstack((np.ones_like(x), x)).T #matrix for linear equations
    C = np.diag(yerr * yerr)
    Y=y.reshape(np.size(y),1)
    dummy = Y-np.dot(A,X)
    chi_sq = np.dot(dummy.T, np.linalg.solve(C, dummy))
    return chi_sq
#Ex10
#print('No outliers Chi^2 = ',Chi_sq(x_1,y_1,yerr_1))
#print('With outliers Chi^2 = ',Chi_sq(x,y,yerr))

#Ex11
def Ex11(x_1,y_1):
    Chi_sq_S = []
    S = np.linspace(1,1500,1500)
    for i in S:
        ones=np.ones_like(yerr_1)
        Chi_sq_S.append(Chi_sq(x_1, y_1, (ones*i)**0.5)[0][0])

    Chi_sq_S=np.array(Chi_sq_S)
    plt.plot(S, Chi_sq_S)
    plt.xlabel('S')
    plt.ylabel('Chi^2')
    plt.ylim((6,24))
    plt.show()
