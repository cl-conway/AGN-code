#Monte Carlo Practise
#Chris Conway 10/11/2016
import numpy as np
import matplotlib.pyplot as plt
import emcee

x = np.array([201, 244, 47 ,287, 203, 58, 210, 202, 198, 158, 165, 201, 157, 131, 166, 160, 186, 125, 218, 146])
y = np.array([592, 401, 583, 402, 495, 173, 479, 504, 510, 416, 393, 442, 317, 311, 400, 337, 423, 334, 533, 344])
yerr = np.array([61, 25, 38, 15, 21, 15, 27, 14, 30, 16, 14, 25, 52, 16, 34, 31, 42, 26, 16, 22])
xerr = np.array([9, 4, 11, 7, 5, 9, 4, 4, 11, 7, 5, 5, 5, 6, 6, 5, 9, 8, 6, 5])
rho = np.array([-0.84, 0.31, 0.64, -0.27, -0.33, 0.67, -0.02, -0.05, -0.84, -0.69, 0.3, -0.46, -0.03, 0.5, 0.73, -0.52, 0.9, 0.4, -0.78, -0.56])

# Excersise 1
x_1= np.delete(x, range(4), None) #array as above with outliers deleted
y_1= np.delete(y, range(4), None)#array as above with outliers deleted
yerr_1= np.delete(yerr, range(4), None)#array as above with outliers deleted
xerr_1= np.delete(xerr, range(4), None)#array as above with outliers deleted


A = np.vstack((np.ones_like(x_1), x_1)).T #matrix for linear equations
C = np.diag(yerr_1 * yerr_1) # diagonal covariance across y values
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A))) # covariance matrix of b and m
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y_1)))

xlin = np.linspace(0, np.max(x), 1000) # x values for plotting

print("standard m uncertainty ", cov[1,1])
"""
#plot data witout ouriters
plt.errorbar(x_1, y_1, yerr_1, fmt='.k')
#overlay plot of least squares fit
plt.plot(xlin, xlin*m_ls + b_ls)
plt.xlabel("x")
plt.ylabel("y")
plt.title("HBL Ex1")
plt.show()
plt.clf()

"""
#Excersise 2
"""
#repeat process of exercise 1 with data that includes outliers
A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))

print("standard m uncertainty ", cov[1,1])
#plot data with outlier and best fit overlayed
plt.errorbar(x, y, yerr, fmt='.k')
xlin = np.linspace(np.min(x), np.max(x), 1000)
plt.plot(xlin, xlin*m_ls + b_ls)
plt.xlabel("x")
plt.ylabel("y")
plt.title("HBL Ex2")
plt.show()
plt.clf()
"""
#Excerise 3
"""
#plot with outliers assuming a quadractic rather than linear relationship
A = np.vstack((np.ones_like(x_1), x_1, x_1**2)).T
C = np.diag(yerr_1 * yerr_1)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls, q_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y_1)))
print("cov matrix", cov)
plt.errorbar(x_1, y_1, yerr_1, fmt='.k')
xlin = np.linspace(0, np.max(x_1), 1000)
plt.plot(xlin, q_ls*xlin**2 + xlin*m_ls + b_ls)
plt.xlabel("x")
plt.ylabel("y")
plt.title("HBL Ex2")
plt.show()
plt.clf()
"""

#Excersise 6

#function for log of likelihood
def lnlike(theta, x, y, yerr):
    m, b, P_b, Y_b, lnV_b = theta
    model = m * x + b
    fg = np.sum(np.log((1-P_b)/(yerr*(2*np.pi)**0.5)))+np.sum( -0.5*(1/yerr**2)*(y-model)**2)
    bg = np.sum(np.log(P_b/(2*np.pi*(yerr**2+np.exp(lnV_b)))**0.5))+np.sum( -0.5*(1/(yerr+np.exp(lnV_b))**2)*(y-Y_b)**2 )
    return np.logaddexp(fg, bg)

#function for log of priors
def lnprior(theta):
    m, b, P_b, Y_b, lnV_b = theta
    if (0.0 < m <5.0 and -100 < b  <500 and 0 < P_b < 1 and 0 <Y_b< 600 and 0 < lnV_b < 10):
        return 0.00
    return -np.inf

# function for log of posterior
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_ls, b_ls, 0.1, 400, 5], args=(x, y, yerr))
m_ml, b_ml, Pb_ml, Yb_ml, lnVb_ml = result["x"]

MCinit = np.array([m_ls, b_ls, 0.1, 400, 5]) # initial estimate for MCMC
ndim, nwalkers = 5, 300
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 100:, :].reshape((-1, ndim)) # remove first 100 (burn in)

plt.figure(1)
plt.subplot(221)
plt.hist2d(samples[:,1],samples[:,0] , bins=100)
plt.colorbar()
plt.xlabel("b")
plt.ylabel("m")

plt.subplot(222)
for m, b, P_b, Y_b, lnV_b in samples[np.random.randint(len(samples), size=10)]:
    plt.plot(xlin, m*xlin+b, color="k", alpha=0.1)
plt.errorbar(x, y, yerr=yerr, fmt=".k")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(223)
plt.hist(samples[:,2], normed = True)
plt.xlabel("Pb")

yerr_2 = yerr/2
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr_2))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

plt.subplot(224)
plt.hist(samples[:,2], normed = True)
plt.xlabel("Pb")

plt.show()
"""
#Ex9
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
"""
