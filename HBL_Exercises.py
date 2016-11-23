
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import emcee
import scipy.optimize as op

def Exercise_One(Data_input_copy):

    #Split the data up into the various components
    Data_input = Data_input_copy[['ID', 'x', 'y', 'yerr']]
    Data_input = Data_input.iloc[4:]

    #Perform the Least Squares Calculations
    A = np.vstack((np.ones_like(Data_input.x), Data_input.x)).T
    C = np.diag(Data_input.yerr * Data_input.yerr)
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
    b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, Data_input.y)))

    #Create model, and calculate the std uncertainty on the slope
    model_x_values = np.linspace(0, np.max(Data_input.x), 1000)
    model_ls = model_x_values * m_ls + b_ls
    sigma_m = math.sqrt(cov[1,1])

    #Print the results
    print('The gradient found is:', m_ls,
            '\nThe y-intercept found is:', b_ls,
            '\nThe standard uncertainty is:', cov[1,1])

    #Plot the figures
    fig, ax = plt.subplots()
    ax.errorbar(Data_input.x, Data_input.y, Data_input.yerr, fmt='.k')
    ax.set(xlabel='x_values',
            ylabel='y_values',
            title='Exercise 1')
    plt.plot(model_x_values, model_ls, color='b', linewidth=1.0)
    plt.show()
    plt.clf()

def Exercise_Two(Data_input_copy):

    #Split the data up into the various components
    Data_input = Data_input_copy[['ID', 'x', 'y', 'yerr']]

    #Perform the Least Squares Calculations
    A = np.vstack((np.ones_like(Data_input.x), Data_input.x)).T
    C = np.diag(Data_input.yerr * Data_input.yerr)
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
    b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, Data_input.y)))

    #Create model, and calculate the std uncertainty on the slope
    model_x_values = np.linspace(0, np.max(Data_input.x), 1000)
    model_ls = model_x_values * m_ls + b_ls
    sigma_m = math.sqrt(cov[1,1])

    #Print the results
    print('The gradient found is:', m_ls,
            '\nThe y-intercept found is:', b_ls,
            '\nThe standard uncertainty is:', cov[1,1])

    #Plot the figures
    fig, ax = plt.subplots()
    ax.errorbar(Data_input.x, Data_input.y, Data_input.yerr, fmt='.k')
    ax.set(xlabel='x_values',
            ylabel='y_values',
            title='Exercise 2')
    plt.plot(model_x_values, model_ls, color='b', linewidth=1.0)
    plt.show()
    plt.clf()

    #Return the LS model parameters
    return [m_ls, b_ls]

def Exercise_Three(Data_input_copy):

    #Split the data up into the various components
    Data_input = Data_input_copy[['ID', 'x', 'y', 'yerr']]
    Data_input = Data_input.iloc[4:]

    #Perform the Least Squares Calculations
    A = np.vstack((np.ones_like(Data_input.x), Data_input.x, (Data_input.x)**2)).T
    C = np.diag(Data_input.yerr * Data_input.yerr)
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
    b_ls, m_ls, q_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, Data_input.y)))

    #Create model, and calculate the std uncertainty on the slope
    model_x_values = np.linspace(0, np.max(Data_input.x), 1000)
    model_ls = (model_x_values)**2 * q_ls + model_x_values * m_ls + b_ls
    sigma_m = math.sqrt(cov[1,1])

    #Print the results
    print('The quadratic co-efficient found is:', q_ls,
            '\nThe x co-efficient found is:', m_ls,
            '\nThe y-intercept found is:', b_ls,
            '\nThe standard uncertainty is:', cov[1,1])

    #Plot the figures
    fig, ax = plt.subplots()
    ax.errorbar(Data_input.x, Data_input.y, Data_input.yerr, fmt='.k')
    ax.set(xlabel='x_values',
            ylabel='y_values',
            title='Exercise 3')
    plt.plot(model_x_values, model_ls, color='b', linewidth=1.0)
    filepath = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/'
    plt.savefig(filepath + 'HBL_ex3', bbox_inches='tight')
    plt.clf()

def Exercise_Six(Data_input_copy, LS_Values):

    #Set the column data values
    Data_input = Data_input_copy[['ID', 'x', 'y', 'yerr']]

    #Define the models x values
    model_x_values = np.linspace(0, np.max(Data_input.x), 1000)

    #Initialie the number of walkers and dimensions and starting postition for MCMC
    ndim, nwalkers = 5, 300
    nll = lambda *args: -lnprob(*args)

    #Minimize the -log probabilty density function i.e. maximize the prob densiy function
    result = op.minimize(nll, [LS_Values[0], LS_Values[1], 0.01, 400, 5], args=(Data_input.x, Data_input.y, Data_input.yerr))
    m_ml, b_ml, Pb_ml, Yb_ml, lnVb_ml = result["x"]

    #Define the starting position for each walker
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    #Read the Data in to seperate variables
    x = np.array(Data_input.x)
    y = np.array(Data_input.y)
    yerr = np.array(Data_input.yerr)

    #Run the EMCEE ensemble sampler and then perform the MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)

    #Remove the first 100 walks in the chain to exclude the burn-in phase
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    #Begin plotting, plots the top left, a colormap of the 2-D parameter space
    plt.figure(1)
    plt.subplot(221)
    plt.hist2d(samples[:,1],samples[:,0] , bins=60)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")

    #Plot the errorbar data with the model fit
    plt.subplot(222)
    for m, b, P_b, Y_b, lnV_b in samples[np.random.randint(len(samples), size=10)]:
        plt.plot(model_x_values, m*model_x_values+b, color="k", alpha=0.1)
    plt.errorbar(Data_input.x, Data_input.y, yerr=Data_input.yerr, fmt=".k")
    plt.xlabel("x")
    plt.ylabel("y")

    #Plot a histogram of Pb (Gaussian probability of bad point)
    plt.subplot(223)
    plt.hist(samples[:,2], normed = True)
    plt.xlabel("Pb")

    #Perform the same analysis with the error bars half as big
    yerr_2 = yerr/2
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr_2))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    #Create the Probability density function histogram plot for the smaller error bars
    plt.subplot(224)
    plt.hist(samples[:,2], normed = True)
    plt.xlabel("Pb")

    #Save the figure
    filepath = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/'
    plt.savefig(filepath + 'HBL_ex6', bbox_inches='tight')
    plt.clf()

def Exercise_Nine(Data_input_copy, LS_Values):

    #Load in the data
    Data_input = Data_input_copy[['ID', 'x', 'y', 'yerr']]
    Data_input = Data_input.iloc[4:]

    #Define the dimensions and walker numbers
    ndim, nwalkers = 5, 300

    #Seperate out the variables within the data
    x = np.array(Data_input.x)
    y = np.array(Data_input.y)
    yerr = np.array(Data_input.yerr)
    yerr_half = yerr / 2

    nll = lambda *args: -lnprob(*args)

    #Minimize the -log probabilty density function i.e. maximize the prob densiy function
    result = op.minimize(nll, [LS_Values[0], LS_Values[1], 0.01, 400, 5], args=(Data_input.x, Data_input.y, Data_input.yerr))
    m_ml, b_ml, Pb_ml, Yb_ml, lnVb_ml = result["x"]

    #Define the starting position for each walker
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    #Run the EMCEE ensemble sampler and then perform the MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)

    #Remove the first 100 results - excludes the burn-in phase
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    #Plot the results
    plt.figure(1)
    plt.subplot(211)
    plt.hist2d(samples[:,1],samples[:,0] , bins=100)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")

    #Perform the same MCMC iterations, but with error bars half the size
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr_half))
    sampler.run_mcmc(pos, 500)

    #Remove the first 100 results - excludes the burn-in phase
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    #Save the filepath
    filepath = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/'

    #Create a new sub-plot within the plot
    plt.subplot(212)
    plt.hist2d(samples[:,1],samples[:,0] , bins=100)
    plt.colorbar()
    plt.xlabel("b")
    plt.ylabel("m")
    plt.savefig(filepath + 'HBL_ex9', bbox_inches='tight')
    plt.clf()

def lnlike(theta, x, y, yerr):

    #Define the model parameters
    m, b, P_b, Y_b, lnV_b = theta
    model = m*x + b

    #Define the equations
    fg = np.log((1-P_b)/(yerr*(2*np.pi)**0.5)) -0.5*((y-model)/yerr)**2
    bg = np.log(P_b/((2*np.pi*(yerr**2+np.exp(lnV_b)))**0.5)) -0.5*(1/(yerr**2+np.exp(lnV_b)))*(y-Y_b)**2

    #Return the summation of the individual arrays
    return np.sum(np.logaddexp(fg, bg))

def lnprior(theta):

    #Define the parameters
    m, b, P_b, Y_b, lnV_b = theta

    #If statement for the relevant values
    if 0.0 < m <5.0 and -100 < b < 500 and 0.0 < P_b < 1.0 and 0 <Y_b< 600 and 0 < lnV_b < 10:
        return 0.00
    return -np.inf

def lnprob(theta, x, y, yerr):

    #Call the prior funcion
    lp = lnprior(theta)

    #If finite result returned from prior, then call the likelihood
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

def main():

    #Read the data in
    Data_location = "C:/Users/User/Documents/University/Year 4/Project/HBL_Data.txt"
    Data_input_copy = pd.read_table(Data_location, sep=' ',header=0)

    #Call the Ex1 function
    #Exercise_One(Data_input_copy)

    #Call the Ex2 function
    Least_Squares_Values = Exercise_Two(Data_input_copy)

    #Call the Ex3 function, set the output equal to a variable
    #Exercise_Three(Data_input_copy)

    #Call the Ex6 function
    #Exercise_Six(Data_input_copy, Least_Squares_Values)

    #Call the Ex9 function
    Exercise_Nine(Data_input_copy, Least_Squares_Values)

main()
