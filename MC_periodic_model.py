""" --------------------------------------------------------------------------
Description:
Similar to PG1302-102program this program takes data from a list of periodic
candidates found by Graham creates a posterior of a sinusoid model of the data
with a flat priors. The posterior marginalised over all parameters other than
period found with an MCMC and a normalised histogram of this is saved.

Inputs: Raw Data (e.g. for PG1302-102)
C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt
C:/Users/User/Documents/University/Year 4/Project/Graham_Candidates_names.txt
or
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Graham_Candidates_names.txt
[Internet connection Required to run]

Exports to:
C:/Users/User/Documents/University/Year 4/Project/MC_periodic_model
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/MC_periodic_model
[Exports LSF image and Light Curve image for each iteration as well as the
results into a text file]

Date Created: 26/11/2016
Authors: Christopher Conway, Nicholas Kinsey
----------------------------------------------------------------------------"""

#Perform the relevant imports
import pandas as pd
import sys
import csv as csv
import scipy
import random
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op

#Multigaussian Monte Carlo Funcitons
#Function for log of likelihood
def lnlike(theta, t, t_0, mag, mag_err):
    #Define the parameters: Amplitude, phase, period, sinusoid offset and errors scaling
    A, phi, p, mu, nu = theta

    #Create the model (sinusoid, defined from set time)
    model = A*np.sin((1/p)*2*np.pi*(t-t_0) + phi) + mu

    #Perform the log of the gaussian function
    fg = -np.log(nu*mag_err*(2*np.pi**0.5)) -0.5*((mag-model)/(nu*mag_err))**2

    #Return the summation of the log gaussian (the product)
    return(np.sum(fg))

    #End of lnlike function

#Function for log of priors
def lnprior(theta,pmax,pmin,mumax,mumin):
    A, phi, p, mu, nu = theta

    #Chosen from physical intuition
    if 0 < A < 10*(mumax-mumin) and 0 < phi <2*np.pi and pmin < p < pmax and mumin < mu < mumax and 0.1 < nu <10:
        return 0.00

    #Return -infinite if above evaluates to false
    return -np.inf

    #End of lnprior function

#Function for log of posterior
def lnprob(theta,pmax,pmin,mumax,mumin, t, t_0, mag, mag_err):
    #Call the prior
    lp = lnprior(theta,pmax,pmin,mumax,mumin)

    #If infinite return infinite value for lnprob function
    if not np.isfinite(lp):
        return -np.inf

    #Else return the prior value of 0.00 + the evaluation of the likelihood
    return lp + lnlike(theta, t, t_0, mag, mag_err)

    #End of lnprob function

#MCMC function
def MCMC(Data, Output_Lightcurve, Output_Variable_Figure, file_loc, Object):

    #Create numpy arrays of the time and errors to pass to L-S model
    Times = Data[['MJD']].as_matrix().ravel()
    Mag_Err = Data[['Magerr']].as_matrix().ravel()
    Mag = Data[['Mag']].as_matrix().ravel()

    #Define the file paths for saving various outputs
    file_path_period = file_loc + '\Period_pdf_' + Object + '.jpg'
    file_path_variables = file_loc + '\Variable_plots_' + Object + '.jpg'
    file_path_Light_Curve = file_loc + '\Light_Curve_' + Object + '.jpg'

    #Define the max and min periods as whole observation time and smallest time difference between two observations
    pmax = Times[np.size(Times)-1] - Times[0]
    pmin = np.min(np.diff(np.sort(Times)))

    #Mu is standard offset of the sinusoid and can be between the highest and lowest magnitude values
    mumax = np.max(Mag)
    mumin = np.min(Mag)
    muave = np.average(Mag)
    Amp_value = ((mumax - muave) + (muave - mumin)) / 2

    #Define an offset time, the earliest observation time for the data
    t0 = Times[0]

    #Estimate of maximum for initialiser for optimiser
    MCinit = [Amp_value, 0, pmax/2, muave, 1]
    #MCinit = [(mumax-mumin)/2, 0, pmax/2, (mumax+mumin)/2, 1]

    #Optimization over the various parameters
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, MCinit, args=(pmax, pmin, mumax, mumin, Times, t0, Mag, Mag_Err))

    #Save the optimized values into a variable
    A_ml, phi_ml, p_ml, mu_ml, nu_ml = result["x"]

    #Define the dimensions (no. of parameters) and no. of walkers
    ndim, nwalkers = 5, 800

    #Initial position of walkers
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    #Define the sampler and run the MCMC iteration, passing across the relevant arguments
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(pmax, pmin, mumax, mumin, Times, t0, Mag, Mag_Err))
    sampler.run_mcmc(pos, 800)

    #For loop to show plot of walker steps for each dimension
    if Output_Variable_Figure == True:
        #Define a new sample chain
        samplesa = sampler.chain

        #Create a diension names to assign titles
        dim_names = ['A','phi','p','mu','nu']

        #Create figure
        plt.figure()

        #Create a subplot for each dimension
        for j in range(ndim):
            for i in range(nwalkers):
                subplot = str(231 + j)
                plt.subplot(subplot)
                plt.plot(range(800), samplesa[i, : , j])
                plt.title(dim_names[j])

        #Save the figure to relevant location
        plt.savefig(file_path_variables, bbox_inches='tight')

    #Remove first 250 (burn in phase), ask will to explain this line
    samples = sampler.chain[:, 250:, :].reshape((-1, ndim))

    #Create a normalised histogram plot, marginalised over all parameters except for the period
    plt.figure()
    plt.hist(samples[:,2], normed=True)
    plt.xlabel("Period(days)")
    plt.savefig(file_path_period, bbox_inches='tight')

    #Set values for each parameter reagrding the value and counts in each bin
    n_A, bins_A = np.histogram(samples[:,0], normed=True)
    n_phi, bins_phi = np.histogram(samples[:,1], normed=True)
    n_p, bins_p = np.histogram(samples[:,2], normed=True)
    n_mu, bins_mu = np.histogram(samples[:,3], normed=True)
    n_nu, bins_nu = np.histogram(samples[:,4], normed=True)

    #Set the maximum marginalized parameters to be ussed in the model
    A_model = bins_A[np.argmax(n_A)]
    phi_model = bins_phi[np.argmax(n_phi)]
    p_model = bins_p[np.argmax(n_p)]
    mu_model = bins_mu[np.argmax(n_mu)]
    nu_model = bins_nu[np.argmax(n_nu)]

    #Create Light Curve
    if Output_Lightcurve == True:
        plt.figure()

        #Plot the error bars, using the original and scaled errors
        plt.errorbar(Times, Mag, Mag_Err, fmt='.k', ecolor='black')
        plt.errorbar(Times,Mag, nu_ml*Mag_Err, fmt='.k', ecolor='b')
        plt.errorbar(Times,Mag, nu_model*Mag_Err, fmt='.k', ecolor='g')

        #Create a linearly spaced array of x-values for the model
        x_lin = np.linspace(np.min(Times), np.max(Times), 1000)

        #Create the models for the optimized and MCMC likelihood values
        model_optimized = A_ml*np.sin((1/p_ml)*2*np.pi*(x_lin-t0) + phi_ml) + mu_ml
        model_mcmc_likelihood = A_model*np.sin((1/p_model)*2*np.pi*(x_lin-t0) + phi_model) + mu_model

        #Create the plots
        plt.plot(x_lin, model_optimized, linewidth=2.0, color = 'b', label='Optimized params')
        plt.plot(x_lin, model_mcmc_likelihood, linewidth=2.0, color = 'g', label = 'Ind MCMC max likelihood')

        #Insert legend and labels
        plt.legend(loc = 'lower right')
        plt.xlabel('Times')
        plt.ylabel('Magnitude')
        plt.title(Object + ' Light curve')

        #Save the figure
        plt.savefig(file_path_Light_Curve, bbox_inches='tight')
        plt.clf()

    #Close all open plots for that Object
    plt.close('all')

    #Return the error scaling parameters, to save to a text file
    return bins_nu[np.argmax(n_nu)], bins_p[np.argmax(n_p)]

    #End of MCMC function

def main():

    #Set the user of the program here
    User = 'N'

    #Define the location of the saved path
    if(User == 'C'):
        Data_location = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt'
    if(User == 'N'):
        Data_location = 'C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt'

    #Output_location
    if (User == 'N'):
        Output_locations = 'C:/Users/User/Documents/University/Year 4/Project/MC_periodic_model'
    elif (User == 'C'):
        Output_locations = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/MC_periodic_model'

    #Set whether to output lightcurve or the figure for the different parameters
    Output_Lightcurve = True
    Output_Variable_Figure = True

    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #Save each column within the dataframe as a seperate numpy array
    Objects = Raw_Data_input[['Name']].as_matrix()
    URLs = Raw_Data_input[['URL']].as_matrix()

    #Open up text file to write number of periods for each object
    MCMC_output_data = open('MCMC_output_values.txt','w')
    MCMC_output_data.write("Object;Error_Scaling;Period\n")

    #Call the MCMC function for each object in the data
    for k in range(0, 10):
        Object = Objects.item(k,0)
        URL = URLs.item(k,0)
        URL_editted = URL.strip()
        Raw_Data = pd.read_csv(URL_editted, sep=',', header=0)

        #MCMC function runs program
        error_scaling, MCMC_period = MCMC(Raw_Data, Output_Lightcurve, Output_Variable_Figure, Output_locations, Object)

        #Write the outputs from the MCMC function to a text file
        MCMC_output_data.write(Object + ';' + str(error_scaling) + ';' + str(MCMC_period) + '\n')

    #Close the error scaling data file
    MCMC_output_data.close()

#End of main function

main()

#End of File

#-----------------------------------------------------------------------------#
