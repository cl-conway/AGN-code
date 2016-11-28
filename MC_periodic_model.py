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

#Multigaussian Monte Carlo Funcitons
#function for log of likelihood
def lnlike(theta, t, t_0, mag, mag_err):
    A, phi, p, mu, nu = theta
    model = A*np.sin((1/p)*2*np.pi*(t-t_0) + phi) + mu
    fg = -np.log(nu*mag_err*(2*np.pi**0.5)) -0.5*((mag-model)/(nu*mag_err))**2
    return(np.sum(fg))

#function for log of priors
def lnprior(theta,pmax,pmin,mumax,mumin):
    A, phi, p, mu, nu = theta
    if 0 < A < 10*(mumax-mumin) and 0 < phi <2*np.pi and pmin < p < pmax and mumin < mu < mumax and 0.1 < nu <10:

        return 0.00
    return -np.inf
# function for log of posterior
def lnprob(theta,pmax,pmin,mumax,mumin, t, t_0, mag, mag_err):
    lp = lnprior(theta,pmax,pmin,mumax,mumin)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, t, t_0, mag, mag_err)

#Iteration function
def MCMC(Data, Output_Lightcurve, file_loc, Object):

    #Create numpy arrays of the time and errors to pass to L-S model
    Times = Data[['MJD']].as_matrix().ravel()
    Mag_Err = Data[['Magerr']].as_matrix().ravel()
    Mag = Data[['Mag']].as_matrix().ravel()

    #Create Light Curve
    if Output_Lightcurve == True:
        plt.figure()
        plt.errorbar(Times, Mag, Mag_Err, fmt='.k')
        file_path_Light_Curve = file_loc + '\Light_Curve_' + Object + '.jpg'
        #Saving the figure
        plt.savefig(file_path_Light_Curve, bbox_inches='tight')
        plt.clf()

    file_path_period = file_loc + '\Period_pdf_' + Object + '.jpg'

    pmax = Times[np.size(Times)-1] - Times[0] # Maximum period is total obervation time
    pmin = Times[1]-Times[0] #Minimum period is the shortest time interval between two data points
    for i in range(1, np.size(Times)-2):
        difference = Times[i+1]-Times[i]
        if difference<pmin:
            pmin = difference

    #Mu is standard offset of the sinusoid and can be between the highest and lowest magnitude values
    mumax = np.max(Mag)
    mumin = np.min(Mag)

    # t0 is an offset time
    t0 = Times[0]

    MCinit = [(mumax-mumin)/2, 0, pmax/2, (mumax+mumin)/2, 1] # Estimate of maximum for initialiser for optimiser

    import scipy.optimize as op
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, MCinit, args=(pmax,pmin,mumax,mumin, Times, t0, Mag, Mag_Err))
    A_ml, phi_ml, p_ml, mu_ml, nu_ml = result["x"]

    ndim, nwalkers = 5, 500
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)] # Initial position of walkers

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(pmax,pmin,mumax,mumin,Times,t0, Mag, Mag_Err))
    sampler.run_mcmc(pos, 800)
    """
    #for loop to show plot of walker steps for each dimension
    samplesa = sampler.chain
    dim_names = ['A','phi','p','mu','nu']
    for j in range(ndim):
        plt.figure(j)
        for i in range(nwalkers):
            plt.plot(range(800), samplesa[i, : , j])
            plt.title(dim_names[j])
        plt.show()
    """
    samples = sampler.chain[:, 250:, :].reshape((-1, ndim)) # remove first 100 (burn in) #ask will to explain this line

    plt.figure()
    plt.hist(samples[:,2], normed=True) # plot marginalised over all but period
    plt.xlabel("p")
    plt.savefig(file_path_period, bbox_inches='tight')
    plt.clf()

    return np.max(samples[:,4])

def main():
    #Set the user of the program here
    User = 'C'
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

    Output_Lightcurve = False
    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #save each column within the dataframe as a seperate numpy array
    Objects = Raw_Data_input[['Name']].as_matrix()
    URLs = Raw_Data_input[['URL']].as_matrix()

    #Open up text file to write number of periods for each object
    Error_Scaling_data = open('Nu_values.txt','w')
    Error_Scaling_data.write("Object    nu\n")

    #Call the calling function for each object in the data
    for k in range(0, len(Raw_Data_input)):
        Object = Objects.item(k,0)
        URL = URLs.item(k,0)
        URL_editted = URL.strip()
        Raw_Data = pd.read_csv(URL_editted, sep=',', header=0)

        if 'PG 1302' in Object:
            #Calling function runs program
            error_scaling = MCMC(Raw_Data,Output_Lightcurve, Output_locations, Object)
            Error_Scaling_data.write( Object + ', ' + str(error_scaling) + '\n')
            break

        #Write the period data into the text file
        #Period_data.write(Object + ' ' + str(pnumber) + '\n')

    #Close the period data file
    Error_Scaling_data.close()

#End of main function
main()

#End of File
