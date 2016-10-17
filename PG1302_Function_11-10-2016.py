""" --------------------------------------------------------------------------
Description:
File for taking data for specific object, performs Gaussian distribution of
magnitudes among the errors for each individual point. A Lomb-Scargle-Fast
model is created for each iteration. Relevant graphs saved to a specific
location and the best period for each LSF model is saved and printed.
All variables to be changed in 'main' function. Rest of code will run
appropriately.

Inputs: Raw Data (e.g. for PG1302-102)
C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt
[Internet connection Required to run]

Exports to:
C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016
[Exports LSF image and Light Curve image for each iteration as well as the
results into a text file]

Date Created: 12/10/2016
Author: Nicholas Kinsey
----------------------------------------------------------------------------"""

#Perform the relevant imports
import pandas as pd
import sys
import csv as csv
import scipy
import random
import numpy as np
import matplotlib.pyplot as plt

from gatspy.periodic import LombScargleFast
import seaborn; seaborn.set()
from datetime import datetime

#Iteration function
def Iteration_function(Data, number, Output_figs, file_loc, Object):

    #Create numpy arrays of the time and errors to pass to L-S model
    Times = Data[['MJD']].as_matrix()
    Errors = Data[['Magerr']].as_matrix()
    Mag_orgs = Data[['Mag']].as_matrix()

    #Define a new list to place the editted data into
    New_array = np.empty([len(Data),1], dtype=float)

    #The first element is always the original data set
    if number==0:

        #Build LSF model and generate file paths
        model_original = LombScargleFast().fit(Times, Mag_orgs, Errors)
        file_path_org = file_loc + '/Light_Curve_org_' + Object + '.jpg'
        file_path_org2 = file_loc + '/LSF_org_' + Object + '.jpg'

        #Build original model
        periods_orig, power_orig = model_original.periodogram_auto(nyquist_factor=len(Data),oversampling=50)

        #Create Light Curve
        if Output_figs == 1:
            fig, ax = plt.subplots()
            ax.errorbar(Times, Mag_orgs, Errors, fmt='.k', color="r", ecolor='gray')
            ax.set(xlabel='Time', ylabel='Magitude',
                title= Object +' Light Curve Original')
            ax.invert_yaxis();
            fig.savefig(file_path_org, bbox_inches='tight')

            #Plot the original LSF model
            fig, ax = plt.subplots()
            ax.plot(periods_orig, power_orig, color="r")
            ax.set(xlim=(0, 2000),
                    xlabel='Period (Days)',
                    ylabel='Lomb-Scargle Power',
                    title=Object + ' LS Power against period plot');
            fig.savefig(file_path_org2, bbox_inches='tight')


        #Choose the best period, conditional on quality of model
        if power_orig.max() < 0.1:
            periods_orig = 0
        else:
            model_original.optimizer.period_range=(100, 2000)
            periods_orig = model_original.best_period

        #Read best period out
        return periods_orig

    #Select the Mag data and apply random gaussian
    for i in range(Data.Mag.count()):
        newvalue = random.gauss(Data['Mag'].iloc[[i]], Data['Magerr'].iloc[[i]])
        New_array[i,0] = newvalue

    #Build L-S model
    model = LombScargleFast().fit(Times, New_array, Errors)

    #Apply the LS mechanism to the data
    periods, power = model.periodogram_auto(nyquist_factor=len(Data), oversampling=50)
    file_path_Light_Curve = file_loc + '/Light_Curve_' + Object + '_' + str(number) + '.jpg'
    file_path_LS_model = file_loc + '/LSF_' + Object + '_' + str(number) + '.jpg'

    #Plotting the Light Curve using seaborn style
    if Output_figs == 1:
        fig, ax = plt.subplots()
        ax.errorbar(Times, New_array, Errors, fmt='.k', color="r", ecolor='gray', label='Light Curve')
        ax.set(xlabel='Time (MJD)', ylabel='Magitude',
            title=Object + ' Light Curve with Red Noise iteration ' + str(number))
        ax.invert_yaxis();
        fig.savefig(file_path_Light_Curve, bbox_inches='tight')

        #Plotting the LS model
        fig, ax = plt.subplots()
        ax.plot(periods, power, color="r", label='LSF')
        ax.set(xlim=(0, 2000),
                xlabel='Period (Days)',
                ylabel='Lomb-Scargle Power',
                title=Object + ' LS Power against period plot');
        fig.savefig(file_path_LS_model, bbox_inches='tight')

    #Calculating the best period of the model, condition on quality of LSF model
    if power.max() < 0.1:
        period = 0
    else:
        model.optimizer.period_range=(100, 2000)
        period = model.best_period

    return period

    #End of Iteration_function

#Calling function calls the iteration function
def calling_function(Data, Iterations, Outputs, Output_loc, Object):

    #Define a period list, to store results
    Period_List = []

    #Convert string variable to number variable
    if Outputs == "Y": Output_n = 1
    elif Outputs == "N": Output_n = 0
    else: print("Error in Output_figures variable")

    #Perform the iterations
    for j in range(0, Iterations):
        New_Period = Iteration_function(Data, j, Output_n, Output_loc, Object)
        Period_List.append(int(New_Period))

    #Call the LS_hist function in order to create the Histogram
    LS_hist(Period_List, Output_loc, Object)

    #End of calling_function

#Historgram function plots a histogram of the created data
def LS_hist(Results, file_loc, Object):

    #Print the number of zero counts within the sample
    print("Zero counts :", Results.count(0))

    #Remove all the zero counts
    for j in range(Results.count(0)):
        Results.remove(0)

    #Plot the Historgram appropriately
    fig, ax = plt.subplots()
    plt.hist(Results, bins='auto')
    plt.title(Object + " Histogram with 'auto' bins")
    file_path_hist = file_loc + '\Hist ' + Object + '.jpg'
    plt.savefig(file_path_hist, bbox_inches='tight')

    #End of LS_hist function

#Main function to be called
def main():

    #Define the location of the saved path
    Data_location = 'C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt'

    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #save each column within the dataframe as a seperate numpy array
    Objects = Raw_Data_input[['Name']].as_matrix()
    URLs = Raw_Data_input[['URL']].as_matrix()

    #Set the number of iterations here
    No_of_Iterations = 1000

    #Whether to output Light curve and L-S figures or not (Y/N)
    Output_figures = 'N'

    #Output_location
    Output_locations = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016'

    #Call the calling function for each object in the data
    for k in range(0, len(Raw_Data_input)):
        Object = Objects.item(k,0)
        URL = URLs.item(k,0)
        URL_editted = URL.strip()
        Raw_Data = pd.read_csv(URL_editted, sep=',', header=0)

        #Calling function runs program
        calling_function(Raw_Data, No_of_Iterations, Output_figures, Output_locations, Object)

    #End of main function

main()

#End of File
