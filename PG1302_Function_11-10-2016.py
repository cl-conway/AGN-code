""" --------------------------------------------------------------------------
Description:
File for taking data for specific object, performs Gaussian distribution of
magnitudes among the errors for each individual point. A Lomb-Scargle-Fast
model is created for each iteration. Relevant graphs saved to a specific
location and the best period for each LSF model is saved and printed.

Inputs: Raw Data (e.g. for PG1302-102)
http://nesssi.cacr.caltech.edu/DataRelease/upload/result_web_file5sP2FZ.csv

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
def Iteration_function(Data, number, Output_figs):

    #Create numpy arrays of the time and errors to pass to L-S model
    Times = Data[['MJD']].as_matrix()
    Errors = Data[['Magerr']].as_matrix()
    Mag_orgs = Data[['Mag']].as_matrix()

    #Define a new list to place the editted data into
    New_list = []
    New_array = np.empty([290,1], dtype=float)

    if number==0:
        #Build LSF model and generate file paths
        model_original = LombScargleFast().fit(Times, Mag_orgs, Errors)
        file_path_org = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/Light_Curve_org.jpg'
        file_path_org2 = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/LSF_org.jpg'

        #Create Light Curve
        if Output_figs == 1:
            fig, ax = plt.subplots()
            ax.errorbar(Times, Mag_orgs, Errors, fmt='.k', color="r", ecolor='gray')
            ax.set(xlabel='Time', ylabel='Magitude',
                title='PG1302-102 Light Curve Original')
            ax.invert_yaxis();
            fig.savefig(file_path_org, bbox_inches='tight')

        #Build original model
        periods_orig, power_orig = model_original.periodogram_auto(nyquist_factor=10)

        print(Mag_orgs.shape)

        #Plot the original LSF model
        if Output_figs == 1:
            fig, ax = plt.subplots()
            ax.plot(periods_orig, power_orig, color="r")
            ax.set(xlim=(0, 2000),
                    xlabel='Period (Days)',
                    ylabel='Lomb-Scargle Power',
                    title='LS Power against period plot');
            fig.savefig(file_path_org2, bbox_inches='tight')


        #Choose the best period, conditional on quality of model
        if power_orig.max() < 0.1:
            periods_orig = 0
        else:
            model_original.optimizer.period_range=(500, 2000)
            periods_orig = model_original.best_period

        #Read best period out
        return periods_orig

    #Select the Mag data and apply random gaussian
    for i in range(Data.Mag.count()):
        newvalue = random.gauss(Data['Mag'].iloc[[i]], Data['Magerr'].iloc[[i]])
        New_array[i,0] = newvalue
        #New_list.append(float(newvalue))

    #Create numpy array from the list, for L-S model input
    print(New_array)

    #Build L-S model
    model = LombScargleFast().fit(Times, New_array, Errors)

    #Apply the LS mechanism to the data
    periods, power = model.periodogram_auto(nyquist_factor=10)
    file_path = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/LSF_' + str(number) + '.jpg'
    file_path2 = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/Light_Curve_' + str(number) + '.jpg'
    iteration_title = 'PG1302-102 Light Curve with Red Noise iteration ' + str(number)

    #Plotting the Light Curve using seaborn style
    if Output_figs == 1:
        fig, ax = plt.subplots()
        ax.errorbar(Times, New_array, Errors, fmt='.k', color="r", ecolor='gray', label='Light Curve')
        ax.set(xlabel='Time (MJD)', ylabel='Magitude',
            title=iteration_title)
        ax.invert_yaxis();
        fig.savefig(file_path2, bbox_inches='tight')

        #Plotting the LS model
        fig, ax = plt.subplots()
        ax.plot(periods, power, color="r", label='LSF')
        ax.set(xlim=(0, 2000),
                xlabel='Period (Days)',
                ylabel='Lomb-Scargle Power',
                title='LS Power against period plot');
        fig.savefig(file_path, bbox_inches='tight')

    #Calculating the best period of the model, condition on quality of LSF model
    if power.max() < 0.1:
        period = 0
    else:
        model.optimizer.period_range=(500, 2000)
        period = model.best_period

    #print(period)

    return period

#Calling function calls the iteration function
def calling_function(Data, Iterations, Outputs):

    #Define a period list, to store results
    Period_List = []

    #Convert string variable to number variable
    if Outputs == "Y": Output_n = 1
    elif Outputs == "N": Output_n = 0
    else: print("Error in Output_figures variable")

    #Perform the iterations
    for j in range(0, Iterations):
        New_Period = Iteration_function(Data, j, Output_n)
        Period_List.append(int(New_Period))

    #Print the list
    print(Period_List)

    #Write the results to a text file
    results_file = open("C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016/Results.txt","w")
    results_file.write("\n".join(map(lambda x: str(x), Period_List)) + "\n")
    results_file.close()

#Read in the Raw Data from the relevant location
Raw_Data = pd.read_csv("http://nesssi.cacr.caltech.edu/DataRelease/upload/result_web_file5sP2FZ.csv", sep=',', header=0)

#Set the number of iterations here
No_of_Iterations = 4

#Whether to output figures or not (Y/N; note for large no. of iterations put N)
Output_figures = 'Y'

#Calling function runs program
calling_function(Raw_Data, No_of_Iterations, Output_figures)

"""
It appears that the random number generation works and indeed that the whole
process works. However, the L-S model, when presented with data that does not
contain a clear period, seems to default with a very small little peak at
a value of the period of 364.5 days. To overcome this, we could limit the data
to above 1000 days. Same problem still exists, but has now changed to 2155.1
this is not good.
"""
