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
or C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt
[Internet connection Required to run]

Exports to:
C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Iteration_figures_11-10-2016
[Exports LSF image and Light Curve image for each iteration as well as the
results into a text file]

Date Created: 12/10/2016
Authors: Nicholas Kinsey, Christopher Conway
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
from selenium import webdriver
from selenium.webdriver.common.keys import Keys

#Iteration function
def Iteration_function(Data, number, Output_figs, file_loc, Object):

    #Create numpy arrays of the time and errors to pass to L-S model
    Times = Data[['MJD']].as_matrix().ravel()
    Errors = Data[['Magerr']].as_matrix().ravel()
    Mag_orgs = Data[['Mag']].as_matrix().ravel()

    #The first element is always the original data set
    if number==0:

        #Build LSF model and generate file paths
        model_original = LombScargleFast().fit(Times, Mag_orgs, Errors)
        file_path_org = file_loc + '/Light_Curve_org_' + Object + '.jpg'
        file_path_org2 = file_loc + '/LSF_org_' + Object + '.jpg'

        #Build original model
        periods_orig, power_orig = model_original.periodogram_auto(nyquist_factor=0.1,oversampling=50)

        # find maximum of periodigram
        max_periods_orig = periods_orig[np.argmax(power_orig)]

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
            plt.clf()
        #Choose the best period, conditional on quality of model
        if power_orig.max() < 0.1:
            best_periods_orig = 0
        else:
            model_original.optimizer.period_range=(500, 2000)
            best_periods_orig = model_original.best_period

        #Read best period out
        return [max_periods_orig, best_periods_orig]

    #Select the Mag data and apply random gaussian
    Mags_mod = Mag_orgs + Errors * np.random.randn(Mag_orgs.size)
    #Build L-S model
    model = LombScargleFast().fit(Times, Mags_mod, Errors)

    #Apply the LS mechanism to the data
    periods, power = model.periodogram_auto(nyquist_factor=0.1, oversampling=50)
    max_period = periods[np.argmax(power)]
    file_path_Light_Curve = file_loc + '/Light_Curve_' + Object + '_' + str(number) + '.jpg'
    file_path_LS_model = file_loc + '/LSF_' + Object + '_' + str(number) + '.jpg'

    #Plotting the Light Curve using seaborn style
    if Output_figs == 1:
        fig, ax = plt.subplots()
        ax.errorbar(Times, Mags_mod, Errors, fmt='.k', color="r", ecolor='gray', label='Light Curve')
        ax.set(xlabel='Time (MJD)',
                ylabel='Magitude',
                title=Object + ' Light Curve with Gaussian Noise iteration ' + str(number))
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
        plt.clf()
    #Calculating the best period of the model, condition on quality of LSF model
    model.optimizer.period_range=(500, 2000)
    best_period = model.best_period

    return [max_period,best_period]

    #End of Iteration_function

#Calling function calls the iteration function
def calling_function(Data, Iterations, Outputs, Output_loc, Object):

    #Define a period list, to store results
    Period_List = []
    Best_Period_List = []
    Times = Data[['MJD']].as_matrix().ravel()

    #Convert string variable to number variable
    if Outputs == "Y": Output_n = 1
    elif Outputs == "N": Output_n = 0
    else: print("Error in Output_figures variable")

    #Perform the iterations
    for j in range(0, Iterations):
        New_Period = Iteration_function(Data, j, Output_n, Output_loc, Object)
        Period_List.append(New_Period[0])
        Best_Period_List.append(New_Period[1])

    #Call the LS_hist function in order to create the Histogram
    LS_hist(Best_Period_List, Output_loc, Object, Iterations, 'Best')
    return (np.max(Times)-np.min(Times))/LS_hist(Period_List, Output_loc, Object, Iterations, 'Max')

    #End of calling_function

#Historgram function plots a histogram of the created data
def LS_hist(Results, file_loc, Object, Iterations, string):
    #Plot the Historgram appropriately
    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(Results, bins='auto' )
    ax.hist(Results, bins='auto')
    ax.set(xlabel = 'Period',
            ylabel = 'No. of Counts',
            title = Object + string + " Histogram " + str(Iterations) + ' Iterations')

    file_path_hist = file_loc + '\Hist ' + Object + ' ' + string + '.jpg'
    fig.savefig(file_path_hist, bbox_inches='tight')
    plt.clf()
    return(bins[np.argmax(n)])
    #End of LS_hist function
"""
def Graham_Candidate_Data_Grab(User):

    #Define the path to the Graham Periodic Candidates
    if User == 'N':
        Graham_Candidate_File = 'C:/Users/User/Documents/University/Year 4/Project/Graham_Candidate_names.txt'
    elif User == 'C':
        Graham_Candidate_File = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Graham_Candidate_names.txt'

    #Bringing in the Graham Candidate input Data
    Graham_Input_Data = pd.read_table(Graham_Candidate_File, sep=';', header=0)

    #Read in the individual column Data
    Names = Graham_Input_Data[['Name']].as_matrix()
    RA_Values = Graham_Input_Data[['RAJ2000']].as_matrix()
    Dec_Values = Graham_Input_Data[['DEJ2000']].as_matrix()

    #Need to re-shape the maticies to arrrays

    #For loop to remove the white spaces from each element
    for i in range(0, len(Names)):
        Names[i,0] = str.strip(Names[i,0])
        RA_Values[i,0] = str.strip(RA_Values[i,0])
        Dec_Values[i,0] = str.strip(Dec_Values[i,0])
        if 'BZQJ1305-1033' in Names[i,0]:
            Names[i,0] = 'PG 1302-102'

    #Define the chrome driver, with an enforced waiting time
    if User == 'N':
        driver = webdriver.Chrome('C:/Users/User/chromedriver')
    if User == 'C':
        driver = webdriver.Chrome('C:/Users/Christopher/Documents/chromedriver')

    driver.implicitly_wait(10)

    #Opens the output file and inserts the headers
    if User == 'N':
        writing_file = open('C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt', 'w')
    if User == 'C':
        writing_file = open('C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt', 'w')

    writing_file.write('Name,URL\n')

    #For loop for each Periodic candidate object
    for j in range(0, len(Names)):

        #Load the webpage with the RA and DEC inputs
        driver.get("http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi")

        #Selects the RA and DEC input boxes appropriately
        RA_text_box = driver.find_element_by_name('RA')
        Dec_text_box = driver.find_element_by_name('Dec')

        #Inserts the relevant text into each input box
        RA_text_box.send_keys(RA_Values[j,0])
        Dec_text_box.send_keys(Dec_Values[j,0])

        #Click the submit button
        driver.find_element_by_name('.submit').click()

        #Grabs the URL link to the data from the webpage
        link = driver.find_element_by_link_text('download')
        URL = link.get_attribute("href")

        #Writes the name of the object and URL in Comma-Delimited format to file
        writing_file.write(Names[j,0] + ',' + URL + '\n')

    #Closes the webpage
    driver.quit()

    #Closes the output file
    writing_file.close()

    #End of Graham_Candidate_Data_Grab function
"""

#Main function to be called
def main():

    User = 'C'

    #Define the location of the saved path
    if(User == 'C'):
        Data_location = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt'
    if(User == 'N'):
        Data_location = 'C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt'

    #Graham_Candidate_Data_Grab(User)

    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #save each column within the dataframe as a seperate numpy array
    Objects = Raw_Data_input[['Name']].as_matrix()
    URLs = Raw_Data_input[['URL']].as_matrix()

    #Set the number of iterations here
    No_of_Iterations = 2

    #Whether to output Light curve and L-S figures or not (Y/N)
    Output_figures = 'N'

    #Output_location
    if (User == 'N'):
        Output_locations = 'C:/Users/User/Documents/University/Year 4/Project/Iteration_figures_11-10-2016'
    elif (User == 'C'):
        Output_locations = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Iteration_figures_11-10-2016'

    f = open("No. of periods.txt", 'w')
    #Call the calling function for each object in the data
    for k in range(0, len(Raw_Data_input)):
        Object = Objects.item(k,0)
        URL = URLs.item(k,0)
        URL_editted = URL.strip()
        Raw_Data = pd.read_csv(URL_editted, sep=',', header=0)

        #Calling function runs program
        pnumber = calling_function(Raw_Data, No_of_Iterations, Output_figures, Output_locations, Object)

        f.write(Object + ' ' + str(pnumber) + '\n')

    f.close()

#End of main function
main()

#End of File
