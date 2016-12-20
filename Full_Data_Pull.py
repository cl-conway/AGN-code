""" --------------------------------------------------------------------------
Description: This program edits the Million_Quasar_Database into a form that
is appropriate for the data from the CRTS to be collected. It is important to
note that a webdriver technique was explored but such an automated process
is prevented by the safeguards introduced by the CRTS team.

Inputs: Raw Data (e.g. for PG1302-102)
C:/Users/User/Documents/University/Year 4/Project/Million_Quasar_Database_Test.txt
or
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Million_Quasar_Database_Test.txt
[Internet connection Required to run]

Exports to:

Date Created: 09/11/2016
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
import math

def Full_data_Grab(User, Radius_Value):

    #Define the path to the Graham Periodic Candidates
    if User == 'N':
        Quasar_Database_Path = 'C:/Users/User/Documents/University/Year 4/Project/Million_Quasar_Database.txt'
    elif User == 'C':
        Quasar_Database_Path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Million_Quasar_Database.txt'

    #Specify the size of the columns from million quasars database
    delim = (11,12,27,5,5,5,4,2,2,7,7,7,4,23,23,23,23)

    #Read in the data, convert to pandas dataframe and name columns
    Test_Data = np.genfromtxt(Quasar_Database_Path, dtype=str, autostrip=True, delimiter=delim, usecols=[0,1,2,3,12])
    Upload_Data = pd.DataFrame(Test_Data)
    Upload_Data.columns = ['RA', 'Dec', 'Name', 'Descrip', 'Qpct']
    #Upload_Data['Radius'] = Radius_Value

    #Strip all the strings for the columns in the dataset
    Upload_Data.Name.str.strip()
    Upload_Data.RA.str.strip()
    Upload_Data.Dec.str.strip()
    Upload_Data.Descrip.str.strip()
    Upload_Data.Qpct.str.strip()

    #Remove the spaces in the names column
    Upload_Data['Name'] = Upload_Data['Name'].str.replace(' ', '')

    #Replace the empty string in Qpct column with NaN and drop these rows
    Upload_Data['Qpct'].replace('', np.nan, inplace=True)
    Upload_Data.dropna(subset=['Qpct'], inplace=True)

    #Only select entries with Qusar percentage greater than 95
    Upload_Data = Upload_Data[Upload_Data.Qpct >= '95']
    Upload_Data = Upload_Data.reset_index(drop=True)

    #Select only the columns to send to CRTS team
    Upload_Data_mod = Upload_Data[['Name', 'RA', 'Dec']]

    #For loop to split the data into chuncks of 100
    #for j in range(math.ceil(len(Upload_Data)/100)):

        #Split the data into chuncks of 100
        #Final_Quasar_Data = Upload_Data_mod[(j*100):((j+1)*100)-1]

    Final_Quasar_Data = Upload_Data_mod

    #Write the Data Upload dataframe to a text document, space delimitted
    if User == 'N':
        Final_Quasar_Data.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Upload_File_Final.txt', header=None, index=None, sep=' ', mode='w')
    elif User == 'C':
        Final_Quasar_Data.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_File_Final.txt', header=None, index=None, sep=' ', mode='w')

    #End of DataGrab function

def main():

    #Define the User
    User = 'N'

    #Define a radius to be used
    Radius = 0.0008333

    #Create a string version of the Radius variable, to pass to function
    Radius_str = str(Radius)

    #Call the full data Grab function
    Full_data_Grab(User, Radius_str)

    #End of main function

main()
