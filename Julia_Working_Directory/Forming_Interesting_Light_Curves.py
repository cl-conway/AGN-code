""" --------------------------------------------------------------------------
Description: This program is designed to pull the data for random objects from
the Million Quasar Database and also to select the data for the 111 Graham
candidates. This data will be saved in seperate text files in the
Julia_Working_Directory, so that a CARMA process may be run on them.

Inputs: Raw Data (e.g. for PG1302-102)
C:/Users/User/Documents/University/Year 4/Project/Upload_File_Final.txt
C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt
[Internet connection Required to run]

Outputs:
C:/Users/User/Documents/University/Year 4/Project/Randoms_Output_Data_...txt
Above file must be manually uploaded to CRTS site, the data may then be attained
and saved as C:/Users/User/Documents/University/Year 4/Project/Randoms_Output_Data_...txt

This will then create a Data file for each random AGN object at location:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/
And a Data file for each Graham Candidate object at location:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/

Date Created: 09/11/2016
Authors: Nicholas Kinsey, Christopher Conway
----------------------------------------------------------------------------"""

#Perform relevant imports
import math
import pandas as pd
import matplotlib.pyplot as plt

def Random_Obj_Make_Light_Curve(User):

    #Define path to interesting random objects
    if User == 'N':
        Interesting_Objects_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Interesting_Randoms.txt'
    elif User == 'C':
        Interesting_Object_path =

    #Read the interesting object names in as pandas dataframe
    Interesting_Obj_Names = pd.read_table(Interesting_Objects_path, sep=' ', header=None)

    for i in range(len(Interesting_Obj_Names)):
        #Path to data
        if User == 'N':
            Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/Data_' + Interesting_Obj_Names.iloc[i,0] +'.txt'
        elif User == 'C':
            Path_to_Data =

        #Read the data
        Data_for_obj = pd.read_table(Path_to_Data, sep=' ', header=None)
        Data_for_obj.columns = ['MJD', 'Mag', 'Magerr']

        #Set the x,y,yerr values
        Times = Data_for_obj[['MJD']].as_matrix().ravel()
        Errors = Data_for_obj[['Magerr']].as_matrix().ravel()
        Mag_orgs = Data_for_obj[['Mag']].as_matrix().ravel()

        #Make the plot and format
        plt.figure()
        plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
        plt.xlabel('Time(MJD)')
        plt.ylabel('Magitude')
        Title = 'Light Curve for ' + Interesting_Obj_Names.iloc[i,0] + ' (a random object)'
        plt.title(Title)
        plt.gca().invert_yaxis()

        #Define the save file path
        if User == 'N':
            file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/' + Interesting_Obj_Names.iloc[i,0] +'.jpg'
        elif User == 'C':
            file_path_Light_Curve =

        #Save the plot
        plt.savefig(file_path_Light_Curve, bbox_inches='tight')


def Graham_Obj_Make_Light_Curve(User):

    #Define path to interesting random objects
    if User == 'N':
        Interesting_Objects_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Interesting_Grahams.txt'
    elif User == 'C':
        Interesting_Object_path = 'C:/Users/Christopher/Documents/UNI/Year 4/AGN-code/Project/Julia_Working_Directory/Interesting_Grahams.txt'

    #Read the interesting object names in as pandas dataframe
    Interesting_Obj_Names = pd.read_table(Interesting_Objects_path, sep=' ', header=None)

    for i in range(len(Interesting_Obj_Names)):
        #Path to data
        if User == 'N':
            Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + Interesting_Obj_Names.iloc[i,0] +'.txt'
        elif User == 'C':
            Path_to_Data = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Grahams_Data/Data_' + Interesting_Obj_Names.iloc[i,0] +'.txt'

        #Read the data
        Data_for_obj = pd.read_table(Path_to_Data, sep=' ', header=None)
        Data_for_obj.columns = ['MJD', 'Mag', 'Magerr']

        #Set the x,y,yerr values
        Times = Data_for_obj[['MJD']].as_matrix().ravel()
        Errors = Data_for_obj[['Magerr']].as_matrix().ravel()
        Mag_orgs = Data_for_obj[['Mag']].as_matrix().ravel()

        #Make the plot and format
        plt.figure()
        plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
        plt.xlabel('Time(MJD)')
        plt.ylabel('Magitude')
        Title = 'Light Curve for ' + Interesting_Obj_Names.iloc[i,0] + ' (a Graham candidate)'
        plt.title(Title)
        plt.gca().invert_yaxis()

        #Define the save file path
        if User == 'N':
            file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Object_Light_Curves/' + Interesting_Obj_Names.iloc[i,0] +'.jpg'
        elif User == 'C':
            file_path_Light_Curve = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_Object_Light_Curves/' + Interesting_Obj_Names.iloc[i,0] +'.jpg'

        #Save the plot
        plt.savefig(file_path_Light_Curve, bbox_inches='tight')

def main():

    User = 'N'

    Random_Obj_Make_Light_Curve(User)

    Graham_Obj_Make_Light_Curve(User)

main()
