""" --------------------------------------------------------------------------
Description: This program is designed to create the light curves from identified
interesting objects. The data is required for this program to run properly.
It also has a built in feature to output a text file containing a value that
represents the fraction of the median time compared to the maximum time for that
particular object. This value allows us to check whether the data is scewed to
one side of the light curve. Therefore, resulting in a poor nested sampling
output.

Inputs:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Interesting_Randoms.txt
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_ID_Values.txt
[Equivalent for Graham candidates]

Outputs:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Object_Light_Curves/
The figures are to be saved to the above locations
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Random.txt
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Graham.txt

Date Created: 12/01/2017
Authors: Nicholas Kinsey, Christopher Conway
----------------------------------------------------------------------------"""

#Perform relevant imports
import math
import pandas as pd
import matplotlib.pyplot as plt

def Random_Obj_Make_Light_Curve(User, check_scewed_data, form_light_curves):

    #Check value of the check_scewed_data variable
    if check_scewed_data == 'Y':

        #Define path to random objects
        if User == 'N':
            Random_IDs_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt'
        elif User == 'C':
            Random_IDs_path = ''

        #Write a file to check objects with scewed data
        if User == 'N':
            writing_file = open('C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Random.txt', 'w')
        if User == 'C':
            writing_file = open('C:/Users/Christopher/Documents/UNI/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Random.txt', 'w')

        #Write the header of the created file
        writing_file.write('Name;Value\n')

        #Load in all the random object data
        All_Random_Objects = pd.read_table(Random_IDs_path, sep=' ', header=None)

        for j in range(len(All_Random_Objects)):
            #Path to data
            if User == 'N':
                Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/Data_' + All_Random_Objects.iloc[j,0] +'.txt'
            elif User == 'C':
                Path_to_Data = ''

            #Read the data
            Data_for_obj = pd.read_table(Path_to_Data, sep=' ', header=None)
            Data_for_obj.columns = ['MJD', 'Mag', 'Magerr']

            #Form a scewed data value
            Scewed_Data_Value = (Data_for_obj['MJD'].median()) / Data_for_obj['MJD'].max()
            Scewed_Data_Value = round(Scewed_Data_Value, 2)

            #Only output if the data is scewed to one side
            if Scewed_Data_Value < 0.4 or Scewed_Data_Value > 0.6:
                writing_file.write(All_Random_Objects.iloc[j,0] + ';' + str(Scewed_Data_Value) +'\n')

        #Close the writing file after for-loop completes
        writing_file.close()

    #Check value of the form_light_curves variable
    if form_light_curves == 'Y':

        #Define path to interesting random objects
        if User == 'N':
            Interesting_Objects_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Interesting_Randoms.txt'
        elif User == 'C':
            Interesting_Object_path = ''

        #Read the interesting object names in as pandas dataframe
        Interesting_Obj_Names = pd.read_table(Interesting_Objects_path, sep=' ', header=None)

        for i in range(len(Interesting_Obj_Names)):
            #Path to data
            if User == 'N':
                Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/Data_' + Interesting_Obj_Names.iloc[i,0] +'.txt'
            elif User == 'C':
                Path_to_Data = ''

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
            Title = 'Light Curve for ' + Interesting_Obj_Names.iloc[i,0]
            plt.title(Title)
            plt.gca().invert_yaxis()

            #Define the save file path
            if User == 'N':
                file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/' + Interesting_Obj_Names.iloc[i,0] +'.jpg'
            elif User == 'C':
                file_path_Light_Curve = ''

            #Save the plot
            plt.savefig(file_path_Light_Curve, bbox_inches='tight')

            #Close all open plots
            plt.close('all')

def Graham_Obj_Make_Light_Curve(User, check_scewed_data, form_light_curves):

    #Check value of the check_scewed_data variable
    if check_scewed_data == 'Y':

        #Define path to random objects
        if User == 'N':
            Graham_IDs_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt'
        elif User == 'C':
            Graham_IDs_path = ''

        #Write a file to check objects with scewed data
        if User == 'N':
            writing_file = open('C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Graham.txt', 'w')
        if User == 'C':
            writing_file = open('C:/Users/Christopher/Documents/UNI/Year 4/Project/Julia_Working_Directory/Scewed_Data_Objects_Graham.txt', 'w')

        #Write the header of the created file
        writing_file.write('Name;Value\n')

        #Load in all the Graham object data
        All_Graham_Objects = pd.read_table(Graham_IDs_path, sep=' ', header=None)

        for j in range(len(All_Graham_Objects)):
            #Path to data
            if User == 'N':
                Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + All_Graham_Objects.iloc[j,0] +'.txt'
            elif User == 'C':
                Path_to_Data = ''

            #Read the data
            Data_for_obj = pd.read_table(Path_to_Data, sep=' ', header=None)
            Data_for_obj.columns = ['MJD', 'Mag', 'Magerr']

            #Form a scewed data value
            Scewed_Data_Value = (Data_for_obj['MJD'].median()) / Data_for_obj['MJD'].max()
            Scewed_Data_Value = round(Scewed_Data_Value, 2)

            #Only output if the data is scewed to one side
            if Scewed_Data_Value < 0.4 or Scewed_Data_Value > 0.6:
                writing_file.write(All_Graham_Objects.iloc[j,0] + ';' + str(Scewed_Data_Value) +'\n')

        #Close the writing file after for-loop completes
        writing_file.close()

    #Check value of the form_light_curves variable
    if form_light_curves == 'Y':

        #Define path to interesting random objects
        if User == 'N':
            Interesting_Objects_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Interesting_Grahams.txt'
        elif User == 'C':
            Interesting_Object_path = ''

        #Read the interesting object names in as pandas dataframe
        Interesting_Obj_Names = pd.read_table(Interesting_Objects_path, sep=' ', header=None)

        for i in range(len(Interesting_Obj_Names)):
            #Path to data
            if User == 'N':
                Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + Interesting_Obj_Names.iloc[i,0] +'.txt'
            elif User == 'C':
                Path_to_Data = ''

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
            Title = 'Light Curve for ' + Interesting_Obj_Names.iloc[i,0]
            plt.title(Title)
            plt.gca().invert_yaxis()

            #Define the save file path
            if User == 'N':
                file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Object_Light_Curves/' + Interesting_Obj_Names.iloc[i,0] +'.jpg'
            elif User == 'C':
                file_path_Light_Curve = ''

            #Save the plot
            plt.savefig(file_path_Light_Curve, bbox_inches='tight')

            #Close all open plots
            plt.close('all')

def main():

    #Set the User
    User = 'N'

    #Decide what to look at
    Check_Scewed_Data = 'Y'
    Form_Light_Curves = 'Y'

    #Random Objects
    Random_Obj_Make_Light_Curve(User, Check_Scewed_Data, Form_Light_Curves)

    #Graham Objects
    Graham_Obj_Make_Light_Curve(User, Check_Scewed_Data, Form_Light_Curves)

main()
