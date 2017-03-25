""" --------------------------------------------------------------------------
Description: This python script is designed to clean and process all the data
outputted in text format by the analysis performed on Tsunami. It does this
by grabbing the object name from the ID_Values text file and looks for the same
objects text file in the results file.

Inputs:
ID_Values for all object types
Results from Tsunami

Outputs:
A cleaned version of the text outputs and writes all the objects for the object
types in singular text files (containing all objects of that type).

Date Created: 14/02/2017
Authors: Nicholas Kinsey, Christopher Conway
----------------------------------------------------------------------------"""

#Perform relevant imports
import pandas as pd
import numpy as np

#Set the User
User = 'N'

#Define Paths
if User == 'N':
    #Path to the ID Values
    Graham_ID_Value_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt"
    Random_ID_Value_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt"

    #Path to the Data
    Graham_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Graham_Text_Output/"
    Random_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Random_Text_Output/"
    Sinusoid_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Sinusoid_Text_Output/"

    #Output Text Locations
    Graham_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graham_Outputs.txt"
    Random_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Random_Outputs.txt"
    Sinusoid_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Sinusoid_Outputs.txt"

elif User == 'C':
    #Path to the ID Values
    Graham_ID_Value_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt"
    Random_ID_Value_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt"

    #Path to the Data
    Graham_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Graham_Text_Output/"
    Random_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Random_Text_Output/"
    Sinusoid_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Sinusoid_Text_Output/"

    #Output Text Locations
    Graham_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graham_Outputs.txt"
    Random_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Random_Outputs.txt"
    Sinusoid_Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Sinusoid_Outputs.txt"

else:
    print("User Error")

#Read the object names into a dataframe
Graham_ID_Values = pd.read_table(Graham_ID_Value_Location, sep=' ', header=None)
Random_ID_Values = pd.read_table(Random_ID_Value_Location, sep=' ', header=None)
Graham_List_Numbers = range(len(Graham_ID_Values))
Sinusoid_ID_Values = Graham_ID_Values.copy()

for i in range(len(Graham_ID_Values)):
    Sinusoid_ID_Values.ix[i,0] = Sinusoid_ID_Values.ix[i,0] + "_Magerr_True"

#Insert the column names
All_Data_columns = ['Object Name', 'Reasonable Frequencies', 'Total Root Number', 'Total oscillatory Root Number', 'Period', 'Upperbound', 'Lowerbound', 'Nu Values', 'Initial Gradient', '~1 Year Gradient', 'Evidence']

#Create empty DataFrames
All_Graham_Data = pd.DataFrame(np.nan, index=[0,1], columns=[0])
All_Random_Data = pd.DataFrame(np.nan, index=[0,1], columns=[0])
All_Sinusoid_Data = pd.DataFrame(np.nan, index=[0,1], columns=[0])

for i in range(len(Graham_ID_Values)):

    #Find the name, grab the data
    Object_Name = Graham_ID_Values.ix[i,0]
    Path_To_Data = Graham_Data_Location + Object_Name + '.txt'
    Objects_Data = pd.read_table(Path_To_Data, sep=';', header=None, skip_blank_lines=True)

    #Form the text and append
    All_Graham_Data = All_Graham_Data.append(Objects_Data, ignore_index=True)

for i in range(len(Random_ID_Values)):

    Object_Name = Random_ID_Values.ix[i,0]
    Path_To_Data = Random_Data_Location + Object_Name + '.txt'
    Objects_Data = pd.read_table(Path_To_Data, sep=';', header=None, skip_blank_lines=True)

    #Form the text and append
    All_Random_Data = All_Random_Data.append(Objects_Data, ignore_index=True)

for i in range(len(Sinusoid_ID_Values)):

    Object_Name = Sinusoid_ID_Values.ix[i,0]
    Path_To_Data = Sinusoid_Data_Location + Object_Name + '.txt'
    Objects_Data = pd.read_table(Path_To_Data, sep=';', header=None, skip_blank_lines=True)

    #Form the text and append
    All_Sinusoid_Data = All_Sinusoid_Data.append(Objects_Data, ignore_index=True)

#Drop the first two rows
All_Graham_Data = All_Graham_Data.drop(All_Graham_Data.index[[0,1]])
All_Random_Data = All_Random_Data.drop(All_Random_Data.index[[0,1]])
All_Sinusoid_Data = All_Sinusoid_Data.drop(All_Sinusoid_Data.index[[0,1]])

#Save the results to a single text file
All_Graham_Data.to_csv(Graham_Output_Location, header=All_Data_columns, index=None, sep=';', mode='w')
All_Random_Data.to_csv(Random_Output_Location, header=All_Data_columns, index=None, sep=';', mode='w')
All_Sinusoid_Data.to_csv(Sinusoid_Output_Location, header=All_Data_columns, index=None, sep=';', mode='w')
