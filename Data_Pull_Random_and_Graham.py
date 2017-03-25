""" --------------------------------------------------------------------------
Description: This program is designed to pull the data for random objects from
the Million Quasar Database and also to select the data for the 111 Graham
candidates. This data will be saved in seperate text files in the
Julia_Working_Directory, so that a CARMA process may be run on them.
It processes all the Random object data, forming the light curves if desired.
It implements a selection criteria, firstly, clipping data points outside a
5 sigma range for all the objects, secondly, it checks that the error bars are
comparable (within 1 sigma) of the Graham object error bars. It also only
chooses objects with more than 50 data-points and uses a 3 arc second radius
when locating the data. Finally, it also checks to see if the data is scewed
to one side (in terms of time measurements) through a defined leniancy value,
only allowing those whose datapoints are adequately spaced to be examined.


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
import numpy as np
import matplotlib.pyplot as plt

def Random_Obj_Data_Grab(User, Iterations_start, Iterations_end):

    #Define the file path to the random object text file (output of Full_Data-Pull.py)
    Random_Obj_Data_path = 'C:/Users/User/Documents/University/Year 4/Project/Upload_File_Final.txt'

    #Read the data into a pandas dataframe
    Rand_Obj_Data = pd.read_table(Random_Obj_Data_path, sep=' ', header=None)

    #Name the columns with the appropriate variables
    Rand_Obj_Data.columns = ['Name', 'RA', 'Dec']

    #Set RA and Dec values to limit of 7 dp's
    pd.to_numeric(Rand_Obj_Data.RA, errors='coerce')
    Rand_Obj_Data['RA'] = Rand_Obj_Data['RA'].round(7)
    pd.to_numeric(Rand_Obj_Data.Dec, errors='coerce')
    Rand_Obj_Data['Dec'] = Rand_Obj_Data['Dec'].round(7)

    for i in range(Iterations_start, Iterations_end):

        #Select the next 100 objects
        Iteration_value_start = i*100
        Iteration_value_end = ((i+1)*100)-1
        Upload_Data = Rand_Obj_Data.loc[Iteration_value_start:Iteration_value_end,]

        #Write the results to a text file
        if User == 'N':
            Upload_Data.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')
        elif User == 'C':
            Upload_Data.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')

def Charisi_Obj_Data_Grab(User):

    Path_To_Location_Data = "C:/Users/User/Documents/University/Year 4/Project/Charisi_Data_Editted.txt"
    Charisi_Obj_Data = pd.read_table(Path_To_Location_Data, sep=' ', header=0, dtype=str)
    Charisi_Obj_Data["Name"] = Charisi_Obj_Data["Name1"] + Charisi_Obj_Data["Name2"]
    Charisi_Obj_Data = Charisi_Obj_Data[["Name", "Ra", "Dec"]]
    #pd.to_numeric(Charisi_Obj_Data.Dec, errors='coerce')

    #Charisi_Obj_Data = Charisi_Obj_Data[Charisi_Obj_Data["Dec"] < 70]
    #print(Charisi_Obj_Data)

    Charisi_Obj_Data.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Charisi_Data_Upload.txt', header=None, index=None, sep=' ', mode='w')


def Process_Random_Obj_Data(User, Iterations_start, Iterations_end, leniancy_value, sigma_level, Form_Light_Curves, Upper_G_Error_Level):

    #Define a list
    Data_Names = {}
    frames = []

    #Use a for loop to read the data in
    for i in range(Iterations_start, Iterations_end):

        #Define the path to the data
        if User == 'N':
            Data_path = 'C:/Users/User/Documents/University/Year 4/Project/Randoms_Output_Data_' + str(i) + '.txt'
        elif User == 'C':
            Data_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Randoms_Output_Data_' + str(i) + '.txt'

        #Define a naming convention and append this to the frames list
        Data_Names["Data{0}".format(i)] = pd.read_table(Data_path, sep=',', header=0)
        DataFrame_Name = "Data" + str(i)
        frames.append(Data_Names[DataFrame_Name])

    #Cocatenate the outputted dataframes
    Total_Data = pd.concat(frames)

    #Reset the index
    Total_Data = Total_Data.reset_index()

    #Make a copy of the IDs
    IDs = pd.DataFrame(Total_Data.InputID)

    #Length of IDs DF indicates how many objects we have
    IDs = IDs.drop_duplicates()
    IDs = IDs.reset_index()
    IDs.columns = ['No_Datapoints', 'ID']

    #Find the number of observations for the objects
    IDs['Difference'] = IDs['No_Datapoints'].shift(-1) - IDs['No_Datapoints']
    Final_Row_Value = int(len(IDs) - 1)
    IDs.Difference.loc[Final_Row_Value] = len(Total_Data) - IDs.No_Datapoints.loc[Final_Row_Value]

    #Only choose objects with 50 or more observations
    IDs_final = IDs.ix[IDs['Difference'] >= 50.0]
    IDs_final = IDs_final.reset_index()
    Observation_No_Rejects = len(IDs) - len(IDs_final)

    #Select specific columns
    Total_Data_selected = Total_Data[['InputID', 'Mag', 'Magerr', 'MJD']]

    #Define an empty DataFrame to store the final IDs that passed selection criteria
    IDs_only = pd.Series(np.nan)

    #Set the counters to 0
    Reject_Error_Counter = 0
    Leniancy_Value_Rejection = 0
    Leniancy_Value_Rejection_Final = 0
    Sigma_Clipped_Counter = 0
    No_of_Points_Removed = 0
    Total_No_Points = 0

    for j in range(len(IDs_final)):

        #Find the data for the objects with more than 50 observations
        Obj_Name = IDs_final.iloc[j]['ID']
        Object_Data = Total_Data_selected.ix[Total_Data_selected['InputID'] == Obj_Name]

        #Force the MJD values to based from a specific point (i.e. 1st date of observation)
        pd.to_numeric(Object_Data.MJD, errors='coerce')
        Object_Data['MJD'] = Object_Data['MJD'] - Object_Data['MJD'].min()
        Object_Data['MJD'] = Object_Data['MJD'].round(5)

        #Check that the errors are a suitable size
        Mean_Err_Value = Object_Data['Magerr'].mean()
        if Mean_Err_Value > Upper_G_Error_Level:
            Error_Selection_Passed = 'N'
            Reject_Error_Counter = Reject_Error_Counter + 1
        else:
            Error_Selection_Passed = 'Y'

        #Remove the objects with data scewed to one-side
        Scewed_Data_Value = (Object_Data['MJD'].median()) / Object_Data['MJD'].max()
        Scewed_Data_Value = round(Scewed_Data_Value, 2)

        if Scewed_Data_Value > (0.5-leniancy_value) and Scewed_Data_Value < (0.5+leniancy_value) and Error_Selection_Passed == 'Y':

            #Append the Object name onto the empty Pandas DataFrame
            Inserting_value = pd.Series([Obj_Name])
            IDs_only = IDs_only.append(Inserting_value, ignore_index=True)

            #Sort in ascending MJD, then reset index
            Object_Data.sort_values('MJD', ascending=True, inplace=True)
            Object_Data.reset_index(inplace=True)
            del Object_Data['index']

            #Sigma clipping inserted
            No_DataPoints_Before_Clipping = len(Object_Data.index)
            Object_Data = Object_Data[abs(Object_Data['Mag'] - Object_Data['Mag'].median()) < sigma_level*Object_Data['Magerr']]
            No_DataPoints_After_Clipping = len(Object_Data.index)

            #Add to counter for sigma-clipped objects
            if No_DataPoints_Before_Clipping != No_DataPoints_After_Clipping:
                Sigma_Clipped_Counter = Sigma_Clipped_Counter + 1
                No_of_Points_Removed = No_of_Points_Removed + (No_DataPoints_Before_Clipping - No_DataPoints_After_Clipping)
                Total_No_Points = Total_No_Points + No_DataPoints_Before_Clipping

            #Select the data needed for the NS process
            Object_Data_selected = Object_Data[['MJD', 'Mag', 'Magerr']]

            #Write the results to a text file
            if User == 'N':
                Object_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/Data_' + Obj_Name + '.txt', header=None, index=None, sep=' ', mode='w')
            elif User == 'C':
                Object_Data_selected.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')

            if Form_Light_Curves == 'Y':

                #Set the x,y,yerr values
                Times = Object_Data[['MJD']].as_matrix().ravel()
                Errors = Object_Data[['Magerr']].as_matrix().ravel()
                Mag_orgs = Object_Data[['Mag']].as_matrix().ravel()

                #Make the plot and format
                plt.figure()
                plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
                plt.xlabel('Time(MJD)')
                plt.ylabel('Magitude')
                Title = 'Light Curve for ' + Obj_Name
                plt.title(Title)
                plt.gca().invert_yaxis()

                #Define the save file path
                if User == 'N':
                    file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/' + Obj_Name + '.jpg'
                elif User == 'C':
                    file_path_Light_Curve = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/' + Obj_Name + '.jpg'

                #Save the plot
                plt.savefig(file_path_Light_Curve, bbox_inches='tight')

                #Close all open plots
                plt.close('all')

        else:
            Leniancy_Value_Rejection = Leniancy_Value_Rejection + 1

    if j == len(IDs_final) - 1:
        Leniancy_Value_Rejection_Final = Leniancy_Value_Rejection - Reject_Error_Counter

    #Drop the inserted NaN value
    IDs_only = IDs_only.dropna()

    #Create a text file containing the IDs only for the Random objects
    if User == 'N':
        IDs_only.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt', header=None, index=None, sep=' ', mode='w')
    elif User == 'C':
        IDs_only.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Randoms_ID_Values.txt', header=None, index=None, sep=' ', mode='w')

    #Formualte the relevant values
    No_Objects_Submitted = (Iterations_end - Iterations_start) * 100
    No_Objects_Examined = len(IDs)
    No_Objects_Survived = No_Objects_Examined - Observation_No_Rejects - Reject_Error_Counter - Leniancy_Value_Rejection_Final
    Survival_Percentage = (No_Objects_Survived / No_Objects_Examined) * 100
    Survival_Percentage = round(Survival_Percentage, 2)

    #Print the values to the output screen
    print("Number of Objects submitted is:", No_Objects_Submitted)
    print("Number of total objects examined is:", No_Objects_Examined)
    print("Number of Objects with less than 50 observations:", Observation_No_Rejects)
    print("Number of Objects with error bars outside 1-sigma range:", Reject_Error_Counter)
    print("Number of Objects with a scewed data value outside the leniancy_value:", Leniancy_Value_Rejection_Final)
    print("Number of Objects with a data point removed from final sample:", Sigma_Clipped_Counter)
    print("Number of Datapoints before clipping:", Total_No_Points)
    print("Number of Data points removed:", No_of_Points_Removed)
    print("Number of Objects surviving selection criteria is:", No_Objects_Survived)
    print("This gives a Survival_Percentage of:", Survival_Percentage)

def Process_Charisi_Obj_Data(User, leniancy_value, sigma_level, Form_Light_Curves, Upper_G_Error_Level):

    #Define the path to the data
    if User == 'N':
        Data_path = 'C:/Users/User/Documents/University/Year 4/Project/Fresh_Data.txt'
    elif User == 'C':
        Data_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Randoms_Output_Data_' + str(i) + '.txt'

    Total_Data = pd.read_table(Data_path, sep='\t', header=0)

    #Make a copy of the IDs
    IDs = pd.DataFrame(Total_Data.InputID)

    #Length of IDs DF indicates how many objects we have
    IDs = IDs.drop_duplicates()
    IDs = IDs.reset_index()
    IDs.columns = ['No_Datapoints', 'ID']

    #Find the number of observations for the objects
    IDs['Difference'] = IDs['No_Datapoints'].shift(-1) - IDs['No_Datapoints']
    Final_Row_Value = int(len(IDs) - 1)
    IDs.Difference.loc[Final_Row_Value] = len(Total_Data) - IDs.No_Datapoints.loc[Final_Row_Value]

    #Only choose objects with 50 or more observations
    IDs_final = IDs.ix[IDs['Difference'] >= 50.0]
    IDs_final = IDs_final.reset_index()
    Observation_No_Rejects = len(IDs) - len(IDs_final)

    #Select specific columns
    Total_Data_selected = Total_Data[['InputID', 'Mag', 'Magerr', 'MJD']]

    #Define an empty DataFrame to store the final IDs that passed selection criteria
    IDs_only = pd.Series(np.nan)

    #Set the counters to 0
    Reject_Error_Counter = 0
    Leniancy_Value_Rejection = 0
    Leniancy_Value_Rejection_Final = 0
    Sigma_Clipped_Counter = 0
    No_of_Points_Removed = 0
    Total_No_Points = 0

    for j in range(len(IDs_final)):

        #Find the data for the objects with more than 50 observations
        Obj_Name = IDs_final.iloc[j]['ID']
        Object_Data = Total_Data_selected.ix[Total_Data_selected['InputID'] == Obj_Name]

        #Force the MJD values to based from a specific point (i.e. 1st date of observation)
        pd.to_numeric(Object_Data.MJD, errors='coerce')
        Object_Data['MJD'] = Object_Data['MJD'] - Object_Data['MJD'].min()
        Object_Data['MJD'] = Object_Data['MJD'].round(5)

        #Check that the errors are a suitable size
        Mean_Err_Value = Object_Data['Magerr'].mean()
        if Mean_Err_Value > Upper_G_Error_Level:
            Error_Selection_Passed = 'N'
            Reject_Error_Counter = Reject_Error_Counter + 1
        else:
            Error_Selection_Passed = 'Y'

        #Remove the objects with data scewed to one-side
        Scewed_Data_Value = (Object_Data['MJD'].median()) / Object_Data['MJD'].max()
        Scewed_Data_Value = round(Scewed_Data_Value, 2)

        if Scewed_Data_Value > (0.5-leniancy_value) and Scewed_Data_Value < (0.5+leniancy_value) and Error_Selection_Passed == 'Y':

            #Append the Object name onto the empty Pandas DataFrame
            Inserting_value = pd.Series([Obj_Name])
            IDs_only = IDs_only.append(Inserting_value, ignore_index=True)

            #Sort in ascending MJD, then reset index
            Object_Data.sort_values('MJD', ascending=True, inplace=True)
            Object_Data.reset_index(inplace=True)
            del Object_Data['index']

            #Sigma clipping inserted
            No_DataPoints_Before_Clipping = len(Object_Data.index)
            Object_Data = Object_Data[abs(Object_Data['Mag'] - Object_Data['Mag'].median()) < sigma_level*Object_Data['Magerr']]
            No_DataPoints_After_Clipping = len(Object_Data.index)

            #Add to counter for sigma-clipped objects
            if No_DataPoints_Before_Clipping != No_DataPoints_After_Clipping:
                Sigma_Clipped_Counter = Sigma_Clipped_Counter + 1
                No_of_Points_Removed = No_of_Points_Removed + (No_DataPoints_Before_Clipping - No_DataPoints_After_Clipping)
                Total_No_Points = Total_No_Points + No_DataPoints_Before_Clipping

            #Select the data needed for the NS process
            Object_Data_selected = Object_Data[['MJD', 'Mag', 'Magerr']]

            #Write the results to a text file
            if User == 'N':
                Object_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Fresh_Data_' + Obj_Name + '.txt', header=None, index=None, sep=' ', mode='w')
            elif User == 'C':
                Object_Data_selected.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')

            if Form_Light_Curves == 'Y':

                #Set the x,y,yerr values
                Times = Object_Data[['MJD']].as_matrix().ravel()
                Errors = Object_Data[['Magerr']].as_matrix().ravel()
                Mag_orgs = Object_Data[['Mag']].as_matrix().ravel()

                #Make the plot and format
                plt.figure()
                plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
                plt.xlabel('Time(MJD)')
                plt.ylabel('Magitude')
                Title = 'Light Curve for ' + Obj_Name
                plt.title(Title)
                plt.gca().invert_yaxis()

                #Define the save file path
                if User == 'N':
                    file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Charisi_Object_Light_Curves/' + Obj_Name + '.jpg'
                elif User == 'C':
                    file_path_Light_Curve = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/Julia_Working_Directory/Random_Object_Light_Curves/' + Obj_Name + '.jpg'

                #Save the plot
                plt.savefig(file_path_Light_Curve, bbox_inches='tight')

                #Close all open plots
                plt.close('all')

        else:
            Leniancy_Value_Rejection = Leniancy_Value_Rejection + 1

    if j == len(IDs_final) - 1:
        Leniancy_Value_Rejection_Final = Leniancy_Value_Rejection - Reject_Error_Counter

    #Drop the inserted NaN value
    IDs_only = IDs_only.dropna()

    #Create a text file containing the IDs only for the Random objects
    if User == 'N':
        IDs_only.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Fresh_ID_Values.txt', header=None, index=None, sep=' ', mode='w')
    elif User == 'C':
        IDs_only.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Randoms_ID_Values.txt', header=None, index=None, sep=' ', mode='w')

    #Formualte the relevant values
    No_Objects_Submitted = 48
    No_Objects_Examined = len(IDs)
    No_Objects_Survived = No_Objects_Examined - Observation_No_Rejects - Reject_Error_Counter - Leniancy_Value_Rejection_Final
    Survival_Percentage = (No_Objects_Survived / No_Objects_Examined) * 100
    Survival_Percentage = round(Survival_Percentage, 2)

    #Print the values to the output screen
    print("Number of Objects submitted is:", No_Objects_Submitted)
    print("Number of total objects examined is:", No_Objects_Examined)
    print("Number of Objects with less than 50 observations:", Observation_No_Rejects)
    print("Number of Objects with error bars outside 1-sigma range:", Reject_Error_Counter)
    print("Number of Objects with a scewed data value outside the leniancy_value:", Leniancy_Value_Rejection_Final)
    print("Number of Objects with a data point removed from final sample:", Sigma_Clipped_Counter)
    print("Number of Datapoints before clipping:", Total_No_Points)
    print("Number of Data points removed:", No_of_Points_Removed)
    print("Number of Objects surviving selection criteria is:", No_Objects_Survived)
    print("This gives a Survival_Percentage of:", Survival_Percentage)

def Graham_candidate_Data(User, leniancy_value, sigma_level, Form_Light_Curves, Upper_G_Error_Level):

    #Define the location of the saved path
    if(User == 'C'):
        Data_location = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt'
    if(User == 'N'):
        Data_location = 'C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt'

    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #Define an empty DataFrame to store the final IDs that passed selection criteria
    Graham_IDs = pd.Series(np.nan)

    #Set the counters to zero
    No_Objects_Clipped = 0
    No_Datapoints_Clipped = 0
    Reject_Error_Counter = 0
    Total_No_Points = 0

    for i in range(len(Raw_Data_input)):

        #Read the Object Name
        Object = Raw_Data_input.Name.iloc[i]

        #Read the raw data from the URL
        Raw_Data = pd.read_csv(Raw_Data_input.URL.iloc[i], sep=',', header=0)

        #Select the data, convert MJD values to a base from 0.0 and set to 5 decimal places
        Raw_Data_selected = Raw_Data[['MJD', 'Mag', 'Magerr']]
        pd.to_numeric(Raw_Data_selected.MJD, errors='coerce')
        Raw_Data_selected['MJD'] = Raw_Data_selected['MJD'] - Raw_Data_selected['MJD'].min()
        Raw_Data_selected['MJD'] = Raw_Data_selected['MJD'].round(5)

        #Perform the sigma clipping
        Before_clipping_Data = len(Raw_Data_selected.index)
        Raw_Data_selected = Raw_Data_selected[abs(Raw_Data_selected['Mag'] - Raw_Data_selected['Mag'].median()) < sigma_level*Raw_Data_selected['Magerr']]
        After_clipping_Data = len(Raw_Data_selected.index)

        #Check for the size of the error bars
        Mean_Err_Value = Raw_Data_selected['Magerr'].mean()
        if Mean_Err_Value > Upper_G_Error_Level:
            Error_Selection_Passed = 'N'
            Reject_Error_Counter = Reject_Error_Counter + 1
        else:
            Error_Selection_Passed = 'Y'

        if Before_clipping_Data != After_clipping_Data:
            No_Objects_Clipped = No_Objects_Clipped + 1
            No_Datapoints_Clipped = No_Datapoints_Clipped + (Before_clipping_Data - After_clipping_Data)
            Total_No_Points = Total_No_Points + Before_clipping_Data

        #Remove the objects with data scewed to one-side
        Scewed_Data_Value = (Raw_Data_selected['MJD'].median()) / Raw_Data_selected['MJD'].max()
        Scewed_Data_Value = round(Scewed_Data_Value, 2)

        if Scewed_Data_Value > (0.5-leniancy_value) and Scewed_Data_Value < (0.5+leniancy_value):

            #Sort the data and reset the index
            Raw_Data_selected.sort_values('MJD', ascending=True, inplace=True)
            Raw_Data_selected.reset_index(inplace=True)
            del Raw_Data_selected['index']

            #Append the Object name onto the empty Pandas DataFrame
            Inserting_value = pd.Series([Object])
            Graham_IDs = Graham_IDs.append(Inserting_value, ignore_index=True)

            #Write the results to a text file
            if User == 'N':
                Raw_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + Object + '.txt', header=None, index=None, sep=' ', mode='w')
            elif User == 'C':
                Raw_Data_selected.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Data_' + Object + '.txt', header=None, index=None, sep=' ', mode='w')

            if Form_Light_Curves == 'Y':

                #Set the x,y,yerr values
                Times = Raw_Data_selected[['MJD']].as_matrix().ravel()
                Errors = Raw_Data_selected[['Magerr']].as_matrix().ravel()
                Mag_orgs = Raw_Data_selected[['Mag']].as_matrix().ravel()

                #Make the plot and format
                plt.figure()
                plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
                plt.xlabel('Time(MJD)')
                plt.ylabel('Magitude')
                Title = 'Light Curve for ' + Object
                plt.title(Title)
                plt.gca().invert_yaxis()

                #Define the save file path
                if User == 'N':
                    file_path_Light_Curve = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Object_Light_Curves/' + Object + '.jpg'
                elif User == 'C':
                    file_path_Light_Curve = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/Julia_Working_Directory/Graham_Object_Light_Curves/' + Object + '.jpg'

                #Save the plot
                plt.savefig(file_path_Light_Curve, bbox_inches='tight')

                #Close all open plots
                plt.close('all')

    #Drop the inserted NaN value
    Graham_IDs = Graham_IDs.dropna()

    #Write the results to a text file
    #if User == 'N':
        #Graham_IDs.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt', header=None, index=None, sep=' ', mode='w')
    #elif User == 'C':
        #Graham_IDs.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Graham_ID_Values.txt', header=None, index=None, sep=' ', mode='w')

    #Print the relevant values to the screen
    print("Number of objects with a datapoint removed:", No_Objects_Clipped)
    print("Number of datapoints removed:", No_Datapoints_Clipped)
    print("Number of objects removed due to large error bars:", Reject_Error_Counter)
    print("Number of Datapoints before clipping:", Total_No_Points)

def main():

    #Define the User
    User = 'N'

    #To remove assignment warnings
    pd.options.mode.chained_assignment = None

    #Define a variable that tells us how many 100 samples of the data are needed
    Iterations_start = 0
    Iterations_end = 50

    #Set the Data leniancy value for the scewing of times to one-side
    Scewed_Data_leniancy_value = 0.2

    #Set the sigma clipping level
    Sigma_Level = 5

    #Set the size of the error bars
    G_Mean = 0.11747274711043741
    G_Std_Dev = 0.03849447039992238
    No_Of_Std_Devs = 1

    Upper_G_Error = G_Mean + (No_Of_Std_Devs*G_Std_Dev)

    #Whether to form the light curves or not
    Form_Light_Curves = 'N'

    #Call the function which grabs the data for the random objects
    #Random_Obj_Data_Grab(User, Iterations_start, Iterations_end)

    #Call the function which processes the data for each object
    #Process_Random_Obj_Data(User, Iterations_start, Iterations_end, Scewed_Data_leniancy_value, Sigma_Level, Form_Light_Curves, Upper_G_Error)

    #Call the function which processes the data for the Graham Candidates
    Graham_candidate_Data(User, Scewed_Data_leniancy_value, Sigma_Level, Form_Light_Curves, Upper_G_Error)

    #Charisi_Obj_Data_Grab(User)
    #Process_Charisi_Obj_Data(User, Scewed_Data_leniancy_value, Sigma_Level, Form_Light_Curves, Upper_G_Error)
main()
