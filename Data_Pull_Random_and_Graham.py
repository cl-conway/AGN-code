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

def Random_Obj_Data_Grab(User, Iterations_needed):

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

    for i in range(Iterations_needed):

        #Select the next 100 objects
        Iteration_value_start = i*100
        Iteration_value_end = ((i+1)*100)-1
        Upload_Data = Rand_Obj_Data.loc[Iteration_value_start:Iteration_value_end,]

        #Write the results to a text file
        if User == 'N':
            Upload_Data.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')
        elif User == 'C':
            Upload_Data.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')

def Process_Random_Obj_Data(User, Iterations_needed):

    #Define a list
    Data_Names = {}
    frames = []

    #Use a for loop to read the data in
    for i in range(Iterations_needed):

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

    #Create a text file containing the IDs only for the Random objects
    IDs_only = IDs_final['ID']
    if User == 'N':
        IDs_only.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt', header=None, index=None, sep=' ', mode='w')
    elif User == 'C':
        IDs_only.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Randoms_ID_Values.txt', header=None, index=None, sep=' ', mode='w')

    #Select specific columns
    Total_Data_selected = Total_Data[['InputID', 'Mag', 'Magerr', 'MJD']]

    for j in range(len(IDs_final)):

        #Find the data for the objects with more than 50 observations
        Obj_Name = IDs_final.iloc[j]['ID']
        Object_Data = Total_Data_selected.ix[Total_Data_selected['InputID'] == Obj_Name]

        #Force the MJD values to based from a specific point (i.e. 1st date of observation)
        pd.to_numeric(Object_Data.MJD, errors='coerce')
        Object_Data['MJD'] = Object_Data['MJD'] - Object_Data['MJD'].min()
        Object_Data['MJD'] = Object_Data['MJD'].round(5)

        #Sort in ascending MJD, then reset index
        Object_Data.sort_values('MJD', ascending=True, inplace=True)
        Object_Data.reset_index(inplace=True)
        del Object_Data['index']

        #Select the data needed for the NS process
        Object_Data_selected = Object_Data[['MJD', 'Mag', 'Magerr']]

        #Write the results to a text file
        if User == 'N':
            Object_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/Data_' + Obj_Name + '.txt', header=None, index=None, sep=' ', mode='w')
        elif User == 'C':
            Object_Data_selected.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Upload_Random_Objects_' + str(i) + '.txt', header=None, index=None, sep=' ', mode='w')

def Graham_candidate_Data(User):

    #Define the location of the saved path
    if(User == 'C'):
        Data_location = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Raw_Data_Inputs.txt'
    if(User == 'N'):
        Data_location = 'C:/Users/User/Documents/University/Year 4/Project/Raw_Data_Inputs.txt'

    #Read in the objects and URLs to be searched
    Raw_Data_input = pd.read_table(Data_location, sep=',', header=0)

    #Select the ID Values only
    Graham_IDs = Raw_Data_input['Name']

    #Write the results to a text file
    if User == 'N':
        Graham_IDs.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt', header=None, index=None, sep=' ', mode='w')
    elif User == 'C':
        Graham_IDs.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Graham_ID_Values.txt', header=None, index=None, sep=' ', mode='w')

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

        #Sort the data and reset the index
        Raw_Data_selected.sort_values('MJD', ascending=True, inplace=True)
        Raw_Data_selected.reset_index(inplace=True)
        del Raw_Data_selected['index']

        #Write the results to a text file
        if User == 'N':
            Raw_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + Object + '.txt', header=None, index=None, sep=' ', mode='w')
        elif User == 'C':
            Raw_Data_selected.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Data_' + Object + '.txt', header=None, index=None, sep=' ', mode='w')

def main():

    #Define the User
    User = 'N'

    #To remove assignment warnings
    pd.options.mode.chained_assignment = None

    #Define a variable that tells us how many 100 samples of the data are needed
    #Iterations_needed = math.ceil(len(Rand_Obj_Data) / 100)
    Iterations_needed = 4

    #Call the function which grabs the data for the random objects
    #Random_Obj_Data_Grab(User, Iterations_needed)

    #Call the function which processes the data for each object
    #Process_Random_Obj_Data(User,Iterations_needed)

    #Call the function which processes the data for the Graham Candidates
    Graham_candidate_Data(User)


main()
