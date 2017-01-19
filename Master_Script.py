""" --------------------------------------------------------------------------
Description: This program is designed to locate the relevant data, open a cmd
window and perform the nested sampling, with p=3; q=2. It will then open another
cmd window and call the julia script which performs the analysis. This is done
for each object, both Random AGN and Graham candidate objects.
Inputs:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/...
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/...
This will then create a Data file:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/state-3-2.dat.txt
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/state-3-2.dat.txt

Date Created: 20/12/2016
Authors: Nicholas Kinsey, Christopher Conway
----------------------------------------------------------------------------"""

#Perform relevant imports
import subprocess
import pandas as pd

def main():

    #Set the User
    User = 'C'

    #Choose whether to anlyse Graham or Random objects
    Analyse_Random = 'N'
    Analyse_Graham = 'Y'

    #Set which objects to examine
    R_iter_start = 5
    R_iter_end = 15
    G_iter_start = 96
    G_iter_end = 106

    Analyse_Random = 'Y'
    Analyse_Graham = 'N'

    #Set which objects to examine
    R_iter_start = 120
    R_iter_end = 130
    G_iter_start = 5
    G_iter_end = 12


    #Set the file paths
    if User == 'N':
        Random_Obj_Name_Path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt'
        Graham_Obj_Name_Path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt'
    elif User =='C':
        Random_Obj_Name_Path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Randoms_ID_Values.txt'
        Graham_Obj_Name_Path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_ID_Values.txt'

    #Read in the object names
    Random_Obj_Names = pd.read_table(Random_Obj_Name_Path, sep=' ', header=None)
    Graham_Obj_Names = pd.read_table(Graham_Obj_Name_Path, sep=' ', header=None)

    if Analyse_Random == 'Y':
        for i in range(R_iter_start, R_iter_end):

            #Read object name, locate relevant data
            Object_Name = Random_Obj_Names.loc[i,0]
            Data_text = 'Data_' + Object_Name + '.txt'

            #Replace spaces in names, to pass obj name to analysis script
            Object_Name_pass = Object_Name.replace(' ', '_')

            #Open command window and perform nested sampling. Re-open command and run julia analysis script
            if User == 'N':
                subprocess.call(r'julia C:/Users/User/.julia/v0.5/CARMA/bin/run_carma.jl ' + Data_text + ' 3 2' , cwd=r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data')
                subprocess.call(r'julia Randoms_Data/Run_NS_Randoms_analysis.jl ' + Object_Name_pass, cwd=r'C:/Users/User/Documents/University/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_ID_Values.txt')
            elif User =='C':
                subprocess.call(r'julia C:/Users/User/.julia/v0.5/CARMA/bin/run_carma.jl ' + Data_text + ' 3 2' , cwd=r'C:/Users/Christopher/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data')
                subprocess.call(r'julia Randoms_Data/Run_NS_Randoms_analysis.jl ' + Object_Name_pass, cwd=r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_ID_Values.txt')

            #Print an output for each object completion
            print('Iteration completed for Random Object', i+1)

    if Analyse_Graham == 'Y':
        for i in range(G_iter_start, G_iter_end):

            #Read object name, locate relevant data
            Object_Name = Graham_Obj_Names.loc[i,0]
            Data_text = 'Data_' + Object_Name + '.txt'

            #Replace spaces in names, to pass obj name to analysis script
            Object_Name_pass = Object_Name.replace(' ', '_')

            #Open command window and perform nested sampling. Re-open command and run julia analysis script
            if User == 'N':
                subprocess.call(r'julia C:/Users/User/.julia/v0.5/CARMA/bin/run_carma.jl "' + Data_text + '" 3 2' , cwd=r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data')
                subprocess.call(r'julia Grahams_Data/Run_NS_Grahams_analysis.jl ' + Object_Name_pass, cwd=r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory')
            elif User == 'C':
                subprocess.call(r'julia C:/Users/Christopher/.julia/v0.5/CARMA/bin/run_carma.jl "' + Data_text + '" 3 2' , cwd=r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Grahams_Data')
                subprocess.call(r'julia Grahams_Data/Run_NS_Grahams_analysis.jl ' + Object_Name_pass, cwd=r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory')

            #Print an output for each object completion
            print('Iteration completed for Graham Object', i+1)

    print("***Tasks Completed***")

main()
