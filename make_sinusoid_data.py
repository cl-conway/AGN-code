"""
Description:
File to produce data files of sinusoids that have the same time stamps of PG1302-102
with different error bars. This is to be used be julia CARMA in order to test whether
periodic can be found.
"""

import math
import pandas as pd
import numpy as np

User= 'C'
using_fake_errors= False

if User == 'N':
    Graham_IDs_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Periods_Medians.txt'
elif User == 'C':
    Graham_IDs_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_Periods_Medians.txt'

#Load in all the Graham object data
All_Graham_Objects = pd.read_table(Graham_IDs_path, sep='\t', header=None, index_col=0)
Object_Names_Text =open('Sinusoids_Obj_Names.txt','w')

MCMC_outputs_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/MCMC_output_value.txt'


for j in range(len(All_Graham_Objects)):
        Object_Name =  All_Graham_Objects.iloc[j,0]
        if User == 'N':
            Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Clipped_Data/Clipped_Data_' + Object_Name +'.txt'
        elif User == 'C':
            Path_to_Data = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Grahams_Clipped_Data/Clipped_Data_' + Object_Name +'.txt'

        #Read the data
        MCMC_Output_Path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/MCMC_output_values.txt'

        MCMCdf = pd.read_table(MCMC_Output_Path, sep =';', header=0)
        Data_for_obj = pd.read_table(Path_to_Data, sep=',', header=0)
        Times = Data_for_obj[['MJD']].as_matrix().ravel()
        Mag_True = Data_for_obj[['Mag']].as_matrix().ravel()
        #Period = (Times.max()-Times.min())/3
        Period = All_Graham_Objects.iloc[j,2]
        if 'PG 1302' in Object_Name:
            Amplitude = MCMCdf.Amplitude.iloc[j]
        else:
            Amplitude = (np.max(Mag_True)-np.min(Mag_True))/2
        Median = All_Graham_Objects.iloc[j,1]
        Phase = np.pi
        Mag = Median + Amplitude*np.sin(2*np.pi*Times/Period + Phase)

        #using Fake Errors
        if using_fake_errors == True:
            for k in range(20):
                Magerr = np.ones_like(Times)*(k+1)/10
                sigma = str((k+1)/20)
                Mag = Mag + Magerr * np.random.randn(Mag.size)
                Object_Names_Text.write('Data_Magerr_' + sigma + '\n')

                if User == 'N':
                    sinusoids_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Sinusoids_Data/' + Object_Name + '_Data_Magerr_' + sigma + '.txt'
                elif User == 'C':
                    sinusoids_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Sinusoids_Data/' + Object_Name + '_Data_Magerr_' + sigma + '.txt'

                d = {'MJD' : pd.Series(Times),
                    'Mag': pd.Series(Mag),
                    'Magerr': pd.Series(Magerr)}
                df = pd.DataFrame(d)
                df.to_csv(sinusoids_path, index=False, sep = ' ',header = None)
        else:
            Found_Phase = False
            Magerr_True= Data_for_obj[['Magerr']].as_matrix().ravel()
            for l in range(1,721):
                Phase = l*2*np.pi/720
                if abs(Amplitude*np.sin(Phase)+Median-Mag_True[0]) < np.median(Magerr_True):
                    Found_Phase = True
                    break
                if l == 720 and Found_Phase == False:
                    print(j," No phase for object ", Object_Name)
            Mag = Median + Amplitude*np.sin(2*np.pi*Times/Period + Phase) + Magerr_True*np.random.randn(Mag.size)
            Object_Names_Text.write('Data_'+ Object_Name+ '_Magerr_True \n')

            if User == 'N':
                sinusoids_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Sinusoids_Data/Data_' + Object_Name + '_Magerr_True.txt'
            elif User == 'C':
                sinusoids_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Sinusoids_Data/Data_' + Object_Name + '_Magerr_True.txt'

            d = {'MJD' : pd.Series(Times),
                'Mag': pd.Series(Mag),
                'Magerr': pd.Series(Magerr_True)}
            df = pd.DataFrame(d)
            df.to_csv(sinusoids_path, index=False, sep = ' ',header = None)

Object_Names_Text.close()
