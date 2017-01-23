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


if User == 'N':
    Graham_IDs_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt'
elif User == 'C':
    Graham_IDs_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_ID_Values.txt'

#Load in all the Graham object data
All_Graham_Objects = pd.read_table(Graham_IDs_path, sep=' ', header=None)

for j in range(len(All_Graham_Objects)):
    if 'PG 1302' in All_Graham_Objects.iloc[j,0]:

        #Path to data
        if User == 'N':
            Path_to_Data = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_' + All_Graham_Objects.iloc[j,0] +'.txt'
        elif User == 'C':
            Path_to_Data = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Grahams_Data/Data_' + All_Graham_Objects.iloc[j,0] +'.txt'


        #Read the data
        Data_for_obj = pd.read_table(Path_to_Data, sep=' ', header=None)
        Data_for_obj.columns = ['MJD', 'Mag', 'Magerr']
        Times = Data_for_obj[['MJD']].as_matrix().ravel()
        Period = (Times.max()-Times.min())/3
        Mag = np.sin(2*np.pi*Times/Period)
        for k in range(20):
            Magerr = np.ones_like(Times)*(k+1)/10
            sigma = str((k+1)/20)

            if User == 'N':
                sinusoids_path = 'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Sinusoids_Data/Data_Magerr_' + sigma + '.txt'
            elif User == 'C':
                sinusoids_path = 'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Sinusoids_Data/Data_Magerr_' + sigma + '.txt'

            d = {'MJD' : pd.Series(Times),
                'Mag': pd.Series(Period),
                'Magerr': pd.Series(Magerr)}
            df = pd.DataFrame(d)
            df.to_csv(sinusoids_path, index=False, sep = ' ')
