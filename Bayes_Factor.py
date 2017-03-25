"""-----------------------------------------------------------------------------
This script examines the Bayes Factors and calculates the various Percentages
regarding the spread of the outputted objects Bayes factors across the Harold
Jeffreys scale.

Inputs:
The all outputs file from processing the data outputted from Tsunami.

Outputs:
The Bayes Factor percentages table. 
-----------------------------------------------------------------------------"""
import pandas as pd
import numpy as np
import math

User = 'N'

if User == 'C':
    All_Outputs_path = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/All_Outputs.txt"
if User == 'N':
    All_Outputs_path = "C:/Users/User/Documents/University/Year 4/Project/All_Outputs.txt"

All_Objects = pd.read_table(All_Outputs_path, sep=';', header=0, index_col=0)
#Oscillatory_Objects = All_Objects[All_Objects["Oscillatory_Stat"] > 0.5

#Oscillatory_Objects = Oscillatory_Objects.append(All_Objects[All_Objects["Object Name"] == "PG 1302-102"], ignore_index = True)
CAR1_Evidence = np.zeros(len(All_Objects))
#print(Oscillatory_Objects)

for j in range(len(All_Objects)):
    Object_Name = All_Objects["Object Name"].iloc[j]
    Identity = All_Objects["Identifier"].iloc[j]

    if User == 'C':
        Evidence_location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Oscillatory_CAR1_Text_Output/ata_" + Object_Name + ".txt"
    if User == 'N':
        if Identity == 0:
            Evidence_location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Oscillatory_CAR1_Text_Output/Randoms/" + Object_Name + ".txt"
        elif Identity == 1:
            Evidence_location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Oscillatory_CAR1_Text_Output/Grahams/" + Object_Name + ".txt"
        else:
            print("Identity Error")
            exit(1)

    Object = pd.read_table(Evidence_location, sep=';', header=None, index_col=None)
    Object.columns = ["wrong_name","CAR1_Evidence"]
    CAR1_Evidence[j] = Object.CAR1_Evidence.iloc[0].round(2)

CAR1_Data_Frame = pd.DataFrame(data=CAR1_Evidence)
CAR1_Data_Frame.columns = ['CAR1_Evidence']

All_Objects['Log_Bayes_Factor'] = All_Objects.Evidence - CAR1_Data_Frame.CAR1_Evidence
All_Objects.Log_Bayes_Factor = All_Objects.Log_Bayes_Factor.round(2)

All_Objects['Bayes_Factor'] = np.exp(All_Objects.Log_Bayes_Factor)
All_Objects['Oscillatory_Stat'] = round(All_Objects.Oscillatory_Stat, 2)
All_Objects['Bayes_Factor'] = round(All_Objects.Bayes_Factor, 2)

Nice_Table = All_Objects[['Object Name', 'Identifier' ,'Bayes_Factor']]

for i in range(0,len(Nice_Table)):
    if Nice_Table.Identifier.iloc[i] == 1:
        Nice_Table.Identifier.iloc[i] = "G"
    else:
        Nice_Table.Identifier.iloc[i] = "R"

Decisive_Data = All_Objects.ix[All_Objects["Bayes_Factor"] >= 10**2]
Very_Strong_Data = All_Objects.ix[(All_Objects["Bayes_Factor"] >= 10**1.5) & (All_Objects["Bayes_Factor"] < 10**2)]
Strong_Data = All_Objects.ix[(All_Objects["Bayes_Factor"] >= 10**1) & (All_Objects["Bayes_Factor"] < 10**1.5)]
Substantial_Data = All_Objects.ix[(All_Objects["Bayes_Factor"] >= 10**0.5) & (All_Objects["Bayes_Factor"] < 10**1)]
Poor_Data = All_Objects.ix[All_Objects["Bayes_Factor"] < 10**0.5]

G_Decisive = len(Decisive_Data[Decisive_Data["Identifier"] == 1])
G_Very_Strong = len(Very_Strong_Data[Very_Strong_Data["Identifier"] == 1])
G_Strong = len(Strong_Data[Strong_Data["Identifier"] == 1])
G_Substantial = len(Substantial_Data[Substantial_Data["Identifier"] == 1])
G_Poor = len(Poor_Data[Poor_Data["Identifier"] == 1])

R_Decisive = len(Decisive_Data[Decisive_Data["Identifier"] == 0])
R_Very_Strong = len(Very_Strong_Data[Very_Strong_Data["Identifier"] == 0])
R_Strong = len(Strong_Data[Strong_Data["Identifier"] == 0])
R_Substantial = len(Substantial_Data[Substantial_Data["Identifier"] == 0])
R_Poor = len(Poor_Data[Poor_Data["Identifier"] == 0])

All_Objects['BF_Interpretation2'] = np.where(All_Objects['Bayes_Factor'] < 1, 'Unfavoured', 'Favoured')
G_Favoured = len(All_Objects[(All_Objects['BF_Interpretation2'] == "Favoured") & (All_Objects['Identifier'] == 1)])
R_Favoured = len(All_Objects[(All_Objects['BF_Interpretation2'] == "Favoured") & (All_Objects['Identifier'] == 0)])

Final_Table = pd.DataFrame(np.array([G_Decisive, G_Very_Strong, G_Strong, G_Substantial, G_Poor, G_Favoured, R_Decisive, R_Very_Strong, R_Strong, R_Substantial, R_Poor, R_Favoured]))
for j in range(0,5):
    Final_Table.iloc[j] = round((Final_Table.iloc[j] / 111) * 100, 2)
for j in range(6,12):
    Final_Table.iloc[j] = round((Final_Table.iloc[j] / 312) * 100, 2)

Final_Table.columns = ["Percentage_of_Objects"]
Final_Table["Names"] = ["G_Decisive", "G_Very_Strong", "G_Strong", "G_Substantial", "G_Poor", "G_Favoured", "R_Decisive", "R_Very_Strong", "R_Strong", "R_Substantial", "R_Poor", "R_Favoured"]

print(Final_Table)

Output_Location_Percentages_Table = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Bayes_Factor_Percentages.txt"
Final_Table.to_csv(Output_Location_Percentages_Table, index=False)

del All_Objects["BF_Interpretation2"]

if User == 'C':
    Output_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Oscillatory_Objects.txt"
if User == 'N':
    Output_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Oscillatory_Objects.txt"

if User == 'C':
    Output_Location_Table = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Oscillatory_Objects_BF_Table.txt"
if User == 'N':
    Output_Location_Table = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Oscillatory_Objects_BF_Table.txt"

All_Objects.to_csv(Output_Location)
Nice_Table.to_csv(Output_Location_Table)
