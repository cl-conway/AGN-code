"""----------------------------------------------------------------------------
Description: Form the new data file for PG 1302-102 for the new data from Charisi
Essentially takes the data from Graham et al. and Charisi et al. and combines
the output. Note this was also hastily done for just the Charisi data
for comparison

Inputs:
Data for Graham and Charisi PG 1302-102

Outputs:
Combined data for PG 1302-102

Date created: 23/03/2017
"""
import pandas as pd
import matplotlib.pyplot as plt

sigma_level = 5
Upper_G_Error_Level = 0.11747274711043741 + 0.03849447039992238

Original_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Grahams_Data/Data_PG 1302-102.txt"
Extra_Data_Location = "C:/Users/User/Documents/University/Year 4/Project/PG_1302_102_Extra_Data.txt"
file_path_Light_Curve = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/PG_1302_LC_Extra_Data.jpg"
file_path_Light_Curve2 = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/PG_1302_LC_Extra_Data_No_Clipping.jpg"

Original_Data = pd.read_table(Original_Data_Location, sep=' ', header=None)
New_Data = pd.read_table(Extra_Data_Location, sep=' ', header=None)

Original_Data.columns = ["Times", "Mag", "Magerr"]
New_Data.columns = ["Times", "Mag", "Magerr"]

New_Data["Magerr"] = New_Data["Magerr"].round(3)

New_Data["Times"] = New_Data["Times"] - New_Data["Times"].min()

Reject_Error_Counter = 0
No_Objects_Clipped = 0
No_Datapoints_Clipped = 0
Total_No_Points = 0

Original_Data["Times"] = Original_Data["Times"] + 53496.20498

Data = [New_Data, Original_Data]

Total_Data = pd.concat(Data)

Total_Data.sort_values("Times", ascending=True, inplace=True)
Total_Data["Times"] = Total_Data["Times"] - Total_Data["Times"].min()

Total_Data["Times"] = Total_Data["Times"].round(5)
Total_Data["Mag"] = Total_Data["Mag"].round(3)
Total_Data["Magerr"] = Total_Data["Magerr"].round(3)

Total_Data = Total_Data.reset_index()
del Total_Data['index']

Reject_Error_Counter = 0
No_Objects_Clipped = 0
No_Datapoints_Clipped = 0
Total_No_Points = 0

New_Data["Times"] = New_Data["Times"].round(5)
New_Data["Mag"] = New_Data["Mag"].round(3)
New_Data["Magerr"] = New_Data["Magerr"].round(3)

New_Data.sort_values("Times", ascending=True, inplace=True)
New_Data = New_Data.reset_index()
del New_Data['index']

New_Data_selected = New_Data.copy()

Before_clipping_Data = len(New_Data.index)
New_Data_selected = New_Data_selected[abs(New_Data_selected['Mag'] - New_Data_selected['Mag'].median()) < sigma_level*New_Data_selected['Magerr']]
After_clipping_Data = len(New_Data_selected.index)

Mean_Err_Value = New_Data_selected['Magerr'].mean()
if Mean_Err_Value > Upper_G_Error_Level:
    Error_Selection_Passed = 'N'
    Reject_Error_Counter = Reject_Error_Counter + 1
else:
    Error_Selection_Passed = 'Y'

if Before_clipping_Data != After_clipping_Data:
    No_Objects_Clipped = No_Objects_Clipped + 1
    No_Datapoints_Clipped = No_Datapoints_Clipped + (Before_clipping_Data - After_clipping_Data)
    Total_No_Points = Total_No_Points + Before_clipping_Data

Scewed_Data_Value = (New_Data_selected['Times'].median()) / New_Data_selected['Times'].max()
Scewed_Data_Value = round(Scewed_Data_Value, 2)
print(Scewed_Data_Value)
print(No_Datapoints_Clipped)

Times = New_Data_selected[['Times']].as_matrix().ravel()
Errors = New_Data_selected[['Magerr']].as_matrix().ravel()
Mag_orgs = New_Data_selected[['Mag']].as_matrix().ravel()

#Make the plot and format
plt.figure()
plt.errorbar(Times, Mag_orgs, Errors, fmt='.k')
plt.xlabel('Time(days)')
plt.ylabel('Magitude')
Title = 'Light Curve for extra PG 1302-102 data'
plt.title(Title)
plt.gca().invert_yaxis()
plt.savefig(file_path_Light_Curve, bbox_inches='tight')
plt.close("all")

Times2 = New_Data[['Times']].as_matrix().ravel()
Errors2 = New_Data[['Magerr']].as_matrix().ravel()
Mag_orgs2 = New_Data[['Mag']].as_matrix().ravel()

#Make the plot and format
plt.figure()
plt.errorbar(Times2, Mag_orgs2, Errors2, fmt='.k')
plt.xlabel('Time(days)')
plt.ylabel('Magitude')
Title = 'Light Curve for extra PG 1302-102 data, no clipping'
plt.title(Title)
plt.gca().invert_yaxis()
plt.savefig(file_path_Light_Curve2, bbox_inches='tight')
plt.close("all")

New_Data_selected.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Data_PG 1302-102_CO.txt', header=None, index=None, sep=' ', mode='w')
New_Data.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Data_PG 1302_CO_NC.txt', header=None, index=None, sep=' ', mode='w')
