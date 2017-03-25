"""-----------------------------------------------------------------------------
This script has the ability to calculate the oscillatory statistic, damping Ratio
and quality factor for the various objects. It outputs the gains charts for
all the objects as a whole for each of those criterion.

Inputs:
The processed results for the objects, processed through Tsunami_Data_Process.
We can then access all the data that the file contains (all objects of that type)

Outputs:
Gains Chart

-----------------------------------------------------------------------------"""
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def Oscillatory_Stat_Function():

    Data_G_location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graham_Outputs.txt"
    Period_Statistic_Hist = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/G_Ocillatory_Stat_Hist.jpg"
    Carma_Nu_Vales_Hist = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/G_Nu_Values_Hist.jpg"
    Data_G = pd.read_table(Data_G_location, sep=';', header=0)

    Reasonable_fqs = Data_G[["Reasonable Frequencies"]].as_matrix().ravel()
    Total_root_no = Data_G[["Total Root Number"]].as_matrix().ravel()
    Oscillatory_Stat = Reasonable_fqs/Total_root_no
    Nu_values = Data_G[['Nu Values']].as_matrix().ravel()

    Data_G["Oscillatory_Stat"] = Data_G["Reasonable Frequencies"]/Data_G["Total Root Number"]
    #Data_G = Data_G.sort_values("Oscillatory_Stat")
    Data_G["Identifier"] = 1

    plt.figure()
    plt.hist(Oscillatory_Stat)
    plt.xlabel("Oscillatory_Statistic")
    plt.ylabel("Counts")
    plt.title("Graham Oscillatory Statistic")
    plt.savefig(Period_Statistic_Hist)
    plt.clf()

    plt.figure()
    plt.hist(Nu_values, normed=True)
    plt.xlabel("Nu Value")
    plt.ylabel("Prob. density")
    plt.title("Graham Nu values from CARMA")
    plt.savefig(Carma_Nu_Vales_Hist)
    plt.clf()

    Data_R_location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Random_Outputs.txt"
    Period_Statistic_Hist = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/R_Ocillatory_Stat_Hist.jpg"
    Carma_Nu_Vales_Hist = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/R_Nu_Values_Hist.jpg"
    Data_R = pd.read_table(Data_R_location, sep=';', header=0)

    Reasonable_fqs = Data_R[["Reasonable Frequencies"]].as_matrix().ravel()
    Total_root_no = Data_R[["Total Root Number"]].as_matrix().ravel()
    Oscillatory_Stat = Reasonable_fqs/Total_root_no
    Nu_values = Data_R[['Nu Values']].as_matrix().ravel()

    Data_R["Oscillatory_Stat"] = Data_R["Reasonable Frequencies"]/Data_R["Total Root Number"]
    #Data_R = Data_R.sort_values("Oscillatory_Stat")
    Data_R["Identifier"] = 0

    plt.figure()
    plt.hist(Oscillatory_Stat)
    plt.xlabel("Oscillatory_Statistic")
    plt.ylabel("Counts")
    plt.title("Random Oscillatory Statistic")
    plt.savefig(Period_Statistic_Hist)
    plt.clf()

    plt.figure()
    plt.hist(Nu_values, normed=True)
    plt.xlabel("Nu Value")
    plt.ylabel("Prob. density")
    plt.title("Random Nu values from CARMA")
    plt.savefig(Carma_Nu_Vales_Hist)
    plt.clf()

    frames = [Data_G,Data_R]
    Data_All = pd.concat(frames)
    Data_All.sort_values("Oscillatory_Stat", inplace=True, ascending = False)
    Data_All['Oscillatory_Stat'] = Data_All['Oscillatory_Stat'].round(3)
    Data_All.reset_index(inplace=True)
    del Data_All['index']
    #Data_All.to_csv(r'C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/All_Outputs.txt', index_col = 0, mode='w', header=list(Data_All), sep=';')

    Data_All['cum_sum'] = Data_All.Identifier.cumsum()
    Data_All['position'] = Data_All.index + 1
    plt.figure()
    Data_All.plot(x='position', y='cum_sum', label="Oscillatory Statistic Ordering")
    y = Data_All.index*(Data_All["cum_sum"].iloc[-1]/Data_All["position"].iloc[-1])
    plt.plot(Data_All[["position"]],y, label="Random line")
    plt.xlabel("Total Number of objects searched")
    plt.ylabel("Number of Graham periodic candidates encountered")
    #plt.title("Gains Chart of Grahams using Oscillatory Statistic")
    plt.savefig("C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Gains_plot.jpg")
    plt.close('all')

def QF_and_DR_Values():

    Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/QF_Text.txt"
    Data = pd.read_table(Data_Location, sep=';', header=None)

    Data.columns = ['Object', 'QF', 'QF_U', 'QF_L', 'Identifier', 'DR']
    Data_DR_sort = Data.sort_values("DR", ascending = True)
    Data_DR_sort['Graham_Objects_ordered_by_Damping_Ratio'] = Data_DR_sort.Identifier.cumsum()
    Data_DR_sort.reset_index(inplace=True)
    del Data_DR_sort['index']
    Data_DR_sort['position'] = Data_DR_sort.index + 1
    """
    plt.figure()
    Data_DR_sort.plot(x='position', y='Graham_Objects_ordered_by_Damping_Ratio')
    y = Data_DR_sort.index*(Data_DR_sort["Graham_Objects_ordered_by_Damping_Ratio"].iloc[-1]/Data_DR_sort["position"].iloc[-1])
    plt.plot(Data_DR_sort[["position"]],y)
    plt.xlabel("Objects Searched")
    plt.ylabel("Grahams Encountered")
    plt.title("Gains Chart of Grahams using Damping Ratio")
    plt.savefig("C:/Users/User/Documents/Universityversity/Year 4/Project/Gains_plot_DR.jpg")
    plt.close('all')
    """
    plt.figure()
    Data_G = Data[Data["Identifier"] == 1]
    Data_G = Data_G[Data_G["QF"] <= 10]
    Data_R = Data[Data["Identifier"] == 0]
    Data_R = Data_R[Data_R["QF"] <= 10]

    #n, bins = np.histogram(Data_G.QF, bins='auto' ) # np arrays of counts and period bin
    plt.hist(Data_G.QF, bins = 'auto', normed = 1)
    plt.show()

    #n, bins = np.histogram(Data_R.QF, bins='auto' ) # np arrays of counts and period bin
    plt.hist(Data_R.QF, bins = 'auto', normed = 1)
    plt.show()
    exit(1)


    Data_QF_sort = Data.sort_values("QF", ascending = False)
    Data_QF_sort['Graham_Objects_ordered_by_Quality_Factor'] = Data_QF_sort.Identifier.cumsum()
    Data_QF_sort.reset_index(inplace=True)
    del Data_QF_sort['index']
    Data_QF_sort['position'] = Data_QF_sort.index + 1

    plt.figure()
    Data_QF_sort.plot(x='position', y='Graham_Objects_ordered_by_Quality_Factor')
    y = Data_QF_sort.index*(Data_QF_sort["Graham_Objects_ordered_by_Quality_Factor"].iloc[-1]/Data_QF_sort["position"].iloc[-1])
    plt.plot(Data_QF_sort[["position"]],y)
    plt.xlabel("Objects Searched")
    plt.ylabel("Grahams Encountered")
    plt.title("Gains Chart of Grahams using Quality_Factor")
    plt.savefig("C:/Users/User/Documents/Universityversity/Year 4/Project/Gains_plot_QF.jpg")
    plt.close('all')

#QF_and_DR_Values()
Oscillatory_Stat_Function()
