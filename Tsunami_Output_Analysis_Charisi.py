import math
import pandas as pd
import matplotlib.pyplot as plt

def Oscillatory_Stat_Function():

    Data_location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Charisi_Outputs.txt"
    Period_Statistic_Hist = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/C_Ocillatory_Stat_Hist.jpg"
    Carma_Nu_Vales_Hist = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/C_Nu_Values_Hist.jpg"
    Data = pd.read_table(Data_location, sep=';', header=0)

    Reasonable_fqs = Data[["Reasonable Frequencies"]].as_matrix().ravel()
    Total_root_no = Data[["Total Root Number"]].as_matrix().ravel()
    Oscillatory_Stat = Reasonable_fqs/Total_root_no
    Nu_values = Data[['Nu Values']].as_matrix().ravel()

    Data["Oscillatory_Stat"] = Data["Reasonable Frequencies"]/Data["Total Root Number"]
    #Data_G = Data_G.sort_values("Oscillatory_Stat")
    Data["Identifier"] = 1

    plt.figure()
    plt.hist(Oscillatory_Stat)
    plt.xlabel("Oscillatory_Statistic")
    plt.ylabel("Counts")
    plt.title("Charisi Oscillatory Statistic")
    plt.savefig(Period_Statistic_Hist)
    plt.clf()


    plt.figure()
    plt.hist(Nu_values)
    plt.xlabel("Nu Value")
    plt.ylabel("Counts")
    plt.title("Graham Nu values from CARMA")
    plt.savefig(Carma_Nu_Vales_Hist)
    plt.clf()

    Data.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Charisi_Outputs.txt', index_col = 0, mode='w', header=list(Data), sep=';')

    Data_R_location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Random_Outputs.txt"
    Period_Statistic_Hist = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/R_Ocillatory_Stat_Hist.jpg"
    Carma_Nu_Vales_Hist = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/R_Nu_Values_Hist.jpg"
    Data_R = pd.read_table(Data_R_location, sep=';', header=0)

    Reasonable_fqs = Data_R[["Reasonable Frequencies"]].as_matrix().ravel()
    Total_root_no = Data_R[["Total Root Number"]].as_matrix().ravel()
    Oscillatory_Stat = Reasonable_fqs/Total_root_no
    Nu_values = Data_R[['Nu Values']].as_matrix().ravel()
    del Data_R["Initial Gradient"]
    del Data_R["~1 Year Gradient"]

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
    plt.hist(Nu_values)
    plt.xlabel("Nu Value")
    plt.ylabel("Counts")
    plt.title("Random Nu values from CARMA")
    plt.savefig(Carma_Nu_Vales_Hist)
    plt.clf()

    frames = [Data,Data_R]
    Data_All = pd.concat(frames)
    Data_All.sort_values("Oscillatory_Stat", inplace=True, ascending = False)
    Data_All['Oscillatory_Stat'] = Data_All['Oscillatory_Stat'].round(3)
    Data_All.reset_index(inplace=True)
    del Data_All['index']
    Data_All.to_csv(r'C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/All_Outputs_CR.txt', index_col = 0, mode='w', header=list(Data_All), sep=';')

    Data_All['Charisi Objects ordered by Oscillatory Statistic'] = Data_All.Identifier.cumsum()
    Data_All['position'] = Data_All.index + 1
    plt.figure()
    Data_All.plot(x='position', y='Charisi Objects ordered by Oscillatory Statistic')
    y = Data_All.index*(Data_All["Charisi Objects ordered by Oscillatory Statistic"].iloc[-1]/Data_All["position"].iloc[-1])
    plt.plot(Data_All[["position"]],y)
    plt.xlabel("Objects Searched")
    plt.ylabel("Grahams Encountered")
    plt.title("Gains Chart of Charisi using Oscillatory Statistic")
    plt.savefig("C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/C_Gains_plot.jpg")
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

    plt.figure()
    Data_DR_sort.plot(x='position', y='Graham_Objects_ordered_by_Damping_Ratio')
    y = Data_DR_sort.index*(Data_DR_sort["Graham_Objects_ordered_by_Damping_Ratio"].iloc[-1]/Data_DR_sort["position"].iloc[-1])
    plt.plot(Data_DR_sort[["position"]],y)
    plt.xlabel("Objects Searched")
    plt.ylabel("Grahams Encountered")
    plt.title("Gains Chart of Grahams using Damping Ratio")
    plt.savefig("C:/Users/User/Documents/University/Year 4/Project/Gains_plot_DR.jpg")
    plt.close('all')

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
    plt.savefig("C:/Users/User/Documents/University/Year 4/Project/Gains_plot_QF.jpg")
    plt.close('all')

Oscillatory_Stat_Function()
