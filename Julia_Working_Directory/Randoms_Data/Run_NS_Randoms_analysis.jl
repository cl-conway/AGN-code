using PyCall
using CARMA
using Ensemble
using PyPlot

@pyimport seaborn as sns

"""-----------------------------------------------------------------------------
Description: This julia script performs the analysis to output a PSD plot for
all the Random AGN objects.

Inputs:
state-3-2.dat file in the same directory as this script
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_ID_Values.txt

This will then create an image file (with its name being the object name):
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Output_Images/

Date Created: 20/12/2016
Authors: Nicholas Kinsey, Christopher Conway
-----------------------------------------------------------------------------"""

usage = "julia Run_NS_Randoms_analysis.jl Object_Name"

#Take the argument as the object name
Object = ARGS[1]
Object_Name = replace(Object, "_", " ")

#Define the location to the output of the nested sampling procedure
Data_Location = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Randoms_Data/state-3-2.dat"

#Open the data file, created by the nested sampling prcedure
post, state = open(deserialize, Data_Location)

postsamples, lnprobs = EnsembleNest.postsample(state)
freqs = Kalman.frequencies(post, postsamples)

fs = Kalman.psdfreq(post)

#Begin a figure plot
fig, ax =subplots()

#Set the number of iterations for the plotting procedure
iterations = 10

for i in 1:(iterations+1)

    p = postsamples[:,rand(1:size(postsamples,2))]
    psd = Kalman.psd(post, p, fs)

    #Create the plot
    ax[:plot](fs, psd, color="k", alpha=0.1, linewidth=2)

end

#Customize the plot
yscale("log")
xscale("log")
xlabel("Frequency")
ylabel("PSD")
title("A PSD plot $iterations iterations")

#Create a figure path, which includes the objects name and save the figure
Figure_path = "C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Random_Output_Images/"
Figure_name = "PSD_Plot_" * Object_Name * ".jpg"
Total_Figure_path = Figure_path * Figure_name
savefig(Total_Figure_path)
