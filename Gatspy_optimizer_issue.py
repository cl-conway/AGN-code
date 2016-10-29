""" --------------------------------------------------------------------------
Description:
Highlighting the optimizer issue within the Gatspy package. Optimizer appears
to fall on a specific value, when the maximum height of the L-S fit is not
at that value.

Inputs:
Raw Data: Gatspy_issue_Data.txt

Outputs:
Creates a historgram (which is shown but not saved)

Date Created: 29/10/2016
Authors: Nicholas Kinsey, Christopher Conway
----------------------------------------------------------------------------"""

#Perform imports
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargleFast

#Grab the raw data
Raw_Data = pd.read_table('Gatspy_issue_Data.txt', sep=',', header=0)

#Seperate out the individual pieces of data
Times = Raw_Data[['MJD']].as_matrix().ravel()
Mag = Raw_Data[['Mag']].as_matrix().ravel()
Errors = Raw_Data[['Magerr']].as_matrix().ravel()

#Set the number of iterations here
iterations = 100

#Define a period list to save all the best periods found by optimizer into
period_best_list = []
period_max_list = []

#Define two variables to keep track of a specific value
count_best = 0
count_max = 0

for i in range(iterations):

    #Perfrom Gaussian distribution among error bars
    Mags_mod = Mag + Errors * np.random.randn(Mag.size)

    #Build the L-S model
    model = LombScargleFast().fit(Times, Mags_mod, Errors)

    #Set the periods and power, note increasing oversampling reduces the percentage
    #of the max values that are the same, as expected. But does not change the
    #percentage of the best (optimizer) periods that are the same. Hence, we
    #believe the optimizer is broken.
    periods, power = model.periodogram_auto(nyquist_factor=0.1, oversampling=500)

    #Find best period using the optimizer, looking at specific range
    model.optimizer.period_range=(500, 2000)
    best_period = model.best_period

    #Instead look at maximum for the L-S fit
    max_period = periods[np.argmax(power)]

    #Keep a count of exact matching values
    if best_period in period_best_list: count_best = count_best + 1
    if max_period in period_max_list: count_max = count_max + 1

    #Save the best_period for that iteration into the period list
    period_best_list.append(best_period)
    period_max_list.append(max_period)

#Percentage of iterations with the exact same best/max period
Percentage_value_best = int((count_best / iterations)*100)
Percentage_value_max = int((count_max / iterations)*100)
print('The percentage of iterations with exact best period twice is:', str(Percentage_value_best), '%')
print('The percentage of iterations with exact max period twice is:', str(Percentage_value_max), '%')

#Create the plot and show
fig, ax = plt.subplots()
n, bins, patches = plt.hist(period_best_list, bins='auto' )
ax.hist(period_best_list, bins='auto')
ax.set(xlabel = 'Period',
        ylabel = 'No. of Counts',
        title = "Histogram, BEST optimizer " + str(iterations) + ' Iterations')

#Create the plot and show
fig, ax = plt.subplots()
n, bins, patches = plt.hist(period_max_list, bins='auto' )
ax.hist(period_max_list, bins='auto')
ax.set(xlabel = 'Period',
        ylabel = 'No. of Counts',
        title = "Histogram, MAX optimizer " + str(iterations) + ' Iterations')

plt.show()

#End of file
