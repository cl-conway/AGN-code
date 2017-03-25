using PyCall
using CARMA
using Ensemble
using PyPlot

@pyimport seaborn as sns

"""-----------------------------------------------------------------------------
Description: This julia script performs the analysis to output a PSD plot for
all the Graham candidate objects. It also performs an estimate of the attained
period for the nested sampling procedure on that object, outputting this
information as a text file. Further, it forms the Quality Factor plot for the
object. This script is designed to run on the Tsunami cluster and is initiated
by the Run_Julia.sub submit file.

Inputs:
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_ID_Values.txt
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_ID_Values.txt

This will then create an image file (with its name being the object name):
C:/Users/User/Documents/University/Year 4/Project/Julia_Working_Directory/Graham_Output_Images/
C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Graham_Output_Images_Values.txt

Date Created: 02/02/2017
Authors: Nicholas Kinsey, Christopher Conway
-----------------------------------------------------------------------------"""

#Read the datafile argument and set the Object name
datafile = ARGS[1]

#Read and interpret the identifier
Identifier = ARGS[2]
if Identifier == "G"
    Object_Name = datafile[37:end-4]
    Output_Text_Location = "/home/clc304/Output/Graham_Text_Output/" * Object_Name * ".txt"
    Output_PSD_Location = "/home/clc304/Output/Graham_PSD_Output/" * Object_Name * ".jpg"
    Output_QF_Location = "/home/clc304/Output/Graham_QF_Output/" * Object_Name * ".jpg"
	  Output_LC_Location = "/home/clc304/Output/Graham_LC_Output/" * Object_Name * ".jpg"
    Output_Scatter_Location = "/home/clc304/Output/Graham_Scatter_Output/" * Object_Name * ".jpg"
	   Output_Post_Location = "/home/clc304/Output/Graham_Post_Output/" * Object_Name * ".dat"
    QF_Error_Location = "/home/clc304/Output/Graham_QF_Output/Error_" * Object_Name * ".txt"

elseif Identifier == "R"
    Object_Name = datafile[37:end-4]
    Output_Text_Location = "/home/clc304/Output/Random_Text_Output/" * Object_Name * ".txt"
    Output_PSD_Location = "/home/clc304/Output/Random_PSD_Output/" * Object_Name * ".jpg"
    Output_QF_Location = "/home/clc304/Output/Random_QF_Output/" * Object_Name * ".jpg"
    Output_LC_Location = "/home/clc304/Output/Random_LC_Output/" * Object_Name * ".jpg"
    Output_Scatter_Location = "/home/clc304/Output/Random_Scatter_Output/" * Object_Name * ".jpg"
	   Output_Post_Location = "/home/clc304/Output/Random_Post_Output/" * Object_Name * ".dat"
    QF_Error_Location = "/home/clc304/Output/Random_QF_Output/Error_" * Object_Name * ".txt"

elseif Identifier == "S"
    Object_Name = datafile[39:end-4]
    Output_Text_Location = "/home/clc304/Output/Sinusoid_Text_Output/" * Object_Name * ".txt"
    Output_PSD_Location = "/home/clc304/Output/Sinusoid_PSD_Output/" * Object_Name * ".jpg"
    Output_QF_Location = "/home/clc304/Output/Sinusoid_QF_Output/" * Object_Name * ".jpg"
    Output_LC_Location = "/home/clc304/Output/Sinusoid_LC_Output/" * Object_Name * ".jpg"
    Output_Scatter_Location = "/home/clc304/Output/Sinusoid_Scatter_Output/" * Object_Name * ".jpg"
	Output_Post_Location = "/home/clc304/Output/Sinusoid_Post_Output/" * Object_Name * ".dat"
    QF_Error_Location = "/home/clc304/Output/Sinusoid_QF_Output/Error_" * Object_Name * ".txt"

elseif Identifier == "C"
      Object_Name = datafile[37:end-4]
      Output_Text_Location = "/home/clc304/Output/Charisi_Text_Output/" * Object_Name * ".txt"
      Output_PSD_Location = "/home/clc304/Output/Charisi_PSD_Output/" * Object_Name * ".jpg"
      Output_QF_Location = "/home/clc304/Output/Charisi_QF_Output/" * Object_Name * ".jpg"
      Output_LC_Location = "/home/clc304/Output/Charisi_LC_Output/" * Object_Name * ".jpg"
      Output_Scatter_Location = "/home/clc304/Output/Charisi_Scatter_Output/" * Object_Name * ".jpg"
  	  Output_Post_Location = "/home/clc304/Output/Charisi_Post_Output/" * Object_Name * ".dat"
      QF_Error_Location = "/home/clc304/Output/Charisi_QF_Output/Error_" * Object_Name * ".txt"


#Kill the program with recognised output
else
    println("R/G/S Argument Error")
    exit(1)
end

#Set the settings for the MCMC
p = 3
q = 2
nlive = 1024
nmcmc = 128

#Read the data file
data = readdlm(datafile)
data = data[sortperm(data[:,1]),:]

#Set the times, magnitude and magnitude error values
ts = Float64[data[1,1]]
ys = Float64[data[1,2]]
dys = Float64[data[1,3]]

#Inserted from run.carma.jl file
for i in 2:size(data,1)
    t = data[i,1]
    y = data[i,2]
    dy = data[i,3]

    if t == ts[end]
        dy2 = dy*dy
        dys2 = dys[end]*dys[end]

        yy = (y*dys2 + ys[end]*dy2)/(dy2 + dys2)
        dyy = dy*dys[end]/sqrt(dys2 + dy2)

        ys[end] = yy
        dys[end] = dyy
    else
        push!(ts, t)
        push!(ys, y)
        push!(dys, dy)
    end
end

#Fix any zeros in dy---set to minimum positive dy
dys[dys.<=0] = minimum(dys[dys.>0])

#Create the output Light Curve Plot
fig = figure()
errorbar(ts, ys, dys, fmt=".k")

#Format the image
xlabel(L"Times(days)")
ylabel(L"Magnitude")
title(L"Object Light Curve")

#Save the figure
savefig(Output_LC_Location)
close(fig)

#Define the posterior
post = Kalman.CARMAKalmanPosterior(ts, ys, dys, p, q)

#Create the nested state
state = EnsembleNest.NestState(x -> Kalman.log_likelihood(post, x), x -> Kalman.log_prior(post, x), Kalman.init(post, nlive), nmcmc)

#Run the nested sampling procedure
EnsembleNest.run!(state, 0.1, false)

#Define a posterior sample
postsamples, lnprobs = EnsembleNest.postsample(state)

#Save the posterior sample
open(stream -> serialize(stream, (post, postsamples)), Output_Post_Location, "w")

#Create the frequency array
T = maximum(post.ts) - minimum(post.ts)
dt_med = median(diff(post.ts))

#Set the maximum frequency
fmax = 1.0
df = 1.0/(2.0*T)

#Create the array
fs = collect(df:df:fmax)

#Inserted from the CARMA tutorial
psds = zeros(size(fs, 1), 1000)
for i in 1:1000
    #Choose a random posterior sample
    p = postsamples[:, rand(1:size(postsamples,2))]
    psds[:,i] = Kalman.psd(post, p, fs)
end

#Define arrays to hold period information for each sample
pmin = zeros(size(fs, 1))
pmax = zeros(size(fs, 1))
pmed = zeros(size(fs, 1))

#Populate the defined arrays
for i in 1:size(fs, 1)
    pmin[i] = quantile(vec(psds[i,:]), 0.16)
    pmed[i] = median(vec(psds[i,:]))
    pmax[i] = quantile(vec(psds[i,:]), 0.84)
end

#Define arrays to hold the Quality factor values, freqs values, Heat maps and nu values, initialize counter
freqs = Float64[]
Qs = Float64[]
nu_values = Float64[]
Tau = Float64[]
P = Float64[]
Unreasonable_Frequencies_Counter = 0

#Populate the frequency value array
for i in 1:size(postsamples,2)
    p = Kalman.to_params(post, postsamples[:,i])
	push!(nu_values, p.nu)
    arroots = p.arroots
    selector = imag(arroots) .> 0
    P_all = (2*pi)./(imag(p.arroots))
    Tau_all = -1./real(p.arroots)
    for j in 1:3
      if P_all[j] == Inf ||P_all[j] == -Inf
        P_all[j] = 0
      end
    end
    for j in 1:3
      if Tau_all[j] == Inf ||Tau_all[j] == -Inf
        Tau_all[j] = 0
      end
    end
    Squares = (Tau_all.^2 + P_all.^2).^0.5
    Scale = maximum(Squares)
    selector_HM = Squares.==Scale
    if length(Tau_all[selector_HM]) != 3
      push!(Tau, Tau_all[selector_HM][1]./(Scale))
      push!(P,-P_all[selector_HM][1]./(Scale))
    end
    if any(selector)

        if imag(arroots[selector][1]/(2.0*pi)) > fmax
            Unreasonable_Frequencies_Counter = Unreasonable_Frequencies_Counter + 1
        else
            Quality_Factor_Value = abs(imag(arroots[selector][1]/(2.0*pi))/real(arroots[selector][1]))
            Frequency_Value = imag(arroots[selector][1]/(2.0*pi))
            if isnan(Frequency_Value) == true
                println("NaN frequency for $(i)th iteration")
            else
                push!(freqs, Frequency_Value)
                push!(Qs, Quality_Factor_Value)
            end
        end
    end
end

#Calculate the period with errors
if isempty(freqs) == false
    Execute_All = "Y"
else
    Execute_All = "N"
end


if Execute_All == "Y"

    Period = median(1.0./freqs)
    Upper_Period_error = quantile(1.0./freqs, 0.84) - median(1.0./freqs)
    Lower_Period_error = quantile(1.0./freqs, 0.16) - median(1.0./freqs)

    #Round all the period values to 2 d.p.
    Period = round(Period, 2)
    Upper_Period_error = round(Upper_Period_error, 2)
    Lower_Period_error = round(Lower_Period_error, 2)

end

#End of CARMA tutorial insert

#Find the mean value of nu parameter from samples
mean_nu = mean(nu_values)
mean_nu = round(mean_nu, 2)

#Begin a figure plot
fig, ax = subplots()

Chosen_fs_values_i0 = Float64[]
Chosen_fs_values_i1 = Float64[]
Chosen_psd_values_i0 = Float64[]
Chosen_psd_values_i1 = Float64[]
Starting_fs_values_i0 = Float64[]
Starting_psd_values_i0 = Float64[]
Starting_fs_values_i1 = Float64[]
Starting_psd_values_i1 = Float64[]
push!(Chosen_fs_values_i0, 10.0 ^ -2.5)
push!(Chosen_fs_values_i1, 10.0^-2.5 + df )
push!(Starting_fs_values_i0, minimum(fs))
push!(Starting_fs_values_i1, minimum(fs)+df)

#Set the number of iterations for the plotting procedure
iterations = 10

for i in 1:(iterations+1)

    #Create a posterior sample, hence form the psd
    p = postsamples[:,rand(1:size(postsamples,2))]
    psd = Kalman.psd(post, p, fs)

    #Obtain value of the psd at the chosen fs values
    Chosen_psd_i0 = Kalman.psd(post, p, Chosen_fs_values_i0)
    Chosen_psd_i1 = Kalman.psd(post, p, Chosen_fs_values_i1)
    Chosen_psd_valuation_i0 = mean(Chosen_psd_i0)
    Chosen_psd_valuation_i1 = mean(Chosen_psd_i1)
    push!(Chosen_psd_values_i0, Chosen_psd_valuation_i0)
    push!(Chosen_psd_values_i1, Chosen_psd_valuation_i1)

    Starting_psd_i0 = Kalman.psd(post, p, Starting_fs_values_i0)
    Starting_psd_i1 = Kalman.psd(post, p, Starting_fs_values_i1)
    Starting_psd_valuation_i0 = mean(Starting_psd_i0)
    Starting_psd_valuation_i1 = mean(Starting_psd_i1)
    push!(Starting_psd_values_i0, Starting_psd_valuation_i0)
    push!(Starting_psd_values_i1, Starting_psd_valuation_i1)

    #Create the plot
    ax[:plot](fs, psd, color="k", alpha=0.1, linewidth=2)
end

#Customize the plot
yscale("log")
xscale("log")
xlabel(L"Frequency")
ylabel(L"PSD")
title("A PSD plot $iterations iterations")

#Save the figure
savefig(Output_PSD_Location)
close(fig)

#Create a Quality Factor Plot, plot the figure using a distplot function
if Execute_All == "Y"
    sns.distplot(Qs)
    xlabel(L"Q")
    ylabel(L"p(Q)")
    savefig(Output_QF_Location)
    close(fig)
end

#Create the heat map plot
fig = figure()
scatter(P, Tau)
xlabel("Scaled Period")
ylabel("Scaled Tau")
savefig(Output_Scatter_Location)
close(fig)

#Take the means of all the associated psd values for the iterations
Starting_psd_value_i0 = mean(Starting_psd_values_i0)
Starting_psd_value_i1 = mean(Starting_psd_values_i1)
Chosen_psd_value_i0 = mean(Chosen_psd_values_i0)
Chosen_psd_value_i1 = mean(Chosen_psd_values_i1)
Starting_fs_value_i0 = mean(Starting_fs_values_i0)
Starting_fs_value_i1 = mean(Starting_fs_values_i1)
Chosen_fs_value_i0 = mean(Chosen_fs_values_i0)
Chosen_fs_value_i1 = mean(Chosen_fs_values_i1)

#Find the Differences for the psd and frequencies for the chosen values
Starting_psd_log_Difference = log(Starting_psd_value_i1) - log(Starting_psd_value_i0)
Starting_log_fs_Difference = log(Starting_fs_value_i1) - log(Starting_fs_value_i0)
Chosen_psd_log_Difference = log(Chosen_psd_value_i1) - log(Chosen_psd_value_i0)
Chosen_log_fs_Difference = log(Chosen_fs_value_i1) - log(Chosen_fs_value_i0)

#Calculate the gradients
Starting_gradient = Starting_psd_log_Difference / Starting_log_fs_Difference
Starting_gradient = round(Starting_gradient, 2)
Chosen_gradient = Chosen_psd_log_Difference / Chosen_log_fs_Difference
Chosen_gradient = round(Chosen_gradient, 2)

#Call the logZ function
loggedZ = EnsembleNest.logZ(state)
loggedZ = round(loggedZ, 2)

#Writing text
if Execute_All == "Y"
    Writing_text_A = "$(Object_Name);$(size(freqs,1));$(size(postsamples, 2));$(Unreasonable_Frequencies_Counter);$(Period);$(Upper_Period_error);$(Lower_Period_error);$(mean_nu);$(Starting_gradient);$(Chosen_gradient);$(loggedZ)"
else
    Writing_text_A = "$(Object_Name);$(size(freqs,1));$(size(postsamples, 2));$(Unreasonable_Frequencies_Counter);Err;Err;Err;$(mean_nu);$(Starting_gradient);$(Chosen_gradient);$(loggedZ)"
end

#Write the text to the files
open(Output_Text_Location, "w") do f
    write(f, Writing_text_A)
end

#Form a print text to the output file
Print_Text = "Completed task for object " * Object_Name * " of category " * Identifier
println(Print_Text)
