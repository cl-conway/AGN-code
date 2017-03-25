using PyCall
using PyPlot;
plt = PyPlot;
using CARMA
using Ensemble
using DataFrames

@pyimport seaborn as sns
sns.set_context("notebook")
sns.set_style("ticks")
sns.set_palette("colorblind")
"""
This script implements the predict function as discussed on 21/02/2017, onto
objects with a Osc. ratio greater than 0.5.

It allows the fitting of the CARMA model to the lightcurve. It also, plots
the binned lightcurve, plots a histogram of nu values, plots the lightcurve
folded onto the period and a shows a plot of the individual datapoint residuals

It also plots the CAR(1) prediction onto the lightcurve to compare with the CARMA
model. 
"""
function bin_folded_lightcurve(ts, ys, dys, bins, P)
    ts_fold = ts % P
    fluxes = zeros(size(bins,1)-1)
    dfluxes = zeros(size(bins,1)-1)

    for i in 1:size(bins, 1)-1
        low = bins[i]
        high = bins[i+1]

        sel = (ts_fold .>= low) & (ts_fold .< high)
        y = ys[sel]
        dy = dys[sel]

        wts = 1.0./(dy.*dy)

        fluxes[i] = sum(wts.*y)./sum(wts)
        dfluxes[i] = sqrt(sum(wts.*wts.*dy.*dy)./(sum(wts).^2))
    end

    fluxes, dfluxes
end

function Call_Predict(Path_To_Post, Path_To_CAR1_Post, Object, Nu_Dist_Save_Loc, Predict_Plot_Save_Loc, Folded_On_Period_Plot_Save_Loc, max_bins, Binned_LightCurve_Save_Loc, Residuals_Save_Loc, Period)

    #Open the data
    post, postsamples = open(deserialize, Path_To_Post)
    post_CAR1, postsamples_CAR1 = open(deserialize, Path_To_CAR1_Post)

    #Define the paths to the save locations
    Nu_Dist_Save_Location = Nu_Dist_Save_Loc * Object * ".jpg"
    Predict_Plot_Save_Location = Predict_Plot_Save_Loc * Object * ".jpg"
    Folded_On_Period_Plot_Save_Location = Folded_On_Period_Plot_Save_Loc * Object * ".jpg"
    Binned_LightCurve_Save_Location = Binned_LightCurve_Save_Loc * Object * ".jpg"
    Residuals_Save_Location = Residuals_Save_Loc * Object * ".jpg"

    #Nu values distribution plot
    fig = sns.distplot(vec([Kalman.to_params(post, postsamples[:,i]).nu for i in 1:size(postsamples,2)]))
    xlabel(L"Nu Value")
    ylabel(L"Counts")
    title("$(Object) Nu Plot")
    savefig(Nu_Dist_Save_Location)
    close("all")

    #Plot the error bars
    errorbar(post.ts, post.ys, post.dys, fmt=".", color="k", alpha=0.3)
    ax = plt.gca()
    ax[:invert_yaxis]()
    ax[:tick_params](axis="x",top="off")
    ax[:tick_params](axis="y",right="off")
    xlabel("Times(days)")
    ylabel("Magnitude")
    title("$(Object) Light Curve, with Predict fit")

    #Call the predict function
    ts_predict = sort(vcat(post.ts, collect(linspace(post.ts[1], post.ts[end], 256))))
    prandom = postsamples[:,rand(1:size(postsamples,2))]
    ys, vys = Kalman.predict(post, prandom, ts_predict)
    dys = sqrt(vys)

    ts_CAR1_predict = sort(vcat(post_CAR1.ts, collect(linspace(post_CAR1.ts[1], post_CAR1.ts[end], 256))))
    prandom_CAR1 = postsamples_CAR1[:,rand(1:size(postsamples_CAR1,2))]
    #Filter = Kalman.make_filter(post_CAR1, prandom_CAR1)
    #ys_CAR1, vys_CAR1 = Kalman.predict(Filter)

    ys_CAR1, vys_CAR1 = Kalman.predict(post_CAR1, prandom_CAR1, ts_CAR1_predict)
    dys_CAR1 = sqrt(vys_CAR1)

    #Save the figure and fill in the errors for the predict function
    plot(ts_predict, ys, color=sns.color_palette()[1], label="CARMA model")
    fill_between(ts_predict, ys+dys, ys-dys, color=sns.color_palette()[1], alpha=0.3)
    plot(ts_CAR1_predict, ys_CAR1, color=sns.color_palette()[2], label="CAR model")
    fill_between(ts_CAR1_predict, ys_CAR1+dys_CAR1, ys_CAR1-dys_CAR1, color=sns.color_palette()[2], alpha=0.3)
    legend(loc="upper left")
    savefig(Predict_Plot_Save_Location)
    close("all")


    #Form the light curve, folded onto the period
    P = Period
    Counter = 0

    while isfinite(P) == false
        prandom = postsamples[:,rand(1:size(postsamples,2))]
        P = 1.0./Kalman.frequencies(post, prandom)[2]
        Counter = Counter + 1

        if Counter > 10
            println("P still Inf")
            exit(1)
        end
    end
    errorbar(post.ts % P, post.ys, post.dys, fmt=".")
    ax = plt.gca()
    ax[:invert_yaxis]()
    ax[:tick_params](axis="x",top="off")
    ax[:tick_params](axis="y",right="off")
    xlabel("Times, folded onto Period (days)")
    ylabel("Magnitude")
    title("$(Object) Data folded onto attained Period")
    savefig(Folded_On_Period_Plot_Save_Location)
    close("all")

    #Define the bins
    bins = collect(linspace(0, P, max_bins))
    bcent = 0.5*(bins[1:end-1] + bins[2:end])

    #Form the folded lightcurve from the data
    f, df = bin_folded_lightcurve(post.ts, post.ys, post.dys, bins, P)

    #Form the folded lightcurve from the prediction
    pf, pdf = bin_folded_lightcurve(ts_predict, ys, dys, bins, P)

    #Plotting
    errorbar(bcent, f, df, fmt=".", label="Original Data")
    errorbar(bcent, pf, pdf, fmt=".", label="Predict fitting")
    ax = plt.gca()
    ax[:invert_yaxis]()
    ax[:tick_params](axis="x",top="off")
    ax[:tick_params](axis="y",right="off")
    legend(loc="upper left")
    xlabel("Times, folded onto Period, binned with $(max_bins) bins. (days)")
    ylabel("Magnitude")
    title("$(Object) Data folded onto attained Period, binned with $(max_bins) bins.")
    savefig(Binned_LightCurve_Save_Location)
    close("all")

    #Form the residuals and plot
    rs, drs = Kalman.residuals(post, prandom)
    errorbar(post.ts, rs, drs, fmt=".")
    ax = plt.gca()
    ax[:tick_params](axis="x",top="off")
    ax[:tick_params](axis="y",right="off")
    xlabel("Times (days)")
    ylabel("Magnitude")
    title("$(Object) Residuals for the data")
    savefig(Residuals_Save_Location)
    close("all")

end


function Plot_Periods(Path_To_Posts, Object, Grahams_Data_Period, Period_Plot_Save_Loc, Damping_Ratio_Save_Loc)

    #Open the data
    post, postsamples = open(deserialize, Path_To_Posts)

    #Define the paths to the save locations
    Period_Plot_Save_Location = Period_Plot_Save_Loc * Object * ".jpg"

    T = maximum(post.ts) - minimum(post.ts)
    Periods = Float64[]
    fmax = 1.0
    fmin = 1.0/T

    #Populate the frequency value array
    for i in 1:size(postsamples,2)
        p = Kalman.to_params(post, postsamples[:,i])
        arroots = p.arroots
        selector = imag(arroots) .> 0

        if any(selector)
            if imag(arroots[selector][1]/(2.0*pi)) < fmax && imag(arroots[selector][1]/(2.0*pi)) > fmin
                push!(Periods, 1.0./imag(arroots[selector][1]/(2.0*pi)))
            end
        end
    end



    #Dist_plot = sns.distplot(Periods)

    if Grahams_Data_Period != -99
        plot([Grahams_Data_Period, Grahams_Data_Period], [0, 0.0025], linewidth=2.0, label="Grahams Value")
        legend(loc="upper right")
    end

    xlabel("Period (days)")
    ylabel("Normalized Distribution")
    savefig(Period_Plot_Save_Location)
    close("all")

end

function Calc_Q_Values(Path_To_Post, Object)

    post, postsamples = open(deserialize, Path_To_Post)

    Qs = Float64[]
    fmax = 1.0

    for i in 1:size(postsamples,2)
        p = Kalman.to_params(post, postsamples[:,i])
        arroots = p.arroots
        selector = imag(arroots) .> 0

        if any(selector)

            if imag(arroots[selector][1]/(2.0*pi)) < fmax
                Quality_Factor_Value = abs(imag(arroots[selector][1]/(2.0*pi))/real(arroots[selector][1]))
                push!(Qs, Quality_Factor_Value)
            end
        end
    end

    Quality_Factor = median(Qs)
    Upper_QF_error = quantile(Qs, 0.84) - median(Qs)
    Lower_QF_error = quantile(Qs, 0.16) - median(Qs)

    Quality_Factor = round(Quality_Factor, 2)
    Upper_QF_error = round(Upper_QF_error, 2)
    Lower_QF_error = round(Lower_QF_error, 2)

    Object, Quality_Factor, Upper_QF_error, Lower_QF_error

end
#-----------------------------------END OF FUNCTIONS------------------------------------------------

#Set the User
User = "N"

#Input data Paths
if User == "N"
    Path_To_Outputs = "C:/Users/User/Documents/University/Year 4/Project/All_Outputs.txt"
    Path_To_Graham_Data = "C:/Users/User/Documents/University/Year 4/SAFE_File/Graham_Period_Data.txt"
    Q_Text_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/QF_Text"
elseif User == "C"
    Path_To_Outputs = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/All_Outputs.txt"
else
    println("User error")
    exit(1)
end

#Read the data
All_Data_Output = readtable(Path_To_Outputs, separator=';', header=true);
Grahams_Period_Values = readtable(Path_To_Graham_Data, separator='\t', header=true);

#Filter data so only those with Osc Stat > 0.5 remain
Filtered_Data = All_Data_Output[All_Data_Output[:Oscillatory_Stat] .> 0.5,:]

for i in 1:nrow(All_Data_Output)
    Object = All_Data_Output[i, :Object_Name]
    if Object == "PG 1302-102"
        append!(Filtered_Data, All_Data_Output[i,:])
    end
end

Quality_Factor_Data = Any[]

#nrow(All_Data_Output)
for i in 1:nrow(All_Data_Output)

    Object_Name = All_Data_Output[i, :Object_Name]
    Period = All_Data_Output[i, :Period]
    max_bins = 15
    Grahams_Data_Period = 0

    if All_Data_Output[i, :Identifier] == 1
        for j in 1:nrow(Grahams_Period_Values)
            if contains(Grahams_Period_Values[j, :Object], Object_Name) == true
                Grahams_Data_Period = Grahams_Period_Values[j, :Period]
                break;
            end
        end
    else
        Grahams_Data_Period = -99
    end

    #If statement to identify nature of object
    if All_Data_Output[i, :Identifier] == 1

        #If statement regarding the user
        if User == "N"
            Path_To_Posts = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Graham_Post_Output/"
            Path_To_CAR1_Posts = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Oscillatory_CAR1_Post_Output/Grahams/"
            Nu_Dist_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Nu_Values/Grahams/[$(i)]"
            Predict_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Predict_Plots/Grahams/[$(i)]"
            Folded_On_Period_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Folded_On_Period/Grahams/[$(i)]"
            Binned_LightCurve_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Binned_Plots/Grahams/[$(i)]"
            Residuals_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Residual_Plots/Grahams/[$(i)]"
            Period_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Period_Plots/Grahams/[$(i)]"
        else
            Path_To_Posts = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Graham_Post_Output/"
            Path_To_CAR1_Posts = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Oscillatory_CAR1_Post_Output/Grahams/"
            Nu_Dist_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Nu_Values/Grahams/[$(i)]"
            Predict_Plot_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Predict_Plots/Grahams/[$(i)]"
            Folded_On_Period_Plot_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Folded_On_Period/Grahams/[$(i)]"
            Binned_LightCurve_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Binned_Plots/Grahams/[$(i)]"
            Residuals_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Residual_Plots/Grahams/[$(i)]"
        end

    elseif All_Data_Output[i, :Identifier] == 0
        if User == "N"
            Path_To_Posts = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Random_Post_Output/"
            Path_To_CAR1_Posts = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Oscillatory_CAR1_Post_Output/Randoms/"
            Nu_Dist_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Nu_Values/Randoms/[$(i)]"
            Predict_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Predict_Plots/Randoms/[$(i)]"
            Folded_On_Period_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Folded_On_Period/Randoms/[$(i)]"
            Binned_LightCurve_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Binned_Plots/Randoms/[$(i)]"
            Residuals_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Residual_Plots/Randoms/[$(i)]"
            Period_Plot_Save_Loc = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Processed_Results/Graphical_Results/Period_Plots/Randoms/[$(i)]"
        else
            Path_To_Posts = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Random_Post_Output/"
            Nu_Dist_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Nu_Values/Randoms/[$(i)]"
            Predict_Plot_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Predict_Plots/Randoms/[$(i)]"
            Folded_On_Period_Plot_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Folded_On_Period/Randoms/[$(i)]"
            Binned_LightCurve_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Binned_Plots/Randoms/[$(i)]"
            Residuals_Save_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Processed_Results/Graphical_Results/Residual_Plots/Randoms/[$(i)]"
        end
    else
        println("Identifier error")
        exit(1)
    end

    #Set the path to the post data
    Path_To_Post_Data = Path_To_Posts * Object_Name * ".dat"
    Path_To_CAR1_Post_Data = Path_To_CAR1_Posts * Object_Name * ".dat"

    #Call predict function
    Call_Predict(Path_To_Post_Data, Path_To_CAR1_Post_Data, Object_Name, Nu_Dist_Save_Loc, Predict_Plot_Save_Loc, Folded_On_Period_Plot_Save_Loc, max_bins, Binned_LightCurve_Save_Loc, Residuals_Save_Loc, Period)

    #Object, QF, QF_U, QF_L = Calc_Q_Values(Path_To_Post_Data, Object_Name)
    #New_Array = [Object, QF, QF_U, QF_L, All_Data_Output[i, :Identifier], Damping_Ratio]
    #append!(Quality_Factor_Data, New_Array)
end
