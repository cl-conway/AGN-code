using PyCall
using CARMA
using Ensemble
using PyPlot

#rc("text", usetex=true)

@pyimport seaborn as sns
"""
This is a Julia Script to open all of the saved postsamples from the carma model
fits to a group of objects (sinusoids, Grahams or randoms) and plot all of the
sinusoids roots for each object on an argand diagram.
"""
identifier = ARGS[1]
fmax = 1


function PSD_Slope(f, df, post, postsamples, iterations)
    p = postsamples[:,rand(1:size(postsamples,2))]
    psd_values_i0 = Float64[]
    psd_values_i1 = Float64[]

    for i in iterations
        psd_i0 = Kalman.psd(post, p, f)
        psd_i1 = Kalman.psd(post, p, f+df)
        psd_valuation_i0 = mean(psd_i0)
        psd_valuation_i1 = mean(psd_i1)
        push!(psd_values_i0, psd_valuation_i0)
        push!(psd_values_i1, psd_valuation_i1)
    end

    psd_value_i0 = mean(psd_values_i0)
    psd_value_i1 = mean(psd_values_i1)
    psd_log_difference = log(psd_value_i1) - log(psd_value_i0)
    fs_log_difference = log(f+df) - log(f)
    slope = psd_log_difference/fs_log_difference[1]
    #slope
end

if identifier == "G"
    Object_Names_Loc = "/home/clc304/Data/Graham_ID_Values.txt"
    Data_Location = "/home/clc304/Output/Graham_Post_Output/"
    Output_Heatmap_Location = "/home/clc304/Output/Graham_HM_Tau_P.jpg"
    Output_Heatmap_Location_Sc = "/home/clc304/Output/Graham_HM_Scat_Tau_P.jpg"
    Output_Violinplot_Location = "/home/clc304/Output/Graham_Violin_PSDSlope2.jpg"
    Output_PSD_Slope_Hist_Location = "/home/clc304/Output/Graham_PSDSlope.jpg"
elseif identifier == "R"

    Object_Names_Loc = "/home/clc304/Data/Randoms_ID_Values.txt"
    Data_Location = "/home/clc304/Output/Random_Post_Output/"
    Output_Heatmap_Location = "/home/clc304/Output/Random_HM.jpg"
    Output_Heatmap_Location_Sc = "/home/clc304/Output/Random_HM_Scat.jpg"
    Output_Violinplot_Location = "/home/clc304/Output/Random_Violin_PSDSlope2.jpg"
    Output_PSD_Slope_Hist_Location = "/home/clc304/Output/Random_PSDSlope.jpg"
    """
    Object_Names_Loc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Randoms_ID_Values.txt"
    Data_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Random_Post_Output/"
    Output_Heatmap_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Random_HM_Tau_P_16th.jpg"
    Output_Heatmap_Location_Sc = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Random_HM_Scat_Tau_P_16th.jpg"
    Output_Violinplot_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Random_Violin_PSDSlope_16th.jpg"
    Output_PSD_Slope_Hist_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Julia_Working_Directory/Random_PSDSlope_16th.jpg"
    """
elseif identifier == "S"
    Object_Names_Loc = "/home/clc304/Data/Graham_ID_Values.txt"
    Data_Location = "/home/clc304/Output/Sinusoid_Post_Output/"
    Output_Heatmap_Location = "/home/clc304/Output/Sinusoid_HM.jpg"
    Output_Heatmap_Location_Sc = "/home/clc304/Output/Sinusoid_HM_Scat.jpg"
    Output_Violinplot_Location = "/home/clc304/Output/Sinusoid_Violin_PSDSlope.jpg"
    Output_PSD_Slope_Hist_Location = "/home/clc304/Output/Sinusoid_PSDSlope.jpg"

elseif identifier == "C"

    Object_Names_Loc = "/home/clc304/Data/Charisi_ID_Values.txt"
    Data_Location = "/home/clc304/Output/Charisi_Post_Output/"
    Output_Heatmap_Location = "/home/clc304/Output/Charisi_HM.jpg"
    Output_Heatmap_Location_Sc = "/home/clc304/Output/Charisi_HM_Scat.jpg"
    Output_Violinplot_Location = "/home/clc304/Output/Charisi_Violin_PSDSlope2.jpg"
    Output_PSD_Slope_Hist_Location = "/home/clc304/Output/Charisi_PSDSlope.jpg"

end


Tau = Float64[]
P = Float64[]
PSD_slopes = Float64[]
psd_1000day_slopes = Float64[]
Object_Names = open(Object_Names_Loc)
df_biggest = 0
f_1000days = Float64[1]
f_1000days[1] = 1/1000
for Object in eachline(Object_Names)
  Object = strip(Object)
  if identifier == "G"
    Object = Object[2:end-1]
  elseif identifier == "S"
    Object = Object[2:end-1]
    Object = Object * "_Magerr_True"
  end
  post, postsamples = open(deserialize, Data_Location * Object * ".dat")
  T = maximum(post.ts) - minimum(post.ts)
  df = 1.0/(2.0*T)
  if df > df_biggest
    df_biggest = df
  end
end
close(Object_Names)
Object_Names = open(Object_Names_Loc)

log_fs = Array{Float64}(1)
log_fs = vcat(log_fs, linspace(log(df_biggest),log(1),10))
deleteat!(log_fs,5:10)

PSD_log_frequencies = Float64[]
println(df_biggest)
df=Float64[]
push!(df, df_biggest)

println("hello")

for Object1 in eachline(Object_Names)
  Object1 = strip(Object1)
  if identifier == "G"
    Object1 = Object1[2:end-1]
  elseif identifier == "S"
    Object1 = Object1[2:end-1]
    Object1 = Object1 * "_Magerr_True"
  end
  post, postsamples = open(deserialize, Data_Location * Object1 * ".dat")
  println("hello")
  #for i in 1:4
  """
  for i in 1:size(postsamples,2)
    println(i)
    p = Kalman.to_params(post, postsamples[:,i])
    arroots = p.arroots
  """
  #  println("imag(p.arroots) $(imag(p.arroots))")
  #  println("real roots $(real(p.arroots))")
  """
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
    selector = Squares.==Scale
    """

    #println("Tau_all $(Tau_all)")
    #println("P_all $(P_all)")
    #println("Selector $(selector)")
    #if length(Tau_all[selector]) != 3
    #  push!(Tau, Tau_all[selector][1]./(Scale))
    #  push!(P,-P_all[selector][1]./(Scale))
    #  println("Taus $(Tau)")
    #  println("Ps $(P)")
    #end
  #end
  """
  for j in log_fs
    iterations = 100
    f = Float64[]
    push!(f, exp(j))
    alpha = PSD_Slope(f, df, post, postsamples, iterations)

    if abs(alpha)<20
      push!(PSD_slopes, alpha )
      push!(PSD_log_frequencies, round(j,2))
    end
  end
  """
  iterations = 500
  push!(psd_1000day_slopes, PSD_Slope(f_1000days,df,post, postsamples,iterations))
end

"""
if length(Tau) == 0
  println("Error Tau size")
  exit(1)
end

println("Making the CHart")

fig = figure()
scatter(Tau, P)
xlabel("Tau")
ylabel("Period")

if identifier == "G"
  title("Scatterplot of all roots in all Graham Objects")
elseif identifier == "R"
  title("Scatterplot of all roots in all Random Objects")
end

savefig(Output_Heatmap_Location_Sc)
close("all")
println("Made the CHart")


PSD_Frequencies = round(exp(PSD_log_frequencies),2)
fig = sns.violinplot(PSD_log_frequencies,PSD_slopes)
xlabel("Frequency")
ylabel("Power law exponent of PSD")


if identifier == "G"
  title("PSD slopes for all Graham Objects")
elseif identifier == "R"
  title("PSD slopes for all Random Objects")
end
savefig(Output_Violinplot_Location)
close("all")
"""
ML_alpha = median(psd_1000day_slopes)
alpha_upperbound = quantile(psd_1000day_slopes, 0.88)
alpha_lowerbound = quantile(psd_1000day_slopes, 0.12)

alpha_new_array = psd_1000day_slopes[(psd_1000day_slopes.>alpha_lowerbound) & (psd_1000day_slopes.<alpha_upperbound)]
fig=sns.distplot(alpha_new_array)
if identifier == "G"
  title("PSD slopes at frequency of 0.01 for all Graham Objects")
elseif identifier == "R"
  title("PSD slopes at frequency of 0.01 for all Random Objects")
end
xlabel("Log-log PSD slope at f=0.01")
ylabel("Probabilty density")
plot( [ML_alpha, ML_alpha], [0,1], linewidth=2.0, label="Most likely alpha value at $(round(ML_alpha,2))" )
legend(loc = "upper left")
savefig(Output_PSD_Slope_Hist_Location)
close("all")
