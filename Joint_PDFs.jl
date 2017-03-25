using PyCall
using CARMA
using Ensemble
using PyPlot
using DataArrays, DataFrames
@pyimport seaborn as sns

function make_Joint_Dist(Periods, Quality_Factors, Object, Joint_Dist_Save_Location)
  #println("Quality factor size $(size(Quality_Factors))")
  Med_Period = median(Periods)
  Period_Upper_limit = quantile(Periods, 0.86) - Med_Period
  Period_Upper_limit_plot = maximum(Periods) +Period_Upper_limit
  Quality_Factor_Upper_Limit_Plot = quantile(Quality_Factors, 0.96)
  #Quality_Factors_new = Quality_Factors[selector_1]
  #Period_new = Periods[selector_1]
  fig = sns.jointplot(Periods, Quality_Factors, kind = "kde", xlim = [0,1500], ylim= [0,5])
  xlabel("Period (days)")
  ylabel("Quality Factor")
  savefig(Joint_Dist_Save_Location)
  close("all")
end

function make_Ps_and_QF_Array(post, postsamples)
  fmax = 1
  Unreasonable_Frequencies_Counter = 0
  Periods = Float64[]
  Quality_Factors = Float64[]
  freqs = Float64[]
  for i in 1:size(postsamples,2)
      p = Kalman.to_params(post, postsamples[:,i])
      arroots = p.arroots
      selector = imag(arroots) .> 0

      if any(selector)
          if imag(arroots[selector][1]/(2.0*pi)) > fmax
              Unreasonable_Frequencies_Counter = Unreasonable_Frequencies_Counter + 1
          else
              Quality_Factor_Value = abs(imag(arroots[selector][1]/(2.0*pi))/real(arroots[selector][1]))
              Period_Value = (2*pi)./(imag(arroots[selector][1]))
              Frequency_Value = imag(arroots[selector][1]/(2.0*pi))
              if isnan(Period_Value) == true
                  println("NaN frequency for $(i)th iteration")
              else
                  push!(Periods, Period_Value)
                  push!(Quality_Factors, Quality_Factor_Value)
                  push!(freqs, Frequency_Value)
              end
          end
      end
  end
  Periods, Quality_Factors, Unreasonable_Frequencies_Counter, freqs
end

function make_Heat_Map(post, postsample_Outputs, identifier, Output_Heatmap_Location_Sc, Output_Heatmap_Location)
  Tau = Float64[]
  P = Float[]
  for i in 1:size(postsamples,2)
    p = Kalman.to_params(post, postsamples[:,i])
    arroots = p.arroots
    Tau_all = -1./real(p.arroots)
    P_all = 2*pi./imag(p.arroots)
    Squares = Tau_all.^2 + P_all.^2
    Scale = maximum(Squares)
    selector = Squares.==Scale
    if length(Tau_all[selector])== 1
      push!(Tau, Tau_all[selector][1]./(Scale^0.5))
      push!(P,P_all[selector][1]./(Scale^0.5))
    end
  end

  fig = figure()
  sns.jointplot(Tau,P, kind = "kde", color = "r")
  xlabel("Tau")
  ylabel("Period")
  if identifier == "G"
    title("Plot of Period against Tau for all Graham Objects")
  elseif identifier == "R"
    title("Plot of Period against Tau for all Random Objects")
  end
  savefig(Output_Heatmap_Location)
  close(fig)

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
end

All_Data_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/All_Outputs.txt"
df = readtable(All_Data_Location, separator = ';', header = true)
df = df[(df[:Oscillatory_Stat] .> 0.5), :]

for i in 1:size(df,1)
  Object_Name = df[i,:Object_Name]
  println(Object_Name)
  Identifier = df[i,:Identifier]
  if Identifier == 1
    Posterior_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Tsunami_Results/Graham_Post_Output/" * Object_Name * ".dat"
    Joint_Dist_Save_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Oscillatory_Joint_PDFs/Grahams/" *Object_Name *".jpg"
  elseif Identifier == 0
    Posterior_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Tsunami_Results/Tsunami_Results/Random_Post_Output/" * Object_Name * ".dat"
    Joint_Dist_Save_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Oscillatory_Joint_PDFs/Randoms/" *Object_Name *".jpg"
  end
  post, postsamples = open(deserialize, Posterior_Location)
  Periods, Quality_Factors, Unreasonable_Frequencies_Counter, freqs = make_Ps_and_QF_Array(post, postsamples)
  make_Joint_Dist(Periods, Quality_Factors, Object_Name, Joint_Dist_Save_Location)
  if Object_Name == "SDSS J144755.57+100040.0" || Object_Name =="SDSSJ144755.57+100040.0"
    println("$(median(Quality_Factors))")
    println("$(quantile(Quality_Factors,0.84))")
    println("$(quantile(Quality_Factors, 0.16))")
    break
  end
end
