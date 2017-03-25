#Script to fit a CAR(1) model to all the objects found to be periodic
#(Oscillatory_Stat >0.5)

using PyCall
using CARMA
using Ensemble
using PyPlot

@pyimport seaborn as sns

#Read the datafile argument and set the Object name
datafile = ARGS[1]
Identity = ARGS[2]

if Identity == "C"
	Object_Name = datafile[38:end-4]
	Output_Text_Location = "/home/clc304/Output/Oscillatory_CAR1_Text_Output/Charisis/" * Object_Name * ".txt"
	Output_Post_Location = "/home/clc304/Output/Oscillatory_CAR1_Post_Output/Charisis/" * Object_Name * ".dat"
else
	Object_Name = datafile[37:end-4]
	if Identity == "R"
		Output_Text_Location = "/home/clc304/Output/Oscillatory_CAR1_Text_Output/Randoms/" * Object_Name * ".txt"
		Output_Post_Location = "/home/clc304/Output/Oscillatory_CAR1_Post_Output/Randoms/" * Object_Name * ".dat"
	elseif Identity == "G"
		Output_Text_Location = "/home/clc304/Output/Oscillatory_CAR1_Text_Output/Grahams/" * Object_Name * ".txt"
		Output_Post_Location = "/home/clc304/Output/Oscillatory_CAR1_Post_Output/Grahams/" * Object_Name * ".dat"
	end

end
#Set the settings for the MCMC
p = 1
q = 0
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

#Define the posterior
post = Kalman.AR1KalmanPosterior(ts, ys, dys)

#Create the nested state
state = EnsembleNest.NestState(x -> Kalman.log_likelihood(post, x), x -> Kalman.log_prior(post, x), Kalman.init(post, nlive), nmcmc)

#Run the nested sampling procedure
EnsembleNest.run!(state, 0.1, false)

#Define a posterior sample
postsamples, lnprobs = EnsembleNest.postsample(state)

#Save the posterior sample
open(stream -> serialize(stream, (post, postsamples)), Output_Post_Location, "w")

#Call the logZ function
loggedZ = EnsembleNest.logZ(state)
loggedZ = round(loggedZ, 2)

#Writing text
Writing_text_A = "$(Object_Name);$(loggedZ)"

#Write the text to the files
open(Output_Text_Location, "w") do f
    write(f, Writing_text_A)
end

#Form a print text to the output file
Print_Text = "Completed task for object " * Object_Name
println(Print_Text)
