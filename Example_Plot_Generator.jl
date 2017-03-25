using PyCall
using PyPlot;
plt = PyPlot;
using CARMA
using Ensemble
using DataFrames
using Distributions

"""-----------------------------------------------------------------------------
This script is designed to create fake plots for the puposes of explaining the
CARMA modelling process. Specifically, highlighting peaks in the PSD for a
sinusoidal process, with this peak changing location as the period of the
oscillations changes. These fake outputs are created using the 'generate'
function, in Kalman.jl

It takes no inputs, however it does generate light curves and PSDs for the
generated data. 
"""-----------------------------------------------------------------------------



@pyimport seaborn as sns

sns.set_context("notebook")
sns.set_style("ticks")
sns.set_palette("colorblind")
#-----------------------------------------------------------------------------------------------------------
function Generate_LC(N, ts, dys, Output_LC_Location, Description, Output_PSD_Location, Period, Q, Real_Root)

  ts = Array(ts)
  dys = Array(dys)

  if Description == "Damped_Random_Walk"
    post = Kalman.AR1KalmanPosterior(ts, zeros(N), dys)

    mu = 20.0
    sigma = 0.7
    nu = 0.66
    tau = 50
    x = [mu, sigma, nu, tau]
    p = Kalman.to_params(post, x)

    T = maximum(post.ts) - minimum(post.ts)
    dt_med = median(diff(post.ts))

    fmax = 1.0/(2.0*dt_med)
    df = 1.0/T

    fs = collect(df:df:fmax)
    Filter = Kalman.make_filter(post, p)

    #Make Psd plot
    loglog(fs, Kalman.psd(Filter, fs))
    y = Kalman.generate(Filter, ts, dys)

  else
    post = Kalman.CARMAKalmanPosterior(ts, zeros(N), dys, 3, 2)
    p = Kalman.to_params(post, zeros(Kalman.nparams(post)))
    p.mu = 20.0
    p.sigma = 0.7
    p.nu = 0.66
    p.maroots = -10*rand(2)

    Real_Of_Complex = -(1/Period)/Q

    if Description == "Sinusoid_Small_Q_Factor"
      p.arroots = [Real_Root, Real_Of_Complex + (2*pi)/Period*1im, Real_Of_Complex - (2*pi)/Period*1im]

    elseif Description == "Sinusoid_Large_Q_Factor"
      p.arroots = [Real_Root, Real_Of_Complex + (2*pi)/Period*1im, Real_Of_Complex - (2*pi)/Period*1im]

    elseif Description == "Test"
      p.arroots = [Real_Root, Real_Of_Complex + (2*pi)/Period*1im, Real_Of_Complex - (2*pi)/Period*1im]

    else
        println("Description Error")
        exit(1)
    end

    fs = Kalman.psdfreq(post)

    #Make Psd plot
    loglog(fs, Kalman.psd(post, p, fs))
    y, dy_true = Kalman.generate(post, p)
  end

  #Customize the plot
  xlabel("Frequency")
  ylabel("PSD")
  title("A Generated PSD Plot of " * replace(Description, "_", " "))
  ax = plt.gca()
  ax[:tick_params](axis="x",top="off")
  ax[:tick_params](axis="y",right="off")

  #Save the figure and close
  savefig(Output_PSD_Location)
  close("all")

  #Make and save Lightcurve
  errorbar(post.ts, y, post.dys, fmt=".")

  #Format the image
  xlabel("Times(days)")
  ylabel("Magnitude")
  title("Generated " * replace(Description, "_", " ") * " Light Curve")
  ax = plt.gca()
  ax[:invert_yaxis]()
  ax[:tick_params](axis="x",top="off")
  ax[:tick_params](axis="y",right="off")

  #Save and close the figure
  savefig(Output_LC_Location)
  close("all")

  post2 = Kalman.CARMAKalmanPosterior(ts, y, dys, 3, 2)
  post2
end

function Use_Predict(post, Output_Predict_Location, Description, Output_Post_Location)

    #Set parameters
    nlive = 1024
    nmcmc = 128

    #Create the nested state
    state = EnsembleNest.NestState(x -> Kalman.log_likelihood(post, x), x -> Kalman.log_prior(post, x), Kalman.init(post, nlive), nmcmc)

    #Run the nested sampling procedure
    EnsembleNest.run!(state, 0.1, false)

    #Define a posterior sample
    postsamples, lnprobs = EnsembleNest.postsample(state)

    open(stream -> serialize(stream, (post, postsamples)), Output_Post_Location, "w")

    #post, postsamples = open(deserialize, Output_Post_Location)

    #Predic procedure
    ts_predict = sort(vcat(post.ts, collect(linspace(post.ts[1], post.ts[end], 256))))
    prandom = postsamples[:,rand(1:size(postsamples,2))]
    ys, vys = Kalman.predict(post, prandom, ts_predict)
    dys = sqrt(vys)

    #Plot the errorbar
    errorbar(post.ts, post.ys, post.dys, fmt=".", color=sns.color_palette()[2])
    ax = plt.gca()
    ax[:invert_yaxis]()
    ax[:tick_params](axis="x",top="off")
    ax[:tick_params](axis="y",right="off")
    xlabel("Times(days)")
    ylabel("Magnitude")
    title("Generated " * replace(Description, "_", " ") * " Light Curve, with Predict fit")

    #Save the figure and fill in the errors for the predict function
    plot(ts_predict, ys, color=sns.color_palette()[3])
    fill_between(ts_predict, ys+dys, ys-dys, color=sns.color_palette()[3], alpha=0.3)
    savefig(Output_Predict_Location)
    close("all")

end
#-----------------------------------------------------------------------------------------------------------
#Set the User
User = "N"

#Input data Paths
if User == "N"
    Output_PSD_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/PSD_"
    Output_LC_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/LC_"
    Output_Predict_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/Predict_"
    Output_Post_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/"
elseif User == "C"
    Output_PSD_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Example_Graphs/Generate_Sinusoid_PDF_"
    Output_LC_Location = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Example_Graphs/Generate_Sinusoid_LC_"
    Output_Predict_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/Predict_"
    Output_Post_Location = "C:/Users/User/Documents/University/Year 4/Project/Tsunami_Results/Fake_Output/"
else
    println("User error")
    exit(1)
end

N = 250
ts = linspace(0,2990,N)
#ts = cumsum(rand(Exponential(5.0), N)) # Exponential waiting times => Poisson process
dys = ones(N) / 10

Period = 500
Real_Root = -10

for i in 2:2

    #Set the Description
    if i == 1
        Description = "Sinusoid_Small_Q_Factor"
        Q = 0.25
    elseif i == 2
        Description = "Sinusoid_Large_Q_Factor"
        Q = 5.0
    elseif i == 3
        Description = "Damped_Random_Walk"
        Q = 5.0
    elseif i == 4
        Description = "Test"
        Q = 1.0
    else
        println("Description Error")
        exit(1)
    end

    println(Description)

    Output_PSD_Loc = Output_PSD_Location * Description * ".jpg"
    Output_LC_Loc = Output_LC_Location * Description * ".jpg"
    Output_Predict_Loc = Output_Predict_Location * Description * ".jpg"
    Output_Post_Loc = Output_Post_Location * Description * ".dat"

    #Set Description and call function
    post = Generate_LC(N, ts, dys, Output_LC_Loc, Description, Output_PSD_Loc, Period, Q, Real_Root)

    #Use_Predict(post, Output_Predict_Loc, Description, Output_Post_Loc)
end
