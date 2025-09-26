#####################################################
#####################################################
#####################################################
####                                            #####   
####  EXAMPLE OF USING SNSENSITIVITY ESTIMATE   #####
####  FOR A SIMPLE 1D USE CASE                  #####
####                                            #####       
#####################################################
#####################################################
#####################################################

# First we load the necessary packages, if they are not installed, do so with:
# using Pkg; Pkg.add("PackageName"); using PackageName
using SNSensitivityEstimate
using BlackBoxOptim
using Distributions
using StatsPlots
using FHist

# first we generate mock data
# we assume gaussian signal process and exponential background
# we assume a dataset of 100 signal events, nS and 1000 background nB 
# since the estimate relies on the efficiency of signal detection 
# (which is calculated as the number of events that pass analysis counts
# np over the total number of simulated events ns; ε = np/ns, we must
# make up some number for ns > nS.

nS = 100_000 # n signal events (after data-cuts)
nB = 1000 # n backgoround events (after data-cuts)
ns = 10_000_000 # mock simulated number of signal events (before datacuts) 
nb = 10_000_000 # mock simulated numbed of background events (before datacuts)

mu_signal = 5 # signal mean value
sigma_signal = 1 # signal sigma
lambda_background = 3 # decay constant of background

n_signal_dist = rand(Normal(mu_signal, sigma_signal), nS)
n_background_dist = rand(Exponential(lambda_background), nB)

# plot results
histogram(n_signal_dist, nbins = 0:0.1:10, label = "signal", st = :stephist, lw = 4, normalize = :pdf)
histogram!(n_background_dist, nbins = 0:0.1:10, label = "background", st = :stephist, lw = 4, normalize = :pdf)

# Next we iniate the DataProcess objects, the core of package
# an DataProcess relies on a number of parameters to be setup, 
# you can look them up with ?DataProcess in REPL:
# The following fields are defined:
#
#    + dataVector::Vector{<:Real} - vector of initial data points
#    + isotopeName::String - isotope name 
#    + signal::Bool - signal or background
#    + activity::Real - activity of the given process in units [Bq]
#    + timeMeas::Real - time of measurement in units [s]
#    + nTotalSim::Real - total number of originally simulated events from Falaise
#    + bins::AbstractRange - binning to be used in format (min:step:max)
#    + amount::Real - mass [kg] or volume [l / m³] of the object where the isotope is present: i.e. source foil
#    + efficiency::Hist2D - the 2D histogram with efficiencies per bin 


signal_isotope_name = "generic_signal"
background_isotope_name = "generic_bkg"

background_activity = 0.01 # mock value
timeMeasYear = 3 # 3 years of measurement
timeMeasSeconds = timeMeasYear*365*24*3600 # 3 years of measurement
bins = 0:0.1:10 # binning
signal_amount = 6 #kg of isotope
background_amount = 100 #kg of background source

signalProcess = DataProcess(
    n_signal_dist,
    signal_isotope_name,
    true,
    1,
    timeMeasSeconds,
    ns,
    bins,
    signal_amount
)

backgroundProcess = DataProcess(
    n_background_dist,
    background_isotope_name,
    false,
    background_activity,
    timeMeasSeconds,
    nb,
    bins,
    background_amount
)

allProcesses = vcat(signalProcess, backgroundProcess)

# Now that we have the signal and background processes, the calculation is super simple
α = 1.65 # ~90% confidence level
W = 0.08192 # molar mass in [kg/mol]
Na = 6.02214e23 # Avogadro's number

# create a 2D map of ROI and tHalf values in each bin
tHalf_map = get_tHalf_map(W, signal_amount, Na, timeMeasYear, α, allProcesses...; approximate="formula")
plot(tHalf_map, title = "tHalf map", xlabel = "min ROI", ylabel = "max ROI", label = "", colorbartitle = "\nsensitivity (yr)",lw = 2, colormap = :viridis, size = (1000, 600), margin = 9Plots.mm, thickness_scaling = 1.5)

# The best sensitivity is found at the point where tHalf_map is maximized
best_tHalf = get_max_bin(tHalf_map)

# We can also extract the efficiency in ROI and background counts
eff = lookup(signalProcess, best_tHalf[:minBinEdge], best_tHalf[:maxBinEdge])
bkg_counts = get_bkg_counts_ROI(best_tHalf, allProcesses...)

println(
    "best ROI found at: ($(best_tHalf[:minBinEdge]) - $(best_tHalf[:maxBinEdge])) a.u.
     with tHalf = $(best_tHalf[:maxBinCount])
     efficiency = $(eff)
     background counts = $(bkg_counts)"
)

