#####################################################
#####################################################
#####################################################
####                                            #####   
####  EXAMPLE OF USING SNSENSITIVITY ESTIMATE   #####
####  FOR A SIMPLE 2D USE CASE                  #####
####                                            #####       
#####################################################
#####################################################
#####################################################

# First we load the necessary packages, if they are not installed, do so with:
# using Pkg; Pkg.add("PackageName"); using PackageName
using SNSensitivityEstimate
using Metaheuristics
using Distributions
using StatsPlots
using FHist


vars = [
    "phi", 
    "sumE", 
    ]

bins = (
    phi = (0,180),
    sumE = (300, 3500),
)

processes = load_ndim_processes("/home/maros/.julia/dev/SNSensitivityEstimate/exampels/data", bins, vars)

signal = get_process("bb0nu_foil_bulk", processes)[1]

background = [
    get_process("bb_foil_bulk", processes)[1],
    get_process("Bi214_PMT_glass_bulk", processes)[1],
]

set_signal!(background[1], false)

set_nTotalSim!( signal, 1e7 )

set_nTotalSim!( background[1], 1e7 )
set_nTotalSim!( background[2], 1e8 )

α= 1.64485362695147


prob(x) = - get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")

x0 = float.([0,100, 2500, 3000])
prob(x0)



function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal)


options = Options(;
    f_tol = 1e-5,
    f_tol_rel = 1e-5,
    f_tol_abs = 1e-5,
    time_limit = 60*3.0,
)


lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float


bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

result = Metaheuristics.optimize(prob, bounds, SA(;options))

@show res = minimizer(result)

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => (best[i], best[i+1]) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best = get_best_ROI_ND(res, signal)

best_sens = get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="table")


h1 = histogram2d(
    getproperty.(signal.data, :sumE),
    getproperty.(signal.data, :phi),
    bins = (0:100:3500,0:10:180),
    xlims = (0,3500),
    ylims = (0,180),
    xlabel = "E (keV)",
    ylabel = "φ (°)", 
    title = "Signal 2D distribution",
    c = :blues,
    margin = 5Plots.mm
)

h2 = histogram2d(
    getproperty.(background[1].data, :sumE),
    getproperty.(background[1].data, :phi),
    bins = (0:100:3500,0:10:180),
    xlims = (0,3500),
    ylims = (0,180),
    xlabel = "E (keV)",
    ylabel = "φ (°)", 
    title = "Background 2D distribution",
    c = :reds,
    margin = 5Plots.mm
)

plot(h1, h2, layout = (1,2), size = (900,400))