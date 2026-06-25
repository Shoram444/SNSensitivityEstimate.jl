import Base: *, maximum, println

"""
    Holds information about sensitivity estimate
    Fields:
    - `roi`: NamedTuple of the ROI
    - `tHalf`: sesnitivity estimate
    - `signalEff`: Signal efficiency
    - `bkgCounts`: Background counts
"""
struct SensitivityEstimateND{T, R, M}
    roi::R
    tHalf::T
    signalEff::T
    bkgCounts::M
end

Base.println(s::SensitivityEstimateND) = println(
    "Sensitivity estimate:
     ROI:   $(s.roi)
     T12 ≥ $(s.tHalf) yr
     ε = $(s.signalEff)
     b = $(s.bkgCounts)
    "
)

"""
    Holds information about ROI efficiency
"""
mutable struct ROIEfficiencyND{T, M}
    roi::Vector{T}
    eff::M
end

*(x::Float64, e::ROIEfficiencyND) = x * e.eff

mutable struct DataProcessND{T, B, NT, S, I} <: AbstractProcess 
    # data::LazyTree
    data::Vector{UnROOT.LazyEvent}
    isotopeName::S
    signal::B
    activity::Union{T, Measurement}
    timeMeas::T
    nTotalSim::T
    bins::NT
    amount::T
    varNames::Vector{S}
    varIdxs::Vector{I}
end



function get_roi_effciencyND(
    # data::LazyTree, 
    data::Vector{UnROOT.LazyEvent}, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys = keys(roi)
)
    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        evt = data[i]
        if passes_roi( evt, roi, roiKeys)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    ROI = Float64[]
    for k in keys(roi)
        push!(ROI, roi[k][1])
        push!(ROI, roi[k][2])
    end

    return ROIEfficiencyND(
        ROI,  
        count[] / nTotalSim  # Extract final atomic count
    )
end

function get_roi_effciencyND(
    data::Vector{UnROOT.LazyEvent}, 
    # data::LazyTree, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys::Vector{String}
)


    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi( data[i], roi, roiKeys)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    ROI = Float64[]
    for k in keys(roi)
        push!(ROI, roi[k][1])
        push!(ROI, roi[k][2])
    end

    return ROIEfficiencyND(
        ROI,  
        count[] / nTotalSim  # Extract final atomic count
    )
end

function get_roi_effciencyND(
    data::Vector{UnROOT.LazyEvent}, 
    # data::LazyTree, 
    roi::Vector{<:Real}, 
    nTotalSim::Real,
    varIdxs::Vector{Int}
)
    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi(data[i], roi, varIdxs)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end
    eff = count[] / nTotalSim  # Extract final atomic count
    return ROIEfficiencyND(
        roi,  
        eff  # Extract final atomic count
    )
end


function get_roi_effciencyND(
    process::DataProcessND,
    roi::NamedTuple
)
    return get_roi_effciencyND(
        process.data,
        roi,
        process.nTotalSim,
        process.varNames
    )
end

function get_roi_effciencyND(
    process::DataProcessND,
    roi::Vector{<:Real}
)
    return get_roi_effciencyND(
        getproperty(process, :data),
        roi,
        getproperty(process, :nTotalSim),
        getproperty(process, :varIdxs)
    )
end

@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::Vector{<:Real},
    varIdxs::Vector{Int}
)
    @inbounds for (r, i) in zip(1:2:length(roi)*2-1, varIdxs)
        if !(roi[r] ≤ abs(data[i]) < roi[r+1])
            return false
        end
    end
    return true
end


@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::Vector{<:Real},
    varNames::Vector{Symbol}
)
    i = 0
    @inbounds for (i,d) in zip(1:2:length(roi), data)
        if !(roi[i] ≤ abs(d) < roi[i+1])
            return false
        end
        i += 1
    end
    i == length(varNames) && return false
    return true
end


@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::NamedTuple,
    varNames::Vector{String}
)
    i = 1
    @inbounds for n in varNames
        range = roi[Symbol(n)]
        d = data[Symbol(n)] 
        if !(range[1] ≤ abs(d) < range[2])
            return false
        end
        i += 1
    end
    i == length(varNames) && return false
    return true
end

function DataProcessND(
    data::LazyTree,
    binsTuple::NamedTuple,
    varNames::Vector{String},
    processDict::Dict
)
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict
    println("creating process: $isotopeName")

    collected_data = collect(data)
    varIdxs = [findfirst(x -> x == n, names(data)) for n in varNames]

    T = promote_type(typeof(activity), typeof(timeMeas), typeof(nTotalSim), typeof(amount))
    I = promote_type(Int64, typeof(varIdxs[1]))

    return DataProcessND{T, Bool, typeof(binsTuple), String, I}(
        collected_data,
        isotopeName,
        signal,
        activity,
        timeMeas,
        nTotalSim,
        binsTuple,
        amount,
        varNames,
        varIdxs
    )
end

function DataProcessND(
    data::LazyTree,
    isotopeName::String, 
    signal::Bool, 
    activity::Real, 
    timeMeas::Real, 
    nTotalSim::Real, 
    binsTuple::NamedTuple, 
    amount::Real,
    varNames::Vector{String}
    )
    println("creating process: $isotopeName")

    collected_data = collect(data)
    varIdxs = [findfirst(x -> x == n, names(data)) for n in varNames]

    T = promote_type(typeof(activity), typeof(timeMeas), typeof(nTotalSim), typeof(amount))
    I = promote_type(Int64, typeof(varIdxs[1]))

    return DataProcessND{T, Bool, typeof(binsTuple), String, I}(
        collected_data,
        isotopeName,
        signal,
        activity,
        timeMeas,
        nTotalSim,
        binsTuple,
        amount,
        varNames,
        varIdxs
    )
end

function get_roi_bkg_counts(
    processes::Vector{<:DataProcessND}, 
    roi::NamedTuple,
)
    bkg = zero(Float64)
    ε = zero(Float64)

    for p in processes
        p.signal && continue
        ε = get_roi_effciencyND(p, roi).eff
        bkg += p.amount * ε * p.activity * p.timeMeas
    end
    return bkg
end

function get_roi_bkg_counts(
    processes::Vector{<:DataProcessND}, 
    roi::Vector{<:Real},
)
    bkg = zero(Float64)

    for p in processes
        p.signal && continue
        bkg += p.amount * get_roi_effciencyND(p, roi).eff * p.activity * p.timeMeas
    end
    return bkg
end

function get_roi_signal_efficiency(
    processes::Vector{<:DataProcessND}, 
    roi::Vector{<:Real},
)
    ε = zero(Float64)

    for p in processes
        p.signal == false && continue
        ε += get_roi_effciencyND(p, roi).eff
    end
    return ε
end

function get_s_to_b(
    SNparams::Dict, 
    α::Float64, 
    processes::Vector{<:DataProcessND}, 
    roi::Vector{Float64};
    approximate="table"
)
    for i in 1:2:length(roi)-1
        if roi[i] >= roi[i+1]
            return -1e20  # penalty for invalid ROI
        end
    end

    signal_id = findfirst(p -> p.signal, processes)
    if signal_id === nothing
        @error "No signal process found!"
        return -1e20 # penalty for no signal process
    end
                    
    ε = get_roi_signal_efficiency(processes, roi)
    ε == 0.0 && return -1e20 # If efficiency is zero, penalize optimization  

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return Measurements.value(ε / S_b)
end

function get_sensitivityND(
    SNparams::Dict, 
    α::Real, 
    processes::Vector{<:DataProcessND}, 
    roi::NamedTuple;
    approximate="table",
    add_mock_bkg =0.0
)

    #check that rois are within the range of the data
    for k in processes[1].varNames
        if (roi[Symbol(k)][1] < processes[1].bins[Symbol(k)][1] || roi[Symbol(k)][2] > processes[1].bins[Symbol(k)][2])
            @error "ROI out of range for $k"
        end
    end

    signal_id = findall([p.signal for p in processes])
    if length(signal_id) > 1 
        @error "Only one signal process allowed! Provided $(length(signal_id))" 
        return SensitivityEstimateND(roi, 0.0, 0.0, 0.0)
    end

    signal_process = processes[first(signal_id)]
    ε = get_roi_effciencyND(signal_process, roi).eff
    ε == 0 && return SensitivityEstimateND(roi, 0.0, 0.0, 0.0) # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi) + add_mock_bkg

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
    S_b = get_FC(b, α; approximate=approximate)

    tHalf = constantTerm * ε / S_b

    # return tHalf
    return SensitivityEstimateND(roi, tHalf, ε, b)
end


function set_activity!(process::DataProcessND, activity::Real)
    process.activity = activity
    return process
end

function set_timeMeas!(process::DataProcessND, timeMeas::Real)
    process.timeMeas = timeMeas
    return process
end

function set_nTotalSim!(process::DataProcessND, nTotalSim::Real)
    process.nTotalSim = nTotalSim
    return process
end


function set_amount!(process::DataProcessND, amount::Real)
    process.amount = amount
    return process
end


function make_stepRange(process::DataProcessND)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


# =====================================================================
#  Stepped (discrete-grid) search support
#  Instead of optimizing physical ROI edges directly, we optimize integer
#  step-indices `j ∈ [0, nsteps]` per variable. A proposed (continuous)
#  index is rounded to the nearest grid node and decoded to the physical
#  edge `min + j*step`. This discretizes the search space and removes the
#  sub-resolution moves that make the objective look flat to the optimizer.
# =====================================================================

"""
    GridSpec

Describes the discrete search grid over the ROI variables of a process.

Fields (all aligned with `keys(process.bins)`):
- `keys`   : variable names (`Vector{Symbol}`)
- `mins`   : lower edge of each variable's range
- `steps`  : step size of each variable
- `nsteps` : number of steps so that `mins + nsteps*steps ≤ max`
"""
struct GridSpec
    keys::Vector{Symbol}
    mins::Vector{Float64}
    steps::Vector{Float64}
    nsteps::Vector{Int}
end

# Resolve the step size for variable `k` from the user-provided `steps`,
# which may be a scalar (same step for all), a NamedTuple, or a Dict.
_get_step(steps::Real, k) = float(steps)
_get_step(steps::NamedTuple, k) = float(steps[k])
_get_step(steps::AbstractDict, k) = float(haskey(steps, k) ? steps[k] : steps[Symbol(k)])

"""
    GridSpec(process::DataProcessND, steps)

Build a `GridSpec` from a process and a `steps` specification. `steps` can be
a single number applied to every variable, or a `NamedTuple`/`Dict` keyed by
variable name giving a per-variable step size, e.g. `(sumE = 100, phi = 10)`.
"""
function GridSpec(process::DataProcessND, steps)
    ks = collect(keys(process.bins))
    mins   = Float64[]
    stepv  = Float64[]
    nsteps = Int[]
    for k in ks
        lo, hi = process.bins[k]
        s = _get_step(steps, k)
        s > 0 || error("Step size for $k must be positive, got $s")
        push!(mins, lo)
        push!(stepv, s)
        push!(nsteps, floor(Int, (hi - lo) / s))  # floor => decoded edge never exceeds hi
    end
    return GridSpec(ks, mins, stepv, nsteps)
end

"""
    make_stepRange(process::DataProcessND, steps) -> Vector{Tuple{Int,Int}}

Grid version: returns the integer index bounds `(0, nsteps)` (lower and upper
edge) for each variable, to be used as optimization bounds together with the
`GridSpec`. Use `GridSpec(process, steps)` to build the decoder.
"""
function make_stepRange(process::DataProcessND, steps)
    grid = GridSpec(process, steps)
    return make_stepRange(grid)
end

function make_stepRange(grid::GridSpec)
    stepRange = Tuple{Int64, Int64}[]
    for n in grid.nsteps
        push!(stepRange, (0, n))  # lower edge index
        push!(stepRange, (0, n))  # upper edge index
    end
    return stepRange
end

"""
    decode_grid_roi(indices, grid::GridSpec) -> Vector{Float64}

Convert a vector of (possibly continuous) step-indices
`[lo₁, hi₁, lo₂, hi₂, …]` into physical ROI edges
`[min₁+lo₁*s₁, min₁+hi₁*s₁, …]`. Indices are rounded to the nearest grid node.
"""
function decode_grid_roi(indices::AbstractVector, grid::GridSpec)
    roi = Vector{Float64}(undef, length(indices))
    @inbounds for (j, vi) in enumerate(1:2:length(indices)-1)
        idx_lo = clamp(round(Int, indices[vi]),   0, grid.nsteps[j])
        idx_hi = clamp(round(Int, indices[vi+1]), 0, grid.nsteps[j])
        roi[vi]   = grid.mins[j] + idx_lo * grid.steps[j]
        roi[vi+1] = grid.mins[j] + idx_hi * grid.steps[j]
    end
    return roi
end

"""
    value_to_index(value, min, step; nsteps=nothing) -> Int

Convert a single physical `value` into the nearest grid step-index given the
variable's lower edge `min` and `step` size, i.e. `round((value - min) / step)`.
If `nsteps` is provided the result is clamped to `[0, nsteps]`.
"""
function value_to_index(value::Real, min::Real, step::Real; nsteps=nothing)
    idx = round(Int, (value - min) / step)
    return nsteps === nothing ? idx : clamp(idx, 0, nsteps)
end

"""
    encode_grid_roi(roi, grid::GridSpec) -> Vector{Int}

Inverse of [`decode_grid_roi`](@ref): convert physical ROI edges into integer
grid step-indices. `roi` may be either

- a flat vector of edges `[lo₁, hi₁, lo₂, hi₂, …]` (same order as `grid.keys`), or
- a `NamedTuple` keyed by variable, e.g. `(sumE = (2700, 3100), phi = (0, 180))`.

Indices are clamped to each variable's `[0, nsteps]` range so the result is a
valid starting point (`x0`) for the optimizer.
"""
function encode_grid_roi(roi::AbstractVector, grid::GridSpec)
    indices = Vector{Int}(undef, length(roi))
    @inbounds for (j, vi) in enumerate(1:2:length(roi)-1)
        indices[vi]   = value_to_index(roi[vi],   grid.mins[j], grid.steps[j]; nsteps=grid.nsteps[j])
        indices[vi+1] = value_to_index(roi[vi+1], grid.mins[j], grid.steps[j]; nsteps=grid.nsteps[j])
    end
    return indices
end

function encode_grid_roi(roi::NamedTuple, grid::GridSpec)
    indices = Int[]
    for (j, k) in enumerate(grid.keys)
        lo, hi = roi[k]
        push!(indices, value_to_index(lo, grid.mins[j], grid.steps[j]; nsteps=grid.nsteps[j]))
        push!(indices, value_to_index(hi, grid.mins[j], grid.steps[j]; nsteps=grid.nsteps[j]))
    end
    return indices
end

"""
    get_s_to_b(SNparams, α, processes, indices, grid::GridSpec; approximate)

Grid-aware objective: decodes integer step-`indices` to a physical ROI via
`grid`, then evaluates the standard `get_s_to_b`. This is the function to wrap
for `Metaheuristics.optimize` when doing a stepped search.
"""
function get_s_to_b(
    SNparams::Dict,
    α::Real,
    processes::Vector{<:DataProcessND},
    indices::AbstractVector,
    grid::GridSpec;
    approximate="formula"
)
    roi = decode_grid_roi(indices, grid)
    return get_s_to_b(SNparams, float(α), processes, roi; approximate=approximate)
end


function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (best[i], best[i+1]) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => (best[i], best[i+1]) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

"""
    get_best_ROI_ND(res::Vector, process, grid::GridSpec) -> NamedTuple

Grid version: decodes the optimizer's best step-indices through `grid` into
physical ROI edges before building the ROI `NamedTuple`. Pass the minimizer
vector, e.g. `get_best_ROI_ND(minimizer(result), process, grid)`.
"""
function get_best_ROI_ND(res::Vector{<:Real}, process, grid::GridSpec)
    best = decode_grid_roi(res, grid)
    return get_best_ROI_ND(best, process)
end

# Obtain the background counts histogram for a given ROI and process and variable 
"""
get_roi_bkg_counts_hist(
    p::DataProcessND, 
    roi::NamedTuple,
    bins,
    var
) -> Hist1D
 
Returns a histogram of background counts for the specified ROI and variable.
    Inputs:
    - `p::DataProcessND`: The data process containing the data and parameters.
    - `roi::NamedTuple`: The region of interest defined as a NamedTuple.
    - `bins`: The bin edges for the histogram.
    - `var`: The variable name for which the histogram is computed.
"""
function get_roi_bkg_counts_hist(
    p::DataProcessND, 
    roi::NamedTuple,
    bins,
    var
)
    d = getproperty.(p.data, var)
    bkg_hist = normalize(Hist1D(d;binedges = bins); width = false)

    ε = get_roi_effciencyND(p, roi).eff
    # @show p.isotopeName, ε

    bkg_hist.bincounts .= bkg_hist.bincounts .* p.amount * ε * p.activity * p.timeMeas 
    return bkg_hist
end


import Base.print
function print(best_roi::SensitivityEstimateND)
    println("tHalf: $(best_roi.tHalf)")
    println("signalEff: $(best_roi.signalEff)")
    println("bkg count: $(best_roi.bkgCounts)")
    println("best ROI:")
    for (k,v) in pairs(best_roi.roi)
        println("  $k : $(v)")
    end
end