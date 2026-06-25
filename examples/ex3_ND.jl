#############################################
#############################################
##
##  LOAD DEPENDENCIES
##
#############################################
#############################################
import Pkg; Pkg.activate("examples")
using SNSensitivityEstimate, UnROOT, Metaheuristics, CSV, DataFrames

const α= 1.64485362695147 # 90% C.L. constant

#############################################
## LOAD SIMULATION Data 
#############################################

# here of course put your own path, maybe absolute rather than relative
# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read("examples/data/data_info.csv", DataFrame)
file_paths = joinpath.("examples/data", data_info.file_name) 
data_files = [UnROOT.ROOTFile(file_path) for file_path in file_paths] 
tree_name = "tree" 
variables = keys(data_files[1][tree_name]) 
data_tables = [UnROOT.LazyTree(data_file, tree_name, variables) for data_file in data_files] 


# STEP 1: define variables and their ranges that we are going to work with

# First we define variables we want to cut
var_names = [ 
    "sumE", 
    "dy", 
    "dz",
    "lPint",
    "lPext",
    ]

# Next we define the boundaries of the variables in the search space. These boundaries will be used to define the search space for the optimization algorithm.
var_bounds = (
    sumE = (300, 3500),
    dy = (0, 200),
    dz = (0, 200),
    lPint = (0., 5.0),
    lPext = (0., 110.0),
)

# STEPPED SEARCH: discrete step size per variable. The optimizer searches over
# integer step-indices `0 .. floor((max-min)/step)` for each ROI edge, which are
# decoded back to physical edges via a `GridSpec`. This discretizes the search
# space and removes the sub-resolution moves that make the objective look flat.
var_steps = (
    sumE  = 100.0,  # keV
    dy    = 5.0,    # mm
    dz    = 5.0,    # mm
    lPint = 0.1,
    lPext = 0.5,
)

selected_tables = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files] 


# STEP 2: Create DataProcessND objects
# We do not filter data here, because that is the role of the optimization process. 
# In the ND example we no longer work with DataProcess object, but rather DataProcessND objects, which are designed to work with N-dimensional data vectors.

processes_nd = [
    SNSensitivityEstimate.DataProcessND(
        selected_tables[i], # Here we pass the whole data-set
        String(data_info.process[i]), 
        data_info.is_signal[i], 
        data_info.activity[i], 
        data_info.time_s[i], 
        data_info.n_sim[i], 
        var_bounds, # bounds for each variable
        data_info.amount[i], 
        var_names, # variable names
    ) for i in 1:length(selected_tables)
]

# Identify the single signal process (used to build the search grid).
signal_process = processes_nd[findfirst(p -> p.signal, processes_nd)]

# Test simple roi, you can play around by hand to get the feeling.
# NOTE: the ROI must specify a window for EVERY variable in `var_names`.
my_roi = (
    sumE = (2000, 3100),
    dy = (0, 200),
    dz = (0, 200),
    lPint = (0.0, 50.0),
    lPext = (0.0, 50.0),
)

# there are 2 approximation methods "table" and "formula", 
# table is more precise but takes longer to calculate, 
# formula is fine for quick tests
my_sensitivity = get_sensitivityND(SNparams, α, processes_nd, my_roi; approximate="table") 


# print your result, it shows the sensitivity, efficiency, background and roi
print(my_sensitivity) 



##############################################
## STEPPED (GRID) OPTIMIZATION BY METAHEURISTICS.jl
##############################################

# Build the grid decoder from the signal process + per-variable step sizes.
# The grid extent comes from `var_bounds`; the resolution from `var_steps`.
grid = GridSpec(signal_process, var_steps)

# The objective: instead of optimizing physical ROI edges, we optimize integer
# step-indices. `get_s_to_b(..., grid)` decodes the proposed indices into a
# physical ROI before evaluating signal-efficiency / FC-background.
# We negate because Metaheuristics.jl MINIMIZES (we want to maximize s/b).
prob(proposed_indices) = - get_s_to_b(SNparams, α, processes_nd, proposed_indices, grid; approximate="formula")


# The search range is now in INDEX units: (0, nsteps) for each ROI edge.
# `make_stepRange(grid)` returns these integer bounds automatically.
searchRange = make_stepRange(grid)


# Lower and upper bounds for the optimizer, derived from the index search range.
lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float


##############################################
# Now we can set up the optimization problem and run it
##############################################

# first the optimization options, basically the hyperparameters of the optimization algorithm
# you can leave them empty - default values will be used,
# or you can try to tune yourself to see what the results give, basically you want to find a good balance between time and quality of results
# For sure set up at least a time_limit and verbosity to true to see progress
# The full list of parameters can be found in https://jmejia8.github.io/Metaheuristics.jl/stable/api/#Metaheuristics.Options
options = Options(;
    # x_tol = 1.0, # tolerance in position space
    # f_tol = 1e-5, # tolerance in function value
    time_limit = 60*20.0, # time limit in seconds, must include .0 to be a Float64 value, doesn't take integers
    verbose = true, # verbosity on
    # iterations = 15, # number of iterations / generations
)

# Set up boxconstraints bounds - just the lower and upper bounds we defined before
bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

# Initial guess for the ROI, IN INDEX UNITS (not physical units).
# Convert a physical edge to an index with: index = (edge - min) / step.
# Order matches `var_names`: sumE, dy, dz, lPint, lPext.
# e.g. sumE in (2700, 3100) with min=300, step=100 -> indices (24, 28).
x0 = float.([
    value_to_index(2700, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]), value_to_index(3100, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]),  # sumE  -> (2700, 3100)
    value_to_index(0, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]), value_to_index(100, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]),   # dy    -> (0, 200)
    value_to_index(0, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]), value_to_index(140, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]),   # dz    -> (0, 200)
    value_to_index(0.0, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]), value_to_index(4.0, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]),   # lPint -> (0.0, 1.0)
    value_to_index(2.0, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]), value_to_index(110.0, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]),   # lPext -> (0.0, 1.0)
])


# Set up optimizer algorithm, I generally use ECA - supports parallel evaluation, PSO - also can use parallel, SA - single core, but pretty good
algo = PSO(;options)

# Set the initial guess for the optimizer
set_user_solutions!(algo, x0, prob)


# Now we run the optimization, depending on how much time you gave it in the options, this can take a while
result = optimize(prob, bounds, algo)

# We can show the best result found
@show minimum(result)

# and store its (index-space) value
@show res = minimizer(result)


# Decode the best step-indices back to a physical ROI (note the extra `grid`
# argument), then calculate the precise sensitivity with the "table" method.
# `res = minimizer(result)` is the best index-vector found above.
best = get_best_ROI_ND(res, signal_process, grid)
best_sens = get_sensitivityND(SNparams, α, processes_nd, best; approximate="table")

# finally print the result
println("\n===== BEST RESULT =====")
print(best_sens)
