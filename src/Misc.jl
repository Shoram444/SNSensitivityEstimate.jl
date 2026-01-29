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


function generate_raw_plots(inDf::DataFrame, isotope; kwargs...)
    phi, e1, e2 = inDf[!, :phi], inDf[!, :reconstructedEnergy1], inDf[!, :reconstructedEnergy2]

    plot_args = Dict(
        :left_margin => 12Plots.mm,
        :right_margin => 12Plots.mm,
        :top_margin => 12Plots.mm,
        :bottom_margin => 12Plots.mm,
        :size => (1200, 900),
        :dpi => 200,
        :titlefontsize => 24,
        :guidefontsize => 22,
        :tickfontsize => 20,
        :legendfontsize => 24,
        :formatter => :auto,
        :legend => :best,
        :lw => 4,
        :thickness_scaling => 1.1
    )

    hPhi = stephist(
        phi;
        nbins=0:5:180,
        xlabel="escape angle " * L"[\degree]",
        ylabel="counts/" * L"5\degree",
        title="angular distribution of $isotope",
        label="",
        plot_args,
        kwargs...
    )

    hEne = stephist(
        e1;
        nbins=0:100:3500,
        xlabel="single-electron energy [keV]",
        ylabel="counts/" * "100 keV",
        title="single electron energy distribution of $isotope",
        label="energy 1",
        plot_args,
    )

    stephist!(
        e2;
        nbins=0:100:3500,
        label="energy 2",
        plot_args,
    )

    stephist!(
        vcat(e1, e2);
        nbins=0:100:3500,
        label="both electrons",
        plot_args,
    )

    safesave(plotsdir("Raw", "Angular", "$(isotope)_raw_angdist.png"), hPhi)
    safesave(plotsdir("Raw", "Single", "$(isotope)_raw_enedist.png"), hEne)

    plot(hPhi, hEne, size=(1800, 600))

end

function fill_from_root_file(inFile::ROOTFile, treeName::String, fieldNames; roi = nothing)
    data = LazyTree(inFile, treeName, [fieldNames...]) 

    if roi !== nothing
        for (field, (minVal, maxVal)) in pairs(roi)
            data = filter(x -> getproperty(x, field) > minVal && getproperty(x, field) < maxVal, data)
        end
    end
    df = DataFrame(data)

    for k in fieldNames
        filter!(string(k) => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df) # filter rows with missing values
    end
    return df
end


function get_2D_vertex_separation(y1, z1, y2, z2)
    return sqrt((y1 - y2)^2 + (z1 - z2)^2)
end

function get_1D_vertex_separation(x1, x2)
    return abs(x1 - x2)
end

function add_vertex_2D_separation_column!(df)
    @rtransform! df :d = get_2D_vertex_separation(
        :y1Escaped,
        :z1Escaped,
        :y2Escaped,
        :z2Escaped
    )
end

function add_vertex_dy_separation_column!(df)
    @rtransform! df :dy = get_1D_vertex_separation(
        :y1Escaped,
        :y2Escaped
    )
end

function add_vertex_dz_separation_column!(df)
    @rtransform! df :dz = get_1D_vertex_separation(
        :z1Escaped,
        :z2Escaped
    )
end


function generate_pseudo_data(h1d::Hist1D)
	n_i = map(x->rand(Poisson(x)), bincounts(h1d))
	b_i = collect(binedges(h1d))
	bin_widths = diff(b_i) # vector of bin widths (in case of non-uniform binning)
	
	data = Vector{Float64}(undef, sum(n_i))
	idx_slice = 1
	for (ni, bi, bw_i) in zip(n_i, b_i, bin_widths)
		if(ni != 0)
			data[idx_slice:idx_slice+ni-1] = rand(Uniform(bi, bi + bw_i ), ni)
			idx_slice+=ni
		end
	end
	return data 
end

get_sigma_MeV(E_MeV, FWHM) = 1/(2*sqrt(2*log(2)))*FWHM*sqrt(1/E_MeV)*E_MeV
get_sigma_keV(E_keV, FWHM) = get_sigma_MeV(E_keV/1000.0, FWHM) * 1000.0



"""
    get_sensitivities_vs_time(signal, background, SNparams; neutron_bkg = 0.0, effFactor = 1.0)
Calculates the sensitivities of a given signal against a background over time, based on the provided SN parameters.

`signal` and `background` are the signal and background processes, respectively.
`SNparams` is a dictionary containing the parameters for the supernova sensitivity calculation.
`neutron_bkg` is the neutron background, defaulting to 0.0.
`effFactor` is a factor to adjust the efficiency, defaulting to 1.0

"""
function get_sensitivities_vs_time(
        signal,
        background,
        SNparams;
        neutron_bkg = 0.0,
        effFactor = 1.0
    )
    t = range(0, 5, 100)
    sensitivities = []
    
    t12(t, e, b) = get_tHalf(
        SNparams["W"],
        SNparams["foilMass"],
        SNparams["Nₐ"],
        t,
        SNparams["a"],
        e*effFactor,
        (b+neutron_bkg)/ SNparams["tYear"] * t,
        α;
        approximate="table"
    )
    t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
    best_t12ESum = get_max_bin(t12MapESum)
    expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
    effbb = lookup(signal, best_t12ESum)
    append!(sensitivities, t12.(t, effbb,expBkgESum))
    return sensitivities
end