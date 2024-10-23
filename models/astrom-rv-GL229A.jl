using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
using DataFrames
using Dates
using OctofitterRadialVelocity

model_deps = [
    Octofitter,
    Distributions,
    CairoMakie,
    PairPlots,
    Pigeons,
    DataFrames,
    Dates,
    OctofitterRadialVelocity,
]

# Old measurements from Brandt 2020 
# The Astronomical Journal, 160:196 (15pp), 2020 October
cd(@__DIR__)
astrom_like_1 = PlanetRelAstromLikelihood(Table(;
    epoch=mjd.(DateTime.([
            "1995 Nov 17",
            "1996 May 25",
            "1996 Nov 09",
            "1999 May 26",
            "2000 May 26",
            "2000 Nov 16",
            "2011 Mar 26",
            "2020 Oct 24",
            "2021 Jan 5",
        ], dateformat"YY u d")),
    sep=[7777.0, 7732.0, 7687.7, 7458.3, 7362.8, 7316.9, 6210, 4922.1, 4890.5],
    σ_sep=[1.7, 2.0, 1.5, 1.6, 1.6, 1.6, 10, 2.3, 2.4],
    pa=deg2rad.([163.224, 163.456, 163.595, 164.796, 165.244, 165.469, 171.2, 179.564, 179.735]),
    σ_pa=deg2rad.([0.015, 0.019, 0.015, 0.015, 0.016, 0.016, 0.1, 0.024, 0.024]),
))
scatter(astrom_like_1.table.sep,astrom_like_1.table.pa)
##
# Brandt 2021 , References: Rosenthal et al. (2021); A. Howard, private communication
# rv_like_1 = StarAbsoluteRVLikelihood(Table(;
#     epoch=jd2mjd.(
#         [2458116.862, 2458117.852, 2458396.142, 2458777.037, 2458794.994, 2458880.798, 2458907.830, 2459101.122, 2459267.794],
#     ),
#     rv=[8.72, 8.71, 6.98, 16.16, 4.30, 15.76, 12.75, 6.79, 19.94],
#     σ_rv=[1.23, 1.24, 1.18, 1.08, 1.14, 1.05, 0.95, 1.15, 1.03],
#     inst_idx=fill(1, 9)
# ))
# rv_like_2 = OctofitterRadialVelocity.HARPS_RVBank_rvs("HD42581")
# rv_like_2.table.inst_idx .= 2

# rvsdf = vcat(
#     DataFrame(rv_like_1.table),
#     DataFrame(rv_like_2.table),
#     cols=:intersect
# )
# rv_like = StarAbsoluteRVLikelihood(rvsdf)
# Plots.plot(rv_like)
##
using CSV
rvsdf = CSV.read("Gl229A_RVs.csv", DataFrame, delim=' ', skipto=2, header=[
    # Epoch[JD] RV[m/s] RV_err[m/s] InstrumentID
    "epoch_jd",
    "rv",
    "σ_rv", 
    "inst_idx"
])
rvsdf.inst_idx .+= 1
rvsdf.epoch .= jd2mjd.(rvsdf.epoch_jd)
rv_like = StarAbsoluteRVLikelihood(rvsdf, instrument_names=[
    "HIRES_k",
    "HIRES_j",
    "HARPSpost",
    "HARPSpre",
    "UVES",
])

##
using StatsBase
function bin_rv(rv_like::MarginalizedStarAbsoluteRVLikelihood{T,jitter}, days) where{T,jitter}
    unique_epochs = unique(round.(rv_like.table.epoch/days).*days)
    bin_epoch = Float64[]
    bin_rv = Float64[]
    bin_σ_rv = Float64[]
    for epoch in unique_epochs
        mask = epoch .== round.(rv_like.table.epoch./days).*days
        if count(mask) > 0
            w = ProbabilityWeights(1 ./ rv_like.table.σ_rv[mask].^2)
            rv = mean(rv_like.table.rv[mask], w)
            # σ_rv = std(rv_like.table.σ_rv[mask], w)
            σ_rv = mean(rv_like.table.σ_rv[mask])/sqrt(count(mask)) # TODO: hack
            push!(bin_epoch, epoch)
            push!(bin_rv, rv)
            push!(bin_σ_rv, σ_rv)
        end
    end
    return MarginalizedStarAbsoluteRVLikelihood(
        Table(rv=bin_rv,epoch=bin_epoch,σ_rv=bin_σ_rv),
        instrument_name=rv_like.instrument_name,
        jitter=jitter
    )
end

rv_ll_1 = OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood(
    [r[1] for r in eachrow(rv_like.table[rv_like.table.inst_idx.==1,:])],
    instrument_name="HIRES_k",
    jitter=:j1
)
rv_ll_1 = bin_rv(rv_ll_1,10)
rv_ll_2 = OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood(
    [r[1] for r in eachrow(rv_like.table[rv_like.table.inst_idx.==2,:])],
    instrument_name="HIRES_j",
    jitter=:j2
)
rv_ll_2 = bin_rv(rv_ll_2,10)
rv_ll_3 = OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood(
    [r[1] for r in eachrow(rv_like.table[rv_like.table.inst_idx.==3,:])],
    instrument_name="HARPSpost",
    jitter=:j3
)
rv_ll_3 = bin_rv(rv_ll_3,10)
rv_ll_4 = OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood(
    [r[1] for r in eachrow(rv_like.table[rv_like.table.inst_idx.==4,:])],
    instrument_name="HARPSpre",
    jitter=:j4
)
rv_ll_4 = bin_rv(rv_ll_4,10)
rv_ll_5 = OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood(
    [r[1] for r in eachrow(rv_like.table[rv_like.table.inst_idx.==5,:])],
    instrument_name="UVES",
    jitter=:j5
)
rv_ll_5 = bin_rv(rv_ll_5,10)

@planet b Visual{KepOrbit} begin
    e ~ Uniform(0.0, 0.999)
    # a ~ LogUniform(10, 1000)
    a ~ LogUniform(1, 500)
    # mass ~ truncated(Normal(70, 5),lower=65,upper=80)
    # mass ~ LogUniform(1,500)
    mass = system.M_sec
    i ~ Sine()
    Ω ~ Uniform(0,2pi)
    ω ~ Uniform(0,2pi)
    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,b,50000.0) # epoch of bunch of astrom measurements 1996
end astrom_like_1

@system GL229A begin
    M_pri ~ LogUniform(0.05, 2)
    M_sec ~ LogUniform(1,500) # mjup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol
    plx ~ gaia_plx(; gaia_id)
    
  
    j1 ~ LogUniform(0.1, 100)
    j2 = system.j1
    j3 = system.j1
    j4 = system.j1
    j5 = system.j1

end rv_ll_1 rv_ll_2 rv_ll_3 rv_ll_4 rv_ll_5 b
model = Octofitter.LogDensityModel(GL229A; autodiff=:ForwardDiff, verbosity=4)

##


##
Octofitter.default_initializer!(model, verbosity=4)

##
chain,pt = octofit_pigeons(
    model,
    n_chains=16,
    n_rounds=8,
    explorer=SliceSampler(),
    n_chains_variational=16,
    variational = GaussianReference(first_tuning_round = 5),
    multithreaded=true,
    extended_traces = true
)

## Increment as needed
increment_n_rounds!(pt,1)
chain,pt = octofit_pigeons(pt)

## Generate corner plots of extended traces
serieses = []
for leg in (:fixed, :variational)
    # leg = :fixed # :fixed
    ii = filter(eachindex(pt.shared.tempering.indexer.i2t)) do i
        nt = pt.shared.tempering.indexer.i2t[i]
        nt.leg == leg
    end
    if leg ==:variational
        temp = pt.shared.tempering.variational_leg.schedule.grids
        temp = reverse(temp)
    elseif leg == :fixed
        temp = pt.shared.tempering.fixed_leg.schedule.grids
    else
        error()
    end
    ii = ii[sortperm(temp,rev=true)]
    temp = sort(temp,rev=true)
    chns_i = map(ii) do i
        chn = Octofitter.result2mcmcchain(  model.arr2nt.(model.invlink.(s[1:model.D] for s in get_sample(pt,i))));
        return i,chn
    end
    n_chns = length(chns_i)
    viz=(
        PairPlots.Scatter(markersize=2),
        # PairPlots.MarginStepHist()
        PairPlots.MarginDensity(
            linestyle= leg == :fixed ? :solid : :dot
        )
    )
    j = 0
    p = map(chns_i) do (i, chn)
        global j
        table_cols = Pair{Symbol,AbstractVector}[
            # :iter=>repeat(1:size(chn,1),outer=size(chn,3)),
        ]
        append!(table_cols, [sym => vec(chn[sym]) for sym in [
            :M_pri
            :M_sec
            :plx
            :j1
            :b_e
            :b_a
            :b_i
            :b_Ω
            :b_ω
            :b_θ
        ]])
        t = FlexTable(;table_cols...)
        # t = FlexTable(;(k=>v[1:10:end] for (k,v) in table_cols)...)
        ser = PairPlots.Series(
            t, label="chain-$i ($(temp[j+1]) )",
            color= Makie.cgrad(:turbo,rev=true)[j/n_chns],# [temp]
            bottomleft= leg == :fixed,
            topright= leg == :variational,
        ) => viz
        j+= 1
        return ser
    end
    append!(serieses, p)
end
figp = pairplot(
    serieses...,
    diagaxis=(
        spinewidth=3,
    )
)