using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
using DataFrames
using Dates
using OctofitterRadialVelocity


# Old measurements from Brandt 2020 
# The Astronomical Journal, 160:196 (15pp), 2020 October
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
# scatter(astrom_like_1.table.sep,astrom_like_1.table.pa)


## Load RVs -- split into 5 instruments with indepdendent zero points.
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

## Bin data to speed up testing the model
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


##
Octofitter.default_initializer!(model,verbosity=4)

## Sample
chain_automala,pt_automala = octofit_pigeons(
    model,
    n_chains=32,
    n_rounds=8,
    explorer=SliceSampler(),
    n_chains_variational=0,
    # variational = GaussianReference(first_tuning_round = 5),
    variational=nothing,
    multithreaded=true,
    extended_traces = true
)


## Run as needed to add rounds
increment_n_rounds!(pt_automala,1)
chain_automala,pt_automala = octofit_pigeons(pt_automala)


## Corner plot of all extended traces
chns_i = map(1:32) do i
    chn = Octofitter.result2mcmcchain(  model.arr2nt.(model.invlink.(s[1:model.D] for s in get_sample(pt_automala,i))));
    return i,chn
end
# temp = pt_automala.shared.tempering.variational_leg.schedule.grids
# temp = pt_automala.shared.tempering.fixed_leg.schedule.grids
temp = pt_automala.shared.tempering.schedule.grids
ii = sortperm(temp)
chns_i = chns_i[ii]
n_chns = length(chns_i)
viz=(
    PairPlots.Scatter(markersize=2),
    # PairPlots.MarginStepHist()
    PairPlots.MarginDensity()
)
j = 0
p = map(chns_i) do (i, chn)
    global j
    table_cols = Pair{Symbol,AbstractVector}[
        :iter=>repeat(1:size(chn,1),outer=size(chn,3)),
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
    ser = PairPlots.Series(
        t, label="chain-$i ($(temp[j+1])",
        color= Makie.cgrad(:turbo,rev=true)[j/n_chns]# [temp]
    ) => viz
    j+= 1
    return ser
end
figp = pairplot(p...,)
