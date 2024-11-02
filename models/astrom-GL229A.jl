using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
using DataFrames
using Dates


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


astrom_like_2 = PlanetRelAstromLikelihood(Table(;
    epoch=mjd.(["2024-10-24"]),
    ra=[-302.41345030000105],
    dec=[-4294.98910568],
    σ_ra=[30],
    σ_dec=[30],
))
##
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0.7, 0.999)
    a ~ LogUniform(1, 500)
    mass = system.M_sec
    i ~ Sine()
    ω_p_Ω ~ UniformCircular()
    ω_m_Ω ~ UniformCircular()
    Ω = (b.ω_p_Ω - b.ω_m_Ω)/2
    ω = (b.ω_p_Ω + b.ω_m_Ω)/2
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50000.0) # epoch of bunch of astrom measurements 1996
end astrom_like_1 astrom_like_2


@system GL229A begin
    M_pri ~ LogUniform(0.05, 2)
    M_sec ~ LogUniform(1,500) # mjup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol
    # M ~ truncated(Normal(0.579, 0.1), lower=0.4, upper=0.7) # (Baines & Armstrong 2011).
    # M ~ Uniform(0.2, 2) # (Baines & Armstrong 2011).
    plx ~ truncated(Normal(173.57398986816406,0.01704956591129303), lower=0.0)
    
end b
model2 = Octofitter.LogDensityModel(GL229A; autodiff=:ForwardDiff, verbosity=4)

##
using Random
Random.seed!(0)
Octofitter.default_initializer!(model2,nruns=1,verbosity=4) # nruns must = 1 or becomes non-deterministic
# Must use Random.seed! and not pass an rng; passing rng to AdvancedHMC is broken.
chain_hmc = octofit(model2,verbosity=4,adaptation=1000,iterations=4100)
##
chain_pt,pt = octofit_pigeons(
    model2,
    # 12 rounds is about the same wall-clock time as the above HMC
    n_rounds=8,
    n_chains=1,
    n_chains_variational=0,
    variational=nothing,#GaussianReference(),
    # explorer=Compose(SliceSampler(),AutoMALA()),
    explorer=AutoMALA(step_size=0.0004),
)
##
octocorner(
    model,
    chain_hmc,
    chain_pt;
    includecols=(:iter,:logpost),
    viz=(PairPlots.Scatter(markersize=2),)
)