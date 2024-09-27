using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
using DataFrames
using Dates
using OctofitterRadialVelocity

##

# Astrometry data
astrom_like = PlanetRelAstromLikelihood(Table(;
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

## Load RV Data
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

## Download HGCA data
gaia_id = 2940856402123426176
hgca_like = HGCALikelihood(; gaia_id)

##

mean_epoch = mean(astrom_like_1.table.epoch)
# @planet b AbsoluteVisual{KepOrbit} begin
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0, 0.999)
    a ~ LogUniform(1, 1000)
    mass ~ LogUniform(1,500)
    i ~ Sine()
    Ω ~ Uniform(0,2pi)
    ω ~ Uniform(0,2pi)
    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,b,$mean_epoch)
end astrom_like_1


@system GL229A_reparam begin
    M ~ truncated(Normal(0.579, 0.1), lower=0) # (Baines & Armstrong 2011).
    plx ~ gaia_plx(; gaia_id)
    
  
    jitter ~ Product([
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
    ])
    # rv0 ~ MvNormal(fill(1500, 5))
    rv0_net ~ Normal(0,1500)
    rv0_each ~ MvNormal(fill(10, 5))
    rv0 = system.rv0_each .+ system.rv0_net


    pmra ~ Normal(-136.42, 1000)
    pmdec ~ Normal(-715.02, 1000)
end rv_like hgca_like b
model = Octofitter.LogDensityModel(GL229A_reparam; autodiff=:ForwardDiff, verbosity=4)


## Initialize it near the mode.
# The posterior appears multi-modal at first, but this period range dominates.
model.starting_points= fill(model.link(
    [-0.4379843968818915, 5.156604091211061, 1.50584325204466, 2.0603093493895033, 1.2023862122831677, 0.959120843046445, 1.506903885206572, 2.381052701959341, -8.069099332206761, 2.294895076733609  …  1.1771693368026892, -145.47959330263868, -705.7990060941062, 1.7681218961823513, -3.406339274322457, -1.8064415810545709, -4.495184238907693, -0.21382190876833038, -0.20192496540804014, -0.13923518661202688
]), 100)
model.ℓπcallback(model.starting_points[1])



## Parallel Tempering
pt = pigeons(
    target=model,
    n_rounds=10,
    n_chains=24,
    explorer = #Compose(
        SliceSampler(p=10),
        #AutoMALA(),
    # ),
    record = [traces; round_trip; record_default();  record_online()],
    multithreaded=true,
    show_report=true,
)
chain_pt = Chains(pt.inputs.target, pt)

## Compare the samplers against eachother in one corner plot
# Wrapper of pairplot to make nice labels and select most the most interesting variables
octocorner(model, chain_pt, small=true)
# Note: pass small=true to just show the most important variables

##
using InferenceReport
InferenceReport.report(pt, writer=InferenceReport.Documenter.LaTeX())