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
mean_epoch = mean(astrom_like.table.epoch)
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0, 0.999)
    a ~ LogUniform(1, 1000)
    mass ~ LogUniform(1,500)
    i ~ Sine()
    # Ω ~ UniformCircular()
    # ω ~ UniformCircular()
    # θ ~ UniformCircular()
    Ω ~ Uniform(0, 2pi)
    ω ~ Uniform(0, 2pi)
    θ ~ Uniform(0, 2pi)
    tp = θ_at_epoch_to_tperi(system,b,$mean_epoch)
end astrom_like

@system GL229A begin
    M ~ truncated(Normal(0.579, 0.1), lower=0) # (Baines & Armstrong 2011).
    plx ~ gaia_plx(; gaia_id)
    jitter ~ Product([
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
        truncated(Normal(0, 10), lower=0),
    ])
    rv0 ~ MvNormal(fill(1500, 5))
    pmra ~ Normal(-136.42, 1000)
    pmdec ~ Normal(-715.02, 1000)
end rv_like hgca_like b

model = Octofitter.LogDensityModel(GL229A; autodiff=:ForwardDiff, verbosity=4)

## Initialize it near the mode.
# The posterior appears multi-modal at first, but this period range dominates.
model.starting_points= fill(model.link([
    0.6666921097218382
    173.5954322176413
      5.460089174960165
      2.9396971616582293
      3.053929939180548
      1.8958506333312857
      4.802605001566408
    189.09099932852274
    197.32536572923578
    195.32333134355542
    198.16288501676306
    199.29659687079726
   -145.1899530896007
   -706.0811075841312
      0.7371378360662447
     40.99933962474703
     71.52220819748763
      0.8210639011591009
    #  -0.9099306940000274
    #   0.1891153547867482
    #  -0.6992656365254634
    #   0.699666397551208
    #  -0.9762798298963133
    #   0.21651257177731173
    rem2pi(atan( 0.1891153547867482, -0.909930694000027),RoundDown)
    rem2pi(atan(0.699666397551208, -0.6992656365254634),RoundDown)
    rem2pi(atan(0.21651257177731173, -0.9762798298963133),RoundDown)
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