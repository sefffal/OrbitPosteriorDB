using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
##

# Note: will need to download a datadep. Can type `y` or set an ENV variable.
planet_name = "hd984b"
gaia_id = 2431157720981843200
star_mass_unc = 1.16, 0.10
# To find the above identifier, search SIMBAD: http://simbad.cds.unistra.fr/simbad/
# with the name of the system (without the b at the end) and copy the 
# ID number next to 'GAIA DR2'

astrom_like_seppa, astrom_like_radec = Octofitter.Whereistheplanet_astrom(planet_name)

# To inspect data: astrom_like.table

# You can select different subsets of the data here (replace 1:end):
# astrom_like = PlanetRelAstromLikelihood(astrom_like_0.table[1:end,1,1])

# "gold standard" results
chain_orbitize = Octofitter.loadhdf5(planet_name)

# Set angle reference epoch to average date of observations
const t_ref = mean(astrom_like_radec.table.epoch)



# I have phrased the model two ways: with "UniformCircular" variables, 
# a shortcut for two standard normal distributions x & y, then calculate
# atan(y,x); and a straightforward way with Uniform priors. The former strategy
# is necessary for NUTS since it allows the sampler to deal with modes that 
# wrap past zero.

##
@planet b Visual{KepOrbit} begin
    a ~ LogUniform(1, 500)
    e ~ Uniform(0.0, 0.99)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,t_ref) 
end astrom_like_radec astrom_like_seppa
@system system_param1 begin
    plx ~ gaia_plx(;gaia_id)
    M   ~ truncated(Normal(star_mass_unc...), lower=0)
end b
model_param1 = Octofitter.LogDensityModel(system_param1)


##
@planet b Visual{KepOrbit} begin
    a ~ LogUniform(1, 500)
    e ~ Uniform(0.0, 0.99)
    i ~ Sine()
    ω ~ Uniform(0,2pi)
    Ω ~ Uniform(0,2pi)
    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,b,t_ref) 
end astrom_like_radec astrom_like_seppa
@system system_param2 begin
    plx ~ gaia_plx(;gaia_id)
    M   ~ truncated(Normal(star_mass_unc...), lower=0)
end b
model_param2 = Octofitter.LogDensityModel(system_param2)

## NUTS with Pathfinder intialization
# chain_nuts = octofit(model_param1,adaptation=4096, iterations=4096)

## Parallel Tempering
pt = pigeons(
    target=model_param2,
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
octocorner(
    model_param2,
    # Compare three different samplers:
    chain_orbitize, 
    # chain_nuts,
    chain_pt,

    viz=(
        PairPlots.Scatter(markersize=2),
        PairPlots.MarginStepHist(linewidth=2),
    ),
    axis=(;
        b_a = (;
            scale=Makie.pseudolog10,
            ticks=2 .^ (0:10)
        )
    ),
    small=true
)
# Note: pass small=true to just show the most important variables
# Note: `tp` might not agree from the "orbitize" chain, but that's not important.


##
using InferenceReport
InferenceReport.report(pt)#, writer=InferenceReport.Documenter.LaTeX())