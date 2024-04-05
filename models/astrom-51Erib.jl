using Octofitter, Distributions
using CairoMakie
using PairPlots
using Pigeons
using DataFrames
##

# Note: will need to download a datadep. Can type `y` or set an ENV variable.
planet_name = "51erib"
gaia_id = 2431157720981843200
star_mass_unc = 1.550, 0.005
# To find the above identifier, search SIMBAD: http://simbad.cds.unistra.fr/simbad/
# with the name of the system (without the b at the end) and copy the 
# ID number next to 'GAIA DR2'
astrom_like = PlanetRelAstromLikelihood(DataFrame([
57009.1  69.3356  -448.917  1.88446  1.81828  -0.00155179
57052.1  78.3783  -444.96   2.0788   2.05287  -0.00117613
57053.1  77.8303  -450.121  2.57     2.38622  -0.0209427
57054.3  76.9638  -455.037  23.8713  24.2003  0.00286876
57266.4  100.052  -443.966  2.22649  2.07075  -0.0110573
57290.0  100.8    -442.0    2.9      3.6      0.0
57291.0  108.7    -440.7    8.8      13.7     0.0
57332.2  108.641  -439.656  5.39868  4.5087   -0.0172164
57374.2  112.918  -441.705  6.18322  4.55306  -0.0635376
57376.2  112.464  -440.892  3.02297  3.39833  0.014624
57403.0  114.3    -442.2    4.5      5.2      0.0
57415.0  110.406  -440.845  6.03121  4.07169  -0.0682188
57649.4  142.053  -432.057  2.01957  2.05944  0.00370176
57652.4  141.521  -428.673  2.66361  2.44022  -0.00564645
57734.0  152.9    -427.1    3.4      4.6      0.0
57739.1  153.258  -422.449  2.15789  2.1202   -0.00251691
58024.0  185.0    -409.2    2.0      2.1      0.0
58068.3  187.509  -406.365  3.01691  3.05686  -0.00638734
58442.2  219.468  -374.674  2.00066  1.76217  0.0759851
58379.0  220.2    -384.8    2.8      3.1      0.0
], [
   :epoch,    :ra, :dec, :σ_ra, :σ_dec, :cor
]))
gaia_id = 3205095125321700480
hgcal   = HGCALikelihood(;gaia_id)

# To inspect data: astrom_like.table

# You can select different subsets of the data here (replace 1:end):
# astrom_like = PlanetRelAstromLikelihood(astrom_like_0.table[1:end,1,1])

# "gold standard" results
# chain_orbitize = Octofitter.loadhdf5(planet_name)

# Set angle reference epoch to average date of observations
const t_ref_epoch = mean(astrom_like.table.epoch)



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
    tp = θ_at_epoch_to_tperi(system,b,t_ref_epoch) 
end astrom_like
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
    tp = θ_at_epoch_to_tperi(system,b,t_ref_epoch) 
end astrom_like
@system system_param2 begin
    plx ~ gaia_plx(;gaia_id)
    M   ~ truncated(Normal(star_mass_unc...), lower=0)
end b
model_param2 = Octofitter.LogDensityModel(system_param2)

## NUTS with Pathfinder intialization
chain_nuts = octofit(model_param1,adaptation=4096, iterations=4096)

## Parallel Tempering
pt = pigeons(
    target=model_param2,
    n_rounds=10,
    n_chains=24,
    explorer = Compose(
        SliceSampler(p=10),
        AutoMALA(),
    ),
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
    # chain_orbitize, 
    chain_nuts,
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
InferenceReport.report(pt, writer=InferenceReport.Documenter.LaTeX())