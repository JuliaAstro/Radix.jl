# The synthetic spectrum of the logo: the emission lines of the Fe K complex of two slabs of gas that Radix computes (a few minutes to half an hour with
# 8 threads), with noise.
#
#     julia -t 8 --project=. docs/logo/spectrum.jl atdb.fits coheat.dat
#
# writes `fek.csv`: the energy (eV), the model and one observation of it (the counts in each bin, with Poisson noise). The gas is two slabs of 10²² cm⁻² and a
# density of 10¹⁰ cm⁻³ that a power law E⁻⁰·⁹ illuminates, at log ξ = 1 (the fluorescence of the cold iron: Kα at 6.40 keV and Kβ at 7.06 keV) and at
# log ξ = 3 (the recombination lines of Fe XXV near 6.7 keV and Fe XXVI at 6.97 keV), seen without the direct continuum, as an obscured nucleus is, over a
# weak scattered continuum.
using Radix, Random

const peak_counts = 400.0          # the counts in the bin of the strongest line (Kα)
const continuum_fraction = 0.03    # the scattered continuum, in units of that peak
const hot_weight = 0.6             # the emission of the log ξ = 3 slab, relative to the log ξ = 1 slab
const window = (5605.0, 7205.0)    # eV: centred on the Fe Kα line at 6.405 keV

# the emission (the lines and the continuum of the zones) that leaves the slab, in units of its incident continuum
emission(spectrum, bins) = (spectrum.diffuse_outward[bins] .+ spectrum.lines_outward[bins]) ./ spectrum.incident[bins]

# the Poisson deviate: the product of uniform deviates for a mean below 30 and the Gaussian limit above
function poisson(rng, λ)
    λ < 30 || return max(0, round(Int, λ + sqrt(λ)*randn(rng)))
    limit, k, p = exp(-λ), 0, 1.0
    while true
        p *= rand(rng)
        p <= limit && return k
        k += 1
    end
end

# the energy, the model and the observation (counts) in the window from the `spectrum` of `tables` of the two slabs
function observation(cold, hot; seed=6400)
    bins = findall(e -> window[1] <= e <= window[2], cold.energy)
    energy = cold.energy[bins]
    lines = emission(cold, bins) .+ hot_weight .* emission(hot, bins)
    model = peak_counts .* (lines ./ maximum(lines) .+ continuum_fraction)
    rng = Xoshiro(seed)
    energy, model, [poisson(rng, m) for m in model]
end

function write_csv(path, (energy, model, observed))
    open(path, "w") do io
        println(io, "energy,model,observed")
        for i in eachindex(energy)
            println(io, round(energy[i]; digits=1), ",", round(model[i]; digits=2), ",", observed[i])
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    atdb, coheat = ARGS
    gas = Model(atdb, coheat)
    spectra = map(((ξ, T),) -> tables(gas, slab(gas; density=1e10, column=1e22, logξ=ξ, luminosity=1e38, α=-0.9, equilibrium=true, T, emult=1.0, steps=3,
        vturb=100.0, ncn2=9999)).spectrum, ((1.0, 100.0), (3.0, 1000.0)))
    write_csv(joinpath(@__DIR__, "fek.csv"), observation(spectra...))
end
