# The spectrum that illuminates the gas: XSTAR's power-law source (ispec4, normalized by ispecgg) and the distance from the source of
# a given ionization parameter.

const luminosity_unit = 1e38            # the spectra of XSTAR are in units of 10³⁸ erg s⁻¹
const spectrum_floor = Float64(1f-24)             # ispec4: the value of the spectrum below its cutoff energy
const spectrum_cutoff = Float64(0.01f0)           # ispec4: ecut (eV): the spectrum is the floor below it
const luminosity_E_min = 13.6                     # ispec4: the energy range of the luminosity, 1 to 1000 Ry (eV), by bins
const luminosity_E_max = 1.36e4
const luminosity_E_min_single = Float64(13.6f0)   # ispecgg compares the energies with 13.6 as a single-precision literal

"""
    power_law(E, α, luminosity)

The spectrum `E^α` on the energy grid `E` (eV) with the luminosity `luminosity` (erg s⁻¹) between 13.6 eV and 13.6 keV, in erg s⁻¹ erg⁻¹ (XSTAR's `ispec4` followed by
`ispecgg`, the `pow` spectrum with `trad = α`): `E^α` above 0.01 eV and 10⁻²⁴ below, normalized to the luminosity in the bins
`nbin(13.6 eV)` to `nbin(13.6 keV)` and then again in the bins of the energies in the range.
"""
function power_law(E::AbstractVector, α, luminosity)
    K = constants()
    n = length(E)
    rad = Radiation(E, zeros(n))
    nb1, nb2 = nbin(rad, luminosity_E_min), nbin(rad, luminosity_E_max)
    shape = [e > spectrum_cutoff ? e^α : spectrum_floor for e in E]
    trapezoid(f, i) = (f[i] + f[i - 1])*(E[i] - E[i - 1])/2
    xlum = luminosity/luminosity_unit
    spectrum = shape .* (xlum/sum(trapezoid(shape, i) for i in nb1:nb2)/K.ergsev)
    total = sum(trapezoid(spectrum, i) for i in 2:n if luminosity_E_min_single <= E[i] <= luminosity_E_max; init=0.0)
    spectrum .*= xlum/total/K.ergsev
    spectrum .* luminosity_unit
end

"""
    source_distance(luminosity, ξ, n)

The distance (cm) from a source of the luminosity (erg s⁻¹) at which the ionization parameter `ξ = L/(n r²)` is `ξ` for the density `n` (cm⁻³).
"""
source_distance(luminosity, ξ, n) = sqrt(luminosity/(ξ*n))
