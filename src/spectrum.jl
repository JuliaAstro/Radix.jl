# The spectrum that illuminates the gas: XSTAR's power-law source (ispec4, normalized by ispecgg) and the distance from the source of
# a given ionization parameter.

const luminosity_unit = 1e38            # the spectra of XSTAR are in units of 10³⁸ erg s⁻¹
const spectrum_floor = Float64(1f-24)             # ispec4: the value of the spectrum below its cutoff energy
const spectrum_cutoff = Float64(0.01f0)           # ispec4: ecut (eV): the spectrum is the floor below it
const luminosity_E_min = 13.6                     # ispec4: the energy range of the luminosity, 1 to 1000 Ry (eV), by bins
const luminosity_E_max = 1.36e4
const luminosity_E_min_single = Float64(13.6f0)   # ispecgg compares the energies with 13.6 as a single-precision literal

# the luminosity (erg s⁻¹) in the bins of a spectrum between 13.6 eV and 13.6 keV, `ispec4`'s sum by the indices of the bins `nbin(13.6 eV)` to `nbin(13.6 keV)` or `ispec`'s by the energies
trapezoid(f, E, i) = (f[i] + f[i - 1])*(E[i] - E[i - 1])/2
band_by_index(f, E) = (rad = Radiation(E, zeros(length(E))); sum(trapezoid(f, E, i) for i in nbin(rad, luminosity_E_min):nbin(rad, luminosity_E_max)))
band_by_energy(f, E) = sum(trapezoid(f, E, i) for i in 2:length(E) if luminosity_E_min <= E[i] <= luminosity_E_max; init=0.0)

# `ispecgg`, applied to every spectrum: the luminosity again, in the bins whose energies are in the range (with 13.6 as a single-precision literal)
function normalize_luminosity(E, spectrum, luminosity)
    total = sum(trapezoid(spectrum, E, i) for i in 2:length(E) if luminosity_E_min_single <= E[i] <= luminosity_E_max; init=0.0)
    spectrum .* (luminosity/luminosity_unit/total/constants().ergsev*luminosity_unit)
end

"""
    power_law(E, α, luminosity)

The spectrum `E^α` on the energy grid `E` (eV) with the luminosity `luminosity` (erg s⁻¹) between 13.6 eV and 13.6 keV, in erg s⁻¹ erg⁻¹ (XSTAR's `ispec4` followed by
`ispecgg`, the `pow` spectrum with `trad = α`): `E^α` above 0.01 eV and 10⁻²⁴ below, normalized to the luminosity in the bins
`nbin(13.6 eV)` to `nbin(13.6 keV)` and then again in the bins of the energies in the range.
"""
function power_law(E::AbstractVector, α, luminosity)
    shape = [e > spectrum_cutoff ? e^α : spectrum_floor for e in E]
    spectrum = shape .* (luminosity/luminosity_unit/band_by_index(shape, E)/constants().ergsev*luminosity_unit)
    normalize_luminosity(E, spectrum, luminosity)
end

const blackbody_kT_inverse = Float64(1.16f-3)      # starf: 1/kT in eV⁻¹ is this over `trad`, the temperature in 10⁷ K
const blackbody_tiny = Float64(1f-37)
const blackbody_x_small = Float64(1f-3)            # below this x = E/kT the Planck function is its Rayleigh-Jeans limit
const blackbody_x_large = Float64(150f0)           # above it the spectrum is a constant over E (and exactly at it, zero)
const blackbody_unit = Float64(3.1415f22)
const blackbody_large_unit = Float64(1f37)

"""
    blackbody(E, T, luminosity)

The spectrum of a black body of the temperature `T` (10⁷ K, XSTAR's `trad` of `spectrum=bbody`) with the luminosity `luminosity` (erg s⁻¹) between 13.6 eV and 13.6 keV (`starf` and `ispecgg`):
`E³/(e^x - 1)` with `x = E/kT` and `1/kT = 1.16×10⁻³/T` eV⁻¹, its limit `E³/x` below `x = 10⁻³`, and the 'tail' `1/E` above `x = 150` (for the luminosity 1: `starf` is called with it and `ispecgg` sets the luminosity). The bins at exactly `x = 150` are 0.
"""
function blackbody(E::AbstractVector, T, luminosity)
    K = constants()
    xkt = blackbody_kT_inverse/(blackbody_tiny + T)
    xlum = 1.0                                         # (`xstar` calls `starf` with a luminosity of 1 and `ispecgg` sets the luminosity)
    shape = map(E) do e
        x = e*xkt
        value = x < blackbody_x_small ? e^3/x : x < blackbody_x_large ? e^3/(expo(x) - 1) : x > blackbody_x_large ? xlum/e/K.ergsev/blackbody_large_unit : 0.0
        blackbody_unit*value
    end
    band = sum(trapezoid(shape, E, i) for i in 2:length(E) if luminosity_E_min_single <= E[i] <= luminosity_E_max; init=0.0)
    normalize_luminosity(E, shape .* (xlum/band/K.ergsev*luminosity_unit), luminosity)
end

"""
    thermal_bremsstrahlung(E, T, luminosity)

The spectrum `e^{-E/kT}` of thermal bremsstrahlung (`spectrum=brems`, `ispec`) with `kT = 861.707 T` eV, `T` being the temperature in 10⁷ K, and the luminosity `luminosity` (erg s⁻¹) between 13.6 eV and 13.6 keV.
"""
function thermal_bremsstrahlung(E::AbstractVector, T, luminosity)
    K = constants()
    ekt = 1000*K.kT_eV*T
    shape = [expo(-e/ekt) for e in E]
    spectrum = shape .* (luminosity/luminosity_unit/band_by_energy(shape, E)/K.ergsev*luminosity_unit)
    normalize_luminosity(E, spectrum, luminosity)
end

const table_floor = 1e-49

"""
    tabulated(E, table_E, table_F, luminosity; units=:energy)

The spectrum of a table (`spectrum=file`, `ispecg`) on the grid `E` (eV): the table `table_F` at the energies `table_E` (ascending, in the units `units`: `:energy` (erg s⁻¹ erg⁻¹), `:photons` (multiplied by the energy, as `rread1` does) or `:log10`)
is interpolated linearly in the logarithms of both, 0 outside it, and normalized to the luminosity (erg s⁻¹) between 13.6 eV and 13.6 keV and then by `ispecgg`.
"""
function tabulated(E::AbstractVector, table_E::AbstractVector, table_F::AbstractVector, luminosity; units=:energy)
    F = units === :photons ? table_F .* table_E : units === :log10 ? 10.0 .^ table_F : float.(table_F)
    n = length(table_E)
    low, high = minmax(table_E[1], table_E[end])
    shape = map(E) do x
        (low <= x <= high) || return 0.0
        j = clamp(searchsortedlast(table_E, x), 1, n - 1)
        zr1, zr2 = log10(max(F[j + 1], table_floor)), log10(max(F[j], table_floor))
        ep1, ep2 = log10(max(table_E[j + 1], table_floor)), log10(max(table_E[j], table_floor))
        alx = min(max(log10(x), ep2), ep1)
        10.0^((zr1 - zr2)*(alx - ep2)/(ep1 - ep2 + table_floor) + zr2)
    end
    band = sum(trapezoid(shape, E, i) for i in 2:length(E) if luminosity_E_min <= E[i] <= luminosity_E_max; init=0.0)*constants().ergsev
    normalize_luminosity(E, shape .* (luminosity/luminosity_unit/band*luminosity_unit), luminosity)
end

"""
    source_distance(luminosity, ξ, n)

The distance (cm) from a source of the luminosity (erg s⁻¹) at which the ionization parameter `ξ = L/(n r²)` is `ξ` for the density `n` (cm⁻³).
"""
source_distance(luminosity, ξ, n) = sqrt(luminosity/(ξ*n))
