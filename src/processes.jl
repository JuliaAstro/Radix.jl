# The processes of the gas that do not depend on the levels of the ions: Compton scattering, free-free absorption, bremsstrahlung and
# Thomson scattering. Each is a type, with its coefficients as functions named by their standard symbols, α for the absorption and σ for the
# scattering coefficient (j for the emission), of the photon energies and the temperature; the density of what absorbs or scatters, the opacity
# and the heating and cooling (XSTAR's freef, bremem, comp2, cmpfnc and the first term of the opacity of calc_emis_all) derive from them.

const compton_Eph_min = 1e-4             # cmpfnc: below this energy (in electron rest masses) its limit 4 Te - Eph is used
const compton_Te_floor = 1e-10          # comp2: added to kT (eV) in the ratio of the rest energy to it
const ion_density_factor = Float64(1.4f0)            # both: the density of the ions, in electron densities (n_i Z² = 1.4 nₑ)

"""
    AbstractContinuum

A process of the gas, with the coefficients (functions of the photon energies `E` (eV) or of a `Radiation`, and the temperature `T` in 10⁴ K)

    α(process, E, T)       # absorption
    σ(process, E, T)       # scattering
    j(process, E, T)       # emission

per unit of the density `density(process, nₑ)` of what absorbs, emits or scatters; the opacity (cm⁻¹) is `opacity(process, E, T, nₑ)`, the
density times `α + σ`. A process is an `AbstractAbsorption` (free-free), an `AbstractScattering` (Thomson and Compton) or emits (bremsstrahlung); the scatterers are the part of the opacity that is in the continuum opacity `opakcont` as well as the total.
The heating and cooling (erg cm⁻³ s⁻¹) in the radiation field at the electron density `nₑ` derive from the coefficients:
`heating(process, radiation, T, nₑ)`, `cooling(...)` and `heating_cooling(...)` for both. What a process does not have is 0.

- `FreeFree()`: absorption `α` per electron and per ion (cm⁵; density `nₑ nᵢ`, `nᵢ = 1.4 nₑ`): heats the gas and has an opacity (`freef`);
- `Bremsstrahlung()`: the emission `j` of the free-free process, per electron and per ion: cools the gas (`bremem`, integrated by `heatf`);
- `Thomson(cfrac=0.0)`: scattering `σ = σ_T (1 - cfrac)` (cm², per electron) for the covering fraction `cfrac`: an opacity only;
- `Compton(Te, Eph, Es)`: Compton scattering with the energy it exchanges between the radiation and the electrons, from the table of XSTAR's
  `coheat.dat` (`load(Compton, path)`): heating and cooling (`comp2`, `cmpfnc`).

`standard_processes(compton; cfrac)` makes the four of a run of XSTAR.
"""
abstract type AbstractContinuum end
abstract type AbstractAbsorption <: AbstractContinuum end
abstract type AbstractScattering <: AbstractContinuum end

struct FreeFree <: AbstractAbsorption end
struct Bremsstrahlung <: AbstractContinuum end
struct Thomson <: AbstractScattering
    cfrac::Float64
end
Thomson() = Thomson(0.0)

"""
    Compton(Te, Eph, Es)

Compton scattering with the table of the Compton heating (XSTAR's `coheat.dat`, read by `load(Compton, path_or_io)`): `Es[i, j]` at the electron
temperature `Te[i]` (in electron rest masses) and the photon energy `Eph[j]` (in electron rest masses).

- `σ(compton, Eph::Real, Te::Real)` is XSTAR's `cmpfnc`: the Compton function at the photon energy `Eph` and the electron temperature `Te` (both in
  electron rest masses), interpolated linearly in the table (its limit `4 Te - Eph` below an energy of `compton_Eph_min`);
- `σ(compton, E, T)` is the scattering cross section, `σ_T`, and `heating(compton, E, T)` and `cooling(compton, E, T)` those of the energy exchanged
  per electron, `σ_T E/mc²` and `σ_T (cmpfnc + E/mc²)`, whose integrals over the radiation are the Compton heating and cooling;
- `integral(compton, radiation, T)` is XSTAR's `comp2`: `(cmp1, cmp2)`, the heating being `cmp1 nₑ` and the cooling `kT cmp2 nₑ` (eV, times erg per eV).

The scattering opacity is that of `Thomson`: `density(compton, nₑ)` is 0.
"""
struct Compton <: AbstractScattering
    Te::Vector{Float64}             # electron temperatures (electron rest masses)
    Eph::Vector{Float64}            # photon energies (electron rest masses)
    Es::Matrix{Float64}             # the tabulated energy exchanged in scattering (the Compton function), Es[i, j] at Te[i] and Eph[j]
end

function load(::Type{Compton}, io::IO)
    entries = [split(line) for line in eachline(io) if !isempty(strip(line))]
    n = maximum(parse(Int, row[1]) for row in entries)
    m = maximum(parse(Int, row[2]) for row in entries)
    compton = Compton(zeros(n), zeros(m), zeros(n, m))
    for row in entries
        i, j = parse(Int, row[1]), parse(Int, row[2])
        compton.Te[i], compton.Eph[j], compton.Es[i, j] = parse(Float64, row[3]), parse(Float64, row[5]), parse(Float64, row[6])
    end
    compton
end
load(::Type{Compton}, path::AbstractString) = open(io -> load(Compton, io), path)

"""
    standard_processes(compton; cfrac=0.0)

The processes of XSTAR: `Compton` scattering (`compton`), `FreeFree()` absorption, `Bremsstrahlung()` and `Thomson(cfrac)` scattering.
"""
standard_processes(compton::Compton; cfrac=0.0) = (compton, FreeFree(), Bremsstrahlung(), Thomson(cfrac))

# the coefficients: zero by default, and of a Radiation on its grid
α(::AbstractContinuum, E::AbstractVector, T) = zeros(length(E))
σ(::AbstractContinuum, E::AbstractVector, T) = zeros(length(E))
j(::AbstractContinuum, E::AbstractVector, T) = zeros(length(E))
for coefficient in (:α, :σ, :j)
    @eval $coefficient(process::AbstractContinuum, rad::Radiation, T) = $coefficient(process, rad.E, T)
end
density(::AbstractContinuum, nₑ) = 0.0
opacity(process::AbstractContinuum, E, T, nₑ) = density(process, nₑ) .* (α(process, E, T) .+ σ(process, E, T))

heating(::AbstractContinuum, ::Radiation, T, nₑ) = 0.0
cooling(::AbstractContinuum, ::Radiation, T, nₑ) = 0.0
heating_cooling(process::AbstractContinuum, rad::Radiation, T, nₑ) = (heating(process, rad, T, nₑ), cooling(process, rad, T, nₑ))

# the energy that a band of the spectrum gives to what has the coefficient: ∫ F coefficient dE (erg cm⁻³ s⁻¹; the energies are in eV)
function absorbed(F, coefficient, E)
    total = 0.0
    for k in 2:length(E)
        total += (F[k]*coefficient[k] + F[k - 1]*coefficient[k - 1])*(E[k] - E[k - 1])/2
    end
    total*constants().ergsev
end

# free-free: the absorption coefficient (1 - e^{-E/kT} is the stimulated emission) and the heating of the radiation it absorbs
function α(::FreeFree, E::AbstractVector, T)
    ekt = T*constants().kT_eV
    [constants().ff_absorption_coeff/sqrt(T)/e^3*(1 - exp(-e/ekt)) for e in E]
end
density(::FreeFree, nₑ) = nₑ*ion_density_factor*nₑ
heating(process::FreeFree, rad::Radiation, T, nₑ) = absorbed(rad.F, opacity(process, rad.E, T, nₑ), rad.E)

# bremsstrahlung: the emission coefficient per electron and ion, and the cooling, its integral over the energy grid
function j(::Bremsstrahlung, E::AbstractVector, T)
    ekt = T*constants().kT_eV
    [constants().ff_emission_coeff*exp(-e/ekt)/sqrt(T) for e in E]
end
density(::Bremsstrahlung, nₑ) = nₑ*ion_density_factor*nₑ
cooling(process::Bremsstrahlung, rad::Radiation, T, nₑ) =
    absorbed(ones(length(rad.E)), density(process, nₑ) .* j(process, rad.E, T), rad.E)

# Thomson scattering: the same coefficient at all energies, per electron (it is in the continuum opacity as well as the total)
σ(process::Thomson, E::AbstractVector, T) = fill(constants().sigma_thomson*max(0.0, 1 - process.cfrac), length(E))
density(::Thomson, nₑ) = nₑ

# hunt3: the index of the table value just below x, between 1 and the length of the table
hunt_index(xx, x) = clamp(searchsortedlast(xx, x), 1, length(xx))

# Compton scattering: cmpfnc, the coefficients of the energy exchanged, and the heating and cooling that derive from them
function σ(compton::Compton, Eph::Real, Te::Real)
    Eph > compton_Eph_min || return 4Te - Eph
    iE = clamp(hunt_index(compton.Eph, Eph), 2, length(compton.Eph))
    iT = clamp(hunt_index(compton.Te, Te), 2, length(compton.Te))
    Es = compton.Es
    dEs_dTe = (Es[iT, iE] - Es[iT - 1, iE] + Es[iT, iE - 1] - Es[iT - 1, iE - 1])/(2*(compton.Te[iT] - compton.Te[iT - 1]))
    dEs_dEph = (Es[iT, iE] - Es[iT, iE - 1] + Es[iT - 1, iE] - Es[iT - 1, iE - 1])/(2*(compton.Eph[iE] - compton.Eph[iE - 1]))
    dEs_dTe*(Te - compton.Te[iT - 1]) + dEs_dEph*(Eph - compton.Eph[iE - 1]) + Es[iT - 1, iE - 1]
end
σ(::Compton, E::AbstractVector, T) = fill(constants().sigma_thomson, length(E))

rest_masses(E) = E ./ constants().electron_rest_eV
heating(::Compton, E::AbstractVector, T) = constants().sigma_thomson .* rest_masses(E)
function cooling(compton::Compton, E::AbstractVector, T)
    K = constants()
    Te = (T*K.kT_eV + compton_Te_floor)/K.electron_rest_eV       # the temperature in electron rest masses
    Eph = rest_masses(E)
    K.sigma_thomson .* (σ.(Ref(compton), Eph, Te) .+ Eph)
end

function integral(compton::Compton, rad::Radiation, T)
    ekt = T*constants().kT_eV
    (absorbed(rad.F, heating(compton, rad.E, T), rad.E)/constants().ergsev,
        absorbed(rad.F, cooling(compton, rad.E, T), rad.E)/constants().ergsev/ekt)
end

function heating_cooling(compton::Compton, rad::Radiation, T, nₑ)
    cmp1, cmp2 = integral(compton, rad, T)
    (cmp1*nₑ*constants().ergsev, T*constants().kT_eV*cmp2*nₑ*constants().ergsev)
end
heating(compton::Compton, rad::Radiation, T, nₑ) = heating_cooling(compton, rad, T, nₑ)[1]
cooling(compton::Compton, rad::Radiation, T, nₑ) = heating_cooling(compton, rad, T, nₑ)[2]
