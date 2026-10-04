# Heating and cooling of a gas: the energy that the level populations exchange with the radiation and the electrons, and the
# Compton, free-free and bremsstrahlung terms. XSTAR's msolvelud/msolvelucy (the sums over the populations), comp2, cmpfnc,
# freef, bremem and heatf.
#
# The rates of ucalc come with the energies ans3, ans4 (the radiative energy: the photons absorbed and emitted, "total" channel)
# and ans5, ans6 (the energy of the electrons: measured from the threshold or the level, "kinetic" channel). Each record puts
# ans4 at the lower level of its transition and -ans3 at the upper, and ans6 and -ans5 in the same way; a positive entry
# multiplied by the population of its level is a cooling, a negative one a heating.

const compton_ee_min = 1e-4             # cmpfnc: below this energy (in electron rest masses) its limit 4 sx - ee is used
const ff_absorption_coeff = Float64(2.614f-37)       # freef: the free-free absorption coefficient (Gaunt factor 1; T in 10⁴ K)
const ff_emission_coeff = Float64(1.032f-13)         # bremem: the bremsstrahlung emissivity coefficient
const ion_density_factor = Float64(1.4f0)            # both: the density of the ions, in electron densities (n_i Z² = 1.4 nₑ)
const compton_sx_floor = 1e-10                       # comp2: added to kT (eV) in the ratio of the rest energy to it
const heating_floor = 1e-37                          # heatf: added to heating + cooling in the relative imbalance

"""
    ComptonTable(sx, e, d)
    load_compton(path_or_io)

The table of the Compton heating (XSTAR's `coheat.dat`) that `compton_integrals` interpolates: `d[i, j]` at the electron temperature
`sx[i]` (in electron rest masses) and the photon energy `e[j]` (in electron rest masses).
"""
struct ComptonTable
    sx::Vector{Float64}
    e::Vector{Float64}
    d::Matrix{Float64}
end

function load_compton(io::IO)
    entries = [split(line) for line in eachline(io) if !isempty(strip(line))]
    n = maximum(parse(Int, row[1]) for row in entries)
    m = maximum(parse(Int, row[2]) for row in entries)
    table = ComptonTable(zeros(n), zeros(m), zeros(n, m))
    for row in entries
        i, j = parse(Int, row[1]), parse(Int, row[2])
        table.sx[i], table.e[j], table.d[i, j] = parse(Float64, row[3]), parse(Float64, row[5]), parse(Float64, row[6])
    end
    table
end
load_compton(path::AbstractString) = open(load_compton, path)

# hunt3: the index of the table value just below x, between 1 and the length of the table
hunt_index(xx, x) = clamp(searchsortedlast(xx, x), 1, length(xx))

"""
    compton_function(table, ee, sxx)

XSTAR's `cmpfnc`: the Compton function at the photon energy `ee` and the temperature `sxx` (both in electron rest masses),
interpolated linearly in the table (its limit 4 sx - ee below an energy of `compton_ee_min`).
"""
function compton_function(table::ComptonTable, ee, sxx)
    ee > compton_ee_min || return 4sxx - ee
    mm = clamp(hunt_index(table.e, ee), 2, length(table.e))
    ll = clamp(hunt_index(table.sx, sxx), 2, length(table.sx))
    d = table.d
    ddedsx = (d[ll, mm] - d[ll - 1, mm] + d[ll, mm - 1] - d[ll - 1, mm - 1])/(2*(table.sx[ll] - table.sx[ll - 1]))
    ddede = (d[ll, mm] - d[ll, mm - 1] + d[ll - 1, mm] - d[ll - 1, mm - 1])/(2*(table.e[mm] - table.e[mm - 1]))
    ddedsx*(sxx - table.sx[ll - 1]) + ddede*(ee - table.e[mm - 1]) + d[ll - 1, mm - 1]
end

"""
    compton_integrals(table, radiation, T)

The integrals of the radiation field for the Compton heating and cooling at the temperature `T` (10⁴ K): `(cmp1, cmp2)`
of XSTAR's `comp2`. The heating is `cmp1 nₑ` and the cooling `kT cmp2 nₑ` (eV, times erg per eV).
"""
function compton_integrals(table::ComptonTable, rad::Radiation, T)
    K = constants()
    E, F = rad.E, rad.F
    ekt = T*K.kT_eV
    sxx = (ekt + compton_sx_floor)/K.electron_rest_eV       # the temperature in electron rest masses
    ee = E[1]/K.electron_rest_eV
    tmp1 = F[1]*compton_function(table, ee, sxx)
    sum1 = sum3 = 0.0
    for kl in 2:length(E)
        tmp1o, eeo = tmp1, ee
        ee = E[kl]/K.electron_rest_eV
        tmp1 = F[kl]*compton_function(table, ee, sxx)
        dE = E[kl] - E[kl - 1]
        sum1 += (tmp1 + tmp1o)*dE/2
        sum3 += (F[kl]*ee + F[kl - 1]*eeo)*dE/2
    end
    hfake = sum3*K.sigma_thomson
    cohc = -sum1*K.sigma_thomson
    (hfake, (-cohc + hfake)/ekt)
end

"""
    free_free_heating(radiation, T, nₑ)

The heating of the free-free absorption of the radiation (erg cm⁻³ s⁻¹) at the temperature `T` (10⁴ K) and the
electron density `nₑ` (XSTAR's `freef`; the opacity is not added to an array).
"""
function free_free_heating(rad::Radiation, T, nₑ)
    K = constants()
    E, F = rad.E, rad.F
    opaff = free_free_opacity(E, T, nₑ)
    heating = 0.0
    for k in 2:length(E)
        heating += (F[k]*opaff[k] + F[k - 1]*opaff[k - 1])*K.ergsev*(E[k] - E[k - 1])/2
    end
    heating
end

"""
    free_free_opacity(E, T, nₑ)

The free-free opacity (cm⁻¹) on the energy grid `E` (eV) at the temperature `T` (10⁴ K) and electron density `nₑ`: XSTAR's `freef`
adds it to the continuum opacity.
"""
function free_free_opacity(E::AbstractVector, T, nₑ)
    K = constants()
    ekt = T*K.kT_eV
    enz2 = ion_density_factor*nₑ
    [ff_absorption_coeff*nₑ*enz2/sqrt(T)/e^3*(1 - exp(-e/ekt)) for e in E]
end

"""
    bremsstrahlung_cooling(radiation, T, nₑ)

The cooling by bremsstrahlung (erg cm⁻³ s⁻¹): the emissivity `brcems` of XSTAR's `bremem` integrated over the energy grid
of `radiation` (`heatf`).
"""
function bremsstrahlung_cooling(rad::Radiation, T, nₑ)
    K = constants()
    E = rad.E
    ekt = T*K.kT_eV
    enz2 = ion_density_factor*nₑ
    emissivity(e) = ff_emission_coeff*nₑ*enz2*exp(-e/ekt)/sqrt(T)
    cooling = 0.0
    b = emissivity(E[1])
    for k in 2:length(E)
        bo, b = b, emissivity(E[k])
        cooling += (b + bo)*(E[k] - E[k - 1])*K.ergsev/2
    end
    cooling
end

# ---------------------------------------------------------------------------
# the energies of the records of an element and the heating and cooling of its populations

# ucalc's energies (ans3, ans4, ans5, ans6) of a record, from the energies of the result of `rate`. The signs are those of the
# final assignment of ucalc, which differs between the types: the energies of the line rates are returned as
# (fenergy, ienergy) for the decays and the photoexcitations, ucalc's ans3 and ans4 being their negatives; the
# collisions have the energies of the electrons only (ans5, ans6)
ucalc_energies(::AbstractRate, r) = (0.0, 0.0, 0.0, 0.0)
ucalc_energies(::Union{AtomicLine2, RadiativeAPED}, r) = (-r.fenergy, -r.ienergy, 0.0, 0.0)
ucalc_energies(::RadiativeFeDecay, r) = (r.ienergy, r.fenergy, 0.0, 0.0)           # (ucalc does not negate them)
ucalc_energies(::Union{RadiativeProb, RadiativeSuper, TwoPhotonDecay}, r) = (-r.ienergy, -r.fenergy, 0.0, 0.0)
ucalc_energies(::Union{ParPhotoIonize1, ParPhotoIonize2, ParPhotoIonize3, PhotoionizeSuper, PhotoRecombX, PhotoionizeFeKedge,
    PhotoionizeDamp}, r) = (-r.ienergy, -r.fenergy, -get(r, :ienergy2, 0.0), -get(r, :fenergy2, 0.0))
ucalc_energies(::Union{ElectronCollision, ElectronImpact1, ElectronImpact2, EffectiveCharge, CollisionHlike1, CollisionHlike2,
    CollisionProb, CollisionHeFine, CollisionHelike, CollisionLS, CollisionHelikeSat, CollisionSuper, CollisionFe19, CollisionAPED,
    CollisionIonize}, r) = (0.0, 0.0, r.ienergy, r.fenergy)

struct EnergyBuffers
    ans3::Vector{Float64}
    ans4::Vector{Float64}
    ans5::Vector{Float64}
    ans6::Vector{Float64}
end

function energy_buffers(n)
    buffers = get!(() -> EnergyBuffers(Float64[], Float64[], Float64[], Float64[]), task_local_storage(), :RadixEnergyBuffers)::EnergyBuffers
    if length(buffers.ans3) < n
        foreach(a -> resize!(a, n), (buffers.ans3, buffers.ans4, buffers.ans5, buffers.ans6))
    end
    buffers
end

function energies_of(coef::AbstractRate, cell::Cell, given::NamedTuple)
    r = rate(coef, cell; NamedTuple{balance_keywords(coef)}(given)...)
    ucalc_energies(coef, r)
end

function group_energies!(buffers::EnergyBuffers, indices, records::AbstractVector, cell, radiation, lfast, escape)
    for k in eachindex(records)
        j, coef = indices[k], records[k]
        ptmp = something(escape_of(escape, j, coef), optically_thin)
        buffers.ans3[j], buffers.ans4[j], buffers.ans5[j], buffers.ans6[j] =
            energies_of(coef, cell, (; radiation, lfast, pesc=ptmp[1] + ptmp[2], ptmp))
    end
end

"""
    element_heating(elements, cell, populations; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast)

The heating and cooling of the element of `elements` with the level `populations` in `cell`, per unit abundance
(erg cm⁻³ s⁻¹ for an abundance of 1 relative to hydrogen: the density of an ion is `cell.ntot` times its population),
as `(heating, cooling, heating2, cooling2)`. They are the sums over the records of the energy of the photons (`heating`,
`cooling`) and of the electrons (`heating2`, `cooling2`) that each exchanges at the levels of its transition, XSTAR's
`msolvelud`: the energy of a record at a level is `energy × population`, a heating if it is negative, a cooling if positive.
"""
function element_heating(layout::Elements, cell::Cell, populations; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast)
    buffers = energy_buffers(length(layout.rates))
    for (indices, records) in layout.groups
        group_energies!(buffers, indices, records, cell, radiation, lfast, escape)
    end
    ht = cl = ht2 = cl2 = 0.0
    xpx = cell.ntot
    function add(x, value, heating, cooling)
        value > 0 ? (return (heating, cooling + x*value)) : (return (heating - x*value, cooling))
    end
    for j in eachindex(layout.rates)
        lo = layout.lo[j]
        lo == 0 && continue
        up = layout.up[j]
        ht, cl = add(populations[lo], buffers.ans4[j]*xpx, ht, cl)
        ht, cl = add(populations[up], -buffers.ans3[j]*xpx, ht, cl)
        # a photoionization that leaves an excited level of the next ion: ucalc takes the energies of its final level from stale
        # memory and its energy from the threshold of the cross section is a difference of nearly equal terms, so the energy of the
        # electrons is left out of the second channel
        k = layout.ion[j]
        layout.up[j] - layout.offset[k] > layout.nlev[k] && continue
        ht2, cl2 = add(populations[lo], buffers.ans6[j]*xpx, ht2, cl2)
        ht2, cl2 = add(populations[up], -buffers.ans5[j]*xpx, ht2, cl2)
    end
    (ht, cl, ht2, cl2)
end
