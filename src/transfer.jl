# The transfer of the radiation through a slab of zones, as far as the escape of the lines and the recombination continua goes:
# XSTAR's pescl and pescv, the line-centre opacities and edge opacities that calc_emisab_ion computes from the populations, and
# stpcut, which adds `opacity × thickness` to the optical depths of each line (`tau0`) and edge (`tauc`) as the zones go by.
# The probabilities with which the photons of the next zones escape follow from those depths (calc_hmc_ion).

const pescl_pi = Float64(3.1415927f0)   # pescl's single-precision pi
const pescl_tau_wide = 1e5              # pescl: the damping-wing scale of the optical depth
const pescl_tau_thin = Float64(1f-5)    # pescl: below this optical depth the escape probability is 1/2 (of 1)
const pescl_log_factor = 0.5
const pescl_gaussian_factor = Float64(1.2f0)
const pescv_floor = Float64(1f-12)      # pescv: smallest escape probability before the factor 1/2
const line_rate_types = (4, 9)          # the lines (and two-photon decays) have a tau0 and a centre opacity
const edge_rate_type = 7                # the photoionization records have a tauc and an edge opacity

"""
    pescl(tau)

The probability that a photon of a line escapes in one direction when the optical depth to the edge of the slab in that direction
is `tau` (XSTAR's `pescl`; 1/2 for a thin line, with the complete redistribution of a Doppler core and the damping wings added
for `tau ≥ 1`).
"""
function pescl(tau)
    p = if tau < 1
        tau < pescl_tau_thin ? one(tau) : (a = 2tau; (1 - exp(-a))/a)
    else
        b = pescl_log_factor*sqrt(log(tau))/(1 + tau/pescl_tau_wide)
        1/(tau*sqrt(pescl_pi)*(pescl_gaussian_factor + b))
    end
    p/2
end

"""
    pescv(tau)

The probability that a photon of a recombination continuum escapes in one direction when the optical depth of the edge is `tau`
(XSTAR's `pescv`: `exp(-tau)/2`, not below 10⁻¹²/2).
"""
pescv(tau) = max(exp(-tau), pescv_floor)/2

"""
    OpticalDepths(mixture)

The optical depths of the lines (rate types 4 and 9, at their centres, XSTAR's `tau0`) and of the recombination edges (rate
type 7, `tauc`) of every record of every element of `mixture`, in the inward (`inward`, XSTAR's first direction) and the
outward direction, in vectors aligned with `mixture.elements[k].rates`. They start at 0 and `add_zone!` accumulates them.
"""
struct OpticalDepths
    inward::Vector{Vector{Float64}}
    outward::Vector{Vector{Float64}}
end
OpticalDepths(mixture::Mixture) = OpticalDepths([zeros(length(el.rates)) for el in mixture.elements],
    [zeros(length(el.rates)) for el in mixture.elements])

# the opacity (cm⁻¹) of record j: the centre of a line (its cross section times the density of its lower level) or the edge of a
# recombination continuum (ucalc's opakab with the populations of its two levels); 0 for the other records
function record_opacity(coef::AbstractRate, cell::Cell, density, abund, radiation, lfast)
    0.0
end
function record_opacity(coef::Union{AtomicLine2, RadiativeAPED, RadiativeFeDecay}, cell::Cell, density, abund, radiation, lfast)
    coef.rtype in line_rate_types || return 0.0
    r = rate(coef, cell; radiation, pesc=1.0)
    r.opacity*density
end
function record_opacity(coef::Union{ParPhotoIonize1, ParPhotoIonize2, ParPhotoIonize3, PhotoRecombX, PhotoionizeSuper}, cell::Cell,
        density, abund, radiation, lfast)
    r = rate(coef, cell; radiation, abund, lfast, ptmp=optically_thin)
    r.opacity
end

"""
    record_opacities(mixture, balance, T, ntot; radiation=NO_RADIATION, lfast=photoionization_lfast)

The centre opacities of the lines and the edge opacities of the recombination continua of every record of every element in the zone
whose `balance` (from `ionization_balance`) is at the temperature `T` (10⁴ K) and hydrogen density `ntot`, in vectors aligned with
the records: XSTAR's `oplin` and `opakab`, as `calc_emisab_ion` makes them. A line has the cross section at its centre times the
density `ntot × abundance × population` of its lower level; an edge the opacity of `photoionization_integrals` for the populations
(times the abundance) of the two levels of its transition.
"""
function record_opacities(mixture::Mixture, balance, T, ntot; radiation=NO_RADIATION, lfast=photoionization_lfast)
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    map(eachindex(mixture.elements)) do k
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        opacity = zeros(length(el.rates))
        for (indices, records) in el.groups
            for i in eachindex(records)
                j = indices[i]
                lo, up = el.lo[j], el.up[j]
                density = x[lo]*ntot*abundance
                opacity[j] = record_opacity(records[i], cell, density, (x[lo]*abundance, x[up]*abundance), radiation, lfast)
            end
        end
        opacity
    end
end

"""
    add_zone!(depths, opacities, Δr; direction=:inward)

Adds the opacities of a zone (`record_opacities`) times its thickness `Δr` (cm) to the optical depths in the given direction,
XSTAR's `stpcut`.
"""
function add_zone!(depths::OpticalDepths, opacities, Δr; direction::Symbol=:inward)
    tau = direction === :inward ? depths.inward : direction === :outward ? depths.outward :
        throw(ArgumentError("direction must be :inward or :outward"))
    for k in eachindex(opacities)
        tau[k] .+= opacities[k] .* Δr
    end
    depths
end

"""
    escape_probabilities(mixture, depths; cfrac=0.0)

The escape probabilities `(ptmp1, ptmp2)` of every record of every element, as the vectors that `escape=` takes (one vector of pairs
per element, in the order of `mixture.Z`): for a line (rate types 4 and 9) `pescl` of its two optical depths with the covering
fraction `cfrac`, for a recombination edge (rate type 7) `pescv` of those of its edge, and 1/2 each for the others. XSTAR takes the
depths of the zones that the radiation has crossed already (`calc_hmc_ion`).
"""
function escape_probabilities(mixture::Mixture, depths::OpticalDepths; cfrac=0.0)
    map(eachindex(mixture.elements)) do k
        el = mixture.elements[k]
        map(eachindex(el.rates)) do j
            t1, t2 = depths.inward[k][j], depths.outward[k][j]
            rtype = el.rates[j].rtype
            if rtype in line_rate_types
                (pescl(t1)*(1 - cfrac), pescl(t2)*(1 - cfrac) + 2pescl(t1 + t2)*cfrac)
            elseif rtype == edge_rate_type
                (pescv(t1)*(1 - cfrac), pescv(t2)*(1 - cfrac) + 2pescv(t1 + t2)*cfrac)
            else
                optically_thin
            end
        end
    end
end

"""
    march_zones(mixture, ntot, processes, E, L, zones; T=100.0, equilibrium=false, xee=1.0, cfrac=0.0, attenuate=true, kw...)

The zones of a slab, one after the other, for a point source of the spectrum `L` (erg s⁻¹ erg⁻¹, 10³⁸ erg/s × the XSTAR file) on the
energy grid `E`: `zones` is a vector of `(r, Δr)`, the distance of the zone from the source and its thickness (cm), each solved at
its own radius with the radiation of the incident spectrum attenuated by the continuum depth of the zones before it,
`exp(-dpthc)` (XSTAR's `trnfrc`) and mapped (`map_spectrum`), and with the escape probabilities of the lines and recombination
edges that the optical depths of those zones give. With `equilibrium=true` the temperature of each is the thermal equilibrium
(`thermal_equilibrium`, from `T` and the `xee` of the previous zone), otherwise `T` is kept and the electron fraction iterated
(`ionization_balance`). `attenuate=false` leaves the spectrum as it is. `processes` is a collection of `AbstractContinuum`es (`standard_processes(compton)`) whose heating, cooling and
opacity are those of the zones; `cfrac` is the covering fraction of the escape probabilities (give `Thomson` the same). Further keywords go to those functions.

After each zone the opacities are added to the depths (`stpcut`): the lines and edges to the `OpticalDepths`, the continuum
opacity (`Continuum` of the `processes`, the strongest photoionization edges of each bin) to `dpthc`.

Returns the named tuples of the zones (the result of `thermal_equilibrium` or `ionization_balance` with `heating_cooling`, and the
fields `r`, `Δr`, `radiation`, `escape`, `opacity`), the final `OpticalDepths` and the continuum depth `dpthc`. The diffuse emission is
not added to the radiation and the lines do not enter the continuum opacity (`linopac`).
"""
function march_zones(mixture::Mixture, ntot, processes, E, L, zones; T=100.0, equilibrium=false, xee=1.0, cfrac=0.0,
        attenuate=true, lfast=photoionization_lfast, kw...)
    depths = OpticalDepths(mixture)
    continuum_gas = Continuum(mixture, processes; lfast)
    dpthc = zeros(length(E))
    results = NamedTuple[]
    for (r, Δr) in zones
        radiation = map_spectrum(point_source(E, attenuate ? L .* exp.(-dpthc) : L, r))
        escape = escape_probabilities(mixture, depths; cfrac)
        zone = if equilibrium
            thermal_equilibrium(mixture, ntot, processes; radiation, escape, T, xee, lfast, kw...)
        else
            balance = ionization_balance(mixture, T, ntot; radiation, escape, xee, lfast, kw...)
            (; heating_cooling(mixture, balance, T, ntot, processes; radiation, escape, lfast)..., T=Float64(T), balance...)
        end
        T, xee = zone.T, zone.xee
        edges = record_opacities(mixture, zone, zone.T, ntot; radiation, lfast)
        continuum = opacity(continuum_gas, zone, zone.T, ntot, radiation, edges)
        push!(results, (; zone..., r, Δr, radiation, escape, opacity=continuum))
        add_zone!(depths, edges, Δr)
        dpthc .+= continuum.total .* Δr
    end
    (; zones=results, depths, dpthc)
end
