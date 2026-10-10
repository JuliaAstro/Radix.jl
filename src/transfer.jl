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

# the results of `f(k)` for every element `k` of `mixture`, one task for each (the elements are independent and of similar sizes)
function elementwise(f, mixture::Mixture)
    results = Vector{Base.promote_op(f, Int)}(undef, length(mixture.elements))
    Threads.@threads for k in eachindex(results)
        results[k] = f(k)
    end
    results
end

# the opacity (cm⁻¹) of record j: the centre of a line (its cross section times the density of its lower level) or the edge of a
# recombination continuum (ucalc's opakab with the populations of its two levels); 0 for the other records
function record_opacity(coef::AbstractRate, cell::Cell, density, abund, radiation, lfast; vturb=default_turbulence)
    0.0
end
function record_opacity(coef::Union{AtomicLine2, RadiativeAPED, RadiativeFeDecay}, cell::Cell, density, abund, radiation, lfast;
        vturb=default_turbulence)
    coef.rtype in line_rate_types || return 0.0
    r = rate(coef, cell; radiation, pesc=1.0, vturb)
    r.opacity*density
end
function record_opacity(coef::Union{ParPhotoIonize1, ParPhotoIonize2, ParPhotoIonize3, PhotoRecombX, PhotoionizeSuper}, cell::Cell,
        density, abund, radiation, lfast; vturb=default_turbulence)
    r = rate(coef, cell; radiation, abund, lfast, ptmp=optically_thin)
    r.opacity
end

"""
    record_opacities(mixture, balance, T, ntot; radiation=NO_RADIATION, lfast=photoionization_lfast, vturb=default_turbulence)

The centre opacities of the lines and the edge opacities of the recombination continua of every record of every element in the zone
whose `balance` (from `ionization_balance`) is at the temperature `T` (10⁴ K) and hydrogen density `ntot`, in vectors aligned with
the records: XSTAR's `oplin` and `opakab`, as `calc_emisab_ion` makes them. A line has the cross section at its centre times the
density `ntot × abundance × population` of its lower level; an edge the opacity of `photoionization_integrals` for the populations
(times the abundance) of the two levels of its transition. `vturb` is the turbulent speed (km/s) of the lines.
"""
function record_opacities(mixture::Mixture, balance, T, ntot; radiation=NO_RADIATION, lfast=photoionization_lfast, vturb=default_turbulence)
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    elementwise(mixture) do k
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        opacity = zeros(length(el.rates))
        for (indices, records) in el.groups
            for i in eachindex(records)
                j = indices[i]
                lo, up = el.lo[j], el.up[j]
                density = x[lo]*ntot*abundance
                opacity[j] = record_opacity(records[i], cell, density, (x[lo]*abundance, x[up]*abundance), radiation, lfast; vturb)
            end
        end
        opacity
    end
end

"""
    line_emissivities(mixture, balance, T, ntot; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast, vturb=default_turbulence)

The emissivity of the lines (the records of type 4, `line_data_type`) of every element in the zone, as the pairs `(rcem1, rcem2)` (erg s⁻¹ cm⁻³ outward and inward) of
`calc_emisab_ion`, in vectors aligned with the records: the energy that the population of the upper level of the record radiates in the decay,
`-abund₂ ans3` split by the escape probabilities `escape` (`escape_probabilities`; the lines are optically thin without it). The other records have 0.
"""
function line_emissivities(mixture::Mixture, balance, T, ntot; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast,
        vturb=default_turbulence)
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    elementwise(mixture) do k
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        emissivity = fill((0.0, 0.0), length(el.rates))
        for (indices, records) in el.groups
            for i in eachindex(records)
                coef, j = records[i], indices[i]
                coef.rtype == line_data_type || continue
                lo, up = el.lo[j], el.up[j]
                abund1, abund2 = x[lo]*ntot*abundance, x[up]*ntot*abundance
                (abund1 > line_abundance_min || abund2 > line_abundance_min) || continue
                ptmp = something(escape_of(element_escape(escape, k), j, coef), optically_thin)
                ans3 = energies_of(coef, cell, (; radiation, lfast, pesc=ptmp[1] + ptmp[2], ptmp))[1]
                emissivity[j] = (-abund2*ans3*ptmp[1]/(ptmp[1] + ptmp[2]), -abund2*ans3*ptmp[2]/(ptmp[1] + ptmp[2]))
            end
        end
        emissivity
    end
end

"""
    edge_emissivities(mixture, balance, T, ntot; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast)

The emission of the recombination continua (the records of rate type 7, `edge_rate_type`) of every element in the zone as the pairs `(cemab1, cemab2)` (erg s⁻¹ cm⁻³ in the two directions) of
`calc_emisab_ion`, in vectors aligned with the records: the energy `|ans3|` that the population of the upper level radiates in the recombination, split by the escape probabilities `escape`
(`pescv` of the depths of the edges, `escape_probabilities`; without it the edges are optically thin). The other records have 0.
"""
function edge_emissivities(mixture::Mixture, balance, T, ntot; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast)
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    elementwise(mixture) do k
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        emissivity = fill((0.0, 0.0), length(el.rates))
        for (indices, records) in el.groups
            for i in eachindex(records)
                coef, j = records[i], indices[i]
                (coef.rtype == edge_rate_type && coef isa opacity_edge_types) || continue
                lo, up = el.lo[j], el.up[j]
                abund1, abund2 = x[lo]*abundance, x[up]*abundance
                (abund1 > line_abundance_min || abund2 > line_abundance_min) || continue
                ptmp = something(escape_of(element_escape(escape, k), j, coef), optically_thin)
                ans3 = energies_of(coef, cell, (; radiation, lfast, pesc=ptmp[1] + ptmp[2], ptmp))[1]
                emissivity[j] = (ptmp[1]*abs(ans3)/(ptmp[1] + ptmp[2])*abund2*ntot, ptmp[2]*abs(ans3)/(ptmp[1] + ptmp[2])*abund2*ntot)
            end
        end
        emissivity
    end
end

"""
    Luminosities(mixture)

The luminosity that the lines (rate type 4) and the recombination edges (rate type 7) of every record of every element add to the spectrum along the slab (XSTAR's `elum` and `elumab`, the `emit_inward` and
`emit_outward` of `xout_lines1.fits` and `xout_rrc1.fits`, in units of 10³⁸ erg s⁻¹), as vectors of pairs aligned with `mixture.elements[k].rates`. `add_luminosities!` adds the zones.
"""
struct Luminosities
    inward::Vector{Vector{Float64}}
    outward::Vector{Vector{Float64}}
end
Luminosities(mixture::Mixture) = Luminosities([zeros(length(el.rates)) for el in mixture.elements], [zeros(length(el.rates)) for el in mixture.elements])

"""
    add_luminosities!(luminosities, line_emissivities, edge_emissivities, r, Δr)

Adds a zone at the distance `r` and of thickness `Δr` (cm) to the `Luminosities`, as `heatt` does: a line adds `rcem Δr 4π r²` in each direction, an edge adds half of `cemab Δr 4π r²` to each (its
two emissivities summed), in units of 10³⁸ erg s⁻¹.
"""
function add_luminosities!(luminosities::Luminosities, lines, edges, r, Δr)
    scale = constants().fourpi*r^2*Δr/luminosity_unit
    for k in eachindex(lines)
        for j in eachindex(lines[k])
            line, edge = lines[k][j], edges[k][j]
            luminosities.inward[k][j] += line[1]*scale + (edge[1] + edge[2])*scale/2
            luminosities.outward[k][j] += line[2]*scale + (edge[1] + edge[2])*scale/2
        end
    end
    luminosities
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

const step_flux_min = Float64(1f-12)    # step: the bins whose spectrum (in these units) is below this do not limit the thickness of the zone
const step_opacity_floor = 1e-49
const step_default_emult = 0.75         # emult: the optical depth that a zone is allowed to have in the bin of the largest opacity
const step_default_taumax = 5.0         # taumax: the bins that are more opaque than this before the zone do not limit it
const step_energy_min = 1.0             # ectt: the lowest energy (eV) that limits the zone
const step_progress_min = 1e-12         # the column that a zone adds, as a fraction of the column of the slab, below which the column is taken as out of reach
const step_default_steps = 2            # numrec: nsteps of XSTAR, its largest zone is `r/numrec`

# the zones of `march_zones`: a vector of (r, Δr) or a function of the zones so far (the results, the continuum depth and the `Transmitted`) that gives the next (r, Δr) or nothing
function zone_source(zones::AbstractVector)
    index = 0
    (results, dpthc, spectrum) -> (index += 1; index <= length(zones) ? zones[index] : nothing)
end
zone_source(zones) = zones

const transmission_tau_min = 0.01       # heatt: below this optical depth of a zone the step of the flux is linear in it
const transmission_floor = 1e-49        # heatt: the opacity is not below this

"""
    Transmitted(L)

The spectrum that passes through the zones of a slab and the diffuse emission that they add: XSTAR's `zrems`. `L` is the spectrum (erg s⁻¹ erg⁻¹) of the source
of the zones, which the zones attenuate and the emission adds to; `inward` and `outward` (XSTAR's `zrems(2)` and `zrems(3)`: the emission in the direction of the first
and second escape probability of each zone, summed over the zones) and `inward_continuum` and `outward_continuum` (`zrems(4)` and `zrems(5)`, which are the emission columns of
the output of XSTAR) have the same units and start at 0. The attenuation and the emission of a zone are added by `transmit`.
"""
struct Transmitted
    L::Vector{Float64}
    inward::Vector{Float64}
    outward::Vector{Float64}
    inward_continuum::Vector{Float64}
    outward_continuum::Vector{Float64}
end
Transmitted(L::AbstractVector) = Transmitted(Float64.(L), (zeros(length(L)) for _ in 1:4)...)

"""
    transmit(spectrum::Transmitted, radiation, continuum, r, Δr; cfrac=0.0)

The `Transmitted` after a zone at the distance `r` of the source and thickness `Δr` (cm), as XSTAR's `heatt` updates `zrems`: `radiation` is the field at the zone (`point_source`
of `spectrum.L`) and `continuum` the result of `opacity(::Continuum, …)` (the opacity `total` and `continuum`, the `emissivity` of the recombination continua and the
`bremsstrahlung`). In each bin the optical depth of the zone `τ = κ Δr` gives the factor `fac = (1 - e^{-τ})/τ` (1 below 0.01) and

    L = max(0, L - (F κ - 4π (ε₁ + ε₂)) fac Δr 4π r²),   inward += 4π ε₁ fac Δr 4π r²,   outward += 4π ε₂ fac Δr 4π r²

with `ε₁ = rccemis₁ + brems (1 - cfrac)/2` and `ε₂ = rccemis₂ + brems (1 + cfrac)/2`; the last two use the factor of the continuum opacity for `inward_continuum` and
`outward_continuum`. XSTAR adds the bremsstrahlung, which is in all directions, to the emission per steradian of the recombination continua, so it is 4π too large in `L` and the sums.
"""
function transmit(spectrum::Transmitted, radiation::Radiation, continuum, r, Δr; cfrac=0.0)
    K = constants()
    fpr2 = K.fourpi*r^2
    n = length(spectrum.L)
    step(opacity) = (tau = max(transmission_floor, opacity)*Δr; tau > transmission_tau_min ? (1 - exp(-tau))/tau : one(tau))
    L, inward, outward, inward_continuum, outward_continuum = (zeros(n) for _ in 1:5)
    for kl in 1:n
        ε₁ = continuum.emissivity[1, kl] + continuum.bremsstrahlung[kl]*(1 - cfrac)/2/K.fourpi
        ε₂ = continuum.emissivity[2, kl] + continuum.bremsstrahlung[kl]*(1 + cfrac)/2/K.fourpi
        absorbed = radiation.F[kl]*max(transmission_floor, continuum.total[kl])
        fac = step(continuum.total[kl])
        L[kl] = max(0.0, spectrum.L[kl] - (absorbed - K.fourpi*(ε₁ + ε₂))*fac*Δr*fpr2)
        inward[kl] = spectrum.inward[kl] + K.fourpi*ε₁*fac*Δr*fpr2
        outward[kl] = spectrum.outward[kl] + K.fourpi*ε₂*fac*Δr*fpr2
        fac = step(continuum.continuum[kl])
        inward_continuum[kl] = spectrum.inward_continuum[kl] + K.fourpi*ε₁*fac*Δr*fpr2
        outward_continuum[kl] = spectrum.outward_continuum[kl] + K.fourpi*ε₂*fac*Δr*fpr2
    end
    Transmitted(L, inward, outward, inward_continuum, outward_continuum)
end

"""
    march_zones(mixture, ntot, processes, E, L, zones; T=100.0, equilibrium=false, xee=1.0, cfrac=0.0, attenuate=true, diffuse=true, lines=true, vturb=default_turbulence, kw...)

The zones of a slab, one after the other, for a point source of the spectrum `L` (erg s⁻¹ erg⁻¹, 10³⁸ erg/s × the XSTAR file) on the
energy grid `E`: `zones` is a vector of `(r, Δr)`, the distance of the zone from the source and its thickness (cm) (or a function of the results, `dpthc` and `Transmitted` so far that gives the next, `march_slab`), each solved at
its own radius with the radiation of the incident spectrum attenuated by the continuum depth of the zones before it,
`exp(-dpthc)` (XSTAR's `trnfrc`) and mapped (`map_spectrum`), and with the escape probabilities of the lines and recombination
edges that the optical depths of those zones give. With `equilibrium=true` the temperature of each is the thermal equilibrium
(`thermal_equilibrium`, from `T` and the `xee` of the previous zone), otherwise `T` is kept and the electron fraction iterated
(`ionization_balance`). `attenuate=false` leaves the spectrum as it is; `lines=false` leaves the lines out of the continuum opacity (`add_line!`), `luminous=false` does not add up the luminosities of the lines and edges (`Luminosities`), with the turbulent speed `vturb` (km/s). `passes` (odd, XSTAR's `npass`) repeats the march over the same zones: the even passes go from the outer edge inwards at the temperatures of the zones to give the optical depths of the lines and edges beyond each zone, which the odd passes use as the outward depths of the escape probabilities. `processes` is a collection of `AbstractContinuum` processes (`standard_processes(compton)`) whose heating, cooling and
opacity are those of the zones; `cfrac` is the covering fraction of the escape probabilities (give `Thomson` the same). Further keywords go to those functions.

After each zone the opacities are added to the depths (`stpcut`): the lines and edges to the `OpticalDepths`, the continuum
opacity (`Continuum` of the `processes`, the strongest photoionization edges and lines of each bin) to `dpthc`, and the part of it that is `continuum` (not
the lines or free-free: XSTAR's `opakcont`) to `dpthcont`, which is the depth in the transmitted spectrum of XSTAR's output.

Returns the named tuples of the zones (the result of `thermal_equilibrium` or `ionization_balance` with `heating_cooling`, and the
fields `r`, `Δr`, `radiation`, `escape`, `opacity`), the final `OpticalDepths` and the continuum depths `dpthc` and `dpthcont`. The diffuse emission is
not added to the radiation.
"""
function march_zones(mixture::Mixture, ntot, processes, E, L, zones; T=100.0, equilibrium=false, xee=1.0, cfrac=0.0,
        attenuate=true, diffuse=true, lines=true, luminous=true, vturb=default_turbulence, lfast=photoionization_lfast, passes=1, kw...)
    isodd(passes) || throw(ArgumentError("passes must be odd: the last pass goes in the same direction as the first"))
    result = sweep(mixture, ntot, processes, E, L, zones; T, equilibrium, xee, cfrac, attenuate, diffuse, lines, luminous, vturb, lfast, kw...)
    for pass in 2:passes
        grid = [(zone.r, zone.Δr) for zone in result.zones]
        if iseven(pass)
            outer = backward(mixture, result, processes, E, L; cfrac, lfast, vturb, kw...)
            result = (; result..., outer)
        else
            result = sweep(mixture, ntot, processes, E, L, grid; T, equilibrium, xee, cfrac, attenuate, diffuse, lines, luminous, vturb, lfast,
                outer=result.outer, kw...)
        end
    end
    (; result..., passes)
end

# backward: the sweep from the outer zone inwards at the temperatures and densities of the zones of `result` (XSTAR's second pass), with the
# radiation of the incident spectrum attenuated by the continuum depth before each zone, to give the optical depths of the lines and edges
# between it and the outer edge: the `outward` depths of each zone (before it is added), as the vectors of `OpticalDepths`
function backward(mixture::Mixture, result, processes, E, L; cfrac, lfast, vturb, kw...)
    depths = OpticalDepths(mixture)
    dpthc = zeros(length(E))
    before = map(result.zones) do zone
        previous = copy(dpthc)
        dpthc .+= zone.opacity.total .* zone.Δr
        previous
    end
    outer = Vector{Vector{Vector{Float64}}}(undef, length(result.zones))
    for i in reverse(eachindex(result.zones))
        zone = result.zones[i]
        outer[i] = map(copy, depths.outward)
        radiation = map_spectrum(point_source(E, L .* exp.(-before[i]), zone.r))
        escape = escape_probabilities(mixture, OpticalDepths(result.inner[i], depths.outward); cfrac)
        n = zone.ntot
        balance = ionization_balance(mixture, zone.T, n; radiation, escape, xee=zone.xee, lfast, kw...)
        edges = record_opacities(mixture, balance, zone.T, n; radiation, lfast, vturb)
        add_zone!(depths, edges, zone.Δr; direction=:outward)
    end
    outer
end

function sweep(mixture::Mixture, ntot, processes, E, L, zones; T, equilibrium, xee, cfrac, attenuate, diffuse, lines, luminous, vturb, lfast, outer=nothing, kw...)
    depths = OpticalDepths(mixture)
    inner = Vector{Vector{Float64}}[]
    continuum_gas = Continuum(mixture, processes; lfast, vturb)
    dpthc = zeros(length(E))
    dpthcont = zeros(length(E))
    spectrum = Transmitted(L)
    luminosities = Luminosities(mixture)
    results = NamedTuple[]
    next_zone = zone_source(zones)
    while true
        step = next_zone(results, dpthc, spectrum)
        step === nothing && break
        r, Δr = step
        incident = point_source(E, spectrum.L, r)
        radiation = map_spectrum(incident)
        push!(inner, map(copy, depths.inward))
        escape = escape_probabilities(mixture, outer === nothing ? depths : OpticalDepths(depths.inward, outer[length(results)+1]); cfrac)
        zone = if equilibrium
            thermal_equilibrium(mixture, t -> gas_density(ntot, t, r), processes; radiation, escape, T, xee, lfast, kw...)
        else
            n = gas_density(ntot, T, r)
            balance = ionization_balance(mixture, T, n; radiation, escape, xee, lfast, kw...)
            (; heating_cooling(mixture, balance, T, n, processes; radiation, escape, lfast)..., T=Float64(T), ntot=n, balance...)
        end
        T, xee, n = zone.T, zone.xee, zone.ntot
        edges = record_opacities(mixture, zone, zone.T, n; radiation, lfast, vturb)
        emissivities = lines ? line_emissivities(mixture, zone, zone.T, n; radiation, escape, lfast, vturb) : nothing
        if luminous
            lined = emissivities === nothing ? line_emissivities(mixture, zone, zone.T, n; radiation, escape, lfast, vturb) : emissivities
            add_luminosities!(luminosities, lined, edge_emissivities(mixture, zone, zone.T, n; radiation, escape, lfast), r, Δr)
        end
        continuum = opacity(continuum_gas, zone, zone.T, n, radiation, edges; emissivities)
        attenuate && (spectrum = transmit(spectrum, incident, diffuse ? continuum : (; continuum..., emissivity=zero(continuum.emissivity),
            bremsstrahlung=zero(continuum.bremsstrahlung)), r, Δr; cfrac))
        push!(results, (; zone..., r, Δr, radiation, escape, opacity=continuum))
        add_zone!(depths, edges, Δr)
        dpthc .+= continuum.total .* Δr
        dpthcont .+= continuum.continuum .* Δr
    end
    (; zones=results, depths, dpthc, dpthcont, spectrum, luminosities, vturb, inner, outer)
end

"""
    march_slab(mixture, ntot, processes, E, L; r, column, emult=0.75, taumax=5.0, ectt=1.0, steps=2, kw...)

The zones of a slab of the hydrogen column density `column` (cm⁻²) and the constant density `ntot` (cm⁻³) starting at the distance `r` (cm) of the source, with the thickness of the zones chosen
as XSTAR's `step` does, and each solved as in `march_zones` (to which `kw` goes). The first zone has no thickness. The thickness of each of the others is the smallest of `r/steps`,
the column that is left and `emult/κ` for the bins above `ectt` eV whose depth `dpthc` before the zone is not above `taumax` and whose spectrum
is above 10⁻¹² (in 10³⁸ erg s⁻¹ erg⁻¹), with `κ` the opacity of the zone before (lines included). The zones go on while the column of those so far is below `column`.
Returns what `march_zones` does.
"""
function march_slab(mixture::Mixture, ntot, processes, E, L; r, column, emult=step_default_emult, taumax=step_default_taumax,
        ectt=step_energy_min, steps=step_default_steps, kw...)
    position, depth, last = float(r), 0.0, nothing
    source = function (results, dpthc, spectrum)
        added = 0.0
        if last !== nothing
            position += last
            added = gas_density(ntot, results[end].T, position)*last      # (`xstar` updates the density of a power law before it adds the column: the one at the outer edge)
            depth += added
        end
        depth < column || return nothing
        # a density that falls faster than the zones grow (at most r/steps each) never reaches the column: the zones go on to densities of 10⁻³⁰⁰ and numbers that are not
        length(results) > 2 && added < step_progress_min*column && column - depth > step_progress_min*column &&
            throw(ArgumentError("the column of $column cm⁻² is out of reach: the last zone added $added cm⁻², and the zones are at most r/steps = r/$steps thick while the density falls (try a lower column or more steps)"))
        thickness = isempty(results) ? 0.0 :
            step_thickness(results[end].opacity.total, E, spectrum.L, dpthc, position, depth, gas_density(ntot, results[end].T, position), column, emult, taumax, ectt, steps)
        last = thickness
        (position, thickness)
    end
    march_zones(mixture, ntot, processes, E, L, source; kw...)
end

# step: the thickness of the next zone from the opacity of the one before
function step_thickness(opacity, E, L, dpthc, r, depth, ntot, column, emult, taumax, ectt, steps)
    thickness = min(column/ntot, r/steps)
    for kl in eachindex(E)
        limit = emult/max(opacity[kl], step_opacity_floor)
        if E[kl] > ectt && dpthc[kl] <= taumax && L[kl] > step_flux_min*luminosity_unit
            thickness = min(thickness, limit)
        end
    end
    min(thickness, (column - depth)/ntot)
end

"""
    slab_model(mixture, processes; density, column, logξ, luminosity, α=-1.0, spectrum=(E, L) -> power_law(E, α, L), ncn2=999, kw...)
    slab_model(mixture, processes; pressure, column, logξ, luminosity, ...)
    slab_model(mixture, processes; density, radexp, ...)

A slab of gas of the `density` (cm⁻³) at the source distance, constant or `density (r/r₀)^radexp` for `radexp ≠ 0` (XSTAR's `radexp`), or the constant `pressure` (dyn cm⁻², `ConstantPressure`, XSTAR's `lcpres=1`: the density of each zone is `P/(k T)` at its temperature, `logξ` is then the logarithm of the ionization parameter of the pressure `Ξ`) and hydrogen column `column` (cm⁻²) illuminated by the power-law source `power_law(E, α, luminosity)` (erg s⁻¹) on XSTAR's energy grid
of `ncn2` points (or the spectrum `spectrum(E, luminosity)`, e.g. `(E, L) -> blackbody(E, 0.5, L)`), at the distance where the ionization parameter is `10^logξ` (`source_distance`): the run of XSTAR with `spectrum=pow trad=α` that `march_slab` solves (its keywords
`T` and `xee` of the first zone, `equilibrium`, `emult`, `taumax`, `steps` (`nsteps`), `vturb`, `lines`, ... go to it). Returns what `march_slab` does and the energies `E`, the incident spectrum
`incident`, the distance `r` and the `transmitted` spectrum of XSTAR's output, `incident × exp(-dpthcont)`.
"""
function slab_model(mixture::Mixture, processes; density=nothing, pressure=nothing, radexp=0.0, column, logξ, luminosity, α=-1.0,
        spectrum=(E, L) -> power_law(E, α, L), ncn2=999, kw...)
    (density === nothing) == (pressure === nothing) && throw(ArgumentError("give the density (constant density) or the pressure (constant pressure)"))
    E = xstar_energy_grid(ncn2)
    incident = spectrum(E, luminosity)
    ntot, r = density !== nothing ? (density, source_distance(luminosity, 10^logξ, density)) :
        (ConstantPressure(pressure), source_distance_pressure(luminosity, 10^logξ, pressure))
    iszero(radexp) || density === nothing || (ntot = PowerLawDensity(density, r, radexp))
    result = march_slab(mixture, ntot, processes, E, incident; r, column, kw...)
    (; result..., E, incident, r, transmitted=incident .* exp.(-result.dpthcont))
end
