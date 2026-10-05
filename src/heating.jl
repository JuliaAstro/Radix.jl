# Heating and cooling of a gas: the energy that the level populations exchange with the radiation and the electrons, and the
# processes that are not tied to levels (src/processes.jl) and the totals. XSTAR's msolvelud/msolvelucy (the sums over the populations).
#
# The rates of ucalc come with the energies ans3, ans4 (the radiative energy: the photons absorbed and emitted, "total" channel)
# and ans5, ans6 (the energy of the electrons: measured from the threshold or the level, "kinetic" channel). Each record puts
# ans4 at the lower level of its transition and -ans3 at the upper, and ans6 and -ans5 in the same way; a positive entry
# multiplied by the population of its level is a cooling, a negative one a heating.

const heating_floor = 1e-37                          # heatf: added to heating + cooling in the relative imbalance

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
