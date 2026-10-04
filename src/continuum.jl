# The continuum opacity of a zone and the attenuation of the spectrum: XSTAR's calc_emis_all (Thomson scattering, the photoionization
# opacity of the records, free-free) with the ranking of the edges of rlbin and calc_emis_ion, and the accumulation of stpcut.

const rank_per_bin = 10                 # nrank: the edges of each energy bin that enter the continuum opacity (the strongest)
const rank_floor = 1e-37                # rlbin: edges weaker than this are not ranked
const rank_energy_floor = 0.1           # setup: the lowest edge energy (eV) of the records that are not of type 49
const damp_rate_type = 42              # the records of rate type 42 (PhotoionizeDamp) are always called
const type49_edge_factor = Float64(13.598f0)      # setup: the edge of a record of type 49 is its first tabulated energy (Ry) times this

# the edge energy (eV) of a photoionization record that xstarsetup ranks it by
rank_energy(coef::ParPhotoIonize1) = isempty(coef.E_grid) ? 0.0 : Float64(coef.E_grid[1])*type49_edge_factor
function rank_energy(coef::AbstractRate)
    hasproperty(coef, :level) || return 0.0
    lo = get(coef.levels, (coef.ion, Int(coef.level)), nothing)
    lo === nothing ? 0.0 : max(rank_energy_floor, Float64(lo.E_inf) - Float64(lo.E))
end

const opacity_edge_types = Union{ParPhotoIonize1, ParPhotoIonize2, ParPhotoIonize3, PhotoionizeSuper, PhotoRecombX}

"""
    continuum_opacity(mixture, balance, T, ntot, radiation, edges; cfrac=0.0, nrank=rank_per_bin, lfast=photoionization_lfast)

The continuum opacity (cm⁻¹) on the energy grid of `radiation` of the zone whose `balance` (from `ionization_balance`) is at the
temperature `T` (10⁴ K) and hydrogen density `ntot`, as `calc_emis_all` makes it: Thomson scattering `nₑ σ_T (1 - cfrac)`, the
photoionization opacity of the records, and the free-free opacity (the last only in `total`, not in `continuum`). `edges` are the
edge opacities of the records (`record_opacities`), which rank them: as `rlbin` does, only the `nrank` strongest edges of each energy
bin (by the edge energy of `rank_energy`) enter, and every record of rate type 42. The lines (`linopac`) are not added.

Returns `(; total, continuum, edges)`, the two arrays of the opacity and the number of records called.
"""
function continuum_opacity(mixture::Mixture, balance, T, ntot, radiation::Radiation, edges; cfrac=0.0, nrank=rank_per_bin,
        lfast=photoionization_lfast)
    K = constants()
    E = radiation.E
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    opacity = Opacity(length(E))
    thomson = balance.nₑ*K.sigma_thomson*max(0.0, 1 - cfrac)
    fill!(opacity.total, thomson); fill!(opacity.continuum, thomson)

    # the strongest `nrank` edges of each bin, in the order of the records (a later record of equal strength ranks lower)
    ranked = Dict{Int, Vector{Tuple{Float64, Int, Int}}}()
    for k in eachindex(mixture.elements), (indices, records) in mixture.elements[k].groups, i in eachindex(records)
        coef = records[i]
        (coef.rtype == edge_rate_type && coef isa opacity_edge_types) || continue
        j = indices[i]
        strength = edges[k][j]
        strength < rank_floor && continue
        eth = rank_energy(coef)
        (E[1] < eth < E[end]) || continue
        list = get!(ranked, nbin(radiation, eth), Tuple{Float64, Int, Int}[])
        position = findfirst(entry -> strength > entry[1], list)
        position === nothing ? length(list) < nrank && push!(list, (strength, k, j)) : (insert!(list, position, (strength, k, j));
            length(list) > nrank && pop!(list))
    end
    called = 0
    for k in eachindex(mixture.elements)
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        selected = Set(j for list in values(ranked) for (_, kk, j) in list if kk == k)
        for (indices, records) in el.groups, i in eachindex(records)
            coef, j = records[i], indices[i]
            (j in selected || (coef isa PhotoionizeDamp && coef.rtype == damp_rate_type)) || continue
            abund = (x[el.lo[j]]*abundance, x[el.up[j]]*abundance)
            rate(coef, cell; radiation, abund, lfast, ptmp=optically_thin, opacity)
            called += 1
        end
    end
    total = opacity.total .+ free_free_opacity(E, T, balance.nₑ)
    (; total, continuum=opacity.continuum, edges=called)
end
