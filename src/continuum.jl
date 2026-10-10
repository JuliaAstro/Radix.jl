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
    Continuum(mixture, processes; nrank=rank_per_bin, lfast=photoionization_lfast, vturb=default_turbulence)

The continuum opacity of the gas of `mixture` (cm⁻¹), as `calc_emis_all` makes it. Its `opacity`,

    opacity(continuum, balance, T, ntot, radiation, edges; emissivities=nothing)

of the zone whose `balance` (from `ionization_balance`) is at the temperature `T` (10⁴ K) and hydrogen density `ntot`, on the energy grid of
`radiation`: the `opacity` of each of the `processes` (`standard_processes`: Thomson scattering and free-free absorption, the processes
times their densities; the scatterers (`AbstractScattering`) are in `continuum` as well as in `total`) and the photoionization opacity of
the records. `edges` are the edge opacities of the records (`record_opacities`), which rank them: as `rlbin` does, only the `nrank` strongest
edges of each energy bin (by the edge energy of `rank_energy`) enter, and every record of rate type 42. With the `emissivities` of the lines
(`line_emissivities`) the lines are put in the bins too, as `linopac` does (`add_line!`): the `nrank` strongest lines of each bin, by their emissivity, enter.

Returns `(; total, continuum, emissivity, bremsstrahlung, edges)`: the two arrays of the opacity, `emissivity` (`rccemis`, 2 × n, per steradian: the recombination
continua of the records called), `bremsstrahlung` (`brcems`: the `emissivity` of the processes, in all directions) and the number of records called.
"""
struct Continuum{M<:Mixture, P, C<:NamedTuple}
    mixture::M
    processes::P
    nrank::Int
    lfast::Int
    vturb::Float64
    widths::Dict{Tuple{Int, Int}, Float64}      # the Auger widths of the K-vacancy levels (`auger_widths`)
    catalog::C                                   # the records that the zones rank, in the order of the records (`record_catalog`)
end
Continuum(mixture::Mixture, processes; nrank=rank_per_bin, lfast=photoionization_lfast, vturb=default_turbulence) =
    Continuum(mixture, processes, nrank, lfast, Float64(vturb), auger_widths(mixture), record_catalog(mixture))

# the records that every zone ranks, found once (there are a million of them in the records, and a few hundred thousand of these):
# `edges` the photoionization records with their edge energy (k, j, eth), `lines` the lines with their wavelength (k, j, λ), `two_photon` the two-photon decays (k, j)
# and `rated` for each element the (j, always) of the records that can be called, with `always` for the records of rate type 42
function record_catalog(mixture::Mixture)
    edges, lines, two_photon = Tuple{Int, Int, Float64}[], Tuple{Int, Int, Float64}[], Tuple{Int, Int}[]
    rated = [Tuple{Int, Bool}[] for _ in mixture.elements]
    for k in eachindex(mixture.elements), (indices, records) in mixture.elements[k].groups, i in eachindex(records)
        coef, j = records[i], indices[i]
        if coef.rtype == edge_rate_type && coef isa opacity_edge_types
            push!(edges, (k, j, rank_energy(coef)))
            push!(rated[k], (j, false))
        elseif coef isa PhotoionizeDamp && coef.rtype == damp_rate_type
            push!(rated[k], (j, true))
        end
        coef.rtype == line_data_type && hasproperty(coef, :λ) && push!(lines, (k, j, abs(Float64(coef.λ))))
        coef isa Union{AtomicLine2, TwoPhotonDecay} && coef.rtype == two_photon_data_type && push!(two_photon, (k, j))
    end
    (; edges, lines, two_photon, rated)
end

# the lines in the bins: the `nrank` strongest of each bin by emissivity (rlbin with lopak = 0 over the lines of data type 4), each put in with `add_line!`
function add_lines!(arrays, continuum::Continuum, balance, T, ntot, radiation::Radiation, edges, emissivities)
    mixture, E = continuum.mixture, radiation.E
    hc = constants().hc_eVÅ_single
    λmin, λmax = hc/E[end], hc/E[1]
    ranked = Dict{Int, Vector{Tuple{Float64, Int, Int}}}()
    for (k, j, elin) in continuum.catalog.lines
        emissivity = sum(emissivities[k][j])
        (edges[k][j] < line_rank_floor && emissivity < line_rank_floor) && continue
        (λmin <= elin <= λmax) || continue
        list = get!(ranked, nbin(radiation, hc/(line_bin_floor + elin)), Tuple{Float64, Int, Int}[])
        rank!(list, (emissivity, k, j), continuum.nrank)
    end
    for bin in sort!(collect(keys(ranked))), (_, k, j) in ranked[bin]
        coef = mixture.elements[k].rates[j]
        add_line!(arrays, radiation, coef, edges[k][j], emissivities[k][j], T; vturb=continuum.vturb, widths=continuum.widths)
    end
    arrays
end

# the continua of the two-photon decays (the records of type 9), which `calc_emis_ion` calls with the lines optically thin
function add_two_photon_continua!(arrays, continuum::Continuum, balance, ntot, radiation::Radiation)
    mixture = continuum.mixture
    for (k, j) in continuum.catalog.two_photon
        el = mixture.elements[k]
        add_two_photon!(arrays, radiation, el.rates[j], balance.populations[k][el.up[j]]*ntot*mixture.abundance[k], optically_thin)
    end
    arrays
end

function opacity(continuum::Continuum, balance, T, ntot, radiation::Radiation, edges; emissivities=nothing)
    mixture, nrank, lfast, vturb = continuum.mixture, continuum.nrank, continuum.lfast, continuum.vturb
    E = radiation.E
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    arrays = Opacity(length(E))
    for process in continuum.processes
        o = opacity(process, E, T, balance.nₑ)
        arrays.total .+= o
        process isa AbstractScattering && (arrays.continuum .+= o)
    end

    # the strongest `nrank` edges of each bin, in the order of the records (a later record of equal strength ranks lower)
    ranked = Dict{Int, Vector{Tuple{Float64, Int, Int}}}()
    for (k, j, eth) in continuum.catalog.edges
        strength = edges[k][j]
        strength < rank_floor && continue
        (E[1] < eth < E[end]) || continue
        list = get!(ranked, nbin(radiation, eth), Tuple{Float64, Int, Int}[])
        rank!(list, (strength, k, j), nrank)
    end
    chosen = [falses(length(el.rates)) for el in mixture.elements]
    for list in values(ranked), (_, k, j) in list
        chosen[k][j] = true
    end
    emissivities === nothing || add_lines!(arrays, continuum, balance, T, ntot, radiation, edges, emissivities)
    add_two_photon_continua!(arrays, continuum, balance, ntot, radiation)
    called = 0
    for k in eachindex(mixture.elements)
        el, x, abundance = mixture.elements[k], balance.populations[k], mixture.abundance[k]
        for (j, always) in continuum.catalog.rated[k]
            (always || chosen[k][j]) || continue
            abund = (x[el.lo[j]]*abundance, x[el.up[j]]*abundance)
            rate(el.rates[j], cell; radiation, abund, lfast, ptmp=optically_thin, opacity=arrays)
            called += 1
        end
    end
    bremsstrahlung = sum((emissivity(process, E, T, balance.nₑ) for process in continuum.processes); init=zeros(length(E)))
    (; total=arrays.total, continuum=arrays.continuum, emissivity=arrays.emissivity, bremsstrahlung, edges=called)
end
