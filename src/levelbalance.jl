# Level populations of one ion: the rate matrix of XSTAR's calc_hmc_ion and the solution of msolvelud
# (prototype; the iterative solver msolvelucy that XSTAR uses converges to the same populations).
#
# The unknowns are the populations of the levels of the ion, the last level being the continuum
# (the next ion). For every rate record with levels (init, final) the two levels are ordered by
# energy (llo, lup) and ucalc's two rates ans1 = llo→lup and ans2 = lup→llo go into the matrix
#     A[lup, llo] += ans1    A[lup, lup] -= ans2
#     A[llo, lup] += ans2    A[llo, llo] -= ans1
# and the last row is replaced by the normalisation Σ x = 1.

const energy_order_tolerance = 1e-8     # calc_hmc_ion keeps (init, final) as (llo, lup) when init/final - 1 is below this
const energy_order_floor = 1e-24        # added to the energy of `final` in that ratio
const unordered_types = (7, 41)         # rate types (photoionization) whose levels are never reordered
const photoionization_lfast = 2         # lfpi of calc_hmc_element: photoionization and recombination, no opacities
const unconnected_tolerance = 1e-39     # msolvelud drops the levels whose rates are all below this
const probe_temperature = 1e6           # (10⁴ K) hot enough that no record's ΔE/kT cutoff hides its levels from `index=true`
const first_level = 1                   # the ground level
const photoionization_type = 1          # rate type 1 counts the photoionization of the ground level only
const total_rate_types = (8, 15)        # rate types of the totals (recombination, ionization) that XSTAR leaves out
const level_resolved_data_type = ParPhotoIonize2   # its rate type 1 records are left out as well

# ucalc's (ans1, ans2) from the NamedTuple that `rate` returns. The line rates of Radix are named by the
# direction init → final (decay first), while ucalc returns the photoexcitation first.
ucalc_rates(coef::AbstractRate, r) = (r.frate, r.irate)
ucalc_rates(coef::Union{AtomicLine2, RadiativeAPED, RadiativeFeDecay}, r) = (r.irate, r.frate)

# the keywords of `rate` that the level balance supplies to each type of record: the radiation field, the speed switch
# of the photoionization integrals and the escape probabilities (the lines take their sum, `pesc`, the others the
# pair `ptmp`). Listed by type because several rates forward their keywords (`kw...`) and would accept any.
balance_keywords(::AbstractRate) = ()
balance_keywords(::Union{AtomicLine2, RadiativeAPED, RadiativeFeDecay}) = (:radiation, :pesc)
balance_keywords(::Union{ParPhotoIonize1, ParPhotoIonize2, ParPhotoIonize3, PhotoionizeDamp, PhotoionizeFeKedge,
    PhotoRecombX, PhotoionizeSuper}) = (:radiation, :lfast, :ptmp)
balance_keywords(::PhotoionizeDelta) = (:radiation,)
balance_keywords(::RadiativeSuper) = (:ptmp,)

# the rates of the ion that go into the level matrix
in_level_matrix(coef::AbstractRate) = !(coef isa Union{AtomicLevel, AtomicLevelFe}) && !(coef.rtype in total_rate_types) &&
    !(coef.rtype == photoionization_type && coef isa level_resolved_data_type)

"""
    ion_rates(records, ion)

The rate records of `records` (a loaded database) that belong to ion `ion` and enter its level balance.
"""
ion_rates(records, ion) = [r for r in records
    if r isa AbstractRate && hasproperty(r, :ion) && r.ion == ion && in_level_matrix(r)]

"""
    element_ions(records, Z)

The ion indices of the element of atomic number `Z`, from the neutral to the last bound stage (its
`Ion` records, in order of ionization stage).
"""
element_ions(records, Z) = [r.ion for r in sort!(filter(r -> r isa Ion && r.Z == Z, records); by=r -> r.stage)]

const optically_thin = (0.5, 0.5)                          # ucalc's ptmp1, ptmp2 without line trapping

"""
    Elements(rates, levels, ions)
    Elements(records, levels, Z)

The structure of the level balance of an element: what does not depend on the gas, on the radiation or on the
position, and so is built once and shared by all cells of a grid. `ions` are the ion indices of the element from the
neutral up and `rates[k]` the records of `ions[k]` (see `ion_rates`); the second form takes them from the
database `records` (`element_ions`, `ion_rates`) for the element of atomic number `Z`.

The ions are stacked as in XSTAR: the last level of an ion (its continuum) is the first level of the
next, so the unknowns are the `nlevels - 1` levels of each ion and the bare nucleus last (`N` in all), and a
photoionization that leaves the next ion in an excited level has a `final` beyond the last level of its ion.
For each record the layout holds the matrix entries it feeds (`lo` and `up`, the global indices of its lower and
upper level, 0 if it does not enter) and whether it keeps its forward rate (`ground_only`), found from the
connectivity of the record (`index=true`) and the energies of the levels, and the keywords of `rate` that it takes.
"""
struct Elements{L<:Levels}
    levels::L
    ions::Vector{Int}
    nlev::Vector{Int}           # levels of each ion, its continuum included
    offset::Vector{Int}         # the global index of level 1 of each ion, minus 1
    N::Int                      # the number of unknowns
    rates::Vector{AbstractRate} # the records of all the ions
    ion::Vector{Int}            # position in `ions` of the ion of each record
    lo::Vector{Int}
    up::Vector{Int}
    zero_forward::Vector{Bool}  # photoionization (rate type 1) from an excited level counts no forward rate
    keywords::Vector{Tuple{Vararg{Symbol}}}
end

function Elements(rates, levels::Levels, ions)
    ions = collect(Int, ions)
    nlev = [nlevels(levels, ion) for ion in ions]
    offset = cumsum([0; nlev[1:end - 1] .- 1])
    N = sum(nlev .- 1) + 1
    probe = Cell(probe_temperature, 1.0, 1.0, 1.0)
    records, ionpos, lo, up, zero_forward, keywords = AbstractRate[], Int[], Int[], Int[], Bool[], Tuple{Vararg{Symbol}}[]
    for (k, ion) in enumerate(ions), coef in rates[k]
        connection = rate(coef, probe; index=true)
        i1, i2 = connection.init, connection.final
        l, u = 0, 0
        if i1 > 0 && i2 > 0 && i1 != i2 && offset[k] + max(i1, i2) <= N
            l, u = i1, i2
            if !(coef.rtype in unordered_types)
                if max(i1, i2) <= nlev[k]            # (the energies of the next ion's levels are not at hand)
                    e1, e2 = levels[(ion, i1)].E, levels[(ion, i2)].E
                    e1/(energy_order_floor + e2) - 1 < energy_order_tolerance || ((l, u) = (i2, i1))
                else
                    l = u = 0
                end
            end
            l > 0 && ((l, u) = (l + offset[k], u + offset[k]))
        end
        push!(records, coef); push!(ionpos, k); push!(lo, l); push!(up, u)
        push!(zero_forward, coef.rtype == photoionization_type && i1 != first_level)
        push!(keywords, balance_keywords(coef))
    end
    Elements(levels, ions, nlev, offset, N, records, ionpos, lo, up, zero_forward, keywords)
end

Elements(records, levels::Levels, Z::Integer) =
    (ions = element_ions(records, Z); Elements([ion_rates(records, ion) for ion in ions], levels, ions))

Base.size(layout::Elements) = (layout.N, layout.N)

# the escape probabilities (ptmp1, ptmp2) of record j: `escape` is nothing (the rates' own default), a function of the
# record or a vector with one pair per record of the layout
escape_of(::Nothing, j, coef) = nothing
escape_of(escape::Function, j, coef) = escape(coef)
escape_of(escape::AbstractVector, j, coef) = escape[j]

"""
    element_matrix!(A, layout, cell; radiation=NO_RADIATION, lfast=photoionization_lfast, escape=nothing)
    element_matrix(layout, cell; kw...)

The rate matrix `A` (s⁻¹, `dx/dt = A x`) of the element of `layout` for the gas in `cell`, written into `A` or
newly allocated; the unknowns are described at `Elements`. Each record that enters the matrix puts the
two rates of ucalc, `ans1` (lower to upper level) and `ans2` (upper to lower), into the four entries
    A[up, lo] += ans1    A[up, up] -= ans2    A[lo, up] += ans2    A[lo, lo] -= ans1.

- `radiation` is the radiation seen by the photoionization rates and the lines: a `Radiation` (the mean intensity
  `J` at the cell, whatever the geometry that gave it; none by default).
- `lfast` is the speed switch of the photoionization integrals (2: photoionization and recombination, no
  opacities, as XSTAR's `calc_hmc_element`).
- `escape` gives the escape probabilities `(ptmp1, ptmp2)` of the lines and the recombination
  continua: a function of the record or a vector with a pair for each record of `layout.rates`; the rates'
  default (no line trapping) when `nothing`. XSTAR obtains them from the optical depths along its ray; another
  geometry supplies its own.
"""
function element_matrix!(A::AbstractMatrix, layout::Elements, cell::Cell;
        radiation=NO_RADIATION, lfast=photoionization_lfast, escape=nothing)
    size(A) == size(layout) || throw(DimensionMismatch("A must be $(layout.N) × $(layout.N)"))
    fill!(A, 0)
    for (j, coef) in enumerate(layout.rates)
        lo = layout.lo[j]
        lo == 0 && continue
        up = layout.up[j]
        p = escape_of(escape, j, coef)
        given = p === nothing ? (; radiation, lfast) : (; radiation, lfast, pesc=p[1] + p[2], ptmp=p)
        r = rate(coef, cell; NamedTuple{filter(key -> key in keys(given), layout.keywords[j])}(given)...)
        ans1, ans2 = ucalc_rates(coef, r)
        layout.zero_forward[j] && (ans1 = zero(ans1))
        A[up, lo] += ans1;  A[up, up] -= ans2
        A[lo, up] += ans2;  A[lo, lo] -= ans1
    end
    A
end

element_matrix(layout::Elements, cell::Cell; kw...) =
    element_matrix!(zeros(typeof(float(cell.T)), size(layout)), layout, cell; kw...)

"""
    rate_matrix(rates, levels, ion, cell; kw...)

The rate matrix of the levels `1:nlevels(levels, ion)` of a single `ion` (`element_matrix` of its layout).
"""
rate_matrix(rates, levels::Levels, ion, cell::Cell; kw...) = element_matrix(Elements([rates], levels, [ion]), cell; kw...)

"""
    element_populations(layout, cells; radiation=NO_RADIATION, kw...)

The level populations (see `level_populations`) of the element of `layout` in each of the `cells`, as an
`N × length(cells)` matrix (one column per cell; the cells may be any array, the columns follow its linear
indices). The cells do not depend on each other and are solved on the threads of Julia.
`radiation` and `escape` (see `element_matrix!`) are fixed values for all the cells, or functions of the index
of the cell that return that cell's value, so that each cell has its own radiation field and line trapping.
"""
function element_populations(layout::Elements, cells::AbstractArray{<:Cell}; radiation=NO_RADIATION, escape=nothing, kw...)
    at(value::Function, i) = value(i)
    at(value, i) = value
    R = typeof(float(first(cells).T))
    x = zeros(R, layout.N, length(cells))
    Threads.@threads for i in eachindex(cells)
        A = element_matrix(layout, cells[i]; radiation=at(radiation, i), escape=at(escape, i), kw...)
        x[:, i] = level_populations(A)
    end
    x
end

"""
    ion_fractions(x, layout)

The fractions of the ions of `layout` and of the bare nucleus last, from the level populations `x` of
`element_matrix`.
"""
function ion_fractions(x, layout::Elements)
    stop = cumsum(layout.nlev .- 1)
    [sum(x[(k == 1 ? 0 : stop[k - 1]) + 1:stop[k]]) for k in eachindex(layout.ions)] |> f -> push!(f, x[end])
end

"""
    level_populations(A)

The populations `x` (fractions of the element, Σ x = 1) that solve `A x = 0`: the last equation is replaced by
the normalisation, as in XSTAR. Levels that no rate connects (`|A|` below `unconnected_tolerance` in their
row and column) are dropped, as XSTAR's `msolvelud` does, and have population 0.
"""
function level_populations(A::AbstractMatrix)
    n = size(A, 1)
    connected(i) = i == n || any(>(unconnected_tolerance), abs.(A[i, :])) || any(>(unconnected_tolerance), abs.(A[:, i]))
    use = filter(connected, 1:n)
    M = A[use, use]
    m = length(use)
    M[m, :] .= 1
    b = zeros(eltype(A), m)
    b[m] = 1
    x = zeros(eltype(A), n)
    x[use] = M \ b
    x
end

"""
    lte_populations(levels, ion, cell)

The populations of the levels of `ion` (the continuum last) in local thermodynamic equilibrium at the
temperature and electron density of `cell`: XSTAR's `levwk`.
"""
function lte_populations(levels::Levels, ion, cell::Cell)
    K = constants()
    n = nlevels(levels, ion)
    cont = levels[(ion, n)]
    kT = K.kT_eV*cell.T
    rs = K.saha_coeff*cell.nₑ*(cell.T*T_unit)^-1.5/cont.g
    x = [levels[(ion, l)].g*rs*exp(max(Float64(cont.E - levels[(ion, l)].E)/kT, 0.0)) for l in 1:n - 1]
    push!(x, 1.0)
    x ./ sum(x)
end
