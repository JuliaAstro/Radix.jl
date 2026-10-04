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
const first_level = 1                   # the ground level
const photoionization_type = 1          # rate type 1 counts the photoionization of the ground level only
const total_rate_types = (8, 15)        # rate types of the totals (recombination, ionization) that XSTAR leaves out
const level_resolved_data_type = ParPhotoIonize2   # its rate type 1 records are left out as well

# ucalc's (ans1, ans2) from the NamedTuple that `rate` returns. The line rates of Radix are named by the
# direction init → final (decay first), while ucalc returns the photoexcitation first.
ucalc_rates(coef::AbstractRate, r) = (r.frate, r.irate)
ucalc_rates(coef::Union{AtomicLine2, RadiativeAPED, RadiativeFeDecay}, r) = (r.irate, r.frate)

# the keywords of `kw` that `rate` accepts for `coef`
accepted(coef::AbstractRate, cell, kw::NamedTuple) =
    NamedTuple{filter(k -> hasmethod(rate, Tuple{typeof(coef), typeof(cell)}, (k,)), keys(kw))}(kw)

# the rates of the ion that go into the level matrix
in_level_matrix(coef::AbstractRate) = !(coef isa AtomicLevel) && !(coef.rtype in total_rate_types) &&
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

"""
    element_matrix(rates, levels, ions, cell; radiation=NO_RADIATION, lfast=photoionization_lfast)

The rate matrix `A` (s⁻¹, `dx/dt = A x`) of an element for the gas in `cell`: `ions` are its ion indices
from the neutral up (see `element_ions`) and `rates[k]` the records of `ions[k]` (see `ion_rates`).
The ions are stacked as in XSTAR: the last level of an ion (its continuum) is the first level of the
next, so the unknowns are the `nlevels - 1` levels of each ion and the bare nucleus last, and a
photoionization that leaves the next ion in an excited level has a `final` beyond the last level of its ion.
`radiation` is the spectrum seen by the photoionization rates and the lines (none by default) and `lfast`
the speed switch of the photoionization integrals (2: photoionization and recombination, no opacities,
as XSTAR's `calc_hmc_element`).
"""
function element_matrix(rates, levels::LevelTable, ions, cell::Cell; radiation=NO_RADIATION, lfast=photoionization_lfast)
    kw = (; radiation, lfast)
    nlev = [nlevels(levels, ion) for ion in ions]
    offset = cumsum([0; nlev[1:end - 1] .- 1])      # the index of level 1 of each ion, minus 1
    N = sum(nlev .- 1) + 1
    A = zeros(typeof(float(cell.T)), N, N)
    for (k, ion) in enumerate(ions), coef in rates[k]
        r = rate(coef, cell; accepted(coef, cell, kw)...)
        i1, i2 = r.init, r.final
        (i1 > 0 && i2 > 0 && i1 != i2 && offset[k] + max(i1, i2) <= N) || continue
        ans1, ans2 = ucalc_rates(coef, r)
        coef.rtype == photoionization_type && i1 != first_level && (ans1 = zero(ans1))
        lo, up = i1, i2
        if !(coef.rtype in unordered_types)
            max(i1, i2) <= nlev[k] || continue          # (no energies of the levels of the next ion)
            e1, e2 = levels[(ion, i1)].E, levels[(ion, i2)].E
            e1/(energy_order_floor + e2) - 1 < energy_order_tolerance || ((lo, up) = (i2, i1))
        end
        lo += offset[k]; up += offset[k]
        A[up, lo] += ans1;  A[up, up] -= ans2
        A[lo, up] += ans2;  A[lo, lo] -= ans1
    end
    A
end

"""
    rate_matrix(rates, levels, ion, cell; radiation=NO_RADIATION, lfast=photoionization_lfast)

The rate matrix of the levels `1:nlevels(levels, ion)` of a single `ion` (`element_matrix` of one ion).
"""
rate_matrix(rates, levels::LevelTable, ion, cell::Cell; kw...) = element_matrix([rates], levels, [ion], cell; kw...)

"""
    ion_fractions(x, levels, ions)

The fractions of the ions `ions` of an element, and of its bare nucleus last, from the level populations `x`
of `element_matrix`.
"""
function ion_fractions(x, levels::LevelTable, ions)
    nlev = [nlevels(levels, ion) for ion in ions]
    stop = cumsum(nlev .- 1)
    [sum(x[(k == 1 ? 0 : stop[k - 1]) + 1:stop[k]]) for k in eachindex(ions)] |> f -> push!(f, x[end])
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
function lte_populations(levels::LevelTable, ion, cell::Cell)
    K = constants()
    n = nlevels(levels, ion)
    cont = levels[(ion, n)]
    kT = K.kT_eV*cell.T
    rs = K.saha_coeff*cell.nₑ*(cell.T*T_unit)^-1.5/cont.g
    x = [levels[(ion, l)].g*rs*exp(max(Float64(cont.E - levels[(ion, l)].E)/kT, 0.0)) for l in 1:n - 1]
    push!(x, 1.0)
    x ./ sum(x)
end
