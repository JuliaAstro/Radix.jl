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
in_level_matrix(coef::AbstractRate) =
    !(coef.rtype in total_rate_types) && !(coef.rtype == photoionization_type && coef isa level_resolved_data_type)

"""
    ion_rates(records, ion)

The rate records of `records` (a loaded database) that belong to ion `ion` and enter its level balance.
"""
ion_rates(records, ion) = [r for r in records
    if r isa AbstractRate && hasproperty(r, :ion) && r.ion == ion && in_level_matrix(r)]

"""
    rate_matrix(rates, levels, ion, cell; radiation=NO_RADIATION, lfast=photoionization_lfast)

The rate matrix `A` (s⁻¹, `dx/dt = A x`) of the levels `1:nlevels(levels, ion)` of `ion` for the
gas in `cell`, from the records `rates` (see `ion_rates`). `radiation` is the spectrum seen by the
photoionization rates and the lines (none by default) and `lfast` the speed switch of the photoionization
integrals (2: photoionization and recombination, no opacities, as XSTAR's `calc_hmc_element`).
"""
function rate_matrix(rates, levels::LevelTable, ion, cell::Cell; radiation=NO_RADIATION, lfast=photoionization_lfast)
    kw = (; radiation, lfast)
    n = nlevels(levels, ion)
    A = zeros(typeof(float(cell.T)), n, n)
    for coef in rates
        r = rate(coef, cell; accepted(coef, cell, kw)...)
        i1, i2 = r.init, r.final
        (i1 > 0 && i2 > 0 && i1 <= n && i2 <= n && i1 != i2) || continue
        ans1, ans2 = ucalc_rates(coef, r)
        coef.rtype == photoionization_type && i1 != first_level && (ans1 = zero(ans1))
        lo, up = i1, i2
        if !(coef.rtype in unordered_types)
            e1, e2 = levels[(ion, i1)].E, levels[(ion, i2)].E
            e1/(energy_order_floor + e2) - 1 < energy_order_tolerance || ((lo, up) = (i2, i1))
        end
        A[up, lo] += ans1;  A[up, up] -= ans2
        A[lo, up] += ans2;  A[lo, lo] -= ans1
    end
    A
end

"""
    level_populations(A)

The populations `x` (fractions of the ion, Σ x = 1) that solve `A x = 0`: the last equation is replaced by
the normalisation, as in XSTAR.
"""
function level_populations(A::AbstractMatrix)
    n = size(A, 1)
    M = copy(A)
    M[n, :] .= 1
    b = zeros(eltype(A), n)
    b[n] = 1
    M \ b
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
