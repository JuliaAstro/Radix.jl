# Electron-impact effective collision strengths for the k−i transition of the
# N-electron ion (CHIANTI 2016 fit, Burgess & Tully 1992, Dere et al. 1997).
# The integers are [i, k, it, ionN] (lower and upper level, transition type,
# ion) and the reals [∆E (Ryd), gf?, C, x(1..n), Υred(1..n)], as read by
# ucalc type 98 (xstarlib/src/ucalc.f90); the comment of the XSTAR manual has the
# integers in a different order.

# XSTAR data type: 98

const ElectronImpact2Desc = "chianti2016 collisional rates"

struct ElectronImpact2{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    kind::I                     # transition type (1-6)
    ion::I                      # ion index (XSTAR ionN)
    ΔE::R                       # transition energy (Ry)
    gf::R                       # second real, not used by ucalc (likely the oscillator strength)
    C::R                        # scale parameter
    x_grid::Vector{R}           # nodes of the reduced fit, in [0, 1]
    Υ::Vector{R}                # reduced effective collision strengths at the nodes
end

function ElectronImpact2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    n = (length(rvec) - 3) ÷ 2
    ElectronImpact2(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3],
        ivec[4], rvec[1], rvec[2], rvec[3], Vector(rvec[4:3 + n]),
        Vector(rvec[4 + n:end]))
end

"""
    rate(coef::ElectronImpact2, cell; levels, index=false)

Electron-impact excitation and de-excitation from the CHIANTI 2016 fits (XSTAR
ucalc type 98), like `ElectronCollision` but with an n-point spline through the
nodes of the record. `levels` is a `level_table`. `init` is the lower and `final`
the upper level (ordered by energy), `frate` the excitation and `irate` the
de-excitation rate (s⁻¹), related by detailed balance at the record's transition
energy. The fit temperature is floored at `ΔE/50k`. Both rates are 0 if a level is
missing from `levels` or `ΔE ≤ 0`.
"""
function rate(coef::ElectronImpact2, cell::Cell; levels, index=false,
    verbose=false)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    eij = Float64(coef.ΔE)*Ry_eV                 # eV
    elin = hc_eVÅ/eij                            # Å
    (eij <= 0 || elin <= tiny) && return none
    index && return (; none..., init=lo.level, final=up.level)

    T = max(cell.T*T_unit, T_floor_coeff/elin)               # K
    Υ = chianti_upsilon(coef.kind, coef.ΔE, coef.C, coef.x_grid, coef.Υ, T)
    cji = collision_rate_coeff*Υ/sqrt(cell.T)/up.g
    cij = cji*up.g*expo(-eij/(kT_eV*cell.T))/lo.g
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate,
        fenergy=frate*eij*ergsev, ienergy=irate*eij*ergsev)
end
