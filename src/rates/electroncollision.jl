# Electron-impact effective collision strength for the k−i transition of
# N-electron ion (CHIANTI fit [12, 17]): r1 = ∆E (Ryd); r2 = C; r3−r7 =
# (Υ red( j), j= 1,5) (reduced effective collision strength); i1= it
# (transition type); i2= i (lower level); i3= k (upper level); i4= Z; i5= ionN

# XSTAR data type: 51

const ElectronCollisionDesc = "op and chianti line coll rates"

struct ElectronCollision{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    kind::I             # transition type
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    ΔE::R               # Ry
    C::R
    Υ::Vector{R}  # reduced effective collision strengths
    levels::L              # the level data of the database (a LevelTable)
end

function ElectronCollision(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronCollision(Int8(rate), label, ivec[1], Transition(ivec[2], ivec[3]), ivec[4:5]..., rvec[1:2]..., Vector(rvec[3:end]), levels)
end

"""
    rate(coef::ElectronCollision, cell; index=false)

Electron-impact excitation and de-excitation from the CHIANTI fits (XSTAR ucalc
type 51). The levels come from the coefficient's level table (its `levels` field). `init` is the lower and `final` the
upper level (ordered by energy), `frate` the excitation and `irate` the
de-excitation rate (s⁻¹), related by detailed balance with the transition
energy `ΔE` of the record. The fit temperature is floored at `ΔE/50k`. Records
carry 5 knots (`splinem`) or 9 knots on x = 0, 1/8, …, 1 (`upsiln`); others are
skipped. Both rates are 0 if a level is missing from the level table or `ΔE ≤ 0`.

Υ is clamped at 0, which ucalc does not do: for 4 records the 5-point spline
dips slightly below zero (rates of about -1e-7 in XSTAR). XSTAR's matrix
assembly decides which level is lower with a ratio test on the level energies
that can swap them for nearly degenerate levels; that quirk is not reproduced.
"""
function rate(coef::ElectronCollision, cell::Cell; index=false)
    levels = coef.levels
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    eij = Float64(coef.ΔE)*K.Ry_eV                 # eV
    elin = K.hc_eVÅ/eij                            # Å
    (eij <= 0 || elin <= tiny) && return none
    n = length(coef.Υ)
    n in (5, 9) || return none
    index && return (; none..., init=lo.level, final=up.level)

    T = max(cell.T*T_unit, K.T_floor_coeff/elin)               # K
    Υ = n == 5 ?
        chianti_upsilon(coef.kind, coef.ΔE, coef.C, coef.Υ, T) :
        chianti_upsilon(coef.kind, coef.ΔE, coef.C, range(0, 1, length=9), coef.Υ, T)
    Υ = max(0., Υ)
    cji = K.collision_rate_coeff*Υ/sqrt(cell.T)/up.g
    cij = cji*up.g*expo(-eij/(K.kT_eV*cell.T))/lo.g
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate,
        fenergy=frate*eij*K.ergsev, ienergy=irate*eij*K.ergsev)
end
