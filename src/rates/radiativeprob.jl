# Radiative transition probability Aki for the k−i transition of N-electron
# ion computed by quantum defect theory (or hydrogenic): r1 = 0.0E + 0;
# i1= i (lower level); i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 54

const RadiativeProbDesc = "h-like cij, bautista (hlike ion)"

@with_levels struct RadiativeProb{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    A::R      # s⁻¹
end

function RadiativeProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeProb(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1])
end

const radiativeprob_min_dE = 1e-24     # added to the level energy difference before inverting it

"""
    rate(coef::RadiativeProb, cell; index=false)

Hydrogenic radiative decay (XSTAR ucalc type 54) between the two levels of the record, from
Gordon's formula (`anl1`): `irate` is the decay rate `A` (s⁻¹) from the upper to the lower level, from `(n, l+1)`
or `(n, l-1)` depending on the orbital quantum numbers, and `frate` is 0. The stored `A` is not used.
Nothing is returned if the two levels have the same `n`. `ienergy` is the rate times ΔE/kT and the erg per eV,
which is what `ucalc` puts in its energy output (it multiplies by the dimensionless ΔE/kT, not by ΔE);
`fenergy` is 0.
"""
function rate(coef::RadiativeProb, cell::Cell; index=false)
    levels = levels_of(coef)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    dE = abs(Float64(up.E) - Float64(lo.E))
    ni, nf = Int(up.n), Int(lo.n)
    ni == nf && return none
    ni < nf && ((ni, nf) = (nf, ni))
    li, lf = Int(up.L), Int(lo.L)
    alm, alp = anl1(ni, nf, lf, Int(coef.Z))
    A = li < lf ? alm : alp
    (; init=lo.level, final=up.level, frate=0., irate=A, fenergy=0.,
        ienergy=A*dE/(K.kT_eV*cell.T)*K.ergsev)
end
