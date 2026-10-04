# Two-photon radiation rate for (k−i) transition of N-electron ion:
# r1 = A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= 1; i4= ionN;
#  s1= transition identifier

# XSTAR data type: 76

const TwoPhotonDecayDesc = "2 photon decay"

@with_levels struct TwoPhotonDecay{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I    # ion index (XSTAR ionN)
    A::R      # s⁻¹
end

function TwoPhotonDecay(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TwoPhotonDecay(Int8(rate), label::String, Transition(ivec[1], ivec[2]), ivec[4], rvec[1])
end

"""
    rate(coef::TwoPhotonDecay, cell; index=false)

Two-photon decay of a level (XSTAR ucalc type 76): `irate` is the stored `A` (s⁻¹), `frate` is 0, `init` the higher and
`final` the lower level by energy (both in `1:nlev`), and `ienergy` is `A` times the energy difference of the levels.
`ucalc` also adds the two-photon continuum, a spectrum `E²(E_max − E)` normalized to `A`, to its emissivity array;
that is not included.
"""
function rate(coef::TwoPhotonDecay, cell::Cell; index=false)
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)
    A = Float64(coef.A)
    dE = abs(Float64(up.E) - Float64(lo.E))
    (; init=up.level, final=lo.level, frate=0., irate=A, fenergy=0., ienergy=A*K.ergsev*dE)
end
