# 

# XSTAR data type: 9

const ChargeExHeDesc = "charge exch. H0 Kingdon and Ferland"

struct ChargeExHe{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    level::I         # level index
    parent::Parent{I}           # recombined ion (0: not stored) and level
    a::R
    b::R
    c::R
    d::R
    T1::R            # K
    T2::R            # K
    ΔE::R            # ΔE/k (10⁴ K)
end

function ChargeExHe(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    # records with several integers give two levels (own level, level of the
    # recombined ion); a lone ionN means the ground levels, stored as level = 0
    z = zero(eltype(ivec))
    level, parent = length(ivec) > 1 ? (ivec[1], Parent(z, ivec[2])) :
        (z, Parent(z, z))
    ChargeExHe(Int8(rate), label, level, parent, rvec...)
end

"""
    rate(coef::ChargeExHe, cell; index=false, nlev=0)

He charge exchange (XSTAR ucalc type 9). The rate goes into `irate`
(`frate` is 0): `1e-9 a min(T,1000)^b (1 + c e^{dT}) 0.1 nₕ`, divided by 6 for records
with two levels. `init` and `final` are the levels (`final` counts from the first
level of the recombined ion, `nlev`, the number of levels of the ion).
"""
function rate(coef::ChargeExHe, cell::Cell; index=false, verbose=false, nlev=0)
    T = cell.T
    two = coef.level != 0
    init, final = two ? (coef.level, nlev + coef.parent.level - 1) : (1, nlev)
    index && return (; init, final, frate=0., irate=0.)
    res = 1e-9*coef.a*min(T, 1000.0)^coef.b*(1 + coef.c*expo(coef.d*T))
    irate = res*cell.nₕ*0.1/(two ? 6 : 1)
    (; init, final, frate=0., irate)
end
