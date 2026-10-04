# Charge exchange rate coefficient with H⁺ (Kingdon and Ferland): r1 = a; r2 = b; r3 = c; r4 = d;
# r5 = T_min, r6 = T_max (K, the range of the fit); r7 = ΔE/k (10⁴ K); r8 = ?;
# i1 = i (level; the first integer of the record, whatever it is); i2 = ionN (absent in one record)

# XSTAR data type: 10

const ChargeExHpDesc = "charge exchange H+ Kingdon and Ferland"

struct ChargeExHp{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    level::I  # level index (the first integer of the record)
    ion::I    # ion index (XSTAR ionN), 0 if the record has none
    a::R
    b::R
    c::R
    d::R
    T_min::R  # K
    T_max::R  # K
    E_k::R    # ΔE/k (10⁴ K)
    extra::R  # not used by ucalc
end

function ChargeExHp(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ChargeExHp(Int8(rate), label, ivec[1], length(ivec) > 1 ? ivec[2] : zero(ivec[1]), rvec...)
end

"""
    rate(coef::ChargeExHp, cell; nlev=0, index=false)

Charge exchange of an ion with H⁺ (XSTAR ucalc type 10): `frate = nHp 10⁻⁹ a T^b (1 + c e^{dT}) e^{-ΔE/kT}`
with T in 10⁴ K and `nHp = cell.ntot - cell.nₕ` the density of H⁺ (XSTAR's `xh1`: the hydrogen that is not neutral).
`irate` is 0. `init` is the level of the record and `final` the continuum level `nlev`. The range of the fit is not
applied, nor is the sum `1 + c e^{dT}` limited at 0 (unlike type 2).
"""
function rate(coef::ChargeExHp, cell::Cell; nlev=0, index=false, verbose=false)
    init, final = Int(coef.level), nlev
    index && return (; init, final, frate=0., irate=0.)
    T = cell.T
    r = cx_unit*Float64(coef.a)*T^Float64(coef.b)*(1 + Float64(coef.c)*expo(Float64(coef.d)*T))*
        expo(-Float64(coef.E_k)/T)
    (; init, final, frate=r*(cell.ntot - cell.nₕ), irate=0.)
end
