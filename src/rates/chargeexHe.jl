# 

# XSTAR data type: 9

const ChargeExHeDesc = "charge exch. H0 Kingdon and Ferland"

struct ChargeExHe{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    a::R
    b::R
    c::R
    d::R
    T1::R
    T2::R
    ΔE::R
end

function ChargeExHe(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ChargeExH0(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::ChargeExHe, cell::Cell; verbose=false)
    res = 1e-9*coef.a*min(T, 1000.0)^coef.b*(1 + coef.c*expo(coef.d*T))

    init, final, frate = ndit > 1 ? (coef.i, nlev*coef.k-1, res/6) :
        (1, nlev, 0.0)

    (; init=init, final=final, frate=frate, irate=res)
end
