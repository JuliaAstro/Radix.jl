# 

# XSTAR data type: 10

const ChargeExHpDesc = "charge exchange H+ Kingdon and Ferland"

struct ChargeExHp{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    a::R
    b::R
    c::R
    d::R
    e::R
    T1::R
    T2::R
    ΔE::R
end

function ChargeExHp(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    if length(ivec) > 1
        ChargeExHp(Int8(rate), label, ivec..., rvec...)
    else
        ChargeExHp(Int8(rate), label, Int32(0), ivec[1], rvec...)
    end
end

function rate(coef::ChargeExHp)
end