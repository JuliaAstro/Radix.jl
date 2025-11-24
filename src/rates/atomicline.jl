#

# XSTAR data type: 4

const AtomicLineDesc = "line data radiative: mendosa; raymond and smith"

struct AtomicLine{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    j::I
    e::R
    f::R
    a::R
end

function AtomicLine(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine(Int8(rate), label, ivec[1], ivec[2], rvec[1], rvec[2], rvec[5])
end

function rate(coef::AtomicLine, cell::Cell)
    
end
