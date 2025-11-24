# 

# XSTAR data type: 5

const TwoPhotonRadDesc = "2 photon transition collisional"

struct TwoPhotonRad{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    j::I
    r1::R
    r2::R
end

function TwoPhotonRad(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TwoPhoton(Int8(rate), label, ivec[1], ivec[2], rvec[1], rvec[2])
end

function rate(coef::TwoPhotonRad)
end