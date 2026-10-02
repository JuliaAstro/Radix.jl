# 

# XSTAR data type: 5

const TwoPhotonRadDesc = "2 photon transition collisional"

struct TwoPhotonRad{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    r1::R
    r2::R
end

function TwoPhotonRad(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TwoPhotonRad(Int8(rate), label, Transition(ivec[1], ivec[2]), rvec[1], rvec[2])
end
