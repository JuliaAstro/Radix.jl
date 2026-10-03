#

# XSTAR data type: 4

const AtomicLineDesc = "line data radiative: mendosa; raymond and smith"

struct AtomicLine{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    λ::R      # wavelength
    f::R      # oscillator strength
    A::R      # Einstein A (s⁻¹)
end

function AtomicLine(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine(Int8(rate), label, Transition(ivec[1], ivec[2]), rvec[1], rvec[2], rvec[5])
end
