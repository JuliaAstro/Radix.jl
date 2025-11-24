# Radiative transition probability Aki for the k−i transition of N-electron
# ion computed by quantum defect theory (or hydrogenic): r1 = 0.0E + 0;
# i1= i (lower level); i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 54

const RadiativeProbDesc = "h-like cij, bautista (hlike ion)"

struct RadiativeProb{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    A::R
end

function RadiativeProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeProb(Int8(rate), label, ivec..., rvec[1])
end

function rate(coef::RadiativeProb)

end
