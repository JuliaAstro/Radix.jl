# Effective ion charge for i-th level of N-electron ion: r1 = Zeff; i1= n;
# i2= L; i3= 2J; i4= Z; i5= i; i6= ionN

# XSTAR data type: 57

const EffectiveChargeDesc = "effective charge to be used in coll. ion."

struct EffectiveCharge{I, R} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    J2::I
    Z::I
    i::I
    N::I
    Zeff::R
end

function EffectiveCharge(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 6 ? (ivec..., Int32(0)) : ivec
    EffectiveCharge(Int8(rate), label, iv..., rvec[1])
end

function rate(coef::EffectiveCharge)

end
