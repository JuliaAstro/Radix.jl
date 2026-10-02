# Effective ion charge for i-th level of N-electron ion: r1 = Zeff; i1= n;
# i2= L; i3= 2J; i4= Z; i5= i; i6= ionN

# XSTAR data type: 57

const EffectiveChargeDesc = "effective charge to be used in coll. ion."

struct EffectiveCharge{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I      # principal quantum number
    L::I      # orbital angular momentum
    twoJ::I   # 2J
    Z::I      # atomic number
    level::I  # level index
    ion::I    # ion index (XSTAR ionN)
    Zeff::R   # effective charge
end

function EffectiveCharge(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 6 ? (ivec..., Int32(0)) : ivec
    EffectiveCharge(Int8(rate), label, iv..., rvec[1])
end
