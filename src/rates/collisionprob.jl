# Collisional transition probability Cik for N-electron ion computed by
# quantum defect theory (or hydrogenic): i1= 1; i2= i (lower level);
# i3= k (upper level); i4= Z; i5= ionN

# XSTAR data type: 63

const CollisionProbDesc = "h-like cij, bautista (hlike ion)"

struct CollisionProb{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
end

function CollisionProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionProb(Int8(rate), label, ivec[2:end]...)
end

function rate(coef::CollisionProb)

end
