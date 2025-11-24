# Collision strengths for Fe XIX [10]: r1 = Υ(k,i); i1= i (lower level);
# i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 81

const CollisionFe19Desc = "Bhatia Fe XIX collision strengths"

struct CollisionFe19{I,R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    Υ::R
end

function CollisionFe19(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionFe19(Int8(rate), label, ivec..., rvec[1])    
end

function rate(coef::CollisionFe19)
    
end
