# Fit to effective collision strengths for satellite levels of He-like ions
# [47]: r1−rj7 = fit coefficients; i1= i (lower level); i2= j (upper level);
# i3= Z; i4= ionN

# XSTAR data type: 73

const CollisionHelikeSatDesc =  "Fit to coll. strengths satellite lvls Helike ion"

struct CollisionHelikeSat{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    c::Tuple
end

function CollisionHelikeSat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelike(Int8(rate), label, ivec..., Tuple(rvec))   
end

function rate(coef::CollisionHelikeSat)
    
end
