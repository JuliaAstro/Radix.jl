# Fits to LS collision strengths for He-like ions [32]: r1−rjmax = 
# coefficients; i1= i (lower level); i2= k (upper level); i3= Z; i8= ionN

# XSTAR data type: 69

const CollisionLSDesc = "Kato & Nakazaki (1996) fit to Helike coll. strgt"

struct CollisionLS{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    c::Tuple
end

function CollisionLS(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionLS(Int8(rate), label, ivec..., Tuple(rvec))
end

function rate(coef::CollisionLS)

end
