# Analytic fits for effective collision strengths in He-like ions [53]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= Z;
# i8= ionN

# XSTAR data type: 68

const CollisionHelikeDesc = ""

struct CollisionHelike{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    c::Tuple
end

function CollisionHelike(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelike(Int8(rate), label, ivec..., Tuple(rvec))
end

function rate(coef::CollisionHelike)

end
