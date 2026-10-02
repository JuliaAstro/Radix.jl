# Analytic fits for effective collision strengths in He-like ions [53]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= Z;
# i8= ionN

# XSTAR data type: 68

const CollisionHelikeDesc = ""

struct CollisionHelike{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
end

function CollisionHelike(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelike(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec))
end
