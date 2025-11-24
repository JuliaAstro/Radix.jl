# Analytic fits for effective collision strengths in H-like ions [14]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= 1;
# i8= ionN ; s1= Transition

# XSTAR data type: 60

const CollisionHlike1Desc = ""

struct CollisionHlike1{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    N::I
    c::Tuple
end

function CollisionHlike1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHlike1(Int8(rate), label, ivec[1], ivec[2], ivec[4], Tuple(rvec))
end

function rate(coef::CollisionHlike1)

end
