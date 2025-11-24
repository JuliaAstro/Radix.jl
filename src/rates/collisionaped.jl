# APED collision strengths [20]: r1−rjmax = (Te( j), j= 1, jmax) (K);
# rjmax+1−rj2*max = (Υ( j), j= 1, jmax); i1= 1; i2= i (lower level); i3= k
# (upper level); i4= Z; i5= ionN

# XSTAR data type: 92

const CollisionAPEDDesc = "aped collision strengths"

struct CollisionAPED{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    T::Tuple
    Υ::Tuple
end

function CollisionAPED(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionAPED(Int8(rate), label, ivec[2:end]..., Tuple(rvec[1:end÷2]),
        Tuple(rvec[end÷2+1:end]))
end

function rate(coef::CollisionAPED)
    
end
