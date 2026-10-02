# APED collision strengths [20]: r1−rjmax = (Te( j), j= 1, jmax) (K);
# rjmax+1−rj2*max = (Υ( j), j= 1, jmax); i1= 1; i2= i (lower level); i3= k
# (upper level); i4= Z; i5= ionN

# XSTAR data type: 92

const CollisionAPEDDesc = "aped collision strengths"

struct CollisionAPED{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    T_grid::Vector{R}   # temperatures (K)
    Υ::Vector{R}  # effective collision strengths
end

function CollisionAPED(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionAPED(Int8(rate), label, Transition(ivec[2], ivec[3]), ivec[4:5]..., Vector(rvec[1:end÷2]),
        Vector(rvec[end÷2+1:end]))
end
