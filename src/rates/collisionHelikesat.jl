# Fit to effective collision strengths for satellite levels of He-like ions
# [47]: r1−rj7 = fit coefficients; i1= i (lower level); i2= j (upper level);
# i3= Z; i4= ionN

# XSTAR data type: 73

const CollisionHelikeSatDesc =  "Fit to coll. strengths satellite lvls Helike ion"

struct CollisionHelikeSat{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
end

function CollisionHelikeSat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelikeSat(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec))   
end
