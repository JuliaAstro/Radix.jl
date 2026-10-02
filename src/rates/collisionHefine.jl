# Fits to fine-structure collision strengths for He-like ions [32]: r1−rjmax =
# coefficients; i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 66

const CollisionHeFineDesc = "Like type 69 but, data in fine structure"

struct CollisionHeFine{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    lower::I           # lower level
    upper::I           # upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
end

function CollisionHeFine(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHeFine(Int8(rate), label, ivec..., Vector(rvec))
end
