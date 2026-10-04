# Analytic fits for effective collision strengths in H-like ions [14]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= 1;
# i8= ionN ; s1= Transition

# XSTAR data type: 62

const CollisionHlike2Desc = ""

struct CollisionHlike2{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
    levels::L              # the level data of the database (a LevelTable)
end

function CollisionHlike2(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHlike2(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[4], Vector(rvec), levels)
end

"""
    rate(coef::CollisionHlike2, cell; index=false)

As `rate(::CollisionHlike1, ...)` (XSTAR ucalc type 62 runs the same code) with the second form of the fit:
the polynomial stops three coefficients before the end and `c[m-2] ln(c[m-1] τ) e^{-c[m] τ}` is added.
"""
rate(coef::CollisionHlike2, cell::Cell; index=false) =
    hlike_collision(coef, cell, index, true)
